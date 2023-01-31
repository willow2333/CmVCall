#!/anaconda3/envs/FEALPy/bin python3.8
# -*- coding: utf-8 -*-
# ---
# @Software: PyCharm
# @File: pileupSNV_v4.py
# @Author: willow
# @Site: 
# @Time: 7月 13, 2022
# ---

import argparse
import os

import pandas as pd
import re
import pysam
import time
from pathlib import Path
import warnings
from multiprocessing import Pool

warnings.filterwarnings('ignore')


def DNA_complement2(sequence):
    trantab = sequence.maketrans('acgt', 'ACGT')  # trantab = str.maketrans(intab, outtab)
    string = sequence.translate(trantab)  # str.translate(trantab)
    return string


cigar_oper = {0: 'M', 1: 'I', 2: 'D', 3: 'N', 4: 'S', 5: 'H', 6: 'P', 7: '=', 8: 'X', 9: 'B'}

ref = '/data/qinliu/Work/SCU_Forensic_mtDNA/reference/rCRS_NC_012920_flank_dloop.fasta'
with open(ref, 'r') as f:
    lines = list(f.readlines()[1].strip())
    pos_ref = {k + 1: v for k, v in enumerate(lines)}
compose = r'/data/qinliu/Work/SCU_Forensic_mtDNA/reference/compose.fasta'
compose_fasta = {}
with open(compose, 'r') as h:
    lines1 = h.readlines()
    for i in range(0, len(lines1), 2):
        id = lines1[i].strip().replace('>', '')
        fa = lines1[i + 1].strip()
        compose_fasta[id] = DNA_complement2(fa)


def allelcount(alignment, base):
    finalallele1 = {}
    finalallele2 = {}
    newalignment = alignment.upper()
    match = newalignment.count('.') + newalignment.count(',')
    sitedel = newalignment.count('*')
    insertpattern = re.compile(r'(\+[0-9]+[ACGT]+)')
    delpattern = re.compile(r'(-[0-9]+[ACGT]+)')
    Apattern = re.compile(r'([^0-9])(A)')
    Gpattern = re.compile(r'([^0-9])(G)')
    Cpattern = re.compile(r'([^0-9])(C)')
    Tpattern = re.compile(r'([^0-9])(T)')
    insertlist = insertpattern.findall(newalignment)
    newinsertlist, other1 = splitindel(insertlist)
    newalignment = re.sub(insertpattern, '', newalignment)
    dellist = delpattern.findall(newalignment)
    newdellist, other2 = splitindel(dellist)
    newalignment = re.sub(delpattern, '', newalignment)
    Alist = [i[1] for i in Apattern.findall(newalignment)]
    Glist = [i[1] for i in Gpattern.findall(newalignment)]
    Clist = [i[1] for i in Cpattern.findall(newalignment)]
    Tlist = [i[1] for i in Tpattern.findall(newalignment)]
    indelmis1 = other1 + other2 + Alist + Glist + Clist + Tlist
    indelmis_set1 = set(indelmis1)
    for i in indelmis_set1:
        finalallele1[i] = indelmis1.count(i)
    finalallele1[base] = match
    finalallele1['del'] = sitedel

    indelmis2 = newinsertlist + newdellist + other1 + other2 + Alist + Glist + Clist + Tlist
    indelmis_set2 = set(indelmis2)
    for i in indelmis_set2:
        finalallele2[i] = indelmis2.count(i)
    finalallele2[base] = match
    finalallele2['del'] = sitedel

    return (finalallele1, finalallele2)


def splitindel(indelist):
    pattren = re.compile(r'[+-]([0-9]+)([ACGT]+)')
    finallist = []
    other = []
    for i in indelist:
        newlist = pattren.findall(i)
        newpattern = re.compile(r"([+-][0-9]+[ACGT]{" + newlist[0][0] + "})")
        finallist += newpattern.findall(i)
        other.append(re.sub(r"([+-][0-9]+[ACGT]{" + newlist[0][0] + "})", '', i))
    newother = ''.join(other)
    return finallist, list(newother)


def bamextract(pmlist):
    print('bam start: {}'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))))
    bamfile = pmlist[0] / pmlist[1]
    samfile = pysam.AlignmentFile(bamfile, 'r')
    id_seq_q = pd.DataFrame()
    id_seq_q['pos'] = [i for i in range(1, 1223)]
    for x in samfile:
        expandseq = {}
        name = x.query_name
        rawcigarlist = x.cigar
        sequence = x.query_sequence
        rstart = x.reference_start
        expandseq.update(dict([(i, 'na') for i in range(rstart)]))
        n = rstart
        m = 0
        if rawcigarlist[0][0] == 4:
            newsequence = sequence[rawcigarlist[0][1]:]
            cigarlist = rawcigarlist[1:]
        elif rawcigarlist[0][0] == 5:
            newsequence = sequence
            cigarlist = rawcigarlist[1:]
        else:
            newsequence = sequence
            cigarlist = rawcigarlist

        for i in range(len(cigarlist)):
            if n < 1222:
                if cigarlist[i][0] == 0:
                    refpos = range(n, n + cigarlist[i][1])
                    qrpos = range(m, m + cigarlist[i][1])
                    for j in range(len(refpos)):
                        if refpos[j] < 1222:
                            expandseq[refpos[j]] = newsequence[qrpos[j]]
                        else:
                            continue
                    n += cigarlist[i][1]
                    m += cigarlist[i][1]
                elif cigarlist[i][0] == 1:
                    expandseq[n - 1] = [expandseq[n - 1], newsequence[m:m + cigarlist[i][1]]]
                    m += cigarlist[i][1]
                elif cigarlist[i][0] == 2:
                    refpos = range(n, n + cigarlist[i][1])
                    for j in range(len(refpos)):
                        if refpos[j] < 1222:
                            expandseq[refpos[j]] = '*'
                        else:
                            continue
                    n += cigarlist[i][1]
                elif cigarlist[i][0] == 4 or cigarlist[i][0] == 5:
                    refpos = range(n, n + cigarlist[i][1])
                    for j in range(len(refpos)):
                        if refpos[j] < 1222:
                            expandseq[refpos[j]] = 'na'
                        else:
                            continue
                    n += cigarlist[i][1]
                else:
                    print(name, x.cigarstring)
            else:
                continue
        if len(expandseq) < 1222:
            for a in range(n, 1222):
                expandseq[a] = 'na'
            id_seq_q['{}_base'.format(name)] = list(expandseq.values())
        else:
            finalseq = {k: v for k, v in expandseq.items() if k < 1222}
            id_seq_q['{}_base'.format(name)] = list(finalseq.values())

    print('bam end: {}'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))))
    return id_seq_q


def Splitbam(P, bamfilelist):
    pool = Pool(len(bamfilelist))
    result = pool.map(bamextract, [[P, i] for i in bamfilelist])
    pool.close()
    pool.join()
    df = result[0]
    for i in range(1, len(result)):
        df = pd.merge(df, result[i], on='pos')
    df.to_csv(P / 'ExpandSeqQual.txt',
              sep='\t', index=False)


def stat(P):
    newdf = pd.read_csv(P / 'ExpandSeqQual.txt', sep='\t')
    print('Statistical !')
    rcrs_id = P / 'rCRS.id'
    sshs_id = P / 'SSHS.id'
    rcrs_df = pd.DataFrame()
    sshs_df = pd.DataFrame()
    rcrs = []
    sshs = []
    cols = newdf.columns
    for i in cols:
        if i != 'pos':
            id = i
            lines = newdf[i].tolist()[17:]
            seq = ''.join(lines).replace('na', '').replace("'", '').replace('[', '').replace(']', '').replace(' ',
                                                                                                              '').replace(
                '*', '').replace(',', '')
            # print(Edit_Distance(seq, compose_fasta['rCRS']))
            dif = abs(len(seq) - len(compose_fasta['rCRS']))
            if dif <20:
                print(id)
                print(dif)
            if len(seq) == len(compose_fasta['rCRS']):
                rcrs.append(id)
            elif len(seq) == len(compose_fasta['SSHS']):
                sshs.append(id)
    rcrs_df['id'] = rcrs
    sshs_df['id'] = sshs
    # rcrs_df.to_csv(rcrs_id, sep='\t', index=False, header=False)
    # sshs_df.to_csv(sshs_id, sep='\t', index=False, header=False)


import numpy as np


def Edit_Distance(word_a, word_b):
    len_a = len(word_a)
    len_b = len(word_b)
    dp = np.zeros(1222, dtype=int)
    for j in range(1, len_b + 1):
        dp[j] = j
    for i in range(1, len_a + 1):
        t1 = dp[0]
        dp[0] = dp[0] + 1
        for j in range(1, len_b + 1):
            t2 = dp[j]
            if word_a[i - 1] == word_b[j - 1]:
                dp[j] = t1
            else:
                dp[j] = min(t1, dp[j - 1], dp[j]) + 1
            t1 = t2
    # print(dp[len_b])
    return dp[len_b]


if __name__ == '__main__':
    # parser = argparse.ArgumentParser()
    # parser.add_argument("--bamdir")
    # args = parser.parse_args()
    # bamdir = args.bamdir
    # P = Path(bamdir)
    # bamfilelist = [i for i in os.listdir(bamdir) if i.endswith('_mapq40_uniq.sort.bam')]
    # newdf = Splitbam(P, bamfilelist)
    # callsnv(P, newdf)

    bamdir = r'/data/qinliu/Work/SCU_Forensic_mtDNA/compose_new_NTG/Q9/YF1071_strand/test/YF1071_FILTER_ACC_95/for/tmp'
    P = Path(bamdir)
    bamfilelist = [i for i in os.listdir(bamdir) if i.endswith('_mapq40_uniq.sort.bam')]
    # newdf = Splitbam(P, bamfilelist)
    stat(P)
