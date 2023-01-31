#!/anaconda3/envs/FEALPy/bin python3.8
# -*- coding: utf-8 -*-
# ---
# @Software: PyCharm
# @File: pileupSNV_v3.py
# @Author: willow
# @Site: 
# @Time: 6月 14, 2022
# ---

import argparse
import os

import pandas as pd
import re
import pysam
from tqdm import tqdm, trange
import time
from pathlib import Path
import warnings
from multiprocessing import Pool
import numpy as np
from scipy.stats import zscore
from math import log

from pyspark import SparkContext
from pyspark.sql import SQLContext
# from pyspark.pandas import spark
# os.environ["MODIN_ENGINE"] = "ray"

# from distributed import Client
# client = Client(n_workers=15)
# import ray
#
# ray.init()

#import modin.pandas as mpd

warnings.filterwarnings('ignore')

cigar_oper = {0: 'M', 1: 'I', 2: 'D', 3: 'N', 4: 'S', 5: 'H', 6: 'P', 7: '=', 8: 'X', 9: 'B'}

ref = '/data/qinliu/Work/SCU_Forensic_mtDNA/reference/rCRS_NC_012920_flank_dloop.fasta'
with open(ref, 'r') as f:
    lines = list(f.readlines()[1].strip())
    pos_ref = {k + 1: v for k, v in enumerate(lines)}

# dloop_position = [(1,626),(15974,16569)]
dloop_position = [i for i in range(15974,16569+1)] + [i for i in range(1,626+1)]
dloop_position_dict = {i+1:dloop_position[i] for i in range(len(dloop_position))}
# print(dloop_position_dict)



class PileupSnv():
    def allelcount(self,alignment, base):
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
        newinsertlist, other1 = self.splitindel(insertlist)
        newalignment = re.sub(insertpattern, '', newalignment)
        dellist = delpattern.findall(newalignment)
        newdellist, other2 = self.splitindel(dellist)
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


    def splitindel(self,indelist):
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


    def DNA_complement2(self,sequence):
        trantab = sequence.maketrans('acgt', 'ACGT')  # trantab = str.maketrans(intab, outtab)
        string = sequence.translate(trantab)  # str.translate(trantab)
        return string


    def bamextract(self,pmlist):
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


    def Splitbam(self,P, bamfilelist):
        pool = Pool(len(bamfilelist))
        result = pool.map(self.bamextract, [[P, i] for i in bamfilelist])
        pool.close()
        pool.join()
        df = result[0]
        for i in range(1, len(result)):
            df = pd.merge(df, result[i], on='pos')
        df.to_csv(P / 'ExpandSeqQual.txt',
                  sep='\t', index=False)
        # return df


    def Depth(self,P, pileupmatrix):
        print('{}:csv load start!'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))))
        newdf = pd.DataFrame()
        pos = []
        depth = []
        aligns = []
        with open(pileupmatrix, 'r') as f:
            lines = f.readlines()
            for i in tqdm(lines):
                if not i.startswith('pos'):
                    tmp = i.strip().split('\t')
                    align = [tmp[i] for i in range(1, len(tmp), 2) if tmp[i] != 'na']
                    aligns.append(align)
                    depth.append(len(align))
                    pos.append(tmp[0])
        print('{}:csv load successfully!'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))))
        newdf['pos'] = pos
        newdf['depth'] = depth
        newdf['align'] = aligns

        newdf.to_csv(P / 'SiteQualPileup.txt',
                     sep='\t', index=False)
        return newdf


    def callsnv(self,P, newdf):
        print('Start call variant')
        infodf = pd.DataFrame()
        pos = []
        ref_pos = []
        ref = []
        depth = []
        entropy = []
        allele1 = []
        allele1_af = []
        allele2 = []
        allele2_af = []

        for j in tqdm(range(len(newdf))):
            # if j == 88:
            pos.append(newdf['pos'][j])
            ref_pos.append(newdf['ref_pos'][j])
            ref.append(pos_ref[j + 1])

            baselist = newdf['align'][j]
            base_count = self.Count(baselist)
            entropy.append(self.calcShannonEnt(base_count))
            result = self.Zscore_entropy(base_count)
            allele1.append(result[0])
            allele1_af.append(result[1])

            allele2.append(result[2])
            allele2_af.append(result[3])
            depth.append(result[4])

        infodf['pos'] = pos
        infodf['ref_pos'] = ref_pos
        infodf['ref'] = ref
        infodf['depth'] = depth
        infodf['entropy'] = entropy
        infodf['allele1'] = allele1
        infodf['allele1_af'] = allele1_af

        infodf['allele2'] = allele2
        infodf['allele2_af'] = allele2_af

        infodf.to_csv(P / 'RawVcf.txt', sep='\t', index=False)

    def Count(self,baselist):
        # depth =depth
        base_count = {}

        for i in range(len(baselist)):
            b = baselist[i]
            if len(b) == 1:
                if b == '*':
                    if b not in list(base_count.keys()):
                        base_count[b] = 1

                    else:
                        base_count[b] += 1

                else:
                    if b not in list(base_count.keys()):
                        base_count[b] = 1

                    else:
                        base_count[b] += 1

            else:
                newb = ''.join(b.replace("'", '').replace('[', '').replace(']', '').split(',')).replace(' ', '')

                if newb not in list(base_count.keys()):
                    base_count[newb] = 1

                else:
                    base_count[newb] += 1

            # print(base_count)
        newbase_count = {'A': 0, 'T': 0, 'C': 0, 'G': 0, '*': 0}
        for k, v in base_count.items():
            newbase_count[k[0]] += v
        for key in list(newbase_count.keys()):
            if newbase_count[key] == 0:
                del newbase_count[key]
        # print(newbase_count)
        return newbase_count
        # else:
        #     print('Expand false !!!!')

    def callallsnv(self,P, newdf):
        print('Start call all variant')
        infodf = pd.DataFrame()
        pos = []
        ref_pos = []
        ref = []
        depth = []
        allele1 = []
        allele1_af = []
        allele2 = []
        allele2_af = []
        allele3 = []
        allele3_af = []
        allele4 = []
        allele4_af = []
        allele5 = []
        allele5_af = []

        for j in tqdm(range(len(newdf))):
            # if j == 88:
            baselist = newdf['align'][j]
            base_count = self.Count(baselist)
            if len(list(base_count.keys())) > 0:
                pos.append(newdf['pos'][j])
                ref_pos.append(newdf['ref_pos'][j])
                ref.append(pos_ref[j + 1])
                depth.append(newdf['depth'][j])
                if sum(list(base_count.values())) == newdf['depth'][j]:
                    if 'A' in base_count.keys():
                        allele1.append('A')
                        allele1_af.append(base_count['A'] / newdf['depth'][j])
                    else:
                        allele1.append('A')
                        allele1_af.append(0)
                    if 'T' in base_count.keys():
                        allele2.append('T')
                        allele2_af.append(base_count['T'] / newdf['depth'][j])
                    else:
                        allele2.append('T')
                        allele2_af.append(0)
                    if 'C' in base_count.keys():
                        allele3.append('C')
                        allele3_af.append(base_count['C'] / newdf['depth'][j])
                    else:
                        allele3.append('C')
                        allele3_af.append(0)
                    if 'G' in base_count.keys():
                        allele4.append('G')
                        allele4_af.append(base_count['G'] / newdf['depth'][j])
                    else:
                        allele4.append('G')
                        allele4_af.append(0)
                    if '*' in base_count.keys():
                        allele5.append('*')
                        allele5_af.append(base_count['*'] / newdf['depth'][j])
                    else:
                        allele5.append('*')
                        allele5_af.append(0)
                else:
                    print('This site depth is wrong !')
                    print(base_count)

        infodf['pos'] = pos
        infodf['ref_pos'] = ref_pos
        infodf['ref'] = ref
        infodf['depth'] = depth
        infodf['alleleA'] = allele1
        infodf['alleleA_af'] = allele1_af
        infodf['alleleT'] = allele2
        infodf['alleleT_af'] = allele2_af
        infodf['alleleC'] = allele3
        infodf['alleleC_af'] = allele3_af
        infodf['alleleG'] = allele4
        infodf['alleleG_af'] = allele4_af
        infodf['allele*'] = allele5
        infodf['allele*_af'] = allele5_af

        infodf.to_csv(P / 'RawallVcf.txt', sep='\t', index=False)





    def Zscore_entropy(self,base_count):
        sorted_base_count = dict(sorted(base_count.items(), key=lambda dict: dict[1], reverse=True))
        # print(sorted_base_count)
        if len(sorted_base_count) > 1:
            try:
                if sorted_base_count['*'] > 0.8 * sum(list(sorted_base_count.values())):
                    pass
                else:
                    del sorted_base_count['*']

            except:
                print('This sites is not contain deletion !')
            if len(sorted_base_count) > 1:
                allele1 = list(sorted_base_count.keys())[0]
                allele1_af = sorted_base_count[allele1] / sum(list(sorted_base_count.values()))
                allele2 = list(sorted_base_count.keys())[1]
                allele2_af = sorted_base_count[allele2] / sum(list(sorted_base_count.values()))
                return allele1, allele1_af, allele2, allele2_af, sum(list(sorted_base_count.values()))
            else:
                allele1 = list(sorted_base_count.keys())[0]
                allele1_af = sorted_base_count[allele1] / sum(list(sorted_base_count.values()))
                return allele1, allele1_af, ' ', 0, sum(list(sorted_base_count.values()))
        elif len(sorted_base_count) == 1:
            allele1 = list(sorted_base_count.keys())[0]
            allele1_af = sorted_base_count[allele1] / sum(list(sorted_base_count.values()))
            return allele1, allele1_af, ' ', 0, sum(list(sorted_base_count.values()))
        else:
            return ' ', 0, ' ', 0, sum(list(sorted_base_count.values()))


    def calcShannonEnt(self,base_count):
        if sum(list(base_count.values())) != 0:
            labelCounts = base_count
            shannonEnt = 0.0
            for key in labelCounts:
                prob = float(labelCounts[key]) / sum(list(base_count.values()))
                shannonEnt -= prob * log(prob, 2)
            return shannonEnt
        # print((shannonEnt))
        else:
            return 'na'


    def run(self,inputpath):
        P = Path(inputpath)
        bamfilelist = [i for i in os.listdir(inputpath) if i.endswith('_mapq40_uniq.sort.bam')]
        self.Splitbam(P, bamfilelist)
        newdf = self.Depth(P, P / 'ExpandSeqQual.txt')
        newdf['ref_pos'] = newdf['pos'].apply(lambda x: dloop_position_dict[int(x)])
        self.callallsnv(P, newdf)
        self.callsnv(P, newdf)

