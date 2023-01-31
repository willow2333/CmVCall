#!/anaconda3/envs/FEALPy/bin python3.8
# -*- coding: utf-8 -*-
# ---
# @Software: PyCharm
# @File: run.py
# @Author: willow
# @Site: 
# @Time: 1月 19, 2023
# ---

import argparse
import os
from pathlib import Path
from pileupSNV import PileupSnv
from mtDNAcaller import Filter
import time
from tqdm import tqdm

samtools='/data/qinliu/software/Anaconda3/envs/Nano/bin/samtools'
minimap2 = '/data/qinliu/software/Anaconda3/envs/Nano/bin/minimap2'
seqkit = '/data/qinliu/software/Anaconda3/envs/Nano/bin/seqkit'
Consent_correct = '/data/qinliu/software/Anaconda3/envs/Nano/bin/CONSENT-correct'



class Run():

    def ListCreate(self,inputdir):
        P = Path(inputdir)
        pathlist = [P/i[0:-6] for i in os.listdir(inputdir) if i.endswith('.fq.gz')]
        name = [i[0:-6] for i in os.listdir(inputdir) if i.endswith('.fq.gz')]
        name_path = dict(zip(name,pathlist))
        return name_path


    def run(self,inputdir,correct,a1,a2,FC,depth):
        filelist = self.ListCreate(inputdir)
        for k,v in filelist.items():
            if os.path.exists(v):
                pass
            else:
                os.system('mkdir {0}'.format(v))
            os.chdir(v)
            os.system('gunzip -c ../{0}.fq.gz |NanoFilt  -q 10 -l 1100 --maxlength 1500 > {1}_filter1.fastq'.format(k,k))
            os.system('mkdir PAF')
            print('{}:Remove NUMT Satrt!'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))))
            os.system('minimap2 -cx map-ont -t 20 --secondary=no /data/qinliu//software/RtN/Calabrese_Dayama_Smart_Numts_filt_3.fa {0}_filter1.fastq >PAF/{1}_numt.paf'.format(k,k))
            os.system("awk  -F'\t' '{{if ($12 > 30) print $1}}' PAF/{0}_numt.paf |sort|uniq >PAF/numt.id".format(k))
            os.system('minimap2 -cx map-ont -t 20 --secondary=no /data/qinliu/Work/SCU_Forensic_mtDNA/reference/rCRS_NC_012920_flank_dloop.fasta {0}_filter1.fastq >PAF/{1}_dloop.paf'.format(k,k))
            os.system("awk -F'\t' '{{if ($8 < 50 && $9 > 1173) print $1}}' PAF/{0}_dloop.paf |sort|uniq >PAF/dloop.id".format(k))
            os.system('comm -23 PAF/dloop.id PAF/numt.id >PAF/final.id')
            print('{}:Remove NUMT End!'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))))
            try:

                if correct == 'yes':
                    print('{}:Correction Satrt!'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))))
                    os.system('seqkit grep -f PAF/final.id {0}_filter1.fastq >{1}_filter.fastq'.format(k, k))
                    os.system('rm -r PAF')
                    os.system('CONSENT-correct --in {0}_filter.fastq --out {1}.fasta --type ONT -j 40'.format(k,k))
                    print('{}:Correction End!'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))))
                    os.system('less {0}.fasta|gzip > {1}_filter.fq.gz'.format(k,k))
                    os.system('minimap2 -ax map-ont -t 20 --secondary=no /data/qinliu/Work/SCU_Forensic_mtDNA/reference/rCRS_NC_012920_flank_dloop.fasta {0}_filter.fq.gz >{1}.sam'.format(k,k))
                    os.system('samtools view -Sb {0}.sam >{1}.bam'.format(k,k))
                    os.system('samtools sort -@6 -O bam -o {0}.sort.bam {1}.bam'.format(k,k))
                    os.system('samtools index {0}.sort.bam'.format(k))
                    os.system("samtools view {0}.sort.bam|awk -F'\t' '{{if ($2==0) print$1}}'|sort|uniq >for.id".format(k))
                    os.system("samtools view {0}.sort.bam|awk -F'\t' '{{if ($2==16) print$1}}'|sort|uniq >rev.id".format(k))
                    os.system('mkdir for')
                    os.system('mkdir rev')
                    os.system('seqkit grep -f for.id {0}_filter.fq.gz |gzip >for/for.fq.gz'.format(k))
                    os.system('seqkit grep -f rev.id {0}_filter.fq.gz |gzip >rev/rev.fq.gz'.format(k))
                    os.system('minimap2 -ax map-ont -t 20 --secondary=no /data/qinliu/Work/SCU_Forensic_mtDNA/reference/rCRS_NC_012920_flank_dloop.fasta for/for.fq.gz >for/for.sam')
                    os.system('samtools view -Sb for/for.sam >for/for.bam')
                    os.system('samtools sort -@6 -O bam -o for/for.sort.bam for/for.bam')
                    os.system('samtools index for/for.sort.bam')
                    os.system('minimap2 -ax map-ont -t 20 --secondary=no /data/qinliu/Work/SCU_Forensic_mtDNA/reference/rCRS_NC_012920_flank_dloop.fasta rev/rev.fq.gz >rev/rev.sam')
                    os.system('samtools view -Sb rev/rev.sam >rev/rev.bam')
                    os.system('samtools sort -@6 -O bam -o rev/rev.sort.bam rev/rev.bam')
                    os.system('samtools index rev/rev.sort.bam')
                    os.chdir('for')
                    os.system('seqkit split for.fq.gz -p 15 -O tmp')
                    os.chdir('tmp')
                    forpartlist = [i for i in os.listdir('./') if i.endswith('.fq.gz')]
                    for forpart in tqdm(forpartlist):
                        n1 = forpart[0:-6]
                        os.system('minimap2 -ax map-ont -t 20 --secondary=no /data/qinliu/Work/SCU_Forensic_mtDNA/reference/rCRS_NC_012920_flank_dloop.fasta {0} >{1}.sam'.format(forpart,n1))
                        os.system('samtools view -Sb {0}.sam >{1}.bam'.format(n1,n1))
                        os.system('samtools sort -@6 -O bam -o {0}.sort.bam {1}.bam'.format(n1,n1))
                        os.system('samtools index {0}.sort.bam'.format(n1))
                        os.system('samtools view -q 40 -b {0}.sort.bam >{1}_mapq40.bam'.format(n1,n1))
                        os.system('samtools sort -@6 -O bam -o {0}_mapq40.sort.bam {1}_mapq40.bam'.format(n1,n1))
                        os.system('samtools index {0}_mapq40.sort.bam'.format(n1))
                        os.system("less {0}.sam|grep '^@' >{1}_mapq40_uniq.sam".format(n1,n1))
                        os.system("samtools view {0}_mapq40.bam|awk -F'\t' '{{if($2==0) print $0}}' >>{1}_mapq40_uniq.sam".format(n1,n1))
                        os.system("samtools view {0}_mapq40.bam|awk -F'\t' '{{if($2==16) print $0}}' >>{1}_mapq40_uniq.sam".format(n1,n1))
                        os.system("samtools view -Sb {0}_mapq40_uniq.sam >{1}_mapq40_uniq.bam".format(n1,n1))
                        os.system('samtools sort -@6 -O bam -o {0}_mapq40_uniq.sort.bam {1}_mapq40_uniq.bam'.format(n1,n1))
                        os.system('samtools index {0}_mapq40_uniq.sort.bam'.format(n1))
                        pileupfor = PileupSnv()
                        pileupfor.run('./')
                    os.chdir('../../')
                    os.chdir('rev')
                    os.system('seqkit split rev.fq.gz -p 15 -O tmp')
                    os.chdir('tmp')
                    revpartlist = [i for i in os.listdir('./') if i.endswith('.fq.gz')]
                    for revpart in tqdm(revpartlist):
                        n2 = revpart[0:-6]
                        os.system(
                            'minimap2 -ax map-ont -t 20 --secondary=no /data/qinliu/Work/SCU_Forensic_mtDNA/reference/rCRS_NC_012920_flank_dloop.fasta {0} >{1}.sam'.format(
                                revpart, n2))
                        os.system('samtools view -Sb {0}.sam >{1}.bam'.format(n2, n2))
                        os.system('samtools sort -@6 -O bam -o {0}.sort.bam {1}.bam'.format(n2, n2))
                        os.system('samtools index {0}.sort.bam'.format(n2))
                        os.system('samtools view -q 40 -b {0}.sort.bam >{1}_mapq40.bam'.format(n2, n2))
                        os.system('samtools sort -@6 -O bam -o {0}_mapq40.sort.bam {1}_mapq40.bam'.format(n2, n2))
                        os.system('samtools index {0}_mapq40.sort.bam'.format(n2))
                        os.system("less {0}.sam|grep '^@' >{0}_mapq40_uniq.sam".format(n2))
                        os.system(
                            "samtools view {0}_mapq40.bam|awk -F'\t' '{{if($2==0) print $0}}' >>{1}_mapq40_uniq.sam".format(
                                n2, n2))
                        os.system(
                            "samtools view {0}_mapq40.bam|awk -F'\t' '{{if($2==16) print $0}}' >>{1}_mapq40_uniq.sam".format(
                                n2, n2))
                        os.system("samtools view -Sb {0}_mapq40_uniq.sam >{1}_mapq40_uniq.bam".format(n2, n2))
                        os.system(
                            'samtools sort -@6 -O bam -o {0}_mapq40_uniq.sort.bam {1}_mapq40_uniq.bam'.format(n2, n2))
                        os.system('samtools index {0}_mapq40_uniq.sort.bam'.format(n2))
                        pileupfor = PileupSnv()
                        pileupfor.run('./')
                    os.chdir('../../')
                    F = Filter()
                    print(os.getcwd())
                    F.run(Path('./'),depth,a1,a2,FC)
                    os.chdir('../')

                elif correct == 'no':
                    os.system('seqkit grep -f PAF/final.id {0}_filter1.fastq |gzip >{1}_filter.fq.gz'.format(k,k))
                    os.system('rm -r PAF')
                    os.system(
                        'minimap2 -ax map-ont -t 20 --secondary=no /data/qinliu/Work/SCU_Forensic_mtDNA/reference/rCRS_NC_012920_flank_dloop.fasta {0}_filter.fq.gz >{1}.sam'.format(
                            k, k))
                    os.system('samtools view -Sb {0}.sam >{1}.bam'.format(k, k))
                    os.system('samtools sort -@6 -O bam -o {0}.sort.bam {1}.bam'.format(k, k))
                    os.system('samtools index {0}.sort.bam'.format(k))
                    os.system(
                        "samtools view {0}.sort.bam|awk -F'\t' '{{if ($2==0) print$1}}'|sort|uniq >for.id".format(k))
                    os.system(
                        "samtools view {0}.sort.bam|awk -F'\t' '{{if ($2==16) print$1}}'|sort|uniq >rev.id".format(k))
                    os.system('mkdir for')
                    os.system('mkdir rev')
                    # os.system('c1=$(cat for.id | wc -l)')
                    # os.system('c2=$(cat rev.id | wc -l)')
                    rows = min([len(open("for.id", 'rU').readlines()), len(open("rev.id", 'rU').readlines())])
                    os.system('seqkit grep -f for.id {0}_filter.fq.gz|seqkit sample --number {1} |gzip >for/for.fq.gz;seqkit grep -f rev.id {0}_filter.fq.gz|seqkit sample --number {1} |gzip >rev/rev.fq.gz'.format(k,rows))
                    os.system(
                        'minimap2 -ax map-ont -t 20 --secondary=no /data/qinliu/Work/SCU_Forensic_mtDNA/reference/rCRS_NC_012920_flank_dloop.fasta for/for.fq.gz >for/for.sam')
                    os.system('samtools view -Sb for/for.sam >for/for.bam')
                    os.system('samtools sort -@6 -O bam -o for/for.sort.bam for/for.bam')
                    os.system('samtools index for/for.sort.bam')
                    os.system(
                        'minimap2 -ax map-ont -t 20 --secondary=no /data/qinliu/Work/SCU_Forensic_mtDNA/reference/rCRS_NC_012920_flank_dloop.fasta rev/rev.fq.gz >rev/rev.sam')
                    os.system('samtools view -Sb rev/rev.sam >rev/rev.bam')
                    os.system('samtools sort -@6 -O bam -o rev/rev.sort.bam rev/rev.bam')
                    os.system('samtools index rev/rev.sort.bam')
                    os.chdir('for')
                    os.system('seqkit split for.fq.gz -p 15 -O tmp')
                    os.chdir('tmp')
                    forpartlist = [i for i in os.listdir('./') if i.endswith('.fq.gz')]
                    for forpart in tqdm(forpartlist):
                        n1 = forpart[0:-6]
                        os.system(
                            'minimap2 -ax map-ont -t 20 --secondary=no /data/qinliu/Work/SCU_Forensic_mtDNA/reference/rCRS_NC_012920_flank_dloop.fasta {0} >{1}.sam'.format(
                                forpart, n1))
                        os.system('samtools view -Sb {0}.sam >{1}.bam'.format(n1, n1))
                        os.system('samtools sort -@6 -O bam -o {0}.sort.bam {1}.bam'.format(n1, n1))
                        os.system('samtools index {0}.sort.bam'.format(n1))
                        os.system('samtools view -q 40 -b {0}.sort.bam >{1}_mapq40.bam'.format(n1, n1))
                        os.system('samtools sort -@6 -O bam -o {0}_mapq40.sort.bam {1}_mapq40.bam'.format(n1, n1))
                        os.system('samtools index {0}_mapq40.sort.bam'.format(n1))
                        os.system("less {0}.sam|grep '^@' >{0}_mapq40_uniq.sam".format(n1))
                        os.system(
                            "samtools view {0}_mapq40.bam|awk -F'\t' '{{if($2==0) print $0}}' >>{1}_mapq40_uniq.sam".format(
                                n1, n1))
                        os.system(
                            "samtools view {0}_mapq40.bam|awk -F'\t' '{{if($2==16) print $0}}' >>{1}_mapq40_uniq.sam".format(
                                n1, n1))
                        os.system("samtools view -Sb {0}_mapq40_uniq.sam >{1}_mapq40_uniq.bam".format(n1, n1))
                        os.system(
                            'samtools sort -@6 -O bam -o {0}_mapq40_uniq.sort.bam {1}_mapq40_uniq.bam'.format(n1, n1))
                        os.system('samtools index {0}_mapq40_uniq.sort.bam'.format(n1))
                        pileupfor = PileupSnv()
                        pileupfor.run('./')
                    os.chdir('../../')
                    os.chdir('rev')
                    os.system('seqkit split rev.fq.gz -p 15 -O tmp')
                    os.chdir('tmp')
                    revpartlist = [i for i in os.listdir('./') if i.endswith('.fq.gz')]
                    for revpart in tqdm(revpartlist):
                        n2 = revpart[0:-6]
                        os.system(
                            'minimap2 -ax map-ont -t 20 --secondary=no /data/qinliu/Work/SCU_Forensic_mtDNA/reference/rCRS_NC_012920_flank_dloop.fasta {0} >{1}.sam'.format(
                                revpart, n2))
                        os.system('samtools view -Sb {0}.sam >{1}.bam'.format(n2, n2))
                        os.system('samtools sort -@6 -O bam -o {0}.sort.bam {1}.bam'.format(n2, n2))
                        os.system('samtools index {0}.sort.bam'.format(n2))
                        os.system('samtools view -q 40 -b {0}.sort.bam >{1}_mapq40.bam'.format(n2, n2))
                        os.system('samtools sort -@6 -O bam -o {0}_mapq40.sort.bam {1}_mapq40.bam'.format(n2, n2))
                        os.system('samtools index {0}_mapq40.sort.bam'.format(n2))
                        os.system("less {0}.sam|grep '^@' >{0}_mapq40_uniq.sam".format(n2))
                        os.system(
                            "samtools view {0}_mapq40.bam|awk -F'\t' '{{if($2==0) print $0}}' >>{1}_mapq40_uniq.sam".format(
                                n2, n2))
                        os.system(
                            "samtools view {0}_mapq40.bam|awk -F'\t' '{{if($2==16) print $0}}' >>{1}_mapq40_uniq.sam".format(
                                n2, n2))
                        os.system("samtools view -Sb {0}_mapq40_uniq.sam >{1}_mapq40_uniq.bam".format(n2, n2))
                        os.system(
                            'samtools sort -@6 -O bam -o {0}_mapq40_uniq.sort.bam {1}_mapq40_uniq.bam'.format(n2, n2))
                        os.system('samtools index {0}_mapq40_uniq.sort.bam'.format(n2))
                        pileupfor = PileupSnv()
                        pileupfor.run('./')
                    os.chdir('../../')
                    F = Filter()
                    print(os.getcwd())
                    F.run(Path('./'), depth,a1,a2,FC)
                    os.chdir('../')


            except ValueError:
                print("Please choose 'yes' or 'no' in CORRECT parameter")
















if __name__ == '__main__':
    # classrun = Run()
    # classrun.run('/data/qinliu/Work/SCU_Forensic_mtDNA/CmVCallTest/correct', 'yes', 0.85,0.005,32,20)


    classrun = Run()
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help='The path of files like *fq.gz.')
    parser.add_argument("--CORRECT",default='yes', help='IF you want to correct the raw reads please choose "yes" else "no", default yes! ')
    parser.add_argument("--a1", help='The frequency of the first allele')
    parser.add_argument("--a2", default=0.005,help='The frequency of the second allele')
    parser.add_argument("--FC", help='The ratio of a1/a2')
    parser.add_argument("--depth", default=20,help='The depth of the least reads.')
    args = parser.parse_args()
    inputdir = Path(args.input)
    correct = args.CORRECT
    a1 = float(args.a1)
    a2 = float(args.a2)
    FC = float(args.FC)
    depth = float(args.depth)
    classrun.run(inputdir,correct,a1,a2,FC,depth)
