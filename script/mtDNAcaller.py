#!/anaconda3/envs/FEALPy/bin python3.8
# -*- coding: utf-8 -*-
# ---
# @Software: PyCharm
# @File: newsitesFilter.py
# @Author: willow
# @Site:
# @Time: 7月 01, 2022
# ---


import argparse

import pandas as pd
from pathlib import Path
import numpy as np
from scipy import stats
import os

#from sklearn.metrics import roc_curve, auc
#import matplotlib.pyplot as  plt

class Filter():

    def RawSitesFilter2(self,P, depth, a1, a2, FC=5):

        for_file = P / 'for/tmp/RawVcf.txt'
        rev_file = P / 'rev/tmp/RawVcf.txt'
        for_df = pd.read_csv(for_file, sep='\t')
        rev_df = pd.read_csv(rev_file, sep='\t')
        filter_for_df1 = for_df[
            (for_df['allele1'] != for_df['ref']) & (for_df['depth'] > depth) & (for_df['allele1_af'] >= a1)]
        filter_for_df1['allele2'] = filter_for_df1['allele1']
        filter_for_df2 = for_df[
            (for_df['allele1'] == for_df['ref']) & (for_df['allele2_af'] > a2) & (for_df['depth'] > depth)]

        filter_for_df3 = for_df[
            (for_df['allele1'] != for_df['ref']) & (for_df['depth'] > depth) & (for_df['allele1_af'] < a1)]

        filter_for_df = pd.concat([filter_for_df1, filter_for_df2, filter_for_df3], axis=0)

        gentypes1 = []
        for i in range(len(filter_for_df)):
            try:
                if filter_for_df['allele1_af'].tolist()[i] / filter_for_df['allele2_af'].tolist()[i] < FC:
                    if filter_for_df['allele1'].tolist()[i] != filter_for_df['allele2'].tolist()[i]:
                        gentypes1.append(
                            ''.join([filter_for_df['allele1'].tolist()[i], filter_for_df['allele2'].tolist()[i]]))
                    else:
                        gentypes1.append(filter_for_df['allele1'].tolist()[i])
                else:
                    gentypes1.append(filter_for_df['allele1'].tolist()[i])
            except:
                gentypes1.append(filter_for_df['allele1'].tolist()[i])
        filter_for_df['genotypes'] = gentypes1
        filter_for_df['genotypes'] = filter_for_df['genotypes'].apply(lambda x: ''.join(sorted(x)))

        filter_rev_df1 = rev_df[
            (rev_df['allele1'] != rev_df['ref']) & (rev_df['depth'] > depth) & (rev_df['allele1_af'] >= a1)]
        filter_rev_df1['allele2'] = filter_rev_df1['allele1']

        filter_rev_df2 = rev_df[
            (rev_df['allele1'] == rev_df['ref']) & (rev_df['allele2_af'] > a2) & (rev_df['depth'] > depth)]
        filter_rev_df3 = rev_df[
            (rev_df['allele1'] != rev_df['ref']) & (rev_df['depth'] > depth) & (rev_df['allele1_af'] < a1)]
        filter_rev_df = pd.concat([filter_rev_df1, filter_rev_df2, filter_rev_df3], axis=0)

        gentypes2 = []
        print('1 FC {}'.format(FC))
        for i in range(len(filter_rev_df)):
            try:
                if filter_rev_df['allele1_af'].tolist()[i] / filter_rev_df['allele2_af'].tolist()[i] < FC:
                    if filter_rev_df['allele1'].tolist()[i] != filter_rev_df['allele2'].tolist()[i]:
                        gentypes2.append(
                            ''.join([filter_rev_df['allele1'].tolist()[i], filter_rev_df['allele2'].tolist()[i]]))
                    else:
                        gentypes2.append(filter_rev_df['allele1'].tolist()[i])
                else:
                    gentypes2.append(filter_rev_df['allele1'].tolist()[i])
            except:
                gentypes2.append(filter_rev_df['allele1'].tolist()[i])
        filter_rev_df['genotypes'] = gentypes2
        filter_rev_df['genotypes'] = filter_rev_df['genotypes'].apply(lambda x: ''.join(sorted(x)))

        newdf = pd.merge(filter_for_df, filter_rev_df, on='pos')
        # print(newdf)
        newdf = newdf[(newdf['genotypes_x'] == newdf['genotypes_y'])]

        finaldf = pd.DataFrame()
        finaldf['pos'] = newdf['pos']
        finaldf['ref_pos'] = newdf['ref_pos_x']
        finaldf['ref'] = newdf['ref_x']
        finaldf['depth'] = newdf['depth_x'] + newdf['depth_y']
        finaldf['allele1'] = newdf['allele1_x']
        finaldf['allele1_af'] = (newdf['allele1_af_x'] * newdf['depth_x'] + newdf['allele1_af_y'] * newdf['depth_y']) / \
                                finaldf['depth']
        finaldf['allele2'] = newdf['allele2_x']
        finaldf['allele2_af'] = (newdf['allele2_af_x'] * newdf['depth_x'] + newdf['allele2_af_y'] * newdf['depth_y']) / \
                                finaldf['depth']

        gentypes = []
        for i in range(len(finaldf)):
            try:
                if finaldf['allele1_af'].tolist()[i] / finaldf['allele2_af'].tolist()[i] < FC:
                    if finaldf['allele1'].tolist()[i] != finaldf['allele2'].tolist()[i]:
                        gentypes.append(''.join([finaldf['allele1'].tolist()[i], finaldf['allele2'].tolist()[i]]))
                    else:
                        gentypes.append(finaldf['allele1'].tolist()[i])
                else:
                    gentypes.append(finaldf['allele1'].tolist()[i])
            except:
                gentypes.append(finaldf['allele1'].tolist()[i])
        finaldf['genotypes'] = gentypes
        finaldf = finaldf[finaldf['genotypes'] != finaldf['ref']]
        finaldf = finaldf.sort_values(by='pos')
        return finaldf

    def Genotypes(self,finaldf, ratios, a2, FC):
        print('1 FC {}'.format(FC))
        allele2 = []
        finalgenotypes = []
        for i in range(len(finaldf)):
            ratio = ratios[i]
            if ratio > FC:
                finalgenotypes.append(finaldf['allele1'].tolist()[i])
                allele2.append(finaldf['allele1'].tolist()[i])
            else:
                if finaldf['allele1'].tolist()[i] == finaldf['allele2'].tolist()[i]:
                    if finaldf['allele2_af'].tolist()[i] > a2:
                        if finaldf['allele1'].tolist()[i] == finaldf['ref'].tolist()[i]:
                            finalgenotypes.append(
                                ','.join([finaldf['allele1'].tolist()[i], finaldf['ref'].tolist()[i]]))
                            allele2.append(finaldf['ref'].tolist()[i])
                        else:
                            if finaldf['allele1_af'].tolist()[i] < 0.75:
                                finalgenotypes.append(
                                    ','.join([finaldf['allele1'].tolist()[i], finaldf['ref'].tolist()[i]]))
                                allele2.append(finaldf['ref'].tolist()[i])
                            else:
                                finalgenotypes.append(finaldf['allele1'].tolist()[i])
                                allele2.append(finaldf['allele1'].tolist()[i])
                    else:
                        finalgenotypes.append(finaldf['allele1'].tolist()[i])
                        allele2.append(finaldf['allele1'].tolist()[i])
                else:
                    finalgenotypes.append(','.join([finaldf['allele1'].tolist()[i], finaldf['allele2'].tolist()[i]]))
                    allele2.append(finaldf['allele2'].tolist()[i])
        return finalgenotypes, allele2

    def run(self,P, depth, a1, a2, FC):
        finaldf = self.RawSitesFilter2(P, depth, a1, a2, 5)
        finaldf = finaldf.drop_duplicates('pos')
        finaldf = finaldf.sort_values(by='pos')
        finaldf['ratios'] = finaldf['allele1_af'] / finaldf['allele2_af']
        finaldf['allele2'] = self.Genotypes(finaldf, finaldf['ratios'].tolist(), a2, FC)[1]
        finaldf['genotypes*'] = self.Genotypes(finaldf, finaldf['ratios'].tolist(), a2, FC)[0]
        finaldf = finaldf.drop(labels=['genotypes', 'ratios'], axis=1)
        finaldf.to_csv(P / 'filter_variants.txt', sep='\t', index=False)

