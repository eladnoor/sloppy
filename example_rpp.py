#!/usr/bin/python
from python.optknock import OptKnock
import logging
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

CARBON_SOURCES = ['g6p',
                  'f6p',
                  '6pgc',
                  'r5p',
                  'succ',
                  'xu5p_D',
                  '2pg',
                  'ac',
                  'dhap',
                  'electrons']

SINGLE_KOS = ['', 'GLCpts','PGI','PFK','FBP','FBA','TPI','GAPD,PGK','PGM',
              'ENO','PYK','PPS','PDH','PFL','G6PDH2r,PGL','GND',
              'EDD,EDA','RPE','RPI','TKT1,TKT2','TALA','PPC','PPCK',
              'CS','ACONTa,ACONTb','ALCD2x','ACALD']

#SINGLE_KOS = ['RPE','RPI','G6PDH2r,PGL','EDD,EDA']

if __name__ == "__main__":
    logger = logging.getLogger('')
    logger.setLevel(logging.INFO)
    if not os.path.exists('res'):
        os.mkdir('res')
    
    max_ko = 2
    
    title = 'rPP'; target_reaction = 'RBC'; knockins = 'RBC,PRK'
    
    df_fname = 'res/data_%s_%d.csv' % (title, max_ko)
    if not os.path.exists(df_fname):    
        df = OptKnock.analyze_kos(CARBON_SOURCES, SINGLE_KOS,
                                  target_reaction, knockins, max_ko)
        with open(df_fname, 'w') as fp:
            df.round(3).to_csv(fp)
    else:
        df = pd.DataFrame.from_csv(df_fname)

    # create a 2D heatmap of all double knockouts and their slopes
    # in the same order as figure S3 (A) in Antonovsky et al. 2016:
    # doi: 10.1016/j.cell.2016.05.064
    
    df = df[df['carbon source'] != 'electrons']
    
    kos = SINGLE_KOS[1:]
    css = CARBON_SOURCES[0:9]
    index_list = df.knockouts.str.split('|').apply(
            lambda l: list(map(SINGLE_KOS.index, l))).tolist()
    df['ko1_idx'], df['ko2_idx'] = zip(*index_list)
    
    # replace the ko1 = 0 of the single knockouts with the same value of ko2,
    # so they will appear on the diagonal of the matrix
    df.loc[df.ko1_idx == 0, 'ko1_idx'] = df[df.ko1_idx == 0].ko2_idx
    df['cs_idx'] = df['carbon source'].apply(CARBON_SOURCES.index)
    
    df['row'] = (df.ko2_idx-1) * 3 + df.cs_idx // 3
    df['col'] = (df.ko1_idx-1) * 3 + df.cs_idx % 3
    data_mat = np.ones((len(kos)*3, len(kos)*3)) * -10
    
    for i, r in df.iterrows():
        if r['yield'] < 1e-9:
            data_mat[r['row'], r['col']] = -10
        else:
            data_mat[r['row'], r['col']] = r['slope']

    fig, ax = plt.subplots(figsize=(10, 8))
    labels = [kos[i//3] if i%3==1 else '' for i in range(len(kos)*3)]
    sns.heatmap(data_mat, vmin=-10,
                cmap=sns.dark_palette("palegreen"), ax=ax,
                xticklabels=labels, yticklabels=labels)
    #hm = HeatMap(data_mat, x='KO1', y='KO2', values='slope')
    #show(hm)