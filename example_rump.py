#!/usr/bin/python
from python.optknock import OptKnock
import logging
import os
import pandas as pd

CARBON_SOURCES = ['methanol,g6p',
                  'methanol,f6p',
                  'methanol,6pgc',
                  'methanol,r5p',
                  'methanol,succ',
                  'methanol,xu5p_D',
                  'methanol,2pg',
                  'methanol,ac',
                  'methanol,dhap',
                  'methanol']

SINGLE_KOS = ['', 'GLCpts','PGI','PFK','FBP','FBA','TPI','GAPD,PGK','PGM',
              'ENO','PYK','PPS','PDH','PFL','G6PDH2r,PGL','GND',
              'EDD,EDA','RPE','RPI','TKT1,TKT2','TALA','PPC','PPCK',
              'CS','ACONTa,ACONTb','ALCD2x','ACALD']

if __name__ == "__main__":
    logger = logging.getLogger('')
    logger.setLevel(logging.INFO)
    if not os.path.exists('res'):
        os.mkdir('res')
    
    max_ko = 2
    
    title = 'rump'; target_reaction = 'H6PS'; knockins = 'MEDH,H6PS,H6PI,H4MPTP,FDH'
    
    df_fname = 'res/data_%s.csv' % title
    if True:    
        df = OptKnock.analyze_kos(CARBON_SOURCES, SINGLE_KOS,
                                  target_reaction, knockins, max_ko)
        with open(df_fname, 'w') as fp:
            df.round(3).to_csv(fp)
    else:
        df = pd.DataFrame.from_csv(df_fname)

    # find knockouts that can grow without extra carbon
    # and also have a small slope for H6PS for at least one
    # other carbon
    df_yield = df.pivot(index='knockouts', columns='carbon source', values='yield').round(3)
    df_slope = df.pivot(index='knockouts', columns='carbon source', values='slope').round(3)
    
    methanol_growers_slopes = df_slope[df_yield['methanol'] > 0]
    methanol_growers_slopes = methanol_growers_slopes[methanol_growers_slopes.max(axis=1) > 0]

    with open('res/methanol_growers_slopes.csv', 'w') as fp:
        methanol_growers_slopes.to_csv(fp)