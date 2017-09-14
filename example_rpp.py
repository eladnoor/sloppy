#!/usr/bin/python
from python.optknock import OptKnock
import logging
import os
import pandas as pd

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
    
    df_fname = 'res/data_%s.csv' % title
    if False:    
        df = OptKnock.analyze_kos(CARBON_SOURCES, SINGLE_KOS,
                                  target_reaction, knockins, max_ko)
        with open(df_fname, 'w') as fp:
            df.round(3).to_csv(fp)
    else:
        df = pd.DataFrame.from_csv(df_fname)

    # create a 2D heatmap of all double knockouts and their slopes
    
    
