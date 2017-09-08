#!/usr/bin/python
from python.optknock import OptKnock
from python.models import knockin_reactions, clone_model, init_wt_model, \
                          set_exchange_bounds, knockout_reactions
import logging
import sys, os
import pandas as pd
from itertools import combinations

def analyze(target_reaction, knockins="", max_knockouts=2):
    """
        Args:
            target_reaction - the reaction for which the coupling to BM yield is made
            knockins        - extra reactions to add to the model
            max_knockouts   - the maximum number of simultaneous knockouts
    """
    carbon_uptake_rate = 50 # mmol C / (gDW*h)

    #single_ko_list = ['RPE','RPI','G6PDH2r,PGL','EDD,EDA']

    single_ko_list = ['GLCpts','PGI','PFK','FBP','FBA','TPI','GAPD,PGK','PGM',
                      'ENO','PYK','PPS','PDH','PFL','G6PDH2r,PGL','GND',
                      'EDD,EDA','RPE','RPI','TKT1,TKT2','TALA','PPC','PPCK',
                      'CS','ACONTa,ACONTb','ALCD2x','ACALD']
    carbon_sources = ['g6p', 'f6p', '6pgc', 'r5p', 'succ', 'xu5p_D',
                      '2pg', 'ac', 'dhap', '']
    
    core_model = init_wt_model('core', {}, BM_lower_bound=0.1)
    knockin_reactions(core_model, 'EDD,EDA', 0, 1000)
    knockin_reactions(core_model, 'EX_g6p,EX_f6p,EX_xu5p_D,EX_r5p,EX_dhap,'
                                  'EX_2pg,EX_e4p,EX_6pgc', 0, 0)
    
    wt_model = clone_model(core_model)
    if knockins is not None:
        knockin_reactions(wt_model, knockins, 0, 1000)

    set_exchange_bounds(wt_model, 'methanol', lower_bound=-999)
    
    sys.stdout.write("There are %d single knockouts\n" % len(single_ko_list))
    sys.stdout.write("There are %d carbon sources: %s\n" %
                     (len(carbon_sources), ', '.join(carbon_sources)))
    data = []
    
    for i in xrange(max_knockouts+1):
        for kos in combinations(single_ko_list, i):
            sys.stdout.write('KOs = ' + ', '.join(kos) + ': ')
            for carbon_source in carbon_sources:
                sys.stdout.write(carbon_source + ', ')
                temp_model = clone_model(wt_model)

                if carbon_source == '':
                    pass
                elif carbon_source == 'electrons':
                    knockin_reactions(temp_model, 'RED', 0,
                                      carbon_uptake_rate*2)
                else:
                    # find out how many carbon atoms are in the carbon source
                    # and normalize the uptake rate to be in units of mmol 
                    # carbon-source / (gDW*h) 
                    nC = 0
                    for cs in carbon_source.split(','):
                        met = wt_model.metabolites[
                            wt_model.metabolites.index(cs + '_c')]
                        nC += met.elements['C']
                    uptake_rate = carbon_uptake_rate / float(nC)
        
                    for cs in carbon_source.split(','):
                        set_exchange_bounds(temp_model, cs,
                                            lower_bound=-uptake_rate)
                
                for ko in kos:
                    knockout_reactions(temp_model, ko)
                
                yd = OptKnock(temp_model).solve_FBA() or 0
                if (target_reaction is not None) and (yd == 0):
                    slope = 0
                else:
                    slope = OptKnock(temp_model).get_slope(target_reaction)
                data.append([len(kos)] + 
                            list(kos) + 
                            [None] * (max_knockouts - len(kos)) + 
                            [carbon_source, yd, slope])
                
            sys.stdout.write('\n')

    ko_cols = ['knockout %d' % (i+1) for i in xrange(max_knockouts)]

    df = pd.DataFrame(data=data,
                      columns=['number of KOs'] + ko_cols +
                              ['carbon source', 'yield', 'slope'])

    df.set_index()

    return df

if __name__ == "__main__":
    logger = logging.getLogger('')
    logger.setLevel(logging.INFO)
    if not os.path.exists('res'):
        os.mkdir('res')
    
    max_ko = 3
    
    title = 'rump'; target_reaction = 'H6PS'; knockins = 'MEDH,H6PS,H6PI,H4MPTP,FDH'
    
    df_fname = 'res/data_%s.csv' % title
    if True:    
        df = analyze(target_reaction, knockins, max_ko)
        with open(df_fname, 'w') as fp:
            df.to_csv(fp)
    else:
        df = pd.DataFrame.from_csv(df_fname)

    # find knockouts that can grow without extra carbon
    # and also have a small slope for H6PS for at least one
    # other carbon
    ko_cols = ['knockout %d' % (i+1) for i in xrange(max_ko)]
    df[ko_cols] = df[ko_cols].fillna('')
    df['carbon source'] = df['carbon source'].fillna('methanol')
    df['all knockouts'] = df['knockout 1'].str.cat([df[c] for c in ko_cols[1:]], sep='|')
    df_yield = df.pivot(index='all knockouts', columns='carbon source', values='yield').round(3)
    df_slope = df.pivot(index='all knockouts', columns='carbon source', values='slope').round(3)
    
    methanol_growers_slopes = df_slope[df_yield['methanol'] > 0]
    methanol_growers_slopes = methanol_growers_slopes[methanol_growers_slopes.max(axis=1) > 0]

    with open('res/methanol_growers_slopes.csv', 'w') as fp:
        methanol_growers_slopes.to_csv(fp)