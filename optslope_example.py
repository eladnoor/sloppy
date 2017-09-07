#!/usr/bin/python
from python.optknock import OptKnock
from python.models import knockin_reactions, clone_model, init_wt_model, \
                          set_exchange_bounds, knockout_reactions
import logging
import sys, os
import pandas as pd
from itertools import combinations

def main():
    #analyze('wild_type', None, None)
    #analyze('rubisco', 'RBC', 'RBC,PRK', 3)
    analyze('rump', 'H6PS', 'MMO,H6PS,H6PI', 2)
    #analyze('deoxyribose', 'DXS', 'DXS')
    #analyze('POG', 'MCS', 'MCS,MCL')

def analyze(title, target_reaction, knockins="", max_knockouts=2):
    """
        Args:
            title           - prefix to be used for saving the result files
            target_reaction - the reaction for which the coupling to BM yield is made
            knockins        - extra reactions to add to the model
    """
    carbon_uptake_rate = 50 # mmol C / (gDW*h)

    #single_ko_list = ['FBA','PFL','G6PDH2r,PGL']

    single_ko_list = ['GLCpts','PGI','PFK','FBP','FBA','TPI','GAPD,PGK','PGM',
                      'ENO','PYK','PPS','PDH','PFL','G6PDH2r,PGL','GND',
                      'EDD,EDA','RPE','RPI','TKT1,TKT2','TALA','PPC','PPCK',
                      'CS','ACONTa,ACONTb','ALCD2x','ACALD']
    carbon_sources = ['g6p', 'f6p', '6pgc', 'r5p', 'succ', 'xu5p_D',
                      '2pg', 'ac', 'dhap', 'electrons']
    
    core_model = init_wt_model('core', {}, BM_lower_bound=0.1)
    knockin_reactions(core_model, 'EDD,EDA', 0, 1000)
    knockin_reactions(core_model, 'EX_g6p,EX_f6p,EX_xu5p_D,EX_r5p,EX_dhap,EX_2pg,EX_e4p,EX_6pgc', 0, 0)
    
    wt_model = clone_model(core_model)
    if knockins is not None:
        knockin_reactions(wt_model, knockins, 0, 1000)

    set_exchange_bounds(wt_model, 'methanol', lower_bound=-999)
    
    sys.stdout.write("There are %d single knockouts\n" % len(single_ko_list))
    sys.stdout.write("There are %d carbon sources: %s\n" % (len(carbon_sources), ', '.join(carbon_sources)))
    data = []
    
    for i in xrange(max_knockouts+1):
        for kos in combinations(single_ko_list, i):
            sys.stdout.write('KOs = ' + ', '.join(kos) + ': ')
            for carbon_source in carbon_sources:
                sys.stdout.write(carbon_source + ', ')
                temp_model = clone_model(wt_model)

                if carbon_source == 'electrons':
                    knockin_reactions(temp_model, 'RED', 0, carbon_uptake_rate*2)
                else:
                    # find out how many carbon atoms are in the carbon source
                    # and normalize the uptake rate to be in units of mmol carbon-source / (gDW*h) 
                    nC = 0
                    for cs in carbon_source.split(','):
                        met = wt_model.metabolites[wt_model.metabolites.index(cs + '_c')]
                        nC += met.elements['C']
                    uptake_rate = carbon_uptake_rate / float(nC)
        
                    for cs in carbon_source.split(','):
                        set_exchange_bounds(temp_model, cs, lower_bound=-uptake_rate)
                
                for ko in kos:
                    knockout_reactions(temp_model, ko)
                
                yd = OptKnock(temp_model).solve_FBA() or 0
                if (target_reaction is not None) and (yd == 0):
                    slope = 0
                else:
                    slope = OptKnock(temp_model).get_slope(target_reaction)
                data.append(['|'.join(kos), carbon_source, yd, slope])
                
            sys.stdout.write('\n')

    df = pd.DataFrame(data=data,
                      columns=['knockouts', 'carbon source', 'yield', 'slope'])

    with open('res/data2D_%s.csv' % title, 'w') as fp:
        df.to_csv(fp)

if __name__ == "__main__":
    logger = logging.getLogger('')
    logger.setLevel(logging.INFO)
    if not os.path.exists('res'):
        os.mkdir('res')
    main()
