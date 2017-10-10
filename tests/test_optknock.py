#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 17:22:11 2017

@author: noore
"""

import unittest

class TestReactionParsing(unittest.TestCase):
    
    def test_rump_optlang(self):
        try:
            import optlang
        except ImportError:
            self.fail('optlang is not installed, please run "pip install optlang"')

        from python import optknock_optlang
        
        target_reaction = 'H6PS'
        knockins = 'MEDH,H6PS,H6PI,H4MPTP,FDH'
        df = optknock_optlang.OptKnock.analyze_kos(
                ['methanol,succ', 'methanol,xu5p_D'],
                ['FBP', 'RPI'],
                target_reaction, knockins,
                n_knockouts=1, solver='glpk')

        yield_df = df.pivot('knockouts', 'carbon source', 'yield')        
        slope_df = df.pivot('knockouts', 'carbon source', 'slope')        
        self.assertAlmostEqual(yield_df['methanol,succ']['FBP'], 0.682, 3)
        self.assertAlmostEqual(slope_df['methanol,xu5p_D']['RPI'], 9.396, 3)

    def test_rump_pulp(self):
        try:
            import pulp
        except ImportError:
            self.fail('PuLP is not installed, please run "pip install pulp"')
        
        from python import optknock_pulp
        
        target_reaction = 'H6PS'
        knockins = 'MEDH,H6PS,H6PI,H4MPTP,FDH'
        df = optknock_pulp.OptKnock.analyze_kos(
                ['methanol,succ', 'methanol,xu5p_D'],
                ['FBP', 'RPI'],
                target_reaction, knockins,
                n_knockouts=1, solver='glpk')

        yield_df = df.pivot('knockouts', 'carbon source', 'yield')        
        slope_df = df.pivot('knockouts', 'carbon source', 'slope')        
        self.assertAlmostEqual(yield_df['methanol,succ']['FBP'], 0.682, 3)
        self.assertAlmostEqual(slope_df['methanol,xu5p_D']['RPI'], 9.396, 3)
        
if __name__ == '__main__':
    unittest.main()
