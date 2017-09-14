from copy import deepcopy
from cobra.manipulation.modify import convert_to_irreversible
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.core import Reaction, Metabolite
from cobra.core.formula import Formula


class Model(object):

    def __init__(self):
        self.cobra_model = None
    
    @staticmethod
    def initialize(model_name='core',
                   carbon_sources={},
                   irreversible=False,
                   ATP_maintenance=False, BM_lower_bound=0.1):
        
        if model_name == 'core':
            model = create_cobra_model_from_sbml_file(
                'data/ecoli_core.xml', old_sbml=True)
            if irreversible:
                convert_to_irreversible(model)
            # the core model has these annoying '_b' metabolites that are used as
            # 'ghost' metabolites that balance the exchange reactions. they should
            # be ignored in the mass-balance equation and therefore the best way to
            # deal with them is to remove them from all the reactions
            for m in model.metabolites:
                if str(m).endswith('_b'):
                    for r in m.get_reaction():
                        coeff = r.get_coefficient(m)
                        r.add_metabolites({m : -coeff})
    
            rxns = dict([(r.id, r) for r in model.reactions])
            if not ATP_maintenance:
                rxns['ATPM'].lower_bound = 0 # remove the ATP maintenance requirement
            rxns['EX_glc_e'].lower_bound = 0 # remove the default carbon source
        elif model_name == 'full':
            model = create_cobra_model_from_sbml_file('data/iJO1366.xml')
            rxns = dict([(r.id, r) for r in model.reactions])
            if not ATP_maintenance:
                rxns['ATPM'].lower_bound = 0 # remove the ATP maintenance requirement
            rxns['EX_glc_e'].lower_bound = 0 # remove the default carbon source
        elif model_name == 'toy':
            model = create_cobra_model_from_sbml_file('data/toymodel.xml')
            
        for key, val in carbon_sources.items():
            rxns['EX_' + key + '_e'].lower_bound = val
            
        # set BM lower bound
        for r in model.reactions:
            if r.objective_coefficient != 0:
                r.lower_bound = BM_lower_bound
        
        m = Model()
        m.cobra_model = model
        if model_name == 'core':
            # the ED pathway and exchange reactions are missing from the
            # core model. Add them now.
            m.knockin_reactions('EDD,EDA', 0, 1000)
            m.knockin_reactions('EX_g6p,EX_f6p,EX_xu5p_D,EX_r5p,EX_dhap,'
                                'EX_2pg,EX_e4p,EX_6pgc', 0, 0)
        return m
    
    def clone(self):
        new_model = Model()
        new_model.cobra_model = deepcopy(self.cobra_model)
        return new_model
    
    @property
    def metabolites(self):
        return self.cobra_model.metabolites

    @property
    def reactions(self):
        return self.cobra_model.reactions
    
    def add_metabolite(self, cid, formula, name, compartment='C'):
        try:
            self.metabolites.index(cid)
        except ValueError:
            met = Metabolite(id=cid, formula=Formula(formula),
                             name=name, compartment=compartment)
            self.cobra_model.add_metabolites([met])
    
    def knockout_reactions(self, ko_reactions):
        for r in ko_reactions.split(','):
            self.cobra_model.remove_reactions(r)
    
    def add_reaction(self, rid, name, sparse,
                     lower_bound=0, upper_bound=1000):
        """
            Adds a reaction to the model
        """
        # convert the sparse representation using the metabolites in the model
        for key in sparse.keys():
            if key not in self.metabolites:
                raise Exception("cannot find the cytoplasmic metabolite %s in the model" % key)
    
        r = dict([(self.metabolites[self.metabolites.index(key)], val)
                  for key, val in sparse.items()])
        reaction = Reaction(name)
        reaction.id = rid
        reaction.add_metabolites(r)
        reaction.lower_bound = lower_bound
        reaction.upper_bound = upper_bound
        self.cobra_model.add_reactions([reaction])
        return reaction
            
    def knockin_reactions(self, ki_reactions, lower_bound=0, upper_bound=1000):
        for rid in ki_reactions.split(','):
            if rid.startswith('EX_'):
                rid = rid[3:]
                self.add_metabolite_exchange(rid, lower_bound, upper_bound)
            elif rid == 'PRK':
                self.add_metabolite('rubp_D_c', 'C5H12O11P2', 'D-ribulose 1,5-bisphosphate')
                sprs = {'ru5p_D_c' : -1, 'atp_c' : -1, 'rubp_D_c' : 1, 'adp_c' : 1}
                self.add_reaction(rid, 'phosphoribulokinase', sprs, lower_bound, upper_bound)
            elif rid == 'RBC':
                self.add_metabolite('rubp_D_c', 'C5H12O11P2', 'D-ribulose 1,5-bisphosphate')
                sprs = {'rubp_D_c' : -1, 'h2o_c' : -1, 'co2_c' : -1, '3pg_c' : 2, 'h_c' : 3}
                self.add_reaction(rid, 'RuBisCO carboxylation', sprs, lower_bound, upper_bound)
            elif rid == 'RBC':
                self.add_metabolite('rubp_D_c', 'C5H12O11P2', 'D-ribulose 1,5-bisphosphate')
                sprs = {'rubp_D_c' : -1, 'h2o_c' : -1, 'co2_c' : -1, '3pg_c' : 2, 'h_c' : 3}
                self.add_reaction(rid, 'RuBisCO carboxylation', sprs, lower_bound, upper_bound)
            elif rid == 'PRK+RBC':
                sprs = {'ru5p_D_c' : -1, 'atp_c' : -1, 'h2o_c' : -1, 'co2_c' : -1, '3pg_c' : 2, 'h_c' : 3, 'adp_c' : 1}
                self.add_reaction(rid, 'PRK+RuBisCO', sprs, lower_bound, upper_bound)
            elif rid == 'EDD':
                self.add_metabolite('2ddg6p_c', 'C6H8O9P', '2-dehydro-3-deoxy-D-gluconate 6-phosphate')
                sprs = {'6pgc_c' : -1, 'h2o_c' : 1, '2ddg6p_c' : 1}
                self.add_reaction(rid, '6-phosphogluconate dehydratase', sprs, lower_bound, upper_bound)
            elif rid == 'EDA':
                self.add_metabolite('2ddg6p_c', 'C6H8O9P', '2-dehydro-3-deoxy-D-gluconate 6-phosphate')
                sprs = {'2ddg6p_c' : -1, 'g3p_c' : 1, 'pyr_c' : 1}
                self.add_reaction(rid, '2-dehydro-3-deoxy-phosphogluconate aldolase', sprs, lower_bound, upper_bound)
            elif rid == 'PKT':
                sprs = {'f6p_c' : -1, 'pi_c' : -1, 'e4p_c' : 1, 'actp_c' : 1, 'h2o_c' : 1}
                self.add_reaction(rid, 'phosphoketolase', sprs, lower_bound, upper_bound)
            elif rid == 'RED':
                sprs = {'nad_c' : -1, 'nadh_c' : 1}
                self.add_reaction(rid, 'free_e', sprs, lower_bound, upper_bound)            
            elif rid == 'ATP':
                sprs = {'adp_c' : -1, 'atp_c' : 1}
                self.add_reaction(rid, 'free_e', sprs, lower_bound, upper_bound)            
            elif rid == 'DXS':
                sprs = {'3pg_c':-1, 'pyr_c':-1}
                self.add_reaction(rid, 'deoxyribose synthase', sprs, lower_bound, upper_bound)
            elif rid == 'MCS':
                self.add_metabolite('malcoa_c', 'C25H40N7O20P3S', 'Malyl-CoA')
                sprs = {'mal_L_c':-1, 'atp_c':-1, 'coa_c':-1, 'malcoa_c':1, 'adp_c':1, 'pi_c':1}
                self.add_reaction(rid, 'malyl-CoA synthase', sprs, lower_bound, upper_bound)
            elif rid == 'MCL':
                self.add_metabolite('malcoa_c', 'C25H40N7O20P3S', 'Malyl-CoA')
                sprs = {'malcoa_c':-1, 'accoa_c':1, 'glx_c':1}
                self.add_reaction(rid, 'malyl-CoA lyase', sprs, lower_bound, upper_bound)
            elif rid == 'SBP':
                self.add_metabolite('sbp_c', 'C7H16O13P2', 'D-sedoheptulose 1,7-bisphosphate')
                sprs = {'sbp_c':-1, 'h2o_c':-1, 's7p_c':1, 'pi_c':1}
                self.add_reaction(rid, 'sedoheptulose bisphosphate phosphatase', sprs, lower_bound, upper_bound)
            elif rid == 'SBA':
                self.add_metabolite('sbp_c', 'C7H16O13P2', 'D-sedoheptulose 1,7-bisphosphate')
                sprs = {'sbp_c':1, 'g3p_c':-1, 'e4p_c':-1}
                self.add_reaction(rid, 'sedoheptulose bisphosphate aldolase', sprs, lower_bound, upper_bound)
            elif rid == 'MEDH':
                self.add_metabolite('methanol_c', 'CH4O', 'methanol')
                self.add_metabolite('formaldehyde_c', 'CH2O', 'formaldehyde')
                sprs = {'methanol_c' : -1, 'nad_c' : -1, 'formaldehyde_c' : 1, 'nadh_c' : 1}
                self.add_reaction(rid, 'methanol dehydrogenase', sprs, lower_bound, upper_bound)
            elif rid == 'H6PS':
                self.add_metabolite('hexulose6p_c', 'CH4O', 'D-hexulose 6-phosphate')
                sprs = {'ru5p_D_c' : -1, 'formaldehyde_c' : -1, 'hexulose6p_c' : 1}
                self.add_reaction(rid, 'hexulose-6-phosphate synthetase', sprs, lower_bound, upper_bound)
            elif rid == 'H6PI':
                self.add_metabolite('hexulose6p_c', 'CH4O', 'D-hexulose 6-phosphate')
                sprs = {'hexulose6p_c' : -1, 'f6p_c' : 1}
                self.add_reaction(rid, 'hexulose-6-phosphate isomerase', sprs, lower_bound, upper_bound)
            elif rid == 'H4MPTP':
                self.add_metabolite('formaldehyde_c', 'CH2O', 'formaldehyde')
                self.add_metabolite('for_c', 'CH2O2', 'formate')
                sprs = {'formaldehyde_c' : -1, 'for_c' : 1}
                self.add_reaction(rid, 'methylene tetrahydromethanopterin pathway', sprs, lower_bound, upper_bound)
            elif rid == 'FDH':
                sprs = {'for_c' : -1, 'nad_c': -1, 'co2_c' : 1, 'nadh_c': 1}
                self.add_reaction(rid, 'formate dehydrogenase', sprs, lower_bound, upper_bound)
            else:
                raise Exception('unknown knockin reaction: ' + rid)
                
    def add_metabolite_exchange(self, metabolite, lower_bound, upper_bound=0):
        try:
            met = self.metabolites[self.metabolites.index(metabolite + '_c')]
        except AttributeError:
            raise KeyError('Model does not have a metabolite with ID: ' + metabolite)
        
        self.add_metabolite(metabolite + '_e', str(met.formula), met.name, 'E')
        self.add_reaction(metabolite + '_transport', met.name + ' permease',
                     {metabolite + '_c' : -1, metabolite + '_e' : 1}, -1000, 1000)
                     
        self.add_reaction('EX_' + metabolite + '_e', met.name + ' exchange',
                     {metabolite + '_e' : -1}, lower_bound, upper_bound)
    
    def set_exchange_bounds(self, metabolite, lower_bound, upper_bound=0):
        rxns = dict([(r.id, r) for r in self.cobra_model.reactions])
        try:
            rxns['EX_' + metabolite + '_e'].lower_bound = lower_bound
            rxns['EX_' + metabolite + '_e'].upper_bound = upper_bound
        except KeyError:
            self.add_metabolite_exchange(metabolite, lower_bound, upper_bound)
    
    def set_single_precursor_objective(self, metabolite, lower_bound=0, upper_bound=1000):
        try:
            met = self.metabolites[self.metabolites.index(metabolite + '_c')]
        except AttributeError:
            raise KeyError('Model does not have a metabolite with ID: ' + metabolite)
        
        for r in self.cobra_model.reactions:
            r.objective_coefficient = 0
    
        if metabolite == 'accoa':
            r = self.add_reaction('Biomass_accoa', met.name + ' biomass',
                             {'accoa_c' : -1, 'coa_c' : 1}, lower_bound, upper_bound)
        else:
            r = self.add_reaction('Biomass_' + metabolite, met.name + ' biomass',
                             {metabolite + '_c' : -1}, lower_bound, upper_bound)
        r.objective_coefficient = 1
        
    def solve(self):
        return self.cobra_model.solve()