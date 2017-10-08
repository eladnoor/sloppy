import numpy as np
from copy import deepcopy
from cobra.core import Solution
import sys
import pandas as pd
from itertools import combinations
import optlang
#optlang.interface = optlang.scipy_interface

from python import models
from python import draw_flux

M = 1000

class OptKnock(object):

    def __init__(self, model, verbose=False):
        self.model = deepcopy(model)
        self.prob = optlang.Model(name='OptKnock')
        self.verbose = verbose

        # locate the biomass reaction
        biomass_reactions = [r for r in self.model.reactions
                             if r.objective_coefficient != 0]
        if len(biomass_reactions) != 1:
            raise Exception('There should be only one single biomass reaction')
        self.r_biomass = biomass_reactions[0]
        
        self.has_flux_as_variables = False

    def add_primal_variables_and_constraints(self):
        # create the continuous flux variables (can be positive or negative)
        self.var_v = {}
        for r in self.model.reactions:
            self.var_v[r] = optlang.Variable('v_%s' % r.id,
                                             lb=r.lower_bound,
                                             ub=r.upper_bound)

        # this flag will be used later to know if to expect the flux
        # variables to exist
        self.has_flux_as_variables = True
        
        # add the mass-balance constraints to each of the metabolites (S*v = 0)
        for m in self.model.metabolites:
            S_times_v = sum(self.var_v[r] * r.get_coefficient(m)
                            for r in m.reactions)
            self.prob.add(optlang.Constraint(S_times_v, lb=0, ub=0))
    
    def add_dual_variables_and_constraints(self):
        # create dual variables associated with stoichiometric constraints
        self.var_lambda = dict([(m, optlang.Variable('lambda_%s' % m.id, 
                                             lb=-M, un=M))
                                for m in self.model.metabolites])

        # create dual variables associated with the constraints on the primal fluxes
        self.var_w_U = dict([(r, optlang.Variable("w_U_%s" % r.id, lb=0, ub=M))
                             for r in self.model.reactions])
        self.var_w_L = dict([(r, optlang.Variable("w_L_%s" % r.id, lb=0, ub=M))
                             for r in self.model.reactions])

        # add the dual constraints:
        #   S'*lambda + w_U - w_L = c_biomass
        for r in self.model.reactions:
            S_times_lambda = sum(self.var_lambda[m] * coeff
                                 for m, coeff in r._metabolites.iteritems()
                                 if coeff != 0)
            row_sum = S_times_lambda + self.var_w_U[r] - self.var_w_L[r]
            self.prob.add(optlang.Constraint(row_sum, lb=r.objective_coefficient,
                                             ub=r.objective_coefficient))
                                   
    def prepare_FBA_primal(self):
        """
            Run standard FBA (primal)
        """
        self.add_primal_variables_and_constraints()
        self.prob.objective = optlang.Objective(self.var_v[self.r_biomass],
                                                direction='max')

    def prepare_FBA_dual(self, use_glpk=False):
        """
            Run shadow FBA (dual)
        """
        self.add_dual_variables_and_constraints()
        
        w_sum_ub = sum([self.var_w_U[r] * r.upper_bound
                        for r in self.model.reactions if r.upper_bound != 0])
        w_sum_lb = sum([self.var_w_L[r] * -r.lower_bound
                        for r in self.model.reactions if r.lower_bound != 0])
        self.prob.objective = optlang.Objective(w_sum_ub + w_sum_lb, direction='min')
    
    def get_reaction_by_id(self, reaction_id):
        if reaction_id not in self.model.reactions:
            return None
        reaction_ind = self.model.reactions.index(reaction_id)
        reaction = self.model.reactions[reaction_ind]
        return reaction
        
    def add_optknock_variables_and_constraints(self):
        # create the binary variables indicating which reactions knocked out
        self.var_y = dict([(r, optlang.Variable("y_%s" % r.id, type='binary'))
                           for r in self.model.reactions])

        # create dual variables associated with the constraints on the primal fluxes
        self.var_mu = dict([(r, optlang.Variable("mu_%s" % r.id))
                             for r in self.model.reactions])

        # equate the objectives of the primal and the dual of the inner problem
        # to force its optimization:
        #   sum_j mu_j - v_biomass = 0
        constr = optlang.Constraint(sum(self.var_mu.values()) - self.var_v[self.r_biomass],
                                    lb=0, ub=0)
        self.prob.add(constr)

        # add the knockout constraints (when y_j = 0, v_j has to be 0)
        for r in self.model.reactions:
            # L_jj * y_j <= v_j
            self.prob.addConstraint(r.lower_bound * self.var_y[r] <= self.var_v[r], 'v_lower_%s' % r.id)
            # v_j <= U_jj * y_j
            self.prob.addConstraint(self.var_v[r] <= r.upper_bound * self.var_y[r], 'v_upper_%s' % r.id)
            
        # set the constraints on the auxiliary variables (mu):
        #    mu_j == y_j * (U_jj * w_u_j - L_jj * w_l_j)
        for r in self.model.reactions:
            w_sum = self.var_w_U[r] * r.upper_bound - self.var_w_L[r] * r.lower_bound

            # mu_j + M*y_j >= 0
            self.prob.add(optlang.Constraint(self.var_mu[r] + M*self.var_y[r], lb=0))
            # -mu_j + M*y_j >= 0
            self.prob.add(optlang.Constraint(-self.var_mu[r] + M*self.var_y[r], lb=0))
            # mu_j - (U_jj * w_u_j - L_jj * w_l_j) + M*(1-y_j) >= 0
            self.prob.add(optlang.Constraint(self.var_mu[r] - w_sum + M*(1-self.var_y[r]), lb=0))
            # -mu_j + (U_jj * w_u_j - L_jj * w_l_j) + M*(1-y_j) >= 0
            self.prob.add(optlang.Constraint(-self.var_mu[r] + w_sum + M*(1-self.var_y[r]), lb=0))

    def add_knockout_bounds(self, ko_candidates=None, num_deletions=5):
        """ 
            construct the list of KO candidates and add a constraint that
            only K (num_deletians) of them can have a y_j = 0
        """
        ko_candidate_sum_y = []
        
        if ko_candidates is None:
            ko_candidates = [r for r in self.model.reactions if r != self.r_biomass]

        for r in set(self.model.reactions).difference(ko_candidates):
            # if 'r' is not a candidate constrain it to be 'active'
            # i.e.   y_j == 1
            self.prob.addConstraint(self.var_y[r] == 1, 'active_%s' % r.id)

        # set the upper bound on the number of knockouts (K)
        #   sum (1 - y_j) <= K
        ko_candidate_sum_y = sum([self.var_y[r] for r in ko_candidates])
        constr = optlang.Constraint(ko_candidate_sum_y, lb=(len(ko_candidate_sum_y) - num_deletions))
        self.prob.add(constr, 'number_of_deletions')

    def prepare_optknock(self, target_reaction_id, ko_candidates=None, 
                         num_deletions=5, use_glpk=False):
        # find the target reaction
        self.r_target = self.get_reaction_by_id(target_reaction_id)

        self.add_primal_variables_and_constraints()
        self.add_dual_variables_and_constraints()
        self.add_optknock_variables_and_constraints()

        # add the objective of maximizing the flux in the target reaction
        self.prob.objective = optlang.Objective(
            self.var_v[self.r_target], direction='max')

        self.add_knockout_bounds(ko_candidates, num_deletions)

    def prepare_optslope(self, target_reaction_id, ko_candidates=None,
                         num_deletions=5, use_glpk=False):
        # add the objective of maximizing the flux in the target reaction
        self.r_target = self.get_reaction_by_id(target_reaction_id)

        # set biomass maximum to 0
        self.r_biomass.lower_bound = 0
        self.r_biomass.upper_bound = 0
        self.r_target.lower_bound = 0
        self.r_target.upper_bound = 0

        self.add_primal_variables_and_constraints()
        self.add_dual_variables_and_constraints()
        self.add_optknock_variables_and_constraints()

        # set the objective as maximizing the shadow price of v_target upper bound
        self.prob.objective = optlang.Objective(
            self.var_w_U[self.r_target] - self.var_w_L[self.r_target],
            direction='max')

        self.add_knockout_bounds(ko_candidates, num_deletions)

    def write_linear_problem(self, fname):
        self.prob.writeLP(fname)

    def solve(self):
        self.prob.optimize()

        if self.prob.status != 'optimal':
            if self.verbose:
                print("LP was not solved because: " + self.prob.status)
            self.solution = Solution(objective_value=None,
                                     status=self.prob.status,
                                     fluxes=None)
        else:
            if self.has_flux_as_variables:
                x = [self.var_v[r].primal for r in self.model.reactions]
            else:
                x = []
            self.solution = Solution(objective_value=self.prob.objective.value,
                                     status=self.prob.status,
                                     fluxes=x)
        return self.solution
    
    def get_objective_value(self):
        if self.solution.status != 'optimal':
            return None
        else:
            return self.prob.objective.value

    def print_primal_results(self, short=True):
        obj = self.get_objective_value()
        if obj is None:
            return
        print("Objective : %6.3f" % obj)
        if not short:
            print("List of reactions : ")
            for r in self.model.reactions:
                print("%30s (%4g <= v <= %4g) : v = %6.3f" % \
                    (r.name, r.lower_bound, r.upper_bound, self.var_v[r].primal))

    def print_dual_results(self, short=True):
        obj = self.get_objective_value()
        if obj is None:
            return
        print("Objective : %6.3f" % obj)
        if not short:
            print("List of reactions : ")
            for r in self.model.reactions:
                print("%30s (%4g <= v <= %4g) : w_L = %5.3f, w_U = %5.3f" % \
                    (r.id, r.lower_bound, r.upper_bound, 
                     self.var_w_L[r].primal, self.var_w_U[r].primal))
            print("List of metabolites : ")
            for m in self.model.metabolites:
                print("%30s : lambda = %5.3f, " % \
                    (m.id, self.var_lambda[m].primal))
                
    def print_optknock_results(self, short=True):
        if self.solution.status != 'optimal':
            return
        print("Objective : %6.3f" % self.get_objective_value)
        print("Biomass rate : %6.3f" % self.var_v[self.r_biomass].primal)
        print("Sum of mu : %6.3f" % np.sum([mu.primal for mu in self.var_mu.values()]))
        print("Knockouts : ")
        print('   ;   '.join(['"%s" (%s)' % (r.name, r.id) for r, val in self.var_y.iteritems() if val.primal == 0]))
        if not short:
            print("List of reactions : ")
            for r in self.model.reactions:
                print('%25s (%5s) : %4g  <=  v=%5g  <=  %4g ; y = %d ; mu = %g ; w_L = %5g ; w_U = %5g' % \
                    ('"' + r.name + '"', r.id,
                     r.lower_bound, self.var_v[r].primal, r.upper_bound,
                     self.var_y[r].primal, self.var_mu[r].primal,
                     self.var_w_L[r].primal, self.var_w_U[r].primal))
            print("List of metabolites : ")
            for m in self.model.metabolites:
                print("%30s : lambda = %6.3f" % \
                    (m.id, self.var_lambda[m].primal))

    def get_optknock_knockouts(self):
        return ','.join([r.id for r, val in self.var_y.iteritems() if val.primal == 0])
    
    def get_optknock_model(self):
        if self.solution.status != 'optimal':
            raise Exception('OptKnock failed, cannot generate a KO model')
        
        optknock_model = deepcopy(self.model)
        knockout_reactions = [r for r, val in self.var_y.iteritems() if val.primal == 0]
        for r in knockout_reactions:
            new_r = optknock_model.reactions[optknock_model.reactions.index(r.id)]
            new_r.lower_bound = 0
            new_r.upper_bound = 0
        return optknock_model

    def solve_FBA(self):
        self.add_primal_variables_and_constraints()
        self.prob.objective = optlang.Objective(self.var_v[self.r_biomass],
                                                direction='max')
        self.solve()
        max_biomass = self.get_objective_value()
        return max_biomass

    def solve_FVA(self, reaction_id):
        """
            Run Flux Variability Analysis on the provided reaction
        """
        self.add_primal_variables_and_constraints()
        self.prob.setObjective(self.var_v[self.r_biomass])
        self.solve()
        max_biomass = self.get_objective_value()
        if max_biomass is None:
            raise Exception("Cannot run FVA because the model is infeasible")
        self.var_v[self.r_biomass].lowBound = max_biomass - 1e-5

        r_target = self.get_reaction_by_id(reaction_id)
        self.prob.objective = optlang.Objective(self.var_v[r_target],
                                                direction='max')
        self.solve()
        max_v_target = self.get_objective_value()

        self.prob.objective = optlang.Objective(self.var_v[r_target],
                                                direction='min')
        self.solve()
        min_v_target = self.get_objective_value()
        
        return min_v_target, max_v_target

    def get_PPP_data(self, reaction_id, bm_range=None):
        """
            Run FVA on a gradient of biomass lower bounds and generate
            the data needed for creating the Phenotype Phase Plane
        """
        self.add_primal_variables_and_constraints()
        self.prob.objective = optlang.Objective(self.var_v[self.r_biomass],
                                                direction='max')
        r_target = self.get_reaction_by_id(reaction_id)
        if r_target is None:
            return None

        self.solve()

        if bm_range is None:
            max_biomass = self.get_objective_value()
            if max_biomass is None:
                return None
            bm_range = np.linspace(1e-5, max_biomass - 1e-5, 50)

        data = []
        for bm_lb in bm_range:
            self.var_v[self.r_biomass].lb = bm_lb

            self.prob.objective = optlang.Objective(self.var_v[r_target],
                                                    direction='max')
            self.solve()
            max_v_target = self.get_objective_value()

            self.prob.objective = optlang.Objective(self.var_v[r_target],
                                                    direction='min')
            self.solve()
            min_v_target = self.get_objective_value()
            
            data.append((bm_lb, min_v_target, max_v_target))
            
        return np.matrix(data)
        
    def get_slope(self, reaction_id, epsilon_bm=0.01):
        data = self.get_PPP_data(reaction_id, bm_range=[epsilon_bm])
        #print(data)
        if data is None:
            return None
        else:
            return data[0, 1] / epsilon_bm
        
    def model_summary(self, html):
        import analysis_toolbox
        analysis_toolbox.model_summary(self.model, self.solution, html)
        
    def draw_svg(self, html):
        # Parse the SVG file of central metabolism
        drawer = draw_flux.DrawFlux('data/CentralMetabolism.svg')
        #drawer = DrawFlux('data/EcoliMetabolism.svg')
        drawer.ToSVG(self.model, self.solution, html)

    @staticmethod
    def analyze_kos(carbon_sources, single_kos,
                    target_reaction, knockins="", n_knockouts=2, n_threads=2,
                    carbon_uptake_rate=50):
        """
            Args:
                target_reaction - the reaction for which the coupling to BM yield is made
                knockins        - extra reactions to add to the model
                max_knockouts   - the maximum number of simultaneous knockouts
                carbon_uptake_rate - in units of mmol C / (gDW*h)
        """
        wt_model = models.Model.initialize()
        
        if knockins is not None:
            wt_model.knockin_reactions(knockins, 0, 1000)
      
        sys.stdout.write("There are %d single knockouts\n" % len(single_kos))
        sys.stdout.write("There are %d carbon sources: %s\n" %
                         (len(carbon_sources), ', '.join(carbon_sources)))
        data = []
        
        kos_and_cs = [(kos, cs) for kos in combinations(single_kos, n_knockouts)
                                for cs in carbon_sources]
            
        def calculate_yield_and_slope(params):
            kos, carbon_source = params
            sys.stderr.write('KOs = ' + ', '.join(kos) + ': ' + carbon_source + '\n')
            temp_model = wt_model.clone()
            if carbon_source == 'electrons':
                temp_model.knockin_reactions('RED', 0, carbon_uptake_rate*2)
            elif carbon_source != '':
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
                    temp_model.set_exchange_bounds(cs, lower_bound=-uptake_rate)
            
            for ko in kos:
                if ko != '':
                    temp_model.knockout_reactions(ko)
            
            yd = OptKnock(temp_model).solve_FBA() or 0
            if (target_reaction is not None) and (yd == 0):
                slope = 0
            else:
                slope = OptKnock(temp_model).get_slope(target_reaction)
            return ('|'.join(kos), carbon_source, yd, slope)
        
        data = map(calculate_yield_and_slope, kos_and_cs)
        df = pd.DataFrame(data=list(data),
                          columns=['knockouts', 'carbon source', 'yield', 'slope'])
    
        return df