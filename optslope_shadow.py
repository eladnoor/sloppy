#!/usr/bin/python

from optknock import OptKnock
from models import *
import matplotlib.pyplot as plt
from analysis_toolbox import *
from html_writer import HtmlWriter
from copy import deepcopy

def main(core=True):
    main_html = HtmlWriter('res/optknock.html')
    main_html.write('<title>OptKnock</title>')
    main_html.write('<h1>OptKnock</h1>\n')

    if core:
        model = init_wt_model('core', {}, BM_lower_bound=0.1)
        knockin_reactions(model, 'EDD,EDA', 0, 1000)
        set_exchange_bounds(model, 'g6p', lower_bound=-10)
        main_html.write('<ul>\n')
        main_html.write('<li>Model used: E.coli core</li>\n')
        main_html.write('<li>Carbon source: glucose-6P, v <= 10 [mmol/g(DW) h]</li>\n')
        main_html.write('<li>Biomass yield lower bound: v >= 0.1 [mmol/g(DW) h]</li>\n')
    else:
        model = init_wt_model('full', {'rib_D': -10}, BM_lower_bound=0.1)
        main_html.write('<ul>\n')
        main_html.write('<li>Model used: iJO1366 (E. coli full model)</li>\n')
        main_html.write('<li>Carbon source: D-ribose, v <= 10 [mmol/g(DW) h]</li>\n')
        main_html.write('<li>Biomass yield lower bound: v >= 0.1 [mmol/g(DW) h]</li>\n')

    knockin_reactions(model, 'PRK,RBC', 0, 1000)
    main_html.write('<li>Added reactions: phosphoribulokinase, RuBisCo</li>\n')
    optknock_reaction = 'RBC'
    
    #knockin_reactions(model, 'DXS', 0, 1000)
    #main_html.write('<li>Added reactions: deoxyribose synthase</li>\n')
    #optknock_reaction = 'DXS'
    
    main_html.write('</ul>\n')

    print "Running standard FBA..."    
    ok = OptKnock(model)
    ok.prepare_FBA_primal()
    ok.prob.writeLP('res/fba_primal.lp')
    ok.solve()
    ok.print_primal_results(short=True)
    
    print '-' * 50

    print "Running dual FBA..."    
    ok = OptKnock(model)
    ok.prepare_FBA_dual()
    ok.prob.writeLP('res/fba_dual.lp')
    ok.solve()
    ok.print_dual_results(short=True)
    
    print '-' * 50

    if False:
        print "Running OptSlope..."
        ok = OptKnock(model)
        ok.prepare_optslope(optknock_reaction, num_deletions=3)
        ok.write_linear_problem('res/optslope.lp')
        solution = ok.solve()
        if solution.status is not None:
            ok.print_optknock_results(short=True)
            ok_model = ok.get_optknock_model()

            #fig, ax = plt.subplots(1)
            #plot_multi_PPP([model, ok_model],
            #               ['Wild-Type', r'OptKnock'],
            #               optknock_reaction, ax)
            #ax.set_title(title)
            #fig.savefig('res/OS_PPP.pdf')

        print '-' * 50
    
    print "Running OptKnock for maximizing flux in %s..." % optknock_reaction
    ok = OptKnock(model)
    ok.prepare_optknock(optknock_reaction, num_deletions=1)
    ok.write_linear_problem('res/optknock.lp')
    solution = ok.solve()
    if solution.status is not None:
        ok.print_optknock_results(short=True)
        knockouts = ok.get_optknock_knockouts()
        ko_model = clone_model(model)
        knockout_reactions(ko_model, knockouts)

        main_html.write('<h2>Knockouts: %s</h2>\n' % knockouts)
    
        # draw the PPP as embedded SVG
        fig, ax = plt.subplots(1, figsize=(6,6))
        wt_PPP, wt_slope = get_PPP(model, optknock_reaction)
        ax.fill_between(wt_PPP[:,0].flat, wt_PPP[:,1].flat, wt_PPP[:,2].flat,
                        facecolor='#E0E0E0', linewidth=0)

        ko_PPP, slope = get_PPP(ko_model, optknock_reaction)
        if slope is None:
            ko_text = 'Not feasible at any Rubisco flux'
            ko_color = '#FF9073'
        else:
            slope = np.round(slope, 1)
            ko_text = ('Slope = %g' % slope)
            if slope > 0:
                ko_color = '#00B64F'
            else:
                ko_color = '#FF7060'

        ax.fill_between(ko_PPP[:,0].flat, ko_PPP[:,1].flat, ko_PPP[:,2].flat,
                        facecolor=ko_color, linewidth=1)
        main_html.embed_matplotlib_figure(fig, width=400, height=400)

        # draw the flux map as embedded SVG
        ok.draw_svg(main_html)

        # write the flux summary for the knockout model as HTML
        ok.model_summary(main_html)

if __name__ == "__main__":
    main(core=True)
