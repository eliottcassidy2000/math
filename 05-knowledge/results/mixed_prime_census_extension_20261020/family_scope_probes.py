#!/usr/bin/env python3
"""Selected p=7 diagonal/product scope probes beyond the complete census.
These four embeddings do not constitute a subgroup census at n=13 or n=14.
"""
import extension_census as e
import json,sys

def main():
    m=e.m;rows=[]
    for n in (13,14):
        t=m.cyc(n,[list(range(7))]);a=tuple([2*i%7 for i in range(7)]+list(range(7,n)));b=m.cyc(n,[[7,8,9]])
        pr=m.census_case(n,'F21_times_C3',[t,a,b]);dr=m.census_case(n,'F21_diagonal_7plus3',[t,m.compose(a,b)])
        e.need(dr['tau_T']==pr['tau_T'] and dr['tau_E']==pr['tau_E'],'same full predicates')
        e.need(dr['relative_tournament_orbits']==2*pr['relative_tournament_orbits'],'index two relative law')
        e.need(dr['tournament_fusion_defect']==2*pr['tournament_fusion_defect']+pr['tau_T'],'tournament fusion law')
        e.need(dr['even_graph_fusion_defect']==pr['even_graph_fusion_defect'],'even fusion law')
        e.need(dr['sign_escape_relative_orbits']==pr['relative_tournament_orbits']+pr['sign_escape_relative_orbits'],'escape law')
        keys=('n','class','order','generators','normalizer_order','pair_orbits','fixed_objects_each','relative_tournament_orbits','relative_even_graph_orbits','tau_T','tau_E','difference','tournament_multiplicity_histogram','even_graph_multiplicity_histogram','escape_graph_multiplicity_histogram','tournament_fusion_defect','even_graph_fusion_defect','sign_escape_relative_orbits','sign_escape_global_types')
        rows.extend({k:r[k] for k in keys} for r in (pr,dr))
        for r in (pr,dr):print(n,r['class'],'R',r['relative_tournament_orbits'],'T',r['tau_T'],'E',r['tau_E'],'F_T',r['tournament_fusion_defect'],'F_E',r['even_graph_fusion_defect'],'Z',r['sign_escape_relative_orbits'],file=sys.stderr,flush=True)
    print(json.dumps({'status':'PASS','scope':'four selected p=7 product/diagonal embeddings at n=13,14; not a subgroup census','checks':m.CHECKS,'rows':rows},sort_keys=True,indent=2))
if __name__=='__main__':main()
