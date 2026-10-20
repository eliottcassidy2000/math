#!/usr/bin/env python3
"""Cross-check the census by direct all-fixed-object classification and Burnside.
Unlike the main path, full type classification does not first divide by N(H).
Graph sign detection is run on every resulting full graph type, including ones
excluded as N-odd by the main path. Prime-degree automorphisms are brute-force
checked in the main script. No external package is needed.
"""
import json
from collections import Counter
import mixed_symmetry_census as c

def main():
    source=json.load(open(__file__.replace('independent_replay.py','census.json')))
    lookup={(r['n'],r['class']):r for r in source['rows']};out=[]
    for n,name,gens in c.cases():
        H=c.closure(gens,n);pairs,index,orbits,base=c.pair_data(n,H)
        N=c.normalizer(n,H,gens);prep=[c.prep(p,pairs,index) for p in N]
        TV=c.masks(orbits,base);GV=c.masks(orbits)
        trs,tt=c.classify([{'mask':x} for x in TV],n,pairs,True)
        grs,gt=c.classify([{'mask':x} for x in GV],n,pairs,False)
        e=sum(c.full_graph_sign(A)[0] for A in grs)
        bt=sum(sum(c.move(x,r,True)==x for x in TV) for r in prep)
        be=sum(sum(c.sign(x,r) for x in GV if c.move(x,r)==x) for r in prep)
        c.need(bt%len(N)==be%len(N)==0,'Burnside integer')
        ref=lookup[n,name]
        c.need(len(tt)==ref['tau_T'] and e==ref['tau_E'],'direct full-cube count differs')
        c.need(bt//len(N)==ref['relative_tournament_orbits'],'Burnside tournament count differs')
        c.need(be//len(N)==ref['relative_even_graph_orbits'],'signed Burnside graph count differs')
        relabel=tuple(reversed(range(n)));rev=c.inverse(relabel)
        conjugated=[c.compose(c.compose(relabel,g),rev) for g in gens]
        rr=c.census_case(n,name,conjugated)
        for field in ['normalizer_order','pair_orbits','tau_T','tau_E','relative_tournament_orbits','tournament_multiplicity_histogram','even_graph_multiplicity_histogram','sign_escape_relative_orbits']:
            c.need(json.loads(json.dumps(rr[field]))==ref[field],'deterministic relabeling changed invariant '+field)
        out.append({'n':n,'class':name,'direct_full_tournament_types':len(tt),'direct_full_graph_types':len(gt),'direct_full_even_graph_types':e,'ordinary_Burnside_numerator':bt,'signed_Burnside_numerator':be,'normalizer_order':len(N),'relative_orbits_each':bt//len(N),'reversed_vertex_labeling_check':'PASS'})
    print(json.dumps({'status':'PASS','checks':c.CHECKS,'rows':out},indent=2,sort_keys=True))
if __name__=='__main__':main()
