#!/usr/bin/env python3
"""Complete profile audit by permutations, independent of subset-cover DP."""
from collections import defaultdict
import json
from pathlib import Path
import sys
from audit_sources import literal,actual_state
from bounded_sources import search,word,graph
from audit_power_sources import need,path_counts,append_zero
HERE=Path(__file__).resolve().parent

def keyset(bank):return {tuple(x['profile']) for x in bank}

def main():
    bank=json.loads((HERE/'strong_profiles256.json').read_text())['profiles']
    other=json.loads((HERE/'strong_profiles256.ears.json').read_text())['profiles']
    need(keyset(bank)==keyset(other),'two complete strong-ear profile covers')
    for filename,Q in [('strong_profiles.json',64),('strong_profiles128.json',128)]:
        small=json.loads((HERE/filename).read_text())['profiles']
        selected=lambda b:b['profile'][0]<Q//2 and b['profile'][3]<=Q and b['profile'][6]<=2*Q
        need(keyset([b for b in bank if selected(b)])==keyset([b for b in small if selected(b)]),
             'fixed-path versus nauty profile completeness')
    groups=defaultdict(list)
    for b in bank:groups[tuple(b['rows'])].append(b)
    checked=0;covers_count=0
    for j,(rows,entries) in enumerate(groups.items()):
        paths,covers=literal(rows);H=len(paths);E=sum(covers.values());covers_count+=E
        for b in entries:
            mask=b['mask'];inside=lambda v:v>=0 and bool(mask>>v&1)
            p=sum(z for (u,v),z in covers.items() if u>=0 and not inside(u))
            q=sum(z for (u,v),z in covers.items() if inside(v))
            s=sum(inside(p[0]) for p in paths)
            t=sum(not inside(p[-1]) for p in paths)
            F=sum(z for (u,v),z in covers.items() if (u<0 or not inside(u)) and (v<0 or inside(v)))
            need([H,E,p,q,s,t,F]==b['profile'],'independent literal profile')
            need(q>=mask.bit_count()*H and E>=1<<len(rows),'two injection bounds')
            checked+=1
        if j%100==0:print(f'literal parents {j}/{len(groups)}',file=sys.stderr,flush=True)
    prev,transitions=search(bank,256,127,512)
    census=json.loads((HERE/'positive256.json').read_text())
    need((len(prev),transitions)==(census['states'],census['transitions']),'saved closure')
    kernel_checks=0
    for state in prev:
        F,q,g,h=state
        if q:
            rows=graph(bank,word(prev,state))
            need(len(rows)-1<=(h*h).bit_length()-1+q.bit_length()-1+1,'logarithmic kernel bound')
            kernel_checks+=1
    direct=0;minima={}
    for b in census['positive_sources']:
        state=tuple(b['state']);rows=graph(bank,b['word'])
        need(actual_state(rows)==state,'actual dyadic source')
        c=b['c'];need(c==state[1]-state[0]-state[2]>0,'deficit')
        minima[c]=min(minima.get(c,10**9),state[1]);direct+=1
    need(9 not in minima,'all-order q<=256 deficit-nine exclusion')
    # Independently recounted upper witness, with x at vertex0.
    parent=[12,121,2,100,45,69,21];mask=67
    rows=[mask<<1]+[(r<<1)|(0 if mask>>v&1 else 1) for v,r in enumerate(parent)]
    need(actual_state(rows)==(449,512,54,117),'deficit-nine q512 upper witness')
    need(path_counts(append_zero(rows,0))[0]==1015,'actual optimal first entry')
    reroutes={}
    for N in [11,23,47,55,119,247]:
        rec=census['arbitrary_slope_one_zero_targets'][str(N)]
        rows=graph(bank,rec['word']);state=actual_state(rows)
        need(state==tuple(rec['state']),'reroute actual state')
        need(path_counts(append_zero(rows,0))[0]==N,'reroute actual output')
        reroutes[N]={'state':state,'rows':rows,'terminal_vertices':len(rows)+1}
    print(json.dumps({'status':'passed','profile_covers_equal':True,'literal_parent_graphs':len(groups),
          'literal_profiles':checked,'literal_covers':covers_count,'kernel_checks':kernel_checks,
          'closure_states':len(prev),'closure_transitions':transitions,'actual_dyadic_sources':direct,
          'exact_minimum_dyadic_slopes_through256':minima,'deficit9_optimal_q':512,
          'deficit9_optimal_first_output':1015,'arbitrary_slope_reroutes':reroutes,
          'scope':'Complete bounded-cost classification plus the all-order source reduction; direct adjacency audits of every positive dyadic state and six hostile reroutes. Universal spectrum completeness remains open.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
