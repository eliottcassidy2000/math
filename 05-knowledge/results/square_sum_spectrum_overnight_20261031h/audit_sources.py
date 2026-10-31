#!/usr/bin/env python3
"""Independent permutation-cover profiles and literal adjacency state audit."""
from collections import Counter,defaultdict
from itertools import permutations
import json
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261031g'))
from audit_power_sources import need,path_counts,append_zero,validate
from bounded_sources import search,word,graph

def literal(rows):
    n=len(rows);paths=[];covers=Counter()
    for order in permutations(range(n)):
        wrong=[]
        for i,(u,v) in enumerate(zip(order,order[1:])):
            if not(rows[u]>>v&1):wrong.append(i+1)
            if len(wrong)>1:break
        if not wrong:paths.append(order)
        breaks=range(n+1) if not wrong else wrong if len(wrong)==1 else []
        for k in breaks:
            covers[(order[k-1] if k else -1,order[k] if k<n else -1)]+=1
    return paths,covers

def actual_state(rows):
    validate(rows)
    parent=[r>>1 for r in rows[1:]]
    h,ends=path_counts(parent) if parent else (1,[])
    g=sum(z for i,z in enumerate(ends,1) if rows[0]>>i&1)
    F,_=path_counts(rows)
    F1,_=path_counts(append_zero(rows,0))
    return F,F1-F-g,g,h

def main():
    bank=json.loads((HERE/'strong_profiles.json').read_text())['profiles']
    grouped=defaultdict(list)
    for b in bank:grouped[tuple(b['rows'])].append(b)
    count=0;cover_permutations=0
    for rows,entries in grouped.items():
        validate(rows);paths,covers=literal(rows)
        cover_permutations+=sum(covers.values())
        for b in entries:
            mask=b['mask'];inside=lambda v:v>=0 and bool(mask>>v&1)
            H=len(paths);E=sum(covers.values())
            p=sum(z for (u,v),z in covers.items() if u>=0 and not inside(u))
            q=sum(z for (u,v),z in covers.items() if inside(v))
            s=sum(inside(p[0]) for p in paths)
            t=sum(not inside(p[-1]) for p in paths)
            F=sum(z for (u,v),z in covers.items() if (u<0 or not inside(u)) and (v<0 or inside(v)))
            need([H,E,p,q,s,t,F]==b['profile'],'literal profile')
            need(q>=mask.bit_count()*H,'cut-before-v injection bound')
            count+=1
    seen,transitions=search(bank,64,31,128)
    positive=[];observed=[];maximum=0
    for state in sorted(seen):
        rows=graph(bank,word(seen,state));maximum=max(maximum,len(rows))
        need(actual_state(rows)==state,'independent actual adjacency state')
        F,q,g,h=state
        if q and not q&(q-1):
            observed.append(state)
            if q>F+g:positive.append(state)
    expected=[(19,64,8,9),(31,64,6,9),(41,64,6,11)]
    need(positive==expected,'complete positive-deficit states for q<=64')
    hits={}
    for state in observed:
        F,q,g,h=state
        for r in range(1,8):
            N=q*(1<<r)+F+g-q
            if N<=127 and N not in hits:
                rows=graph(bank,word(seen,state))
                for _ in range(r):rows=append_zero(rows,0)
                H,_=path_counts(rows)
                need(H==N,'terminal adjacency count')
                hits[N]={'state':state,'zeros':r,'rows':rows}
    need(47 not in hits,'47 route obstruction')
    need(sorted(hits)==[3,5,9,17,33,41,65,67,69,71,73,91,95,101,111],'complete small terminal list')
    print(json.dumps({'status':'passed','literal_profiles':count,'literal_parent_graphs':len(grouped),
       'literal_ordered_covers':cover_permutations,'actual_state_witnesses':len(seen),
       'actual_maximum_vertices':maximum,'transitions':transitions,'q_bound':64,'h_bound':31,'F_bound':128,
       'positive_deficit_states':positive,'complete_zero_tail_targets_3_to127':hits,
       'scope':'All retained profiles checked by permutation covers. Every finite automaton state is checked on actual adjacency; completeness across parent orders uses the SCC theorem and inherited strong floor.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
