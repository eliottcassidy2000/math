#!/usr/bin/env python3
"""Independent exhaustive Aut(R), vertex orbit and Hamilton-cycle card control."""
from itertools import permutations
from collections import Counter
from pathlib import Path
from audit import matrix,scores,extend,require
import json

def flat(A):return tuple(x for r in A for x in r)
def permuted(A,order):return [[A[x][y] for y in order] for x in order]
def canonical(A):return min(flat(permuted(A,order)) for order in permutations(range(len(A))))
def unflat(code):
    n=int(len(code)**.5);return [list(code[i*n:(i+1)*n]) for i in range(n)]
def automorphisms(A):
    return [q for q in permutations(range(len(A))) if permuted(A,q)==A]
def orbits(perms,n):
    todo=set(range(n));ans=[]
    while todo:
        u=min(todo);orb={q[u] for q in perms};todo-=orb;ans.append(sorted(orb))
    return ans

def run(n):
    p=n//2;words=Counter();inner=Counter();examples={}
    for w in range(1<<((n-1)*(n-2)//2)):
        A=matrix(n,w);d=scores(A)
        if set(d)<={p-1,p}:
            c=canonical(A);words[c]+=1;inner[c]+=d[0]==p and d[-1]==p-1;examples[c]=A
    parents={};cards=[]
    for code in words:
        T=unflat(code);a=len(automorphisms(T));E=extend(T);parent=canonical(E)
        # Independently enumerate actual paths, not the word census.
        h=down=0
        for q in permutations(range(n)):
            if all(T[q[i]][q[i+1]] for i in range(n-1)):
                h+=1;down+=sum(T[q[0]])==p and sum(T[q[-1]])==p-1
        require(h==a*words[code] and down==a*inner[code],'actual H/Aut word formula')
        cards.append({'card_flat':code,'parent_flat':parent,'H':h,'Aut':a,'H_down':down,'word_mass':words[code],'internal_degree':inner[code]})
        if parent in parents:continue
        R=unflat(parent);N=n+1;G=automorphisms(R);vo=orbits(G,N)
        cycles=0
        for q in permutations(range(1,N)):
            c=(0,)+q+(0,)
            cycles+=all(R[x][y] for x,y in zip(c,c[1:]))
        reports=[]
        for O in vo:
            z=O[0];survivors=[v for v in range(N) if v!=z]
            card=canonical(permuted(R,survivors))
            autcard=sum(q[z]==z for q in G)
            require(card in words,'deletion card occurs in fixed path universe')
            require(len(automorphisms(unflat(card)))==autcard,'actual Aut extension equals stabilizer')
            require(inner[card]*len(G)==len(O)*cycles,'orbit proportionality')
            require(inner[card]*autcard==cycles,'card cycles/Aut')
            reports.append({'vertex_orbit':O,'Aut_card':autcard,'card_flat':card,'internal_degree':inner[card]})
        require(len({tuple(r['card_flat']) for r in reports})==len(vo),'distinct vertex orbits give nonisomorphic cards')
        parents[parent]={'parent_flat':parent,'Aut_parent':len(G),'Hamilton_cycles':cycles,'vertex_orbits':reports,'total_internal_degree':sum(x['internal_degree'] for x in reports)}
    require(sum(P['total_internal_degree'] for P in parents.values())==sum(inner.values()),'all parent orbit sums exhaust internal words')
    return {'n':n,'ordinary_cards':len(words),'regular_parents':len(parents),'cards':cards,'parents':list(parents.values())}
if __name__=='__main__':
    import argparse
    pa=argparse.ArgumentParser();pa.add_argument('--output',default='parent-cards.json');args=pa.parse_args()
    result={'status':'PASS','controls':[run(n) for n in (2,4,6)]}
    Path(args.output).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result,indent=2,sort_keys=True))
