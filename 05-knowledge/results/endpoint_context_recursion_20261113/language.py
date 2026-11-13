"""Exact bounded selector and density obstruction for the nine half-return matrices."""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import random
import json,gzip
from typed_glue import HERE,need,save
from port_compiler import mul
from repair23 import profile


def run():
    raw=json.loads(gzip.decompress((HERE/'renewal_catalog.json.gz').read_bytes()));rules=raw['records']
    graphs=((2,4,9,3),(2,8,3,5));zs=tuple(profile(T,3) for T in graphs)
    need(zs==((5,2,3,6),(5,1,1,5)),'small strong seed pair')
    factors=[sum(r['target_matrix'][0]) for r in rules]
    floors=[next(k for k in range(1,1000) if (k+1)**4>L**3) for L in factors]
    need(all(k**4<=L**3 for k,L in zip(floors,factors)),'Kraft rational upper bound')
    upper=sum(Fraction(1,k) for k in floors);need(upper<Fraction(3,4),'zero-density word bound')
    cap=100_000_000;seen={(5,9):dict(word=[],profiles=zs,n=4)};stack=[(5,9)]
    while stack:
        h,E=stack.pop();src=seen[h,E]
        for i,r in enumerate(rules):
            N=r['target_matrix'];H=N[0][0]*h+N[0][1]*E;F=N[1][0]*h+N[1][1]*E
            need(H>=factors[i]*h and F>=H,'growth cone')
            if H>cap or (H,F) in seen:continue
            states=tuple(mul(r['matrix'],z) for z in src['profiles'])
            need(all(z[0]==H and z[0]+z[3]-z[1]==F and z[3]+z[1]==z[0]+z[2] for z in states),'quotient retains observations and plane')
            d=tuple(y-x for x,y in zip(*states));need(d[0]==0 and d[2]==2*d[1] and d[3]==d[1] and d[1]!=0,'nonzero recursive repair')
            need(all(all(v>=0 for v in z) for z in states),'actual profile nonnegativity')
            seen[H,F]=dict(word=src['word']+[i],profiles=states,n=src['n']+r['n']-1);stack.append((H,F))
    targets=Counter(h for h,E in seen);obstruction=53084231
    # Independent target-decreasing membership, compared with forward closure.
    @lru_cache(None)
    def inverse(h,E):
        if (h,E)==(5,9):return True
        if h<5 or h%2==0 or E%2==0 or 2*E<3*h or E>2*h:return False
        for r in rules:
            (a,b),(c,d)=r['target_matrix'];det=a*d-b*c
            x=d*h-b*E;y=-c*h+a*E
            if x%det or y%det:continue
            p,q=x//det,y//det
            if 0<p<h and inverse(p,q):return True
        return False
    for h,E in seen:need(inverse(h,E),'independent decreasing positive selector')
    rng=random.Random(20261113);hostiles=0
    for _ in range(1000):
        h=2*rng.randrange(3,500000)+1;E=2*rng.randrange((3*h+3)//4,h+1)+1
        need(inverse(h,E)==((h,E) in seen),'independent membership control');hostiles+=1
    lengths={(5,9):0};word_counts={(5,9):1}
    for h,E in sorted(seen):
        for r in rules:
            (a,b),(c,d)=r['target_matrix'];child=a*h+b*E,c*h+d*E
            if child in seen:
                lengths[child]=max(lengths.get(child,0),lengths[h,E]+1)
                word_counts[child]=word_counts.get(child,0)+word_counts[h,E]
    examples=[dict(H=h,E=E,**seen[h,E]) for h,E in sorted(seen)[:25]]
    largest=max(seen);longest=max(seen,key=lambda k:len(seen[k]['word']))
    obj=dict(status='FINITE-EXACT complete language selector',scope='Every word in the nine <=8 half-return matrices starting from the specified H5 pair, with H<=100000000. Arbitrary context sizes and other directions excluded. Repeated target-pair states safely merged by proved quotient and renewable guard.',cap=cap,rules=rules,seed=dict(rows=graphs,root=3,profiles=zs),growth_factors=factors,kraft_denominators=floors,kraft_upper=str(upper),target_pairs=len(seen),distinct_H=len(targets),largest_H=max(targets),longest_word=max(lengths.values()),construction_words=sum(word_counts.values()),inverse_positive_controls=len(seen),inverse_random_controls=hostiles,obstruction_target=obstruction,obstruction_in_language=obstruction in targets,head=examples,largest_example=dict(H=largest[0],E=largest[1],**seen[largest]),pair_multiplicity_histogram=dict(Counter(targets.values())))
    save('language',obj)
    print('LANGUAGE',obj['target_pairs'],'pairs',obj['distinct_H'],'H','LONGEST',obj['longest_word'],'53084231',obj['obstruction_in_language'],'KRAFT',upper,flush=True)
    print('HEAD',[(x['H'],x['E']) for x in examples],flush=True)


if __name__=='__main__':run()
