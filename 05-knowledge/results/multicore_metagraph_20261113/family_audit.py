"""Exact word-matrix certificate, protected-interface maps and literal counts."""
from pathlib import Path
from itertools import product
from collections import Counter
import sys,json,gzip
import sympy as sp
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'berggren_target_return_20261111'/'reset'))
import funnel_six_state as f
from rooted_growth import need,save,endpoints,relabel
from clock_probe import birth

def build(recipe,k):
    C=(0,) if recipe['order']==1 else f.a.unrank(recipe['order'],recipe['rank'])
    word=tuple(recipe['left'])+(recipe['middle'],)*k+tuple(recipe['right'])
    P,p=f.a.funnel(C,[0],word)
    T=birth(P,p);r=len(T)-1
    module=tuple(range(len(C)+len(recipe['left']),len(C)+len(recipe['left'])+k))
    return T,module,p,r

def direct(T,r):
    h=sum(endpoints(T));rev=tuple(((1<<len(T))-1)^(1<<i)^row for i,row in enumerate(T))
    start=endpoints(rev)[r]
    deleted=relabel(T,tuple(i for i in range(len(T)) if i!=r))
    return h,sum(endpoints(deleted))-start

def run():
    data=json.loads(gzip.decompress((HERE/'family_search.json.gz').read_bytes()))
    bank={(r['A'],r['B'],r['middle']):r for r in data['families']}
    states=[(a,b) for a in range(3) for b in range(2)];matrices=[]
    for bit in (0,1):
        M=sp.zeros(6)
        for j,(last,started) in enumerate(states):
            M[states.index((1 if bit else 2,started)),j]+=1
            if started or not bit:M[states.index((last,1)),j]+=1
        need(M**3-3*M**2+2*M==sp.zeros(6),'all-length word annihilator')
        need(all(x>=0 for x in M**2-M),'nonnegative exponential coefficient')
        matrices.append([list(map(int,M.row(i))) for i in range(6)])
    checked=0;coefficient_checks=0
    # Later points are controls of the matrix proof, not its justification.
    for r in bank.values():
        for k in (1,2,5,11):
            w=tuple(r['left'])+(r['middle'],)*k+tuple(r['right'])
            state=f.advance(r['state'],w)
            need(state[0]+state[4]-1==r['A']*(1<<k)+r['B'],'full six-state target')
            coefficient_checks+=1
    for a,b,aa,bb,d,bit in data['edges']:
        src,dst=bank[a,b,bit],bank[aa,bb,bit]
        need(aa==3*a and bb==3*b+d,'affine edge')
        for k in (1,2,4):
            T,M,p,r=build(src,k);U,N,q,s=build(dst,k)
            old=M+(p,r);new=N+(q,s)
            need(relabel(T,old)==relabel(U,new),'all protected comparisons agree')
            need([i for i,row in enumerate(T) if row>>r&1]==[p],'source sole predecessor')
            need([i for i,row in enumerate(U) if row>>s&1]==[q],'target sole predecessor')
            # Canonical unique-start recognition is independent of H equality.
            parent=f.a.induced(U,list(range(s)))
            path=f.a.unique_start_path(parent,q)
            need(path is not None,'literal next-family readiness witness')
            checked+=1
    selected=[(462,-5,1),(1386,-13,1),(4158,-37,1),(45,7,1),(135,23,1),(455,59,1),(1365,179,1),(4095,539,1)]
    literal=[]
    for key in selected:
        r=bank[key]
        for k in (1,2):
            T,M,p,z=build(r,k);h,defect=direct(T,z)
            need(h==r['A']*(1<<k)+r['B'] and defect==1,'independent literal whole-graph count and readiness')
            literal.append(dict(key=key,k=k,rows=T,module=M,predecessor=p,root=z,H=h,defect=defect))
    candidates=sorted({r['B'] for r in bank.values() if r['A']==405})
    need(71 not in candidates,'405 target excluded within complete bound')
    # Rank-table enumeration contains every core allowed by G>=2^n+1.
    need(data['core_order']>=8 and data['fixed']>=8,'all-order 405 reduction bounds covered')
    save('family_audit',dict(scope='Word-matrix proof checks, all retained recipe observers, every same-interface edge at k1,2,4, independent selected graph counts',matrices=matrices,coefficient_checks=coefficient_checks,protected_edge_checks=checked,edges=len(data['edges']),plus2_edges=sum(e[4]==2 for e in data['edges']),literal=literal,A405_offsets=candidates))
    print('PASS',coefficient_checks,'observer checks',checked,'protected maps',len(literal),'literal counts; A405 offsets',candidates)
if __name__=='__main__':run()
