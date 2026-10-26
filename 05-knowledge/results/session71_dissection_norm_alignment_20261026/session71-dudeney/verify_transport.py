#!/usr/bin/env python3
"""Exact checks for an abstract length-transfer lemma, NOT a dissection solver."""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import argparse,json,hashlib

def need(p,m):
    if not p: raise RuntimeError(m)

def compile_chain(word):
    # Start unrestricted; require a positive input even for the empty word.
    eps,c,lo,hi=1,Q(0),Q(0),None
    for L in map(Q,word):
        # 0 < eps*x+c < L
        a,b=(-c,L-c) if eps==1 else (c-L,c)
        lo=max(lo,a);hi=b if hi is None else min(hi,b)
        if hi<=lo: return None
        eps,c=-eps,L-c
    return eps,c,lo,hi

def evaluate(state,x):
    if state is None: return None
    eps,c,lo,hi=state
    return eps*x+c if x>lo and (hi is None or x<hi) else None

def direct(word,x):
    if x<=0:return None
    for L in map(Q,word):
        if not 0<x<L:return None
        x=L-x
    return x

def compose(f,g):
    # g after f
    if f is None or g is None:return None
    e,c,lo,hi=f;d,b,jlo,jhi=g
    if e==1:a,k=jlo-c,None if jhi is None else jhi-c
    else:a,k=(None if jhi is None else c-jhi),c-jlo
    if a is not None:lo=max(lo,a)
    if k is not None:hi=k if hi is None else min(hi,k)
    if hi is not None and lo>=hi:return None
    return d*e,d*c+b,lo,hi

def normalize(s):
    return None if s is None else [str(t) if t is not None else None for t in s]

def pairs(left,right):return set(product(left,right))

# Transcribed finite case labels from HTML Sec5.1/5.3, not independently generated topologies.
S={14:['A1'],13:['A2','D1','D2'],12:['A3','A4','B1','D3','D4','D6','D7','F1','F2','F3','F4','F11'],11:['A5','B2','D5','D8','D9','D11','D12','E1','E2','E3','E4','F5','F6','F7','F9'],10:['A6','B3','C1','D10','D13','D14','E5','E6','E7','F8','F10'],9:['D15','E8']}
T={12:['D1'],11:['D2','F1'],10:['D3','E1']}
tri_survive=['B2','D5','E1','E2','E3','E4']
filtered=pairs(T[12],S[12])|pairs(T[11],tri_survive)|pairs(T[10],S[10])
A=pairs(['D1'],['D6','F1','F2','F4','F11','F3'])|pairs(['D2','F1'],['B2','E1','E3'])|pairs(['D3','E1'],['F8'])
B=pairs(['D1'],['A3','D7'])|pairs(['D2','F1'],['E2','E4'])
C=pairs(['D1'],['A4','B1','D3','D4'])|pairs(['D2','F1'],['D5'])|pairs(['D3','E1'],['A6','B3','C1','D10','D13','D14','E5','E6','E7','F10'])

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out');args=ap.parse_args()
    samples=[Q(i,4) for i in range(-4,25)]
    chain_count=sample_count=composition_count=0
    words=[w for k in range(7) for w in product(range(1,5),repeat=k)]
    for w in words:
        f=compile_chain(w);chain_count+=1
        for x in samples:
            need(evaluate(f,x)==direct(w,x),f'evaluation {w} {x}');sample_count+=1
        for j in range(len(w)+1):
            need(compose(compile_chain(w[:j]),compile_chain(w[j:]))==f,f'composition {w} {j}');composition_count+=1
    need(evaluate(compile_chain([3,4]),Q(3,2))==Q(5,2),'order hostile1')
    need(evaluate(compile_chain([4,3]),Q(3,2))==Q(1,2),'order hostile2')
    need(direct([3,4,Q(3,2)],Q(3,2)) is None,'Boolean continuation reject')
    need(direct([4,3,Q(3,2)],Q(3,2)) is not None,'Boolean continuation accept')
    need(compile_chain([2,2])[:2]==compile_chain([3,3])[:2],'domain hostile same affine')
    need(evaluate(compile_chain([2,2]),Q(5,2)) is None,'domain hostile reject')
    need(evaluate(compile_chain([3,3]),Q(5,2))==Q(5,2),'domain hostile accept')
    need(len({s for g in S.values() for s in g})==44,'44 square labels')
    need(len({s for g in T.values() for s in g})==5,'5 triangle labels')
    need(A|B|C==filtered and not(A&B or A&C or B&C),'disjoint coverage')
    # Formal coefficient witness in basis (sigma,tau), exact integer arithmetic.
    symbolic_word=[(1,0),(0,1)]*3+[(0,1),(1,0)]*3
    sign=1;offset=(0,0);prefix=[]
    for total in symbolic_word:
        sign=-sign;offset=tuple(a-b for a,b in zip(total,offset));prefix.append((sign,offset))
    need((sign,offset)==(1,(0,0)),'formal identity despite empty domain')
    need(prefix[4]==(-1,(3,-2)),'fifth output 3sigma-2tau-x')
    need(243<256,'exact 3sigma<4: fourth powers')
    need(sum(x[0] for x in symbolic_word[::2])==sum(x[0] for x in symbolic_word[1::2])==3,'balanced sigma')
    need(sum(x[1] for x in symbolic_word[::2])==sum(x[1] for x in symbolic_word[1::2])==3,'balanced tau')
    # Component charge computation; all n,m in explicit finite universe.
    from math import gcd
    charge_checks=0
    for n in range(3,41):
      for m in range(3,41):
       if gcd(n,m)!=1:continue
       for a in range(n+1):
        for b in range(m+1):
         condition=(2*a*m-2*b*n)%(n*m)==0
         divisor_condition=a%(n//gcd(n,2))==0 and b%(m//gcd(m,2))==0
         need(condition==divisor_condition,'charge equivalence');charge_checks+=1
    r={'status':'PASS','scope':'Abstract positive-length reflection chains and transcribed case-list arithmetic; NOT polygon realization or exhaustive topology/case proof','chains':chain_count,'rational_evaluation_checks':sample_count,'composition_checks':composition_count,'coprime_charge_checks':charge_checks,'square_type_count':44,'triangle_type_count':5,'unfiltered_pair_count':220,'cc_only_pair_count':sum(len(v)*len(S[k]) for k,v in T.items()),'angular_filtered_pair_count':len(filtered),'case_counts':{'A':len(A),'B':len(B),'C':len(C)},'symbolic_empty_domain_hostile':{'word':'(sigma,tau)^3(tau,sigma)^3','formal_final':[sign,list(offset)],'fifth_prefix':[prefix[4][0],list(prefix[4][1])],'strict_bound_certificate':'(3sigma)^4=243<256=4^4'},'hostiles':{'ordered_3_4':normalize(compile_chain([3,4])),'ordered_4_3':normalize(compile_chain([4,3])),'identity_on_0_2':normalize(compile_chain([2,2])),'identity_on_0_3':normalize(compile_chain([3,3]))}}
    out=json.dumps(r,indent=2,sort_keys=True)+'\n'
    if args.out:Path(args.out).write_text(out)
    print(out,end='')
if __name__=='__main__':main()
