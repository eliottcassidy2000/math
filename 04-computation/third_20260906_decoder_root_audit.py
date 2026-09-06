#!/usr/bin/env python3
"""Independent prime-valuation and literal actual-entry audit; no producer import."""
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd, prod
from functools import reduce
from pathlib import Path
import json
from sympy import Matrix, factorint

GATES=0

def need(ok,why):
    global GATES
    GATES+=1
    if not ok:
        raise ArithmeticError(why)

def main():
    roots=0
    for n in range(2,9):
        for alpha in product(range(4),repeat=n):
            eta=min(alpha)
            for a in set(alpha):
                lhs=sum(min(a,b) for b in alpha)-a
                need(lhs<=eta+(n-2)*a,'universal prime-valuation mechanism')
                roots+=1
    Q=91**6
    M=(177,31684,6684150,1694040507,468424857663)
    S=(1,177,35483,8350122,2170260903)
    H=(13520696477,124998734,914513,5397,28)
    for a,m,s,h,G in zip(range(2,7),M,S,H,(2,4,9,30,90)):
        if a>2:
            rational=F(356**(a-1)*(a-2)**(a-2),(a-1)**(a-1))
            need(m==rational.numerator//rational.denominator,'inherited cofactor integer ceiling')
            need(s**(a-1)<m**(a-2)<=(s+1)**(a-1),'strict partner endpoint')
        need(h==min(Q//G,a*Q//(7*(14-a)*s)),'automatic cutoff both gates')
    den=(215,251,257,263,273)
    P=prod(den)
    V=tuple(sorted((P,)+tuple(83*P//q for q in den)))
    U=(2,3,4,8,16,19,28)
    t=3251
    physical=tuple(t*v for v in V)+U
    need(len(set(physical))==13 and reduce(gcd,physical)==1,'primitive actual row')
    need(sum(physical)<Q**2,'physical box')
    edges=[]
    for i,j in combinations(range(13),2):
        d=gcd(physical[i],physical[j]); height=(physical[i]+physical[j])//d
        if height<=356 and all(p%3==2 and e<=2 for p,e in factorint(height).items()):
            need((i<6)==(j<6),'no cross decoder edge')
            row=[0]*13;row[i]=physical[j]//d;row[j]=-physical[i]//d
            edges.append(row)
    need(Matrix(edges).rank()==11,'literal decoder rank11 with two components')
    mixed=0
    for a,b in combinations(V,2):
        for u in U:
            D=t*gcd(a,b)
            need(D//gcd(D,u)>Q,'every outside U coefficient forbidden by divisibility')
            mixed+=1
    for u,v in combinations(U,2):
        for a in V:
            need(t*a>Q*(u+v),'every outside V coefficient forbidden by amplitude')
            mixed+=1
    profiles=json.loads(Path(__file__).with_name('overnight12_20260906_lrc_decoder_descent_inherited_profiles.json').read_text())
    count=0
    for size in range(7,13):
        permitted={(d,tuple(word)) for d,word in profiles['levels'][str(13-size)]['profiles']}
        for subset in combinations(range(13),size):
            d=reduce(gcd,(physical[i] for i in subset))
            word=tuple(sorted(gcd(d,x) for j,x in enumerate(physical) if j not in subset))
            need((d,word) in permitted,'every retained complement word')
            count+=1
    phase=F(1301,6502)
    clearance=min(min((x*phase)%1,(-x*phase)%1) for x in physical)
    need(clearance==F(1289,6502)>F(1,14),'physical positive control')
    root=min(V)
    need(prod(gcd(root,x) for x in V if x!=root)==root**4,'root product equality control')
    partner=min(gcd(root,x) for x in V if x!=root)
    need(56*partner*max(U)<=6*Q,'selected pair native consumer')
    need(gcd(2,4)*gcd(2,6)>2,'common-scale deletion hostile')
    sharp=(177*178,177*179,178*179)
    need(min(gcd(a,b) for a,b in combinations(sharp,2))==177,'strict three-shape endpoint attained')
    print('INDEPENDENT_ROOTED_DIVISOR_AND_ACTUAL_ENTRY_AUDIT')
    print(f'prime_valuation_root_states={roots} mixed_supports={mixed} complete_profiles={count}')
    print(f'actual_decoder_rank=11 phase={phase} clearance={clearance}')
    print(f'PASS gates={GATES}')

if __name__=='__main__':
    main()
