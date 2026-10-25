#!/usr/bin/env python3
"""Portable exhaustive controls for the proved labelled prime-power counts."""
import argparse
import itertools
import json
from pathlib import Path
P=Path(__file__).resolve().parent
LINES=((0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6))
def check(x,msg):
    if not x: raise RuntimeError(msg)
def grid(c,a,b,m):
    return tuple(x%m for x in (a,3*c-a-b,b,c-a+b,c,c+a-b,2*c-b,a+b-c,2*c-a))
def magic(g,m): return len({sum(g[i] for i in line)%m for line in LINES})==1
def ordinary(p,h):
    m=p**h
    if (p==2 and h<2) or (p==3 and h==0): return 0
    return m*(m-4)*(m-8) if p==2 else m*(m-3)*(m-5)
def unit(p,k):
    shift=3 if p==2 else 1
    return ordinary(p,k-shift) if k>=shift else 0
def count(p,k): return sum(unit(p,j) for j in range(k,0,-2))
def census_square(m):
    Q={x*x%m for x in range(m)}
    n=0
    for c,a,b in itertools.product(Q,repeat=3):
        g=grid(c,a,b,m)
        if len(set(g))==9 and all(x in Q for x in g):
            check(magic(g,m),'direct eight-line equations')
            n+=1
    return n
def run():
    # No parametrization assumption: brute force all 3^9 grids and compare sets.
    raw={g for g in itertools.product(range(3),repeat=9) if magic(g,3)}
    param={grid(c,a,b,3) for c,a,b in itertools.product(range(3),repeat=3)}
    check(raw==param,'full raw grid control at modulus 3')
    ordinary_rows=[]
    for p,H in ((2,9),(3,6)):
        for h in range(H+1):
            m=p**h
            pairs=sum(len(set(grid(0,u,v,m)))==9 for u,v in itertools.product(range(m),repeat=2))
            n=m*pairs
            check(n==ordinary(p,h),'ordinary census formula')
            ordinary_rows.append({'p':p,'h':h,'modulus':m,'valid_offset_pairs':pairs,'labelled_ordinary_count':n})
    square_rows=[]
    for p,K in ((2,9),(3,5)):
        for k in range(1,K+1):
            n=census_square(p**k)
            check(n==count(p,k),'direct square census recurrence')
            square_rows.append({'p':p,'k':k,'modulus':p**k,'labelled_square_count':n})
    # Multiply power series by claimed denominators using exact integers.
    gf_checks=[]
    for p,shift,numerator in ((2,7,{7:1536}),(3,3,{3:216,4:5832})):
        coeff=[count(p,k) for k in range(31)]
        product=coeff[:]
        for degree,value in ((2,1),(1,p),(1,p*p),(1,p**3)):
            old=product[:]
            product=[old[k]-(value*old[k-degree] if k>=degree else 0) for k in range(len(old))]
        check(all(x==numerator.get(k,0) for k,x in enumerate(product)),'GF exact coefficient identity')
        gf_checks.append({'p':p,'verified_through_degree':30,'numerator':numerator,'coefficients_k0_to_12':coeff[:13]})
    return {'scope':'Labelled distinct square RESIDUES; zero allowed; no integer witness or external novelty claim.','raw_mod3_magic_count':len(raw),'ordinary_censuses':ordinary_rows,'square_censuses':square_rows,'generating_function_controls':gf_checks,'normalized_limits':{'p2':'1/504','p3':'27/728'}}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    out=json.dumps(run(),indent=2,sort_keys=True)+'\n'
    if args.write: (P/'certificate.json').write_text(out)
    else: check((P/'certificate.json').read_text()==out,'certificate exact match')
    print('PASS: raw magic equations, ordinary and square censuses, exact generating-function controls')
if __name__=='__main__': main()
