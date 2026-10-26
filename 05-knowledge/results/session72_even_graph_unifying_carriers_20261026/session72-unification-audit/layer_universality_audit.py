#!/usr/bin/env python3
"""Independent exact set-family/prime-layer universality controls."""
from itertools import combinations,product
from math import prod
import argparse,json
PRIMES=(3,5,7,11,13,17,19,23)
def need(ok,why):
    if not ok:raise RuntimeError(why)
def connected(S,E):
    if not S:return True
    reached={next(iter(S))}
    while True:
        new=reached|{v for u,v in E if u in reached and v in S}|{u for u,v in E if v in reached and u in S}
        if new==reached:return reached==S
        reached=new

def main():
    p=argparse.ArgumentParser();p.add_argument('--inject',action='store_true');args=p.parse_args()
    families=cases=0
    for r in range(1,4):
        subsets=[{i for i in range(r) if mask>>i&1} for mask in range(2**r)]
        pairs=list(combinations(range(r),2))
        for flags in product(range(2),repeat=len(subsets)):
            family=[S for S,b in zip(subsets,flags) if b];families+=1
            primes=PRIMES[:len(family)]
            ds=[prod(p for p,S in zip(primes,family) if i in S) for i in range(r)]
            m=prod(primes) if family else 3
            need(m>=3 and m%2==1 and all(m%d==0 for d in ds),'odd native ambient period')
            actual=[{i for i,d in enumerate(ds) if d%p==0} for p in primes]
            need(actual==family,'literal layer realization')
            need(all(d%(p*p)!=0 for d in ds for p in primes),'no additional higher prime-power layers')
            for edgeflags in product(range(2),repeat=len(pairs)):
                E=[e for e,b in zip(pairs,edgeflags) if b];cases+=1
                need(all(connected(S,E) for S in family)==all(connected(S,E) for S in actual),'all graph/family constraints transported')
    ds=[15,15,1]
    need([{i for i,d in enumerate(ds) if d%p==0} for p in (3,5)]==[{0,1},{0,1}],'duplicate set realization')
    E=[(0,1),(0,2)];S={1,2}
    need(connected({0,1,2},E) and not connected(S,E),'unweighted graph-connectedness hostile')
    if args.inject:need(connected(S,E),'deliberate ordinary-connectivity conflation')
    need((families,cases)==(276,2084),'bounded exhaustive counts')
    print(json.dumps({'status':'PASS','families':families,'graph_family_cases':cases,'r_range':[1,2,3],'duplicate_sets_verified':True,'scope':'representation universality only; no hardness or polynomial-size native construction claim'},sort_keys=True,indent=2))
if __name__=='__main__':main()
