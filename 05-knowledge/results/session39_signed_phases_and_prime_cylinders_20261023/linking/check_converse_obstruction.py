#!/usr/bin/env python3
"""Exact antiisometry obstruction and converse/reflection type controls."""
import argparse, hashlib, json
from pathlib import Path
HERE=Path(__file__).resolve().parent
CHECKS=0
def need(p,m):
    global CHECKS
    CHECKS+=1
    if not p:raise RuntimeError(m)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    numerators=set()
    for a in range(1,256,2):
        need(a*a%8==1,'every odd square modulo eight')
        for t in range(4):
            z=(a*a+64*t)%256;numerators.add(z)
            need(z%8==1,'all possible order256 self-numerators')
            need(z!=255,'no image can have self-link minus1/256')
    need(numerators==set(range(1,256,8)),'complete scalar necessary numerator set')
    need((-1)%8==7 and 7 not in {z%8 for z in numerators},'antiisometry contradiction')
    need(153*13**2%256==1,'Session38-bank generator normalization')
    need(105*27**2%256==1,'incoming-bank generator normalization')
    need(153*9**2%256==105,'cyclic abstract generator reconciliation only')
    need(1%2==(-1)%2,'positive self-negative order-two pairing control')
    controls=[]
    for n in range(3,9):
        free=[(i,j) for i in range(n) for j in range(i+2,n)];m=len(free);full=(1<<m)-1
        R=[free.index((n-1-j,n-1-i)) for i,j in free];closing=free.index((0,n-1))
        need(sorted(R)==list(range(m)),'reflection is a free-cell permutation')
        need(all(R[R[i]]==i for i in range(m)),'reflection involution')
        need(R[closing]==closing,'closing cell fixed under reflection')
        # A coordinate permutation commutes with flipping every coordinate.
        def reflect(x):return sum(((x>>i)&1)<<R[i] for i in range(m))
        test_masks=range(1<<m) if n<=7 else [0,full]+[1<<i for i in range(m)]+[full^(1<<i) for i in range(m)]
        tested=blue=black=0
        for x in test_masks:
            y=reflect(x);need(reflect(y)==x,'literal involution control')
            need(reflect(x^full)==(y^full),'reflection commutes with complement')
            need(y!=(x^full),'closing bit forbids swapping complement endpoints')
            fixed={y,y^full}=={x,x^full};need(fixed==(y==x),'edge fixed iff blue')
            blue+=y==x;black+=y!=x;tested+=1
        controls.append({'n':n,'free_cells':m,'tested_presentations':tested,'coverage':'all presentations' if n<=7 else 'zero/full and every singleton/co-singleton; full permutation identities proved above','blue_presentations_in_control':blue,'black_presentations_in_control':black})
    out={'status':'PASS','scope':'No antiisometry for the inherited native n8 linking form; converse fixes critical divisors but can swap literal black twins.','native_dependency':'Nondegenerate cyclic summand <1/256> with orthogonal complement exponent at most4, certified in the frozen linking packet.','order256_possible_self_numerators':sorted(numerators),'negative_target_numerator':255,'mod8_native':1,'mod8_reversed':7,'normalizing_multipliers':{'153':13,'105':27},'reflection_controls':controls,'checks':CHECKS,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    p=HERE/'converse_obstruction_receipt.json'
    if args.write:p.write_text(json.dumps(out,indent=2)+'\n')
    else:need(json.loads(p.read_text())==out,'frozen antiisometry receipt')
    print(json.dumps({'status':'PASS','checks':out['checks'],'native_unit_mod8':1,'orientation_reversed_unit_mod8':7,'antiisometry':False},sort_keys=True))
if __name__=='__main__':main()
