#!/usr/bin/env python3
"""Independent partition-formula controls for relative tournament orbit counts."""
from fractions import Fraction
from math import factorial,gcd
from collections import Counter
from pathlib import Path
import json
CHECKS=0

def need(x,msg):
    global CHECKS
    CHECKS+=1
    if not x:raise RuntimeError(msg)

def odd_partitions(n,small=1):
    if n==0:yield ();return
    for a in range(small,n+1,2):
        for tail in odd_partitions(n-a,a):yield (a,)+tail

def colored_tournaments(n,c):
    ans=Fraction(0)
    for part in odd_partitions(n):
        z=1
        for a,m in Counter(part).items():z*=a**m*factorial(m)
        e=sum((a-1)//2 for a in part)+sum(gcd(a,b) for i,a in enumerate(part) for b in part[i+1:])
        ans+=Fraction(c**len(part)*2**e,z)
    need(ans.denominator==1,'Burnside sum not integral')
    return ans.numerator

def main():
    rows=json.loads(Path(__file__).with_name('extension_counts.json').read_text())['rows'];checked=[]
    for r in rows:
        n,name=r['n'],r['class'];expect=None
        if name=='F21_fixed':expect=colored_tournaments(n-7,2)
        if name=='C15_5plus3':expect=2*colored_tournaments(n-8,4)
        if name=='C21_7plus3':expect=4*colored_tournaments(n-10,4)
        if name=='F21_times_C3':expect=2*colored_tournaments(n-10,4)
        if name=='F21_diagonal_7plus3':expect=4*colored_tournaments(n-10,4)
        if name=='F55_fixed':expect=colored_tournaments(n-11,2)
        if expect is not None:
            need(expect==r['relative_tournament_orbits'],f'relative Burnside mismatch {n} {name}')
            checked.append({'n':n,'class':name,'relative_count':expect})
    need(len(checked)==12,'scope count')
    print(json.dumps({'status':'PASS','scope':'12 of the 18 census rows; relative orbit counts only, not full fusion/sign descent','checks':CHECKS,'two_colored_tournament_counts':[colored_tournaments(n,2) for n in range(7)],'four_colored_tournament_counts':[colored_tournaments(n,4) for n in range(5)],'rows':checked},sort_keys=True,indent=2))
if __name__=='__main__':main()
