#!/usr/bin/env python3
"""Universal complete norm2 two-root distribution from directed c3,c4."""
from independent_audit import mat,need
from independent_interval_chains import rootpair,c4
import itertools,random,json
from collections import Counter
from pathlib import Path
OUT=Path(__file__).resolve().parent

def check_matrix(z):
    n=len(z);A=n*(n-1)*(n-2)
    c3=sum(all(sum(z[v][w]==1 for w in s)==1 for v in s) for s in itertools.combinations(range(n),3));four=c4(z)
    expected=Counter()
    for g in(-2,2):expected[g,0]=n*(n-1)
    for g in(-1,1):
        for w in(-3,3):expected[g,w]=6*c3
        for w in(-1,1):expected[g,w]=A-6*c3
    for w in(-4,4):expected[0,w]=4*four
    for w in(-2,2):expected[0,w]=A*(n-3)//6+8*c3*(n-3)-16*four
    expected[0,0]=2*A*(n-3)//3-16*c3*(n-3)+24*four
    need(all(x>=0 for x in expected.values()),('negative predicted rootpair',n))
    need(rootpair(z)==+expected,('universal rootpair formula',n,c3,four))

def main():
    rows=[]
    for n in range(1,6):
        for b in range(1<<(n*(n-1)//2)):check_matrix(mat(n,b))
        rows.append(dict(n=n,tournaments=1<<(n*(n-1)//2),exhaustive=True))
    rng=random.Random(360052)
    for n in range(6,11):
        for _ in range(32):check_matrix(mat(n,rng.getrandbits(n*(n-1)//2)))
        rows.append(dict(n=n,tournaments=32,exhaustive=False))
    report=dict(status='PASS',rows=rows)
    (OUT/'rootpair_formula_receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
