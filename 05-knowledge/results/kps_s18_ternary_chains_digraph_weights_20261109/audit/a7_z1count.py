#!/usr/bin/env python3
"""AUDIT: labelled order-5 digraphs with exactly one missing pair and no reciprocal pair, counted by H (cross-check of
Session234 census/general_census.csv row n=5, H=7, absent=1, reciprocal=0: 600)."""
import itertools
n=5
pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
cnt={}
for z in pairs:
    others=[p for p in pairs if p!=z]
    for bits in range(1<<len(others)):
        A=[[0]*n for _ in range(n)]
        for t,(i,j) in enumerate(others):
            if bits>>t&1: A[i][j]=1
            else: A[j][i]=1
        h=sum(1 for p in itertools.permutations(range(n)) if all(A[p[k]][p[k+1]] for k in range(n-1)))
        cnt[h]=cnt.get(h,0)+1
print("order 5, exactly one missing pair, no reciprocal pair: labelled count with H=7:",cnt.get(7,0),"; with H=21:",cnt.get(21,0))
