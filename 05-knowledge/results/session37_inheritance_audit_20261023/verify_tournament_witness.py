#!/usr/bin/env python3
"""Independent direct permutation check: no subset-DP implementation is shared.
Usage: python3 verify_tournament_witness.py --n 9 --masks 29759988540 18875161993
Coefficient arrays use ascending powers of x. Edge bits are upper-pair lex;
bit 1 means i -> j for i<j. Counts directed cycles up to cyclic rotation,
not reversal, directly from the permutation's disjoint-cycle decomposition.
"""
import argparse,json
from collections import Counter
from itertools import permutations
from pathlib import Path


def verify(n,mask):
    A=[[0]*n for _ in range(n)];b=0
    for i in range(n):
        for j in range(i+1,n):
            A[i][j]=(mask>>b)&1;A[j][i]=1-A[i][j];b+=1
    per=[0]*(n+1);detplus=[0]*(n+1);charminus=[0]*(n+1)
    packing=Counter();H=0
    for p in permutations(range(n)):
        H+=all(A[p[i]][p[i+1]] for i in range(n-1))
        if not all(i==p[i] or A[i][p[i]] for i in range(n)):continue
        seen=set();lengths=[]
        for i in range(n):
            if i in seen:continue
            j=i;k=0
            while j not in seen:seen.add(j);k+=1;j=p[j]
            if k>1:lengths.append(k)
        lengths=tuple(sorted(lengths));packing[lengths]+=1
        fixed=n-sum(lengths);per[fixed]+=1
        detplus[fixed]+=(-1)**(sum(lengths)-len(lengths))
        charminus[fixed]+=(-1)**len(lengths)
    H_ocf=sum((2**len(l))*count for l,count in packing.items() if all(k%2 for k in l))
    assert H==H_ocf
    # O_m / E_m are determined by per(xI+A) and det(xI-A).
    O=[(per[n-m]-charminus[n-m])//2 for m in range(n+1)]
    E=[(per[n-m]+charminus[n-m])//2 for m in range(n+1)]
    def o(m):return O[m] if m<=n else 0
    def e(m):return E[m] if m<=n else 0
    D35=packing[(3,5)];T333=packing[(3,3,3)]
    if n<=7:assert H==1+2*(o(3)+o(5)+o(7))+4*e(6)
    if n==8:assert H==1+2*(o(3)+o(5)+o(7))+4*(e(6)+D35)
    if n==9:assert H==1+2*(o(3)+o(5)+o(7)+o(9))+4*(e(6)+D35)+6*T333
    return {'n':n,'mask':mask,'adjacency':A,'H_direct_paths':H,'H_odd_cycle_formula':H_ocf,
            'per_xI_plus_A_ascending':per,'det_xI_plus_A_ascending':detplus,
            'det_xI_minus_A_ascending':charminus,'O_m':O,'E_m':E,
            'cycle_length_partition_counts':{str(l):v for l,v in sorted(packing.items())}}


def main():
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=9)
    p.add_argument('--masks',type=int,nargs='+',default=[29759988540,18875161993]);p.add_argument('--output',type=Path,default=Path(__file__).with_suffix('.json'));a=p.parse_args()
    out=[verify(a.n,m) for m in a.masks]
    if len(out)==2:
        assert out[0]['per_xI_plus_A_ascending']==out[1]['per_xI_plus_A_ascending']
        assert out[0]['det_xI_plus_A_ascending']==out[1]['det_xI_plus_A_ascending']
        assert out[0]['H_direct_paths']!=out[1]['H_direct_paths']
    a.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
