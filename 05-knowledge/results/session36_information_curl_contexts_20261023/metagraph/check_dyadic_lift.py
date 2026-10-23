#!/usr/bin/env python3
"""Exact 2-adic Smith valuations via unit pivots modulo 2^32.
It fails rather than claiming completion if precision is exhausted.
"""
import json,time
from pathlib import Path
from check_dyadic_network import lap,need,val2
HERE=Path(__file__).resolve().parent

def dyadic_smith(A,bits=32):
    A=[list(map(int,row)) for row in A];m=len(A);shift=0;out=[]
    while m:
        need(bits>0,'insufficient 2-adic precision')
        mod=1<<bits
        A=[[v%mod for v in row] for row in A]
        rank=0
        for k in range(m):
            pivot=next(((i,j) for i in range(k,m) for j in range(k,m) if A[i][j]&1),None)
            if pivot is None:break
            i,j=pivot;A[k],A[i]=A[i],A[k]
            if j!=k:
                for row in A:row[k],row[j]=row[j],row[k]
            inv=pow(A[k][k],-1,mod)
            row=[v*inv%mod for v in A[k][k+1:]]
            for i in range(k+1,m):
                a=A[i][k]
                A[i][k+1:]=[(x-a*y)%mod for x,y in zip(A[i][k+1:],row)]
                A[i][k]=0
            rank+=1
        out.extend([shift]*rank)
        A=[[A[i][j]//2 for j in range(rank,m)] for i in range(rank,m)]
        need(all(0<=v<mod//2 for row in A for v in row))
        m-=rank;shift+=1;bits-=1
    return sorted(out)

if __name__=='__main__':
 out=[]
 for a in json.loads((HERE/'dyadic_network_receipt.json').read_text()):
    n=a['n'];t=time.time();L=lap(a['nodes'],a['blue_edges'])+2*lap(a['nodes'],a['paired_black_edges'])
    vals=dyadic_smith(L[:-1,:-1].tolist());need(sum(v>0 for v in vals)==a['critical_two_rank'])
    if n<=6:need(vals==[val2(v) for v in a['smith_invariants']])
    row=dict(n=n,valuation_histogram={str(k):vals.count(k) for k in sorted(set(vals))},spanning_tree_valuation2=sum(vals),two_primary_exponents=[v for v in vals if v],precision_bits=32)
    out.append(row);print({k:v for k,v in row.items() if k!='two_primary_exponents'},'seconds',round(time.time()-t,3),flush=True)
 (HERE/'dyadic_lift_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
 print('ALL CHECKS PASSED')
