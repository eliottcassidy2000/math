#!/usr/bin/env python3
"""Independent score-parity quotient direct masks and Krawtchouk verification."""
import ast
import itertools
import json
from fractions import Fraction
from math import comb
from pathlib import Path

CHECKS=0

def check(ok,message):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(message)

def choose(n,k):return comb(n,k) if 0<=k<=n else 0

def K(n,a,r):return sum((-1)**j*choose(r,j)*choose(n-r,a-j) for j in range(a+1))

def rank(A):
    A=[[Fraction(x) for x in row] for row in A];r=0
    for j in range(len(A[0]) if A else 0):
        p=next((i for i in range(r,len(A)) if A[i][j]),None)
        if p is None:continue
        A[p],A[r]=A[r],A[p]
        z=A[r][j];A[r]=[x/z for x in A[r]]
        for i in range(len(A)):
            if i!=r:
                z=A[i][j];A[i]=[x-z*y for x,y in zip(A[i],A[r])]
        r+=1
    return r

def main():
    rows=[]
    for n in range(33):
        m=comb(n,2)
        weights=[r for r in range(n+1) if r%2==m%2]
        dual=[a for a in range(n//2+1) if not (2*a==n and a%2)]
        check(len(weights)==len(dual),'cut/score type counts')
        A=[[0]*len(weights) for _ in weights]
        for i,r in enumerate(weights):
            for j,s in enumerate(weights):
                A[i][j]=comb(r,2) if s==r-2 else comb(n-r,2) if s==r+2 else r*(n-r) if s==r else 0
            check(sum(A[i])==m,'weighted degree')
        B=[[K(n,a,r) for a in dual] for r in weights]
        check(rank(B)==len(weights),'full Krawtchouk eigenbasis')
        spectrum=[2*a*(n-a) for a in dual]
        for i,r in enumerate(weights):
            for j,a in enumerate(dual):
                check(sum(A[i][k]*B[k][j] for k in range(len(weights)))==(m-spectrum[j])*B[i][j],'Krawtchouk diagonalization')
        if n>=3:check(min(spectrum[1:])==2*(n-1) and len(set(spectrum))==len(spectrum),'simple gap and full simple spectrum')
        else:check(spectrum==[0],'small boundary')
        if n%2==0:
            middle=[K(n,n//2,r) for r in weights]
            check((not any(middle))==(n//2%2==1),'balanced cut obstruction')
        if n<=5:
            edges=list(itertools.combinations(range(n),2))
            seen=set()
            fibres={}
            p_of={}
            for t in range(1<<m):
                p=0
                for e,(u,v) in enumerate(edges):p^=1<<(v if t>>e&1 else u)
                check(p.bit_count()%2==m%2,'score parity total')
                seen.add(p);p_of[t]=p
                fibres[p]=fibres.get(p,0)+1
                direct={r:0 for r in weights}
                for u,v in edges:direct[(p^(1<<u)^(1<<v)).bit_count()]+=1
                r=p.bit_count();i=weights.index(r)
                check([direct[s] for s in weights]==A[i],'direct endpoint flip multiplicities')
            check(seen=={p for p in range(1<<n) if p.bit_count()%2==m%2},'all allowed parity vectors realized')
            check(len(set(fibres.values()))==1,'uniform cycle-space fibres')
            for t in range(1<<m):
                for e,(u,v) in enumerate(edges):
                    check(p_of[t^(1<<e)]==p_of[t]^(1<<u)^(1<<v),'direct tournament edge reversal')
        rows.append({'n':n,'allowed_score_parity_weights':weights,'admissible_cut_shores':dual,'laplacian_eigenvalues':spectrum})
    tree=ast.parse(Path(__file__).read_text())
    check(not any(isinstance(x,ast.Assert) for x in ast.walk(tree)),'no assert nodes')
    check(not any(isinstance(x,ast.Constant) and isinstance(x.value,float) for x in ast.walk(tree)),'no float constants')
    print(json.dumps({'status':'PASS','checks':CHECKS,'universe':'Direct all labelled tournaments n=0..5; all Krawtchouk eigenbasis entries and exact rank n=0..32','rows':rows},indent=2,sort_keys=True))

if __name__=='__main__':main()
