#!/usr/bin/env python3
"""Triangular homogeneous-context tests and the full-profile H37 lift.
The quotient-word DP and literal tournament path enumeration are independent.
"""
from functools import lru_cache
from math import factorial
from pathlib import Path
import importlib.util,json

_spec=importlib.util.spec_from_file_location('reset_literal',Path(__file__).with_name('metagraph_reset_probe.py'))
m=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(m)


def tester(T,r):
    n=len(T);N=n+2*r;R=list(T)+[0]*(2*r)
    for v in range(n):
        for i in range(r):R[v]|=1<<(n+2*i);R[n+2*i+1]|=1<<v
    for i in range(n,N):
        for j in range(i+1,N):R[i]|=1<<j
    return tuple(R)


def coefficient(r,k):
    N=2*r;full=(1<<N)-1
    @lru_cache(None)
    def go(mask,left,last):
        if mask==full and left==0:return 1
        total=0
        # -2 means before the first letter; -1 means the module letter.
        if left and (last==-2 or last>=0 and last%2==1):total+=go(mask,left-1,-1)
        for j in range(N):
            if mask>>j&1:continue
            if last==-2 or last==-1 and j%2==0 or last>=0 and j>last:
                total+=go(mask|1<<j,left,j)
        return total
    return go(0,k,-2)


def main():
    matrix=[]
    for r in range(5):
        row=[coefficient(r,k) for k in range(1,r+3)]
        m.need(row[-1]==0 and row[-2]==factorial(r),'triangular return-gap coefficient')
        matrix.append(row[:-1])
    tests=0
    for n in range(1,5):
        for code in range(1<<(n*(n-1)//2)):
            T=m.decode(n,code);c=m.cover_profile(T)
            for r in range(3 if n<=3 else 2):
                U=tester(T,r);expected=sum(coefficient(r,k)*c[k-1] for k in range(1,n+1))
                m.need(len(m.paths(U))==expected,'literal context word expansion')
                tests+=1
    A=(40,17,3,6,13,30)
    word=((4,2,0,5,3,1),(2,0,5,3,4,1),(2,0,5,4,1,3),(0,1,5,4,2,3))
    graphs=[A]
    for p in word:graphs.append(m.native(graphs[-1],p))
    c=m.cover_profile(A);lift=[]
    for r in (1,2,3):
        values=[m.H(tester(T,r)) for T in graphs]
        expected=sum(coefficient(r,k)*c[k-1] for k in range(1,7))
        m.need(set(values)=={expected},'all five block-supported steps preserve contextual H')
        if r==1:
            m.need(len(m.paths(tester(graphs[0],r)))==len(m.paths(tester(graphs[-1],r)))==expected,'literal eight-vertex endpoint check')
        lift.append({'return_pairs':r,'order':6+2*r,'H_at_each_step':values})
    actual=(37,lift[0]['H_at_each_step'][0],lift[1]['H_at_each_step'][0])
    guard=lambda h:(h[2]-10*h[1]+21*h[0])%12
    m.need(guard(actual)==0 and guard((actual[0],actual[1],actual[2]+2))!=0,'source-cover factorial lattice guard')
    return {'status':'PROVED triangular context mechanism; FINITE-EXACT coefficient and literal controls',
            'coefficient_rows_r0_to4':matrix,'small_literal_contexts':tests,'H37_uniform_lifts':lift,
            'context_mod12_guard':{'actual_queries':actual,'actual_residue':guard(actual),'fake_last_plus2_residue':guard((actual[0],actual[1],actual[2]+2))},
            'scope':'Uniform module contexts; block-supported resets are not whole-graph native moves.'}


if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
