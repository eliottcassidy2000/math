#!/usr/bin/env python3
"""Standard-library exact controls for sparse square-sum trees and cycles."""
from math import isqrt
from collections import Counter
import json
C=Counter()
def check(b,name):
    C[name]+=1
    if not b:raise RuntimeError(name)
def ceilroot(n):
    r=isqrt(n);return r+(r*r<n)
def square(n):return n>=0 and isqrt(n)**2==n
def append(labels,parent):
    M=max(labels);b=labels[parent];r=ceilroot(M+b+1);c=r*r-b
    check(c>M,'greedy_strict_growth');check(b+c==r*r,'greedy_edge')
    labels.append(c)
def close(labels):
    A,B=labels[0],labels[-1];M=max(labels);r=ceilroot(32*M)
    r+=(r-A-B-1)%2
    s=(A+B+1-r*r)//2;v=r*r-B;w=s*s-A
    check((r-A-B-1)%2==0,'closing_parity')
    check(B+v==r*r and v+w==(s-1)**2 and w+A==s*s,'closing_edges')
    check(M<v<w<900*M*M,'closing_distinctness_and_bound')
    return labels+[v,w]
def close_max_endpoint(labels):
    M=max(labels)
    check(labels[-1]==M,'maximum_at_endpoint')
    A=M;B=labels[0];K=ceilroot(M);r=ceilroot(M+B+16*K+1)
    r+=(r-A-B-1)%2
    delta=r*r-A-B-1;s=-delta//2;v=r*r-B;w=s*s-A
    check(16*K<=delta<40*K,'linear_closure_delta_bound')
    check(B+v==r*r and v+w==(s-1)**2 and w+A==s*s,'linear_closure_edges')
    check(M<v<w<1600*M,'linear_closure_distinctness_and_bound')
    return labels+[w,v]
seed=[482,3362,6242]
check(all(square(seed[i]+seed[(i+1)%3]) for i in range(3)),'seed_triangle')
cycle4=seed+[2222]
check(all(square(cycle4[i]+cycle4[(i+1)%4]) for i in range(4)),'exact_four_cycle')
# Exhaust every endpoint pair for small maxima, even without a path between them.
for M in range(2,51):
    for A in range(1,M+1):
        for B in range(1,M+1):
            if A!=B:close([A,M,B])
# Prefix controls for the exact seed, including substantially longer cycles.
labels=seed[:];samples={}
for N in range(5,1006):
    cycle=close_max_endpoint(labels)
    check(len(cycle)==N and len(set(cycle))==N,'cycle_vertex_count_and_distinctness')
    check(cycle[:3]==seed,'exact_seed_retained')
    check(all(square(cycle[i]+cycle[(i+1)%N]) for i in range(N)),'complete_cycle_edges')
    check(max(cycle)<1600*(80+2*(N-5))**2,'global_quadratic_numeric_bound')
    if N in (5,6,10,100,1000):samples[str(N)]={'max_label':max(cycle),'max_bits':max(cycle).bit_length()}
    append(labels,len(labels)-1)
# A balanced binary tree, a star, and a path use the same greedy rule.
for kind in ('binary','star','path'):
    labels=[1]
    for i in range(1,1000):
        parent={'binary':(i-1)//2,'star':0,'path':i-1}[kind]
        append(labels,parent)
        check(max(labels)<=(1+2*i)**2,'tree_quadratic_numeric_bound')
print(json.dumps({'status':'FINITE-EXACT elementary sparse bounds; standard library only','checks':dict(C),'cycle_samples':samples},indent=2))
