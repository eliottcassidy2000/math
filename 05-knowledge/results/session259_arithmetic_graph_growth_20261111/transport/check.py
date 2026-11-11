#!/usr/bin/env python3
"""Focused exact controls for the new causal sum-product graph, not old canon."""
import json
from math import isqrt
from pathlib import Path


def require(test, message):
    if not test:
        raise RuntimeError(message)


def path_to_one(u):
    if u == 1:
        return [1]
    if u == 2:
        return [2, 1]
    if u == 3:
        return [3, 8, 4, 2, 1]
    if u == 4:
        return [4, 2, 1]
    if u % 2 == 0:
        return [u] + path_to_one(u // 2)
    return [u, 4*(u-1), 2*(u-1)] + path_to_one(u-1)


def adjacent_shifted(u, v):
    s, p = sorted([u+4, v+4])
    delta = s*s-4*p
    if delta <= 0:
        return False
    d = isqrt(delta)
    return d*d == delta and (s-d) % 2 == 0 and (s-d)//2 >= 2


def main():
    bound = 500
    pairs = [(a,b) for a in range(2,bound) for b in range(a+1,bound)
             if a+b <= bound]
    for a,b in pairs:
        s,p = a+b,a*b
        require(p>s, 'causality')
        require(s*s-4*p == (b-a)**2, 'discriminant')
        require(p-4 == a*(s-4)-(a-2)**2, 'affine square defect')
    for n in range(5,bound+1):
        products = [a*(n-a) for a in range(2,(n+1)//2)]
        require(len(products) == (n-1)//2-1, 'outdegree')
        require(min(products) == 2*n-4, 'first completion')
        upper = (n*n-1)//4 if n%2 else n*n//4-1
        require(max(products) == upper, 'last completion')
    divisor_sum = 0
    for N in range(1,bound+1):
        divisor_sum += sum(N%d == 0 for d in range(1,N+1))
        if N>=2:
            pending = sum(a+b<=N<a*b for a,b in pairs)
            formula = (N-1)**2//4-(N-2)-(divisor_sum-2*N+2-isqrt(N))//2
            require(pending==formula, f'backlog at {N}')
    max_length = 0
    for u in range(1,2001):
        path = path_to_one(u)
        require(path[-1] == 1, 'path endpoint')
        require(all(adjacent_shifted(a,b) for a,b in zip(path,path[1:])), 'path legality')
        require(max(path)<=max(8,4*u), 'path excursion')
        require(len(path)-1<=4*u.bit_length()+4, 'log path bound')
        max_length = max(max_length,len(path)-1)
    receipt = {'status':'FINITE-EXACT PASS','parent_sum_bound':bound,
               'pairs':len(pairs),'paths':2000,'max_path_length':max_length,
               'examples':{'birth10_products':[16,21,24],
                           'shifted3_path':path_to_one(3),
                           'shifted13_path':path_to_one(13)},
               'scope':'Exact controls of new transport formulas; no external conjecture or inherited-source audit.'}
    print(json.dumps(receipt,indent=2,sort_keys=True))

if __name__=='__main__':
    main()
