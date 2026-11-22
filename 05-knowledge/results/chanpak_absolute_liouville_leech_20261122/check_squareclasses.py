#!/usr/bin/env python3
"""Exact Pell non-collision check on the inherited literal +2 digit outputs."""
from itertools import product
from math import isqrt

BASE=490032
LOW=360532
HIGH=362956
SLOPE=63618252914778017089536
CONST=47124685853114258307047

def need(ok,message):
    if not ok:
        raise ArithmeticError(message)

def val(word):
    v=0
    for x in word:
        v=BASE*v+x
    return SLOPE*v+CONST

def main():
    need(100*HIGH<101*LOW,"literal narrow interval")
    checked=0
    for n in range(1,9):
        vals=[val(w) for w in product((LOW,HIGH),repeat=n)]
        need(len(vals)==len(set(vals))==2**n,"unique digit word")
        X=[v+1 for v in vals]
        need(100*max(X)<101*min(X),"Pell interval")
        A=[v*(v+2) for v in vals]
        for i,a in enumerate(A):
            for j,b in enumerate(A):
                if j>=i:
                    break
                u=a*b
                t=isqrt(u)
                need(t*t!=u,"paired product squareclass collision")
                checked+=1
    print(f"{checked} exact unequal-squareclass tests on actual +2 counts through depth eight")
if __name__=="__main__":
    main()
