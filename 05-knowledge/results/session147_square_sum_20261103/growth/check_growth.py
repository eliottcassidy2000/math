#!/usr/bin/env python3
"""Exact finite checks supporting the all-n elementary growth proof."""
from math import isqrt

def require(ok,msg):
    if not ok: raise AssertionError(msg)

def older(n):
    if type(n)!=int or n<1: raise ValueError('n must be a positive integer')
    return {k*k-n for k in range(isqrt(n)+1,isqrt(2*n-1)+1)}

def components(n):
    adj={i:set() for i in range(1,n+1)}
    for i in adj:
        for j in older(i): adj[i].add(j);adj[j].add(i)
    todo=set(adj); count=0
    while todo:
        count+=1; stack=[todo.pop()]
        while stack:
            for y in adj[stack.pop()] & todo:
                todo.remove(y);stack.append(y)
    return count

def main():
    for bad in (0,-1,True,4.5):
        try:older(bad)
        except ValueError:pass
        else:raise AssertionError('malformed input accepted')
    require(all(components(n)>1 for n in range(2,14)), 'earlier disconnectedness')
    require(components(14)==1, 'G14 disconnected')
    require(sum(len(older(n)) for n in range(1,15))==13,'G14 tree')
    require(older(19)=={6,17}, 'G19 neighbors')
    rank=0
    for n in range(15,10001):
        d=len(older(n)); require(d>=1,'new vertex lacks older neighbor')
        rank+=d-1
        if n<=18:require(rank==1,'G15-G18 not unicyclic')
        if n>=19:require(rank>=2,'cycle rank decreased')
        if n<=22:print('G%d old neighbors %s, cycle rank %d'%(n,sorted(older(n)),rank))
    for n in range(5,10001):
        k=3 if n<9 else isqrt(n)+1
        require(n<k*k<=2*n-1,'square window failure')
    print('Checked graph growth and square window through n=10000; all checks passed.')
if __name__=='__main__': main()
