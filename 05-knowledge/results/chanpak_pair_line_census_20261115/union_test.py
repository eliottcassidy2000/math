#!/usr/bin/env python3
"""Pointwise compare the exact all-TT bank with enlarged physical hosts."""
import sys
from repair import coeffs
cap=59999
wanted=set(range(69,cap+1,2))
profiles=set()
for line in open(sys.argv[1]):
    n,h,p,c,s=line.split()
    n,h,p=map(int,(n,h,p))
    if q:=p+2*h:
        if q<cap: profiles.add((q,h))
profiles=sorted(profiles)
reader=set()
for line in open(sys.argv[2]):
    n,c=map(int,line.split())
    if n<=6:
        for (r,M),mark in coeffs(n,c):
            assert r%2==1 and M>0 and r>=1-M
            reader.add((r,M))
hits=set()
for r,M in reader:
    for q,h in profiles:
        if M*q>2*cap:break
        if (N:=r*h+M*q) in wanted:hits.add(N)
old=set()
for line in open(sys.argv[3]):
    N,tt,small,anti,tri=map(int,line.split())
    if (tt or anti) and N in wanted:old.add(N)
u=hits|old
print("all orders, host cap six:",len(hits),"of",len(wanted))
print("all orders, union with arbitrary TT and prior anti-readers:",len(u),"of",len(wanted))
missing=sorted(wanted-u)
print("misses:",missing)
