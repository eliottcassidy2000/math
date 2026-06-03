#!/usr/bin/env python3
"""lrc_doubling_sporadic_mod3_s609.py -- the doubling-sporadic law for LRC.

From the n=14 enumeration (S608): the sporadic tight config V* is AP[12->24],
the unique on-wall doubling. This session GENERALIZES: for even n, the doubling
    AP{1..n-1} with (n-2) -> 2(n-2)
is LRC-tight (M = 1/n) for some n and not others. The pattern is exact:

  CONJECTURE (doubling-sporadic mod-3 law) [VERIFIED n=4..32]:
     AP[(n-2)->2(n-2)] is tight  <=>  n = 2 (mod 3)  <=>  3 | (2n-1).
  Hits: n = 8 (->{1..5,7,12}), 14 (V*=AP[12->24]), 20, 26, 32, ...

MECHANISM: the pair-sum pinch modulus is 2n-1 (THM-401). When 3 | (2n-1) the
runners that are multiples of 3 form a shell (S592 'prime 3 at n=14', S593
sporadic shells); n-2 = 0 (mod 3) and 2(n-2) = 0 (mod 3) both lie in that shell,
so doubling (n-2)->2(n-2) preserves the antipodal binding pairs {a, n-a} that
hold the wall M=1/n (S608). When 3 does not divide 2n-1, no such shell exists and
the doubling either loosens (M>1/n) or is unavailable.

Session: claude-2026-06-03-S609 (lrc-doubling-sporadic-mod3).
"""
import sys; sys.stdout.reconfigure(line_buffering=True)
from math import gcd
from functools import reduce
from fractions import Fraction as F
from itertools import combinations
def lcm(a,b): return a*b//gcd(a,b)
def is_tight(V):
    n=len(V); L=reduce(lcm,V); D=(n+1)*L; ivs=[]
    for v in V:
        Lv=L//v
        for j in range(v+1):
            lo=((j*(n+1)-1)*Lv)%D; hi=lo+2*Lv
            if hi<=D: ivs.append((lo,hi))
            else: ivs.append((lo,D)); ivs.append((0,hi-D))
    ivs.sort(); pos=0
    for s,e in ivs:
        if s>pos: return False
        if e>pos: pos=e
    return pos>=D
print("Doubling sporadic AP[(n-2)->2(n-2)] tight?  vs  n mod 3 and 3|(2n-1):")
print(f"  {'n':>3} {'n%3':>4} {'3|(2n-1)?':>10} {'doubling tight?':>16}")
for N in range(4,33,2):
    AP=list(range(1,N)); a=N-2
    V=tuple(sorted([x for x in AP if x!=a]+[2*a]))
    t=is_tight(V)
    print(f"  {N:>3} {N%3:>4} {str((2*N-1)%3==0):>10} {str(t):>16}  {'<-- MATCH' if (t==((2*N-1)%3==0)) else '*** MISMATCH'}")
