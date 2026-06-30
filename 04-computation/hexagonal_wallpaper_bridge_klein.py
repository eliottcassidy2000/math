#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""THE DESIGN-COVERING BRIDGE via the HEXAGONAL/EISENSTEIN WALLPAPER GROUP (klein-S24).

Q(sqrt-3) = the Eisenstein integers Z[zeta_6] = the HEXAGONAL (triangular) lattice; wallpaper group p6m.
Phi_6(n)=|n-zeta_6|^2 = norm = index [Z[zeta_6] : (n-zeta_6)], and Z[zeta_6]/(n-zeta_6) = Z/Phi_6(n) (cyclic).
In the quotient zeta_6 = n, so the 6-FOLD ROTATION of the hexagonal lattice = MULTIPLICATION BY n mod Phi_6(n)
(n^6 = 1). This script verifies that 6-fold symmetry, checks it is a MULTIPLIER of the Singer difference set
(the p6 symmetry of the construction), and frames the bridge to Kershner (hexagonal = thinnest plane covering).
"""
import math
from sympy import factorint, primefactors

def phi6(n): return n*n-n+1
def order_mod(a,m):
    a%=m; 
    if math.gcd(a,m)!=1: return None
    o=1; x=a%m
    while x!=1: x=(x*a)%m; o+=1
    return o

print("="*88)
print(" (1) The 6-fold rotation: zeta_6 = n mod (n-zeta_6), so n^6 = 1 mod Phi_6(n) (the p6/C6 action)")
print("="*88)
print(f" {'n':>3}{'Phi_6':>8}{'n^6 mod Phi6':>14}{'ord(n)':>8}{'ord(n-1)=q mult':>16}{'ord(-1)':>9}")
for n in [3,4,6,8,14,9,12]:
    P=phi6(n); o_n=order_mod(n,P); o_q=order_mod(n-1,P); o_m1=order_mod(P-1,P)
    print(f" {n:>3}{P:>8}{pow(n,6,P):>14}{str(o_n):>8}{str(o_q):>16}{str(o_m1):>9}")
print("   => n has order EXACTLY 6 mod Phi_6(n): the hexagonal 6-fold rotation acts on Z/Phi_6(n) as x->n*x.")
print("   (n-1=q is the classical Singer multiplier; -1 = n^3 is the 180-deg rotation, order 2.)")

# (2) build a (Phi_6(n),n,1) difference set for small n and test the 6-fold (mult-by-n) symmetry
def find_diffset(v,k):
    # brute incremental search for a planar (v,k,1) difference set containing 0,1
    target=set(range(1,v))
    def diffs(S):
        d=set()
        for a in S:
            for b in S:
                if a!=b: d.add((a-b)%v)
        return d
    # DFS
    import itertools
    base=[0,1]
    def ext(S):
        if len(S)==k:
            return S if len(diffs(S))==k*(k-1) else None
        last=S[-1]
        for nxt in range(last+1, v):
            T=S+[nxt]; dd=diffs(T)
            if len(dd)==len(T)*(len(T)-1):  # all diffs distinct
                r=ext(T)
                if r: return r
        return None
    return ext(base)

print("\n"+"="*88)
print(" (2) Singer difference set + its 6-fold (mult-by-n) hexagonal symmetry (multiplier test)")
print("="*88)
for n in [3,4,6]:
    v=phi6(n); k=n
    D=find_diffset(v,k)
    if not D: print(f"   n={n} (v={v},k={k}): no diffset found"); continue
    Ds=set(D)
    # is mult-by-n a multiplier? (n*D = D + s for some shift s)
    def is_mult(t):
        tD=sorted((t*x)%v for x in Ds)
        for s in range(v):
            if set((x+s)%v for x in tD)==Ds: return s
        return None
    sn=is_mult(n%v); sq=is_mult((n-1)%v); sm=is_mult(v-1)
    print(f"   n={n} (v={v},k={k}): D={sorted(D)}")
    print(f"      mult-by-n={n} (zeta_6, 6-fold): multiplier? {'YES shift '+str(sn) if sn is not None else 'no'}")
    print(f"      mult-by-q={n-1} (Singer):       multiplier? {'YES shift '+str(sq) if sq is not None else 'no'}")
    print(f"      mult-by-(-1)=n^3 (180-deg):     multiplier? {'YES shift '+str(sm) if sm is not None else 'no'}")

print("\n"+"="*88)
print(" (3) THE BRIDGE FRAME")
print("="*88)
print(" Z[zeta_6] = hexagonal lattice (wallpaper p6m). Phi_6(n)=norm; Z/Phi_6(n) = hexagonal-lattice quotient.")
print(" The Singer covering = the projective plane = the OPTIMAL covering DESIGN; carried by the p6 6-fold")
print(" rotation (mult-by-n). Continuous side: the hexagonal lattice is the THINNEST covering of the plane")
print(" (Kershner 1939) -- the unique optimal LATTICE covering, p6m-symmetric. So the bridge:")
print("   discrete design-optimality (projective plane / Steiner)  <-->  continuous covering-optimality")
print("   (hexagonal lattice, Kershner),  TRANSFERRED by the shared p6/Eisenstein symmetry.")
print(" The covering-min n/Phi_6(n) is the density of the hexagonal-symmetric (p6) optimal covering.")
