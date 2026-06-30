#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""X_0(14) cusps = the Klein four-group under Atkin-Lehner; cusp<->LRC-core map; genus jump (klein-S10).

Moves 1-2: (a) modular data of X_0(N) for the LRC(2p) family -- #cusps, widths, genus; verify
genus(X_0(6))=0 (LRC(6), solved) and genus(X_0(14))=1 (LRC(14), first hard); (b) the 4 cusps of X_0(14)
<-> divisors {1,2,7,14} <-> subsets of {2,7}; Atkin-Lehner W_2,W_7,W_14 = toggles => W(14)=(Z/2)^2 acts
REGULARLY on the 4 cusps => the cusps ARE the Klein four-group = the n=4 tournament classes (THM-584);
(c) map each cusp to its LRC core type + width.
"""
from math import gcd, prod

def divisors(n): return [d for d in range(1, n+1) if n % d == 0]
def psi(n):
    r = n
    for p in set(pf(n)): r = r*(p+1)//p
    return r
def pf(n):
    f=[];d=2;m=n
    while d*d<=m:
        while m%d==0:f.append(d);m//=d
        d+=1
    if m>1:f.append(m)
    return f
def num_cusps(n):
    from math import gcd
    # Euler phi
    def phi(m):
        r=m
        for p in set(pf(m)): r=r//p*(p-1)
        return r if m>0 else 0
    return sum(phi(gcd(d, n//d)) for d in divisors(n))
def nu2(n):
    if n%4==0: return 0
    return sum(1 for x in range(n) if (x*x+1)%n==0)
def nu3(n):
    if n%9==0: return 0
    return sum(1 for x in range(n) if (x*x+x+1)%n==0)
def genus(n):
    g = 1 + psi(n)/12 - nu2(n)/4 - nu3(n)/3 - num_cusps(n)/2
    return round(g)

print("="*82)
print(" X_0(N) modular data for the LRC(2p) family (klein-S10)")
print("="*82)
print(f" {'N=2p':>5} {'psi':>4} {'#cusps':>7} {'nu2':>4} {'nu3':>4} {'genus':>6}   role")
roles={6:'LRC(6) apex 3 -- SOLVED',10:'(apex 5, not Mers/Heeg)',14:'LRC(14) apex 7 -- FIRST HARD/OPEN',
       22:'(apex 11)',26:'(apex 13)'}
for N in [6,10,14,22,26]:
    print(f" {N:>5} {psi(N):>4} {num_cusps(N):>7} {nu2(N):>4} {nu3(N):>4} {genus(N):>6}   {roles.get(N,'')}")
print(" => genus(X_0(2p)) = 0,0,1,2,? for p=3,5,7,11,13; the genus JUMPS 0->1 exactly at N=14")
print("    (among LRC-relevant apices {3,7}: genus 0 = LRC(6) solved, genus 1 = LRC(14) first hard).")

# ---- N=14: the 4 cusps, widths, Atkin-Lehner = Klein four-group ----
N=14; ds=divisors(N)
print("\n"+"="*82); print(" X_0(14): 4 cusps <-> divisors {1,2,7,14} <-> subsets of {2,7}; widths N/c"); print("="*82)
widths={c: N//c for c in ds}   # cusp of denominator c|N has width N/c (gcd(c,N/c)=1)
core={1:'d=1  cusp at oo   : BULK (unbounded-far / interior)',
      2:'d=2  cusp        : DOUBLING (2-adic descent, W_2)',
      7:'d=7  cusp        : APEX-7 (speed-7, 7a/14=a/2 even/odd coupling) = HARD',
      14:'d=14 cusp at 0   : FULL-DENSE (m_R->0, all mod-7 resonance = covering)'}
for c in ds:
    print(f"   cusp d={c:>2}  width={widths[c]:>2}   {core[c]}")
print(f"   sum of widths = {sum(widths.values())} = psi(14) = {psi(14)}  (index check)")

# Atkin-Lehner: cusp labeled by subset S of primes {2,7} (d=prod S); W_q toggles q in S.
prims=[2,7]
def label(d): return frozenset(p for p in prims if d%p==0)
cusps=[label(d) for d in ds]
def Wq(q, S): return S ^ {q}            # symmetric difference toggle
print("\n Atkin-Lehner action on cusps (W_q toggles prime q):")
for q,name in [(2,'W_2'),(7,'W_7')]:
    perm={d: prod([p for p in (label(d)^{q})] or [1]) for d in ds}
    print(f"   {name}: " + ", ".join(f"d{d}->d{perm[d]}" for d in ds))
# W_14 = W_2 W_7
perm14={}
for d in ds:
    S=label(d); S2=(S^{2})^{7}; perm14[d]=prod([p for p in S2] or [1])
print(f"   W_14=W_2 W_7 (Fricke): " + ", ".join(f"d{d}->d{perm14[d]}" for d in ds))
print(" => W(14)={1,W_2,W_7,W_14} = (Z/2)^2 acts REGULARLY (simply transitively) on the 4 cusps")
print("    => THE 4 CUSPS ARE THE KLEIN FOUR-GROUP = the n=4 tournament classes {T,+,-,S} (THM-584).")
print("\n Klein dictionary (THM-584: complement = coordinate swap; SC diagonal vs NS swapped pair):")
print("   d=1 (bulk)  <-> T (00, source&sink)        [identity]")
print("   d=2 (W_2)   <-> + (10, source killed)       [2-adic / doubling involution]")
print("   d=7 (W_7)   <-> - (01, sink killed)         [apex-7 involution] -- THE HARD CUSP")
print("   d=14 (W_14) <-> S (11, neither = strong)    [Fricke / full complement]")
print("   The 2-adic descent peels W_2; the residual apex is the W_7 cusp; the binding doublet sits there.")
