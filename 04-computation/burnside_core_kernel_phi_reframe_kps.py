#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Burnside core kernel, the phi/GCD-matrix reframe, and the metagraph enumerator family.
kind-pasteur-2026-06-15-S7.  Builds on the user's 1-tail-peeling compression of A000568.

A000568(n) = #tournaments on n nodes = sum_{lambda |- n, ALL PARTS ODD} 2^{e(lambda)}/z_lambda,
  e(lambda) = sum_i (lam_i-1)/2 + sum_{i<j} gcd(lam_i,lam_j)   (# free edge-orbits),
  z_lambda  = prod_k k^{m_k} m_k!                               (centralizer / stabilizer size).
USER COMPRESSION: split lambda = mu U 1^r (mu = odd parts >=3 = the CORE; 1^r = the tail):
  e(lambda) = e(mu) + C(r,2) + r*len(mu),   z_lambda = z_mu * r!.
  => a(n) = sum_{m,t} B[m,t] * 2^{C(n-m,2)+(n-m)t}/(n-m)!,   B[m,t] = sum_{mu:|mu|=m,len=t} 2^{e(mu)}/z_mu.

THIS SESSION adds:
 (1) THE phi-REFRAME (the gcd cross-term, isolated as the residual difficulty, made explicit):
       e(mu) = C(t,2) + (1/2) sum_{d odd >=3} phi(d) * M_d^2,   M_d = #parts of mu divisible by d.
     (Smith/GCD-matrix identity gcd = sum_{d|.} phi(d).)  A positive-definite quadratic form on the
     divisor-multiplicity lattice -> a theta-function exponent.  The residual difficulty (no clean
     (m,t) recurrence) is EXACTLY this cross-term: it needs the divisor profile {M_d}, not just (m,t).
 (2) THE EXACT add-a-part RECURRENCE (state = divisor profile): adding an odd part p>=3,
       Delta e = (p-1)/2 + sum_{d|p, d>=3} phi(d) M_d.
 (3) THE ENUMERATOR FAMILY (answering "what metrics am I missing"): the SAME cycle-index machinery
     with base x gives a 2-variable polynomial; x=2 is the tournament count, the odd/all-parts choice
     toggles orientation vs graph:
       T_n(x) = sum_{odd lambda} x^{e}/z   (x=2 -> A000568 tournaments; x=1 -> A000246/n! odd-cycle perms)
       G_n(x) = sum_{ALL  lambda} x^{e_graph}/z (x=2 -> A000088 graphs),  e_graph=sum floor(lam_i/2)+sum gcd
     Both peel the 1-tail identically; the phi-reframe applies to both.  P_n(x)=n!*T_n(x) is an integer poly.
"""
import sys, math
from fractions import Fraction as Fr
from collections import Counter
from functools import lru_cache
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def phi(d):
    r=d; dd=d; p=2
    while p*p<=dd:
        if dd%p==0:
            while dd%p==0: dd//=p
            r-=r//p
        p+=1
    if dd>1: r-=r//dd
    return r

def parts_gen(n, maxp, minp):
    """partitions of n into ODD parts in [minp,maxp], non-increasing."""
    if n==0:
        yield (); return
    p = maxp if maxp%2==1 else maxp-1
    while p>=minp:
        if p<=n:
            for rest in parts_gen(n-p,p,minp): yield (p,)+rest
        p-=2

def all_parts_gen(n, maxp):
    """partitions of n into ANY parts <= maxp, non-increasing."""
    if n==0:
        yield (); return
    for p in range(min(n,maxp),0,-1):
        for rest in all_parts_gen(n-p,p): yield (p,)+rest

def e_tour(parts):
    e=sum((p-1)//2 for p in parts)
    for i in range(len(parts)):
        for j in range(i+1,len(parts)): e+=math.gcd(parts[i],parts[j])
    return e
def e_graph(parts):
    e=sum(p//2 for p in parts)
    for i in range(len(parts)):
        for j in range(i+1,len(parts)): e+=math.gcd(parts[i],parts[j])
    return e
def e_phi_reframe(parts):
    t=len(parts); Mc=Counter()
    for p in parts:
        for d in range(3,p+1,2):
            if p%d==0: Mc[d]+=1
    s=sum(phi(d)*Mc[d]*Mc[d] for d in Mc)
    return math.comb(t,2)+Fr(s,2)
def zf(parts):
    c=Counter(parts); zz=1
    for p,mult in c.items(): zz*=(p**mult)*math.factorial(mult)
    return zz

A000568={1:1,2:1,3:2,4:4,5:12,6:56,7:456,8:6880,9:191536,10:9733056,11:903753248,12:154108311168}
A000246={1:1,2:1,3:3,4:9,5:45,6:225,7:1575,8:11025,9:99225,10:893025}  # all-odd-cycle perms
A000088={1:1,2:2,3:4,4:11,5:34,6:156,7:1044,8:12346}                   # graphs on n nodes

print("=== (1) phi-reframe e(mu)=C(t,2)+1/2 sum_{d>=3 odd} phi(d) M_d^2  (cores, parts>=3) ===")
fails=0;chk=0
for m in range(0,41):
    for mu in parts_gen(m,m,3):
        if e_phi_reframe(mu)!=e_tour(mu): fails+=1
        chk+=1
print(f"   checked {chk} cores (mass<=40): failures={fails}")

@lru_cache(None)
def Bmt(m,t):
    return sum((Fr(2**e_tour(mu),zf(mu)) for mu in parts_gen(m,m,3) if len(mu)==t), Fr(0))

def a_full(n):  return sum((Fr(2**e_tour(l),zf(l)) for l in parts_gen(n,n,1)), Fr(0))
def a_comp(n):
    s=Fr(0)
    for m in range(0,n+1):
        r=n-m
        for t in range(0,m//3+1):
            b=Bmt(m,t)
            if b: s+=b*Fr(2**(math.comb(r,2)+r*t),math.factorial(r))
    return s
print("\n=== (2) full Burnside == compressed == A000568 ===")
allok=True
for n in range(1,17):
    af,ac=a_full(n),a_comp(n); k=A000568.get(n)
    ok=(af==ac) and (k is None or int(af)==k); allok&=ok and af.denominator==1
    if n<=12 or not ok: print(f"   n={n:2d}: a={int(af):>22d}  full==comp:{af==ac}  known:{'OK' if (k is None or int(af)==k) else 'FAIL'}")
print(f"   ALL n<=16 integer & consistent: {allok}")

print("\n=== (3) exact add-a-part recurrence: Delta e = (p-1)/2 + sum_{d|p} phi(d) M_d  (d>=1; d=1 term = #parts) ===")
def e_incremental(parts):
    """build e by adding parts one at a time using the recurrence; compare to e_tour.
       Delta e = (p-1)/2 [internal] + sum_{d|p} phi(d) M_d [cross-gcd, all odd divisors incl d=1]."""
    e=0; Mc=Counter()
    for p in sorted(parts):
        cross=sum(phi(d)*Mc[d] for d in range(1,p+1) if p%d==0)   # all divisors d|p (p odd -> all odd)
        e += (p-1)//2 + cross
        for d in range(1,p+1):
            if p%d==0: Mc[d]+=1
    return e
rfails=0
for m in range(0,36):
    for mu in parts_gen(m,m,3):
        if e_incremental(mu)!=e_tour(mu): rfails+=1
print(f"   recurrence reproduces e on all cores (mass<=35): failures={rfails}")

print("\n=== (4) (m,t) state count vs odd-partition count (the user's ~533x collapse) ===")
def n_odd_parts(n): return sum(1 for _ in parts_gen(n,n,1))
def n_states(n):
    c=0
    for m in range(0,n+1):
        for t in range(0,m//3+1):
            if m>=3*t and (m-t)%2==0: c+=1   # achievable iff m>=3t and m≡t (mod 2)
    return c
for n in (20,40,60,100):
    nm=n_states(n)
    op=n_odd_parts(n) if n<=60 else None
    extra=f", odd partitions={op}, collapse={op/nm:.0f}x" if op else " (odd-partition count skipped; user reports 444793 -> 533x at n=100)"
    print(f"   n={n}: active (m,t) states={nm}{extra}")

print("\n=== (5) the enumerator family: T_n(x) [tournaments/orientations, odd lambda] & G_n(x) [graphs, all lambda] ===")
def T_poly(n):
    d={}
    for l in parts_gen(n,n,1):
        d[e_tour(l)]=d.get(e_tour(l),Fr(0))+Fr(1,zf(l))
    return d
def G_poly(n):
    d={}
    for l in all_parts_gen(n,n):
        d[e_graph(l)]=d.get(e_graph(l),Fr(0))+Fr(1,zf(l))
    return d
def ev(poly,x): return sum(c*Fr(x)**e for e,c in poly.items())
print("   n | T_n(1)*n!=#odd-cyc-perms(A000246) | T_n(2)=tournaments(A000568) | T_n(3)*? | G_n(2)=graphs(A000088)")
for n in range(1,9):
    Tp=T_poly(n); Gp=G_poly(n)
    t1=ev(Tp,1)*math.factorial(n); t2=ev(Tp,2); t3=ev(Tp,3); g2=ev(Gp,2)
    chk246 = '=A000246 OK' if int(t1)==A000246.get(n,-1) else ''
    chk568 = '=A000568 OK' if int(t2)==A000568.get(n,-1) else ''
    chk088 = '=A000088 OK' if int(g2)==A000088.get(n,-1) else ''
    print(f"   {n}: {int(t1):>10d} {chk246:14s}| {int(t2):>12d} {chk568:14s}| 3-col K_n: {ev(Tp,3)} | graphs {int(g2)} {chk088}")

print("\n=== integer polynomial P_n(x) = n! * T_n(x) = sum_{odd-cycle sigma} x^{#edge-orbits} ===")
for n in range(1,9):
    Tp=T_poly(n)
    coeffs={e:int(c*math.factorial(n)) for e,c in sorted(Tp.items())}
    print(f"   P_{n}(x): "+" + ".join(f"{coeffs[e]}x^{e}" for e in sorted(coeffs)) )
