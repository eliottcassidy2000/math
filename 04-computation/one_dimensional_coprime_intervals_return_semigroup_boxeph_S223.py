#!/usr/bin/env python3
"""one_dimensional_coprime_intervals_return_semigroup_boxeph_S223.py -- boxeph-2026-07-21-S223

CORRECTED BY MISTAKE-234: this script verifies support reachability and the
positive-coefficient case only. It does not prove mixed-sign noncancellation
or replace DvdK.

For f = sum_{k in S} c_k z^k (support S, 0 not in S = pure two-sided / charge != 0):
  * Newton polytope = the INTERVAL [min S, max S]; two-sided <=> 0 in its interior.
  * period d = gcd of the support gaps; coprime = aperiodic (d after reduction = 1).
  * the RETURN set R = {m : 0 is an m-fold sum of S} is a NUMERICAL SEMIGROUP (closed under +).
  * for POSITIVE coefficients CT(f^m) != 0  <=>  m in R. For mixed signs,
    reachability alone says nothing about cancellation (MISTAKE-234).
  Two poles: the bare coprime PAIR {-q,p} is PERIODIC (R = multiples of p+q = THM-1840, Frobenius = inf);
  the FILLED coprime interval (endpoints + interior, gcd gaps = 1) is COFINITE (R = all m >= Frobenius#).
"""
from math import gcd
from functools import reduce

print("CORRECTED BY MISTAKE-234: return semigroups prove reachability only; mixed-sign noncancellation and the DvdK bypass remain open.")

def sep(t): print("\n"+"="*72+"\n"+t+"\n"+"="*72)
def period(S):
    ds=[abs(a-b) for a in S for b in S if a!=b]
    return reduce(gcd, ds) if ds else 0
def two_sided(S): return min(S)<0<max(S)
def reachable(S, M):   # {m<=M : 0 is an m-fold sum of S (m terms from S summing to 0)}
    cur={0} if 0 in S else set()   # m=0 gives sum 0 trivially; we track m>=1
    reach=set()
    sums={0}  # achievable sums with 0 terms
    have=[]
    for m in range(1,M+1):
        sums={s+k for s in sums for k in S}
        if 0 in sums: reach.add(m)
    return reach
def frobenius(reach, M):
    # largest m<=M NOT in reach, assuming cofinite; return None if seems periodic (gaps persist to M)
    misses=[m for m in range(1,M+1) if m not in reach]
    if not misses: return -1
    # cofinite if the tail is all-present
    tail_ok=all(m in reach for m in range(M-20,M+1))
    return (max(misses) if tail_ok else None)

# ==========================================================================
sep("A  the coprime PAIR {-q,p} is PERIODIC: R = multiples of p+q (THM-1840); Frobenius = infinity")
for (q,p) in [(2,3),(1,2),(3,4),(2,5)]:
    S=[-q,p]; R=reachable(S,60); d=period(S)
    mult=sorted(R)[:6]
    print(f"  S={{-{q},{p}}} (coprime {gcd(p,q)==1}): period d={d}=p+q={p+q}; R (returns) start {mult} = multiples of {p+q}? "
          f"{all(m%(p+q)==0 for m in R)} ; Frobenius=inf (periodic)")
print("  => the single-character/bare-pair is fully PERIODIC: the coprime return m0=p+q, and ONLY its multiples.")

# ==========================================================================
sep("B  the FILLED coprime interval (endpoints + interior, gcd gaps = 1) is COFINITE: R = all m >= Frob#")
tests=[[-1,1,2],[-2,3,1],[-2,3,2],[-3,2,-1,1],[-2,3,-1,1,2]]
for S in tests:
    S=sorted(set(S)); R=reachable(S,80); d=period(S); F=frobenius(R,80)
    smallR=sorted(R)[:6]
    print(f"  S={S}: two-sided={two_sided(S)}, period d={d}; R starts {smallR}; Frobenius# (last non-return) = {F}"
          f"  => R = all m >= {F+1 if isinstance(F,int) and F>=0 else '?'}")
print("  => adding INTERIOR exponents (gcd of gaps -> 1) turns the periodic pair into a COFINITE return set:")
print("     all large m are SUPPORT returns, with a finite FROBENIUS NUMBER (the last non-return).")
print("     This is an exact reachability characterization, not coefficient-uniform noncancellation.")

# ==========================================================================
sep("C  POSITIVE coefficients: CT(f^m) != 0  <=>  m in R (reachability). Verified.")
def pmul(a,b):
    c={}
    for e1,c1 in a.items():
        for e2,c2 in b.items(): c[e1+e2]=c.get(e1+e2,0)+c1*c2
    return c
def ppow(a,m):
    r={0:1}
    for _ in range(m): r=pmul(r,a)
    return r
for S in [[-2,3],[-1,1,2],[-2,3,1]]:
    f={k:1 for k in S}  # positive (all coeffs 1)
    R=reachable(S,20)
    ok=all((CT_m!=0)==(m in R) for m in range(1,21) for CT_m in [ppow(f,m).get(0,0)])
    print(f"  S={S} (coeffs=1): CT(f^m) != 0  <=>  m in R  for m=1..20 ? {ok}"
          f"  ; CT sample {[ppow(f,m).get(0,0) for m in range(1,9)]}")
print("  => for positive coefficients the return set R IS the nonvanishing set: the coprime-interval semigroup.")

# ==========================================================================
sep("D  the corrected 1D support picture + the three-distance / LRC analogy")
print("""  ONE-DIMENSIONAL COPRIME INTERVALS -- exact SUPPORT reachability:
    interval [min S, max S]  +  period d=gcd(gaps)  =>  the RETURN SEMIGROUP R = {m : 0 in m-fold sum}.
    * bare coprime PAIR {-q,p}: PERIODIC, R = (p+q)Z_{>0}  (THM-1840 single-character seed).
    * APERIODIC two-sided support: COFINITE R, hence all m beyond a support conductor are reachable.
    * positive coeffs: CT(f^m)!=0 <=> m in R.
    * mixed/complex coeffs: a separate noncancellation theorem is required; infinitely many reachable m may vanish.
  Therefore the conductor is an effective support sidecar, not a DvdK replacement or a coefficient-uniform m0.

  THREE-DISTANCE / LRC analogy: coprimality also governs danger-arc orbit gaps. The gaps of
  {k t mod 1} take exactly THREE values (three-distance/Steinhaus), determined by the CONTINUED FRACTION of t
  = orbit arithmetic. This is a useful comparison, not a map from constant-term noncancellation to LRC.""")
def CT(a): return a.get(0,0)
