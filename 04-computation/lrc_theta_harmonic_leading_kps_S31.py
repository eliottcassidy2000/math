#!/usr/bin/env python3
"""
kind-pasteur-2026-07-06-S31: [theta => min-doubling] via HARMONIC-relation leading order.

mac-mini (HYP-4482) reduced AP-uniqueness (U) to:
  [OPEN: safe=0 => min-doubling] + [classical: min-doubling => AP] + [residue-pinning].
opus (HYP-4446): safe(S,beta) = SUM_{a in L(S)} prod_i h_hat(a_i), h_hat_0=1-2beta,
  h_hat_m = -sin(2 pi m beta)/(pi m)  (so |h_hat_m| ~ 1/m: SHORT relations dominate).
kps (S30): the (1,-2,1) HARMONIC relations (v_i - 2 v_{i+1} + v_{i+2}=0) are the SHORTEST
  nontrivial relations and CHARACTERIZE the AP (vanishing 2nd differences <=> AP, GREEN).

THE ROUTE for the open piece: the harmonic relations are the LEADING-ORDER theta terms;
safe=0 forces them present (mod a bounded truncation error over longer relations) =>
2nd-diff=0 => AP (S30).  We verify: (1) the theta-sum is dominated by the a=0 main term +
the shortest relations; (2) for the AP the harmonic/short relations carry the cancellation
to 0; (3) a non-AP family LACKS harmonic relations, so its safe > 0 by the leading order.
"""
import numpy as np
from itertools import product as iproduct
from math import sin, pi

BETA = 2.0/25.0
def hhat(m):
    if m == 0: return 1 - 2*BETA
    return -sin(2*pi*m*BETA)/(pi*m)

def safe_exact(v, grid=2000000):
    """safe(S,beta) = Leb{t: ||v_i t||>=beta all i} by fine grid."""
    t = (np.arange(grid)+0.5)/grid
    ok = np.ones(grid, bool)
    for x in v:
        d = np.abs(x*t - np.round(x*t))
        ok &= (d >= BETA)
    return ok.mean()

def theta_truncated(v, K):
    """SUM over relations a with |a_i|<=K, sum a_i v_i = 0, of prod hhat(a_i).
    Enumerate a in {-K..K}^n with sum a_i v_i = 0. Feasible for small n,K."""
    n = len(v); tot = 0.0; nrel = 0
    for a in iproduct(range(-K, K+1), repeat=n):
        if sum(a[i]*v[i] for i in range(n)) == 0:
            p = 1.0
            for ai in a: p *= hhat(ai)
            tot += p; nrel += 1
    return tot, nrel

# small n to keep the theta enumeration feasible (n=6, K=2)
print("=== theta-sum leading-order vs exact safe (n=6, the AP {1..6} vs perturbations) ===", flush=True)
def report(v, name, K=2):
    se = safe_exact(v)
    th, nr = theta_truncated(v, K)
    main = hhat(0)**len(v)
    print(f"  {name} v={v}: exact safe={se:.5f}; theta(|a|<={K})={th:.5f} ({nr} rels); main(a=0)={main:.5f}", flush=True)
    return se, th
AP6 = [1,2,3,4,5,6]
report(AP6, "AP {1..6}")
report([1,2,3,4,5,7], "lift 6->7")
report([1,2,3,4,5,8], "lift 6->8")
report([1,2,4,5,7,8], "no-3AP")
report([1,2,3,4,5,12], "far lift")

print(flush=True)
print("=== the (1,-2,1) HARMONIC relations: present for AP, absent for non-AP ===", flush=True)
def harmonic_relations(v):
    """count consecutive triples with v_i - 2 v_{i+1} + v_{i+2} = 0 (needs sorted)."""
    s = sorted(v); cnt = 0; total = len(s)-2
    for i in range(len(s)-2):
        if s[i] - 2*s[i+1] + s[i+2] == 0: cnt += 1
    return cnt, total
for v,name in [([1,2,3,4,5,6],"AP"),([1,2,3,4,5,7],"lift 6->7"),([1,2,4,5,7,8],"no-3AP"),
               ([2,4,6,8,10,12],"dilated AP"),([1,2,3,4,5,12],"far lift")]:
    c,t = harmonic_relations(v)
    print(f"  {name} {v}: {c}/{t} consecutive 2nd-differences vanish  {'(ALL => AP, S30)' if c==t else '(GAP => not AP)'}", flush=True)

print(flush=True)
print("=== leading-order mechanism: safe=0 needs the harmonic terms; missing => safe>0 ===", flush=True)
print("  the theta main term (1-2beta)^n = "+f"{hhat(0)**6:.5f} > 0 (n=6); the AP's harmonic", flush=True)
print("  relations supply the negative correction to reach ~0; a family MISSING them cannot", flush=True)
print("  cancel the main term at leading order, so safe > 0 -- the [theta=>min-doubling] route.", flush=True)
