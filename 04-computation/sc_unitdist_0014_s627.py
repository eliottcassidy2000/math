#!/usr/bin/env python3
"""
S627 — the 0.014 coincidence: SC-tournament shape deficit = (claimed) unit-distance exponent,
and the n-2 (=n+2) recursion that powers both, in the partition-function frame (HYP-2245).

SC side (canon s578): non-SC(n) = 2^{C(n-1,2)-n+3} · (1 - deficit(n)); deficit halves each step,
and the correction 2^{...}-non-SC(n) ~ SC(n-2)  (the n-2 recursion).  0.014 = deficit at n=8.
Unit-distance side (A186705): the increments / empirical exponent; look for the 0.014 / halving.
"""
from math import comb, log

# canon sequences (s578), n=3..
nonSC = [1, 3, 14, 121, 1995, 64648, 4163979, 534849295, 137175056830, 70300582005021]
SC    = [1, 5, 50, 903, 30773, 2032504, 264271477, 68184627441, 35047197032002, 35958496436958947]
# index 0 -> n=3

print("SC-TOURNAMENT SHAPE DEFICIT  non-SC(n)/2^{C(n-1,2)-n+3}  (deficit = 1 - ratio)")
print("="*72)
for i, ns in enumerate(nonSC):
    n = i+3
    e = comb(n-1, 2) - n + 3
    lead = 2**e
    ratio = ns/lead
    deficit = 1 - ratio
    corr = lead - ns                      # the correction
    print(f" n={n:2d}: 2^{e:<2d}={lead:<13d} non-SC={ns:<13d} ratio={ratio:.5f} deficit={deficit:.4f} corr={corr}")

print("\nn-2 RECURSION:  correction(n) = 2^{C(n-1,2)-n+3} - non-SC(n)  ~  SC(n-2)")
for i in range(2, len(nonSC)):
    n = i+3
    e = comb(n-1, 2) - n + 3
    corr = 2**e - nonSC[i]
    sc2 = SC[i-2]                          # SC(n-2)
    print(f"  n={n:2d}: corr={corr}  SC(n-2)={sc2}  corr/SC(n-2)={corr/sc2:.6f}")

print("\nDEFICIT halving + the 0.014 point (n=8):")
defs = []
for i, ns in enumerate(nonSC):
    n=i+3; e=comb(n-1,2)-n+3; defs.append((n, 1-ns/2**e))
for n,d in defs: print(f"  n={n}: deficit={d:.4f}" + ("   <-- 0.014" if abs(d-0.014)<0.002 else ""))
print("  ratios deficit(n)/deficit(n-1):", [round(defs[i][1]/defs[i-1][1],3) for i in range(1,len(defs))])

print("\n" + "="*72)
print("UNIT DISTANCE A186705 u(n): empirical exponent and increments")
u = {2:1,3:3,4:5,5:7,6:9,7:12,8:14,9:18,10:20,11:23,12:27,13:30,14:33,15:37,16:41,17:43,18:46,19:50,20:54,21:57,22:60}
ns = sorted(u)
print(" n  u(n)  log u/log n  (4/3 - that)   local-exp   incr")
for k,n in enumerate(ns):
    a = log(u[n])/log(n)
    le = (log(u[n])-log(u[ns[k-1]]))/(log(n)-log(ns[k-1])) if k>0 else float('nan')
    incr = u[n]-u[ns[k-1]] if k>0 else 0
    print(f" {n:2d}  {u[n]:3d}   {a:.4f}      {4/3-a:+.4f}      {le:.4f}    {incr}")
# look for 0.014 as a gap-to-4/3 or as a deficit
print("\n where does ~0.014 appear on the unit-distance side?")
for k,n in enumerate(ns):
    g = 4/3 - log(u[n])/log(n)
    if abs(g-0.014)<0.004: print(f"   n={n}: 4/3 - log u/log n = {g:.4f}  ~ 0.014")
