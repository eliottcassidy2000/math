#!/usr/bin/env python3
r"""
lrc_adaptive_hunter_kps_S73.py   (kind-pasteur-2026-07-07-S73, HYP-5147 part b)

THE FAMILY-ADAPTIVE HUNTER FLOOR at k=8 (THM-638 C3 extended from '>= 6/49' to its exact
family-dependent value), and its honest ceiling against the MISTAKE-123 bar T_8 = 0.6750.

Floor: mu_{1/7}(E) >= P(W_top) >= 1 - 7*theta + maxST sum m(d_a, d_b),
  d_a = e* - e_a (top endpoint; also bottom endpoint, take max), theta = 1/7,
  m(d_a, d_b) = theta^2 + G+(r_a, r_b)/(49 * q_a * q_b)   [THM-638 same-sign law]
  with (q_a, q_b) = (d_a/g, d_b/g), g = gcd, r_i = q_i mod 7, G+ = min(r)(7 - max(r)).

At k=8 the Bonferroni base 1 - 7*theta = 0, so the floor IS the max-spanning-tree mass.
Degree bookkeeping: each pair mass involves elements {e*, e_a, e_b} = DEGREE 3 -- this is
the leading (degree-3) sector of monad's HYP-5097 degree gap; the G-corrections are the
weight-3 balanced relations m1*(e*-e_a) + m2*(e*-e_b) = 0 of the kps-S59 deficit frame,
anchored at the endpoint.

QUESTIONS: (1) exact adaptive floors across the family zoo (structure meter);
(2) what is the MAX of the adaptive Hunter floor over ALL 8-families? (structured search)
(3) both-endpoint union: W_top and W_bot are disjoint-ish events -- can we add?
    (W_top ∩ W_bot both nonempty: need joint control -- report empirical overlap.)
"""
import math, random
from fractions import Fraction as F

TH = F(1, 7)

def Gplus(r1, r2):
    return min(r1, r2) * (7 - max(r1, r2))

def pair_mass(d1, d2):
    g = math.gcd(d1, d2)
    q1, q2 = d1 // g, d2 // g
    r1, r2 = q1 % 7, q2 % 7
    return TH * TH + F(Gplus(r1, r2), 49 * q1 * q2)

def max_spanning_tree_mass(ds):
    """Prim max-ST on complete graph over difference events with THM-638 weights."""
    n = len(ds)
    in_tree = [False] * n
    in_tree[0] = True
    best = [pair_mass(ds[0], ds[j]) for j in range(n)]
    total = F(0)
    for _ in range(n - 1):
        # pick max-weight edge to a new vertex
        j_best, w_best = -1, F(-1)
        for j in range(n):
            if not in_tree[j] and best[j] > w_best:
                j_best, w_best = j, best[j]
        in_tree[j_best] = True
        total += w_best
        for j in range(n):
            if not in_tree[j]:
                w = pair_mass(ds[j_best], ds[j])
                if w > best[j]:
                    best[j] = w
    return total

def hunter_floor(E):
    """max over top/bottom endpoint of the k=8 Hunter floor (= max-ST mass; base 0)."""
    E = sorted(E)
    d_top = [E[-1] - e for e in E[:-1]]
    d_bot = [e - E[0] for e in E[1:]]
    return max(max_spanning_tree_mass(d_top), max_spanning_tree_mass(d_bot))

def mu_numeric(E, res=20000, thr=1.0/7.0):
    c = 0
    for r in range(res):
        x = (r + 0.5) / res
        ph = sorted((e * x) % 1.0 for e in E)
        mg = ph[0] + 1 - ph[-1]
        for i in range(len(ph) - 1):
            mg = max(mg, ph[i + 1] - ph[i])
        if mg > thr:
            c += 1
    return c / res

T8_HONEST = F(1702763, 2522520)

print("=" * 96)
print(f"FAMILY-ADAPTIVE HUNTER FLOOR, k=8.  Bar: honest T_8 = {float(T8_HONEST):.4f} (MISTAKE-123)")
print("=" * 96)
rng = random.Random(73)
zoo = {
    "AP_8 {1..8}": list(range(1, 9)),
    "spread AP d=2": [1 + 2 * j for j in range(8)],
    "spread AP d=5": [1 + 5 * j for j in range(8)],
    "GW-ish {1..7,9}": [1, 2, 3, 4, 5, 6, 7, 9],
    "prim-sat-ish 2*{1..7}+{7}": [2, 4, 6, 7, 8, 10, 12, 14],
    "nearAP {1..7,15}": [1, 2, 3, 4, 5, 6, 7, 15],
    "geometric 2^j": [2 ** j for j in range(8)],
    "primes": [2, 3, 5, 7, 11, 13, 17, 19],
    "random small": sorted(rng.sample(range(1, 60), 8)),
    "random large": sorted(rng.sample(range(1, 10 ** 6), 8)),
    "multiples of 7 (apex-invisible)": [7 * j for j in range(1, 9)],
}
print(f"  {'family':>34} {'HunterFloor':>12} {'float':>7} {'mu (numeric)':>12} {'floor/T8':>9}")
for nm, E in zoo.items():
    hf = hunter_floor(E)
    mu = mu_numeric(E)
    print(f"  {nm:>34} {str(hf):>12} {float(hf):7.4f} {mu:12.4f} {float(hf/T8_HONEST):9.3f}")

print()
print("  Reference: bare floor (all G+ = 0) = 6*theta^2 =", F(6, 49), "=", float(F(6, 49)))
print()
print("(2) MAXIMIZE the adaptive floor over 8-families (hill climb over difference sets):")
def climb(seed, iters=4000):
    rng2 = random.Random(seed)
    E = sorted(rng2.sample(range(1, 40), 8))
    best = hunter_floor(E)
    for _ in range(iters):
        E2 = sorted(set(E))
        j = rng2.randrange(8)
        delta = rng2.choice([-3, -2, -1, 1, 2, 3])
        cand = sorted(set(E2) - {E2[j]} | {max(1, E2[j] + delta)})
        if len(cand) < 8:
            continue
        h = hunter_floor(cand)
        if h > best:
            best, E = h, cand
    return best, E

overall_best, overall_E = F(0), None
for s in range(8):
    b, E = climb(s)
    if b > overall_best:
        overall_best, overall_E = b, E
print(f"    best found: {overall_best} = {float(overall_best):.4f} at E = {overall_E}")
print(f"    vs T_8 honest = {float(T8_HONEST):.4f}  =>  max adaptive Hunter reaches "
      f"{float(overall_best / T8_HONEST) * 100:.0f}% of the bar")
print()
print("(3) both-endpoint UNION (empirical overlap of W_top, W_bot at the AP):")
def W_events(E, res=20000, thr=1.0/7.0):
    E = sorted(E)
    top_ct = bot_ct = both_ct = 0
    for r in range(res):
        x = (r + 0.5) / res
        ph = [(e * x) % 1.0 for e in E]
        pt, pb = ph[-1], ph[0]
        others_t = [((p - pt) % 1.0) for p in ph[:-1]]
        others_b = [((pb - p) % 1.0) for p in ph[1:]]
        wt = all(o > thr for o in others_t)
        wb = all(o > thr for o in others_b)
        top_ct += wt; bot_ct += wb; both_ct += (wt and wb)
    return top_ct / res, bot_ct / res, both_ct / res
for nm in ("AP_8 {1..8}", "random large"):
    t, b, bo = W_events(zoo[nm])
    print(f"    {nm:>34}: P(Wtop)={t:.4f} P(Wbot)={b:.4f} P(both)={bo:.4f} "
          f"union={t + b - bo:.4f}")
print()
print("(4) THE EXACT DEGREE-3 ENDPOINT CEILING (little theorem, proof = finite enumeration):")
print("    per-edge mass theta^2 + G+(r1,r2)/(49 q1 q2) > 1/14  <=>  G+/(q1q2) > 5/2;")
print("    G+ <= 12 forces q1*q2 < 24/5; coprime distinct pairs with product <= 4:")
for (q1, q2) in [(1, 2), (1, 3), (1, 4), (2, 3)]:
    m = pair_mass(q1, q2)
    print(f"      ({q1},{q2}): mass = {m} = {float(m):.5f}  {'= 1/14 MAX' if m == F(1,14) else ''}")
print("    => max per-edge mass = 1/14 EXACTLY, uniquely at reduced (1,2) (doubling pairs).")
print("    Tree on k-1=7 events has 6 edges => endpoint-Hunter floor <= 6/14 = 3/7 = 0.42857")
print("    for EVERY 8-family; attained by doubling-difference chains:")
dbl = [65 - d for d in [64, 32, 16, 8, 4, 2, 1]] + [65]   # E with differences {1,2,4,...,64}
dbl = sorted(dbl)
hf = hunter_floor(dbl)
mu = mu_numeric(dbl)
print(f"      E = {dbl} (differences from top = 2^j): floor = {hf} = {float(hf):.4f}, mu = {mu:.4f}")
print(f"    3/7 = {float(F(3,7)):.4f} < honest T_8 = {float(T8_HONEST):.4f}  (MISTAKE-123 bar)")
print()
print("(5) k-sweep of the endpoint-Hunter ceilings vs honest bars:")
HONEST = {8: F(1702763, 2522520), 9: F(35456, 63063), 10: F(114041, 252252),
          11: F(83549, 252252), 12: F(50285, 252252), 13: F(14249, 252252)}
print(f"    {'k':>3} {'uniform floor (G=0)':>22} {'adaptive ceiling':>18} {'honest T_k':>12} {'reach':>6}")
for k in range(8, 14):
    base = 1 - (k - 1) * TH
    unif = max(F(0), base + (k - 2) * TH * TH)
    adapt = max(F(0), base + (k - 2) * F(1, 14))
    Tk = HONEST[k]
    print(f"    {k:>3} {str(unif):>22} {str(adapt):>18} {float(Tk):>12.4f} "
          f"{float(adapt / Tk):>6.2f}")
print()
print("READING: the adaptive Hunter floor is a DEGREE-3 harvester (shared-endpoint pair")
print("masses = weight-3 balanced relations anchored at the endpoint -- the leading sector")
print("of monad's HYP-5097 degree gap).  Its EXACT ceiling is 3/7 at k=8 (doubling chains),")
print("uniform value 6/49; both far below the honest bar 0.675.  Palindromic families have")
print("W_top = W_bot EXACTLY (equal difference sets: min_j frac(-mx) over the same m-set),")
print("so the 2-endpoint union gains NOTHING for the AP (the S62 palindrome symmetry in the")
print("Hunter frame); spread families gain (union 0.55) but need no help.  Interior anchors")
print("have mixed-sign differences whose pair masses can VANISH (THM-638 C2) => no degree-3")
print("floor at all.  THE AP'S TAIL MASS IS CARRIED BY INTERIOR/ROTATING ANCHORS INVISIBLE")
print("AT DEGREE 3: the responsible gap sits above the Farey-cell element q, which rotates")
print("with the cell -- no fixed element carries it.  k=8 needs degree >= 4 or k-body tools.")
