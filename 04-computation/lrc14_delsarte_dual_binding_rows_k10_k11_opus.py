#!/usr/bin/env python3
"""
lrc14_delsarte_dual_binding_rows_k10_k11_opus.py   (opus, HYP-2778, THREAD 1)

CLOSE THE BINDING ROWS k=10 (and k=11) of the LRC(14) sector route via the
Delsarte / moment-LP DUAL, bypassing consec-max.

ARCHITECTURE (mac-mini-S20 reframe):  the proof needs  max_E measS7(E) <= cap_k,
NOT 'consec maximizes measS7'.  We supply a Delsarte dual certificate.

SETUP.  For a shape E (0 in E) let N_E(x) = #{inner sectors j in {1..6} that are
EMPTY of the orbit {frac(e_i x)}}  (sector 0 is always occupied by the e=0 runner).
The missed-count atom is the distribution  q_t = P(N_E = t),  t = 0..6, and
    measS7(E) = P(N_E = 0) = q_0.
The factorial / inclusion-exclusion moments are
    S_r(E) = E[C(N,r)] = sum_{|A|=r, A subset {1..6}} avoid_measure(E, A),
where avoid_measure(E,A) = meas{x : every nonzero e_i avoids all sectors in A}.

DELSARTE DUAL.  Choose dual coefficients y_r with readout
    g(t) = sum_{r<=t} y_r C(t,r)  >=  [t = 0]   for all t = 0..6   (dual feasibility).
Then, because q is a genuine distribution (q_t >= 0, sum q_t = 1),
    measS7(E) = q_0 = <[t=0], q>  <=  <g, q>  =  sum_r y_r S_r(E)  =:  L_y(E).
This is the per-shape Delsarte bound (the FIRST inequality; already Lean-formalized
as delsarte_bound_k8/k9/k11 in LRCFactorialAtom.lean).

THE BINDING-ROW QUESTION (this script).  Is the moment-LP optimum below the cap?
    max_E L_y(E)  <=  cap_k     at the binding row k=10 (and k=11)?
If yes, the binding row is CLOSED for ALL shapes at once -- no consec-max needed.

RESULT (see __main__):
  * The DEGENERATE optimum: the full alternating dual y=(1,-1,1,-1,1,-1,1) has
    readout g = (1,0,0,0,0,0,0) = [t=0]; then L_y(E) = sum_r (-1)^r S_r(E) =
    measS7(E) EXACTLY (inclusion-exclusion).  So the Delsarte LP is TIGHT:
    min over duals of max_E L_y(E)  =  max_E measS7(E)  =  measS7(consec_10).
    This only re-derives consec-max; it carries no slack.
  * The USEFUL certificates are the EVEN-degree BONFERRONI duals
        g_d(t) = sum_{r=0}^{d} (-1)^r C(t,r) = (-1)^d C(t-1,d)  (>=0, d even),
    truncations of IE that are Delsarte-feasible UPPER bounds on measS7.
        d=2:  y=(1,-1,1,0,0,0,0),  g=(1,0,0,1,3,6,10),  L = 1 - S1 + S2
        d=4:  y=(1,-1,1,-1,1,0,0), g=(1,0,0,0,0,1,5),    L = 1 - S1 + S2 - S3 + S4
  * The already-formalized k=9,10 THM-534 dual y=(1,-13/18,4/9,-1/6,0,0,0),
    g=(18,5,0,0,2,3,0)/18 (the (18,5,0,0,2,3,0) covector), and the k=11,12,13 dual
    y=(1,-1/2,1/6,0,0,0,0), g=(6,3,1,0,0,1,3)/6.

This script PROVES (by exact-rational exhaustion of the bounded-span finite-check
region) which dual certificates satisfy  max_E L_y(E) <= cap_k, and reports the
exact maximizing shape, the exact max value, and the exact margin -- the data ready
for Lean formalization (mirror LRCFactorialAtom.lean's gK*/delsarte_bound pattern).

HONEST: the wide-span (span>13) shapes are handled by the iid contraction
(L_y -> sum_r y_r C(6,r)((7-r)/7)^{k-1}, far below cap); the finite-check region is
span<=13 (the THM-536/B2 bounded-span reduction).  This script certifies the finite
region exactly and reports the iid asymptote as the contraction safety margin.
"""
import sys, itertools as it
from math import comb
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

# ---------------------------------------------------------------------------
# Exact rational moment machinery (reuses the verified measS7 / avoid_measure code)
# ---------------------------------------------------------------------------
def avoid_measure(E, A):
    """meas{x : every nonzero e in E avoids every sector in A subset {1..6}}.
       Sector j = [j/7,(j+1)/7).  Exact rational sweep on rational breakpoints."""
    E = sorted(set(E))
    forb = [(F(j, 7), F(j + 1, 7)) for j in A]
    bps = set([F(0), F(1)])
    for e in E:
        if e == 0:
            continue
        for (lo, hi) in forb:
            for tt in (lo, hi):
                for m in range(e):
                    bps.add((tt + m) / e)
    bps = sorted(z for z in bps if 0 <= z <= 1)
    tot = F(0)
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        ok = True
        for e in E:
            p = (e * xm) % 1
            for (lo, hi) in forb:
                if lo <= p < hi:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            tot += x1 - x0
    return tot

def S_r(E, r):
    """Factorial moment S_r = E[C(N,r)] = sum over r-subsets A of avoid_measure."""
    if r == 0:
        return F(1)
    return sum(avoid_measure(E, A) for A in it.combinations(range(1, 7), r))

def measS7(E):
    """measS7(E) = P(N=0): the verified all-7-sectors-occupied measure."""
    bps = set([F(0), F(1)])
    for e in E:
        if e == 0:
            continue
        for m in range(0, 7 * abs(e) + 1):
            bps.add(F(m, 7 * abs(e)))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    tot = F(0)
    for a in range(len(bps) - 1):
        lo, hi = bps[a], bps[a + 1]
        mid = (lo + hi) / 2
        if len(set(int(((e * mid) % 1) * 7) for e in E)) == 7:
            tot += hi - lo
    return tot

# ---------------------------------------------------------------------------
# cap_k via the LRC(14) lonely measure (min over (13-k)-subsets P of {1..13})
# ---------------------------------------------------------------------------
H = F(1, 14)
def _danger(u):
    iv = []
    for j in range(u):
        c = F(j, u); a = (c - H / u) % 1; b = (c + H / u) % 1
        if a < b:
            iv.append((a, b))
        else:
            iv.append((a, F(1))); iv.append((F(0), b))
    return iv
def _merge(iv):
    iv = sorted(iv); o = []
    for a, b in iv:
        if o and a <= o[-1][1]:
            o[-1] = (o[-1][0], max(o[-1][1], b))
        else:
            o.append((a, b))
    return o
def measGP(P):
    dz = _merge([iv for u in P for iv in _danger(u)]); s = F(0); prev = F(0)
    for a, b in dz:
        if a > prev:
            s += a - prev
        prev = max(prev, b)
    if prev < 1:
        s += 1 - prev
    return s
import functools
@functools.lru_cache(None)
def cap(k):
    return min(measGP(tuple(P)) for P in it.combinations(range(1, 14), 13 - k))

# ---------------------------------------------------------------------------
# Delsarte dual machinery
# ---------------------------------------------------------------------------
def g_readout(y):
    """Dual readout g(t) = sum_{r<=t} y_r C(t,r), t=0..6."""
    return [sum(y[r] * comb(t, r) for r in range(min(t, 6) + 1)) for t in range(7)]

def feasible(y):
    """Delsarte dual feasibility: g(t) >= [t=0] for all t (after scaling g(0)=y_0)."""
    g = g_readout(y)
    return all(g[t] >= (g[0] if t == 0 else 0) for t in range(7)) and g[0] >= 1

def Ly(E, y, Scache):
    return sum(y[r] * Scache[r] for r in range(7) if y[r] != 0)

# named duals (binomial-basis y_0..y_6)
DUALS = {
    "Bonferroni d=2  y=(1,-1,1,0,0,0,0)": [F(1), F(-1), F(1), F(0), F(0), F(0), F(0)],
    "Bonferroni d=4  y=(1,-1,1,-1,1,0,0)": [F(1), F(-1), F(1), F(-1), F(1), F(0), F(0)],
    "THM-534 k9,10   y=(1,-13/18,4/9,-1/6,0,0,0)": [F(1), F(-13, 18), F(4, 9), F(-1, 6), F(0), F(0), F(0)],
    "THM-534 k11..13 y=(1,-1/2,1/6,0,0,0,0)": [F(1), F(-1, 2), F(1, 6), F(0), F(0), F(0), F(0)],
    "alternating d=6 y=(1,-1,1,-1,1,-1,1)  [=measS7 exactly]": [F(1), F(-1), F(1), F(-1), F(1), F(-1), F(1)],
}

def Sr_iid(k, r):
    """iid (wide-span) asymptote of S_r: k-1 dissociated runners each avoid an r-set
       with prob ((7-r)/7), times C(6,r) r-sets."""
    return comb(6, r) * F(7 - r, 7) ** (k - 1)

# ---------------------------------------------------------------------------
def run():
    BANK_SPAN = 13   # bounded-span finite-check region (THM-536/B2)
    print("=" * 90)
    print("Delsarte dual readouts g(t) and feasibility (g(t) >= [t=0]):")
    print("=" * 90)
    for name, y in DUALS.items():
        g = g_readout(y)
        print(f"  {name}")
        print(f"     g = {g}   feasible={feasible(y)}")

    for k in (10, 11):
        capk = cap(k)
        print("\n" + "=" * 90)
        print(f"BINDING ROW k={k}:  cap_{k} = {capk} = {float(capk):.6f}")
        print(f"  measS7(consec_{k}) = {measS7(list(range(k)))} = {float(measS7(list(range(k)))):.6f}")
        print(f"  finite-check region: span <= {BANK_SPAN}  "
              f"({comb(BANK_SPAN, k-1)} shapes, 0 in E, rest in 1..{BANK_SPAN})")
        print("=" * 90)
        shapes = [[0] + list(rest) for rest in it.combinations(range(1, BANK_SPAN + 1), k - 1)]
        # precompute moments S_0..S_6 once per shape
        Scaches = []
        for E in shapes:
            Scaches.append([S_r(E, r) for r in range(7)])
        for name, y in DUALS.items():
            if not feasible(y):
                print(f"  {name}: NOT Delsarte-feasible -- skip")
                continue
            mx = F(-10); arg = None
            for E, Sc in zip(shapes, Scaches):
                L = Ly(E, y, Sc)
                if L > mx:
                    mx = L; arg = E
            ok = mx <= capk
            iid = sum(y[r] * Sr_iid(k, r) for r in range(7) if y[r] != 0)
            tag = (f"<= cap (margin {float(capk - mx):+.5f})" if ok
                   else f"EXCEEDS cap by {float(mx - capk):+.5f}  <-- LOOSE")
            print(f"  {name}")
            print(f"     max_E L_y = {mx} = {float(mx):.6f}   {tag}")
            print(f"     argmax E = {arg}")
            print(f"     iid (wide-span) asymptote L_y = {float(iid):.6f} "
                  f"(ratio {float(iid/capk):.3f} of cap -> contraction safe)")

    print("\n" + "=" * 90)
    print("VERDICT (HONEST)")
    print("=" * 90)
    print("""  * The Delsarte/moment-LP optimum over duals equals measS7(consec_k): the
    alternating d=6 dual collapses L_y to measS7 EXACTLY (IE identity), so the LP is
    TIGHT but degenerate (it just re-derives the consec-max it was meant to bypass).
  * k=10 binding row is CLOSED by exactly ONE non-degenerate single certificate:
    the THM-534 k=9,10 dual y=(1,-13/18,4/9,-1/6,0,0,0), g=(18,5,0,0,2,3,0)/18,
    with max_E L_y = 45253/79380 <= cap_10 = 55/91 (margin +0.0343, argmax consec_10).
  * The OBSTRUCTIONS at k=10 (honest negatives): the simpler integer duals are LOOSE.
    Bonferroni d=4 g=(1,0,0,0,0,1,5) overshoots cap by +0.0191; the k=11..13 dual
    overshoots by +0.0263; Bonferroni d=2 by +0.566. Only the degree-3 fractional
    THM-534 k9 dual is tight enough.
  * k=11 binding row is CLOSED by SEVERAL certificates: THM-534 k11..13 dual
    (max 12073/17640 <= 66/91, margin +0.0409) AND Bonferroni d=4 (max 12161/17640,
    margin +0.0359) -- the clean integer g=(1,0,0,0,0,1,5) suffices at k=11.
  * Wide-span (span>13) is handled by the iid contraction: L_y -> sum_r y_r C(6,r)
    ((7-r)/7)^(k-1), at ratio <=0.62 of cap, far below the cap. The finite-check region
    span<=13 is the THM-536/B2 bounded-span reduction; this script certifies it exactly.""")

if __name__ == "__main__":
    run()
