"""
tight_locus_census_kps.py  (kind-pasteur, 2026-06-16)

CRUX of LRC(14): census the PRIMITIVE TIGHT LOCUS
  { S : L(S)=0, gcd(S)=1, 13 distinct positive ints }
to decide finiteness / bounded-lcm.

L(S) = meas{ tau in [0,1) : ||v tau|| > 1/14  for all v in S }.
Danger set D_v = union_k [(14k-1)/(14v),(14k+1)/(14v)] ; total meas 1/7 each.
L = 1 - meas(union_v D_v).  L>0 = loose, L=0 = tight.

We use EXACT rational arc-sweep (verbatim from shared context).
"""
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
from itertools import combinations
import sys, json, time

# ---------- exact L ----------
def danger_arcs(v):
    w = Fr(1, 14*v); A = []
    for k in range(v+1):
        c = Fr(k, v); lo = c - w; hi = c + w
        if lo < 0:
            A += [(Fr(0), hi), (1+lo, Fr(1))]
        elif hi > 1:
            A += [(lo, Fr(1)), (Fr(0), hi-1)]
        else:
            A.append((lo, hi))
    return A

def L_exact(S):
    A = []
    for v in S:
        A.extend(danger_arcs(v))
    A = sorted((a, b) for a, b in A if b > a)
    tot = Fr(0); cl = ch = None
    for a, b in A:
        if ch is None:
            cl, ch = a, b
        elif a <= ch:
            ch = max(ch, b)
        else:
            tot += ch - cl; cl, ch = a, b
    if ch is not None:
        tot += ch - cl
    return 1 - tot

def gcd_list(S):
    return reduce(gcd, S)

def is_tight(S):
    return L_exact(tuple(S)) == 0

def lcm2(a, b):
    return a*b // gcd(a, b)

def lcm_list(S):
    return reduce(lcm2, S)

# ---------- max-min computation (exact) ----------
# max-min = max_tau min_v ||v tau||.  This is the largest gap left UNCOVERED by
# the OPEN danger arcs D_v^o = {||v tau|| < 1/14}.  Tight (L=0) means the CLOSED
# arcs cover [0,1].  We compute the supremum of min_v ||v tau||.
# The function tau -> min_v ||v tau|| is piecewise linear; its local maxima occur
# at points equidistant between two "danger centers" k/v.  Equivalent: the
# uncovered region by open arcs [(14k-1)/(14v),(14k+1)/(14v)] open.  We sweep the
# complement of OPEN arcs and find the max of min_v||v tau|| over leftover points
# AND arc-boundary crossings.
# Simpler robust approach: max-min is achieved at a point that is a boundary of
# some open danger arc OR a crossing of two boundary functions. Candidate taus:
# all arc endpoints (closed-arc style) give min_v||v tau|| values; plus midpoints.
# We instead compute it via the OPEN-arc complement sweep: the value at any
# leftover open interval's interior is >1/14 only if the closed arcs don't cover.
# For a TIGHT config the closed arcs cover everything, so max-min <= 1/14.
# We want to know if max-min == 1/14 (touch) or < 1/14 (counterexample needs L=0
# strictly via closed-arc cover with margin -> max-min<1/14).
#
# We compute max-min exactly as follows: the candidate optimal taus are the
# rationals k/v +- 1/(14v) (arc endpoints).  At an endpoint tau0 of v's arc,
# ||v tau0|| = 1/14 exactly.  min over all u of ||u tau0|| <= 1/14.  The global
# max-min for a covering config is sup over the closed-arc boundary structure.
# Robust exact method: collect all "transition" taus = endpoints of all OPEN
# arcs, sort, and on each maximal sub-interval where the active nearest-center
# assignment per v is fixed, min_v||v tau|| is piecewise linear with breakpoints
# only at arc endpoints; its sup on the closed interval is at an endpoint.
# So max-min = max over all arc endpoints tau of min_v ||v tau||... BUT the sup
# could be approached in the open interior between two endpoints where one v's
# distance is increasing and another's decreasing -> crossing.  Include pairwise
# crossings  ||a tau|| == ||b tau||  i.e.  a tau - p = +-(b tau - q).
# Implement: gather candidate taus from endpoints + pairwise crossings, eval exact.

def frac_norm(x):
    # ||x|| distance to nearest integer, exact Fraction in [0,1/2]
    f = x - (x.__floor__())
    if f > Fr(1, 2):
        f = 1 - f
    return f

def max_min_exact(S):
    S = list(S)
    cands = set()
    cands.add(Fr(0))
    # arc endpoints for each v: k/v +- 1/(14 v)
    for v in S:
        w = Fr(1, 14*v)
        for k in range(v+1):
            c = Fr(k, v)
            for t in (c - w, c + w, c):
                if 0 <= t <= 1:
                    cands.add(t)
    # pairwise crossings of ||a tau|| and ||b tau||
    # ||a tau|| = |a tau - round|, linear pieces a tau - p (p integer near).
    # crossing a tau - p = b tau - q  -> tau=(p-q)/(a-b)
    # crossing a tau - p = -(b tau - q) -> tau=(p+q)/(a+b)
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            a, b = S[i], S[j]
            for p in range(a+1):
                for q in range(b+1):
                    if a != b:
                        t = Fr(p-q, a-b)
                        if 0 <= t <= 1: cands.add(t)
                    t = Fr(p+q, a+b)
                    if 0 <= t <= 1: cands.add(t)
    best = Fr(0); arg = Fr(0)
    for t in cands:
        m = min(frac_norm(v*t) for v in S)
        if m > best:
            best = m; arg = t
    return best, arg

# ---------- enumeration ----------
def replacements_k(base, k, hi):
    """k-element replacements of `base` (a sorted list) with new entries up to hi.
    Yields candidate 13-sets (primitive check done by caller)."""
    base = list(base)
    n = len(base)
    pool = [x for x in range(1, hi+1)]
    for drop_idx in combinations(range(n), k):
        kept = [base[i] for i in range(n) if i not in drop_idx]
        kept_set = set(kept)
        # choose k new entries not in kept
        choices = [x for x in pool if x not in kept_set]
        for newvals in combinations(choices, k):
            S = sorted(kept + list(newvals))
            if len(set(S)) == 13:
                yield tuple(S)

def primitive(S):
    return gcd_list(S) == 1

def main():
    t0 = time.time()
    AP = list(range(1, 14))  # {1..13}
    tight = []   # list of dicts
    seen = set()

    def consider(S, source):
        S = tuple(sorted(S))
        if S in seen:
            return
        if len(set(S)) != 13:
            return
        if any(x <= 0 for x in S):
            return
        seen.add(S)
        if not primitive(S):
            return
        if is_tight(S):
            tight.append({"S": S, "lcm": lcm_list(S), "max": max(S),
                          "source": source})

    # sanity: the AP itself
    consider(AP, "AP {1..13}")
    print("AP {1..13} tight?", is_tight(AP), flush=True)
    # known sporadic
    spor = list(range(1,12)) + [13, 24]
    consider(spor, "known sporadic {1..11,13,24}")
    print("sporadic {1..11,13,24} tight?", is_tight(spor), flush=True)
    # claimed inf config
    inf_cfg = list(range(1,12)) + [13, 36]
    print("{1..11,13,36} L =", L_exact(tuple(inf_cfg)), flush=True)

    # (a) k=1 replacements, entries up to 120
    print("\n=== k=1 replacements, hi=200 ===", flush=True)
    cnt = 0
    for S in replacements_k(AP, 1, 200):
        consider(S, "k=1 repl hi200")
        cnt += 1
    print(f"  scanned {cnt} k=1 candidates; tight so far {len(tight)}", flush=True)

    # (b) k=2 replacements, entries up to 90
    print("\n=== k=2 replacements, hi=90 ===", flush=True)
    cnt = 0
    for S in replacements_k(AP, 2, 90):
        consider(S, "k=2 repl hi90")
        cnt += 1
    print(f"  scanned {cnt} k=2 candidates; tight so far {len(tight)}", flush=True)

    # (c) k=3 replacements, entries up to 40 (combinatorial blowup control)
    print("\n=== k=3 replacements, hi=40 ===", flush=True)
    cnt = 0
    for S in replacements_k(AP, 3, 40):
        consider(S, "k=3 repl hi40")
        cnt += 1
    print(f"  scanned {cnt} k=3 candidates; tight so far {len(tight)}", flush=True)

    # (d) scaled cores  d*({1..13}\{j}) ∪ {w}, small d
    print("\n=== scaled-perturbed cores ===", flush=True)
    for d in range(1, 7):
        for j in range(1, 14):
            core = [d*x for x in range(1,14) if x != j]
            for w in range(1, d*40+1):
                S = sorted(core + [w])
                if len(set(S)) == 13:
                    consider(S, f"scaled d={d} drop {j}*d add {w}")
    print(f"  tight so far {len(tight)}", flush=True)

    # (e) structured: Sidon, geometric-ish, mixed - quick probes
    print("\n=== structured probes ===", flush=True)
    structured = [
        list(range(1,14)),                       # AP
        [2*i+1 for i in range(13)],              # odd AP
        [i*i for i in range(1,14)],              # squares
        [1,2,4,8,16,32,64,128,256,512,1024,2048,4096],  # powers of 2
        list(range(2,15)),                       # 2..14 (gcd!) -> not primitive
        list(range(1,14)),
    ]
    for S in structured:
        consider(S, "structured")
    print(f"  tight so far {len(tight)}", flush=True)

    # report
    print("\n\n===== TIGHT CONFIGS FOUND =====", flush=True)
    tight_sorted = sorted(tight, key=lambda d: (d["lcm"], d["max"]))
    for d in tight_sorted:
        print(f"  lcm={d['lcm']:>8} max={d['max']:>4}  S={list(d['S'])}  [{d['source']}]", flush=True)
    print(f"\nTotal primitive tight configs found: {len(tight)}", flush=True)
    if tight:
        print(f"Max lcm among tight: {max(d['lcm'] for d in tight)}", flush=True)
        print(f"Max entry among tight: {max(d['max'] for d in tight)}", flush=True)

    # max-min for each tight config (look for counterexample < 1/14)
    print("\n===== MAX-MIN (exact) per tight config =====", flush=True)
    threshold = Fr(1,14)
    counterexamples = []
    for d in tight_sorted:
        mm, arg = max_min_exact(d["S"])
        flag = ""
        if mm < threshold:
            flag = "  <<<< COUNTEREXAMPLE max-min < 1/14 !!!"
            counterexamples.append((d["S"], mm, arg))
        print(f"  S={list(d['S'])}  max-min={mm}={float(mm):.6f}  (1/14={float(threshold):.6f}) arg={arg}{flag}", flush=True)

    # save
    out = {
        "tight": [{"S": list(d["S"]), "lcm": d["lcm"], "max": d["max"], "source": d["source"]} for d in tight_sorted],
        "n_tight": len(tight),
        "max_lcm": max((d["lcm"] for d in tight), default=None),
        "max_entry": max((d["max"] for d in tight), default=None),
        "counterexamples": [{"S": list(s), "maxmin_num": mm.numerator, "maxmin_den": mm.denominator} for s, mm, arg in counterexamples],
        "elapsed_sec": time.time()-t0,
    }
    with open("05-knowledge/results/tight_locus_census_kps.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nElapsed {time.time()-t0:.1f}s. Saved JSON.", flush=True)

if __name__ == "__main__":
    main()
