#!/usr/bin/env python3
"""
lrc14_bounded_diameter_census_boxeph_S20.py  (HYP-6150, boxeph-2026-07-12-S20)

THE BOUNDED-DIAMETER FINITE CHECK, executed to Vmax <= VMAX (default 30).

Exhaustive census of primitive divisor-complete 13-subsets of [1..VMAX], one pass,
bucketed by Vmax = max(v):
  (a) min-clear q per family (kps cont.42 convention: smallest q >= 15, 14 nmid q,
      admitting p with all v_i*p mod q in the band [ceil(q/14), q-ceil(q/14)]);
      per-stratum worst min-clear => the window growth law vs the blocking bound.
  (b) M-floor detector: clearing at q <= 24 forces M >= ceil(q/14)/q >= 1/12
      (band-edge, opus-S235 PROVED), so only the min-clear >= 25 tail can carry
      M < 1/12 -- run the capped exact evaluator (cap 1/12) on that tail only.
      Extends kps cont.41's floor census (1/12 at Vmax <= 22) by 8 strata.
  (c) hard tail (min-clear >= 25) family list + blocked-moduli profile.

Worker mode: process a slice of (first, second)-element chunks.
  python3 ... --worker K NWORKERS VMAX   (writes JSON to scratchpad)
  python3 ... --combine NWORKERS VMAX    (merges + reports)
Pure Python; integer arithmetic only in hot paths.
"""
import sys, json, os
from math import gcd, comb
from functools import reduce
from itertools import combinations

SCRATCH = "/tmp/claude-1000/-home-claude-ephrepos-math/4f3f31d5-9a2c-4d2f-850f-3b8f9f9ee557/scratchpad"

FULL = (1 << 13) - 1
def cover_mask(x):
    m = 0
    for d in range(2, 15):
        if x % d == 0: m |= 1 << (d - 2)
    return m

def make_qplan(qmax=60):
    return [q for q in range(15, qmax + 1) if q % 14]

def min_clear(v, qplan):
    for q in qplan:
        band = -(-q // 14); hi = q - band
        for p in range(1, q):
            ok = True
            for x in v:
                r = (x * p) % q
                if r < band or r > hi: ok = False; break
            if ok: return q
    return None

def M_capped_below(v, capn=1, capd=12):
    """Return exact (m, q) if M(v) < capn/capd, else None. Small q first."""
    qs = list(range(2, 61))
    big = set()
    n = len(v)
    for i in range(n):
        for j in range(i, n):
            big.add(v[i] + v[j])
            if v[i] != v[j]: big.add(abs(v[i] - v[j]))
    qs += sorted(q for q in big if q > 60)
    mb, qb = 0, 1
    for q in qs:
        for p in range(1, q // 2 + 1):
            m = q
            for x in v:
                r = (x * p) % q
                if r > q - r: r = q - r
                if r < m:
                    m = r
                    if m * qb <= mb * q: break
            if m * qb > mb * q:
                mb, qb = m, q
                if mb * capd >= capn * qb:
                    return None
    return (mb, qb)

def chunks_for(vmax):
    """(a,b) first-two-element chunks with sizes, largest first."""
    out = []
    for a in range(1, vmax - 11):
        for b in range(a + 1, vmax - 10):
            out.append((comb(vmax - b, 11), a, b))
    out.sort(reverse=True)
    return out

def worker(k, nworkers, vmax):
    masks = [0] * (vmax + 1)
    for x in range(1, vmax + 1): masks[x] = cover_mask(x)
    qplan = make_qplan()
    mych = [c for i, c in enumerate(chunks_for(vmax)) if i % nworkers == k]
    # stats per Vmax stratum: counts + worst min-clear + min-clear histogram
    stats = {V: {"n": 0, "hist": {}, "worst": 0, "worstfam": None} for V in range(14, vmax + 1)}
    hard = []   # (family, minclear or None)
    below = []  # families with M < 1/12: (family, m, q)
    rng = range
    for size, a, b in mych:
        base_mask = masks[a] | masks[b]
        rest = rng(b + 1, vmax + 1)
        for tail in combinations(rest, 11):
            m = base_mask
            for x in tail: m |= masks[x]
            if m != FULL: continue
            v = (a, b) + tail
            g = 0
            for x in v:
                g = gcd(g, x)
                if g == 1: break
            if g != 1: continue
            V = v[-1]
            q = min_clear(v, qplan)
            st = stats[V]
            st["n"] += 1
            key = str(q)
            st["hist"][key] = st["hist"].get(key, 0) + 1
            qq = q if q is not None else 999
            if qq > st["worst"]:
                st["worst"] = qq; st["worstfam"] = list(v)
            if q is None or q >= 25:
                hard.append((list(v), q))
                r = M_capped_below(list(v))
                if r is not None:
                    below.append((list(v), r[0], r[1]))
    out = {"stats": {str(V): s for V, s in stats.items() if s["n"]},
           "hard": hard, "below": below}
    with open(f"{SCRATCH}/census_w{k}_{vmax}.json", "w") as f:
        json.dump(out, f)
    print(f"worker {k}/{nworkers} done: {sum(s['n'] for s in stats.values())} DC-primitive families, "
          f"{len(hard)} hard, {len(below)} below-1/12")

def combine(nworkers, vmax):
    agg = {}; hard = []; below = []
    for k in range(nworkers):
        with open(f"{SCRATCH}/census_w{k}_{vmax}.json") as f:
            d = json.load(f)
        for V, s in d["stats"].items():
            V = int(V)
            if V not in agg: agg[V] = {"n": 0, "hist": {}, "worst": 0, "worstfam": None}
            agg[V]["n"] += s["n"]
            for q, c in s["hist"].items():
                agg[V]["hist"][q] = agg[V]["hist"].get(q, 0) + c
            if s["worst"] > agg[V]["worst"]:
                agg[V]["worst"] = s["worst"]; agg[V]["worstfam"] = s["worstfam"]
        hard += d["hard"]; below += d["below"]
    print(f"BOUNDED-DIAMETER FINITE CHECK, Vmax <= {vmax} (HYP-6150)")
    print(f"{'Vmax':>4s} {'#prim-DC':>10s} {'cum#':>10s} {'worst min-clear (this stratum)':>32s} {'cum-worst':>9s}")
    cum = 0; cworst = 0
    for V in sorted(agg):
        s = agg[V]; cum += s["n"]; cworst = max(cworst, s["worst"])
        wf = s["worstfam"] if s["worst"] >= 25 else ""
        print(f"{V:4d} {s['n']:10d} {cum:10d} {s['worst']:32d} {cworst:9d}  {wf}")
    print(f"\nhard tail (min-clear >= 25 or none <= 60): {len(hard)} families")
    nn = [h for h in hard if h[1] is None]
    print(f"  min-clear NONE <= 60: {len(nn)}")
    for f_, q_ in nn[:10]: print(f"    {f_}")
    print(f"\nfamilies with M < 1/12: {len(below)}")
    for f_, m_, q_ in sorted(below, key=lambda t: t[1]/t[2])[:15]:
        print(f"    M = {m_}/{q_} = {m_/q_:.6f}  at {f_}")
    if not below:
        print("    NONE -- the DC M-floor 1/12 (kps cont.41, {1,2,3,4,10..18}) STANDS through Vmax <=", vmax)
    print("\nband-edge consequence: every family clearing at q <= 24 has M >= 1/12; the floor census")
    print("above is therefore complete (only the hard tail could dip, and it was checked exactly).")

if __name__ == '__main__':
    if sys.argv[1] == '--worker':
        worker(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
    elif sys.argv[1] == '--combine':
        combine(int(sys.argv[2]), int(sys.argv[3]))
