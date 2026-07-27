#!/usr/bin/env python3
"""Round 7: beam search over irregular hybrids.

State = 13-set; moves = replace one element with a value from the candidate
pool (duty-stacked values in [61,300] + [15,60] fillers); objective =
lexicographic (#bad moduli in [14,60], total excess depth, then -v_max).
Seeds: the 59 survivors + best round-1 shapes.  Any state with 0 bad moduli
in [14,60], v_max >= 65, primitive, covering -> fold check + exact M.
"""
import sys, time, heapq
from math import gcd
from itertools import combinations
sys.path.insert(0, ".")
from lrc_engine import exact_M, covering_check, d_gap, W_of_q, THR, FLOOR

def badness(V):
    """(count of bad q in [14,60], total excess) using exact per-q max-min."""
    nbad = 0
    excess = 0
    for q in range(14, 61):
        res = [v % q for v in V]
        if 0 in res:
            continue
        d = d_gap(q)
        worst = 0
        for a in range(1, q):
            if gcd(a, q) != 1:
                continue
            m = min(min((a * r) % q, q - (a * r) % q) for r in res)
            if m > worst:
                worst = m
        if worst > d:
            nbad += 1
            excess += worst - d
    return nbad, excess

# duty-stacked pool
POOL = sorted(set(
    [v for v in range(15, 301)
     if any(v % q == 0 for q in range(7, 55))]  # carries at least one duty
))

def neighbors(V):
    Vs = set(V)
    for i in range(13):
        for x in POOL:
            if x in Vs:
                continue
            W = sorted(Vs - {V[i]} | {x})
            yield tuple(W)

def main():
    SURV = [tuple(int(x) for x in l.split(":")[1].split())
            for l in open("survivors.txt") if ":" in l]
    seeds = set(SURV)
    seen = {}
    beam = []
    for S in seeds:
        b = badness(S)
        seen[S] = b
        beam.append((b, S))
    beam.sort()
    BW = 400
    t0 = time.time()
    hits = []
    for depth in range(4):
        newbeam = []
        for (b, S) in beam[:BW]:
            for W in neighbors(S):
                if W in seen:
                    continue
                bb = badness(W)
                seen[W] = bb
                if bb[0] == 0:
                    V = list(W)
                    g = 0
                    for x in V:
                        g = gcd(g, x)
                    if (g == 1 and max(V) >= 65 and not covering_check(V)
                            and len(set(V)) == 13):
                        hits.append(V)
                        print(f"  ZERO-BAD state: {V}", flush=True)
                newbeam.append((bb, W))
        beam = sorted(set(beam[:BW]) | set(newbeam))[:BW * 3]
        bcounts = {}
        for (b, S) in beam[:BW]:
            bcounts[b[0]] = bcounts.get(b[0], 0) + 1
        print(f"depth {depth+1}: beam best={beam[0][0]} "
              f"bad-count histogram(top beam)={dict(sorted(bcounts.items()))} "
              f"explored={len(seen)} [{time.time()-t0:.0f}s]", flush=True)
    print(f"\nzero-bad assembled candidates: {len(hits)}")
    for V in hits:
        folds = set()
        for i, a in enumerate(V):
            for b2 in V[i+1:]:
                if a + b2 >= 61:
                    folds.add(a + b2)
                if b2 - a >= 61:
                    folds.add(b2 - a)
            folds.add(2 * a)
        kill = None
        for q in sorted(folds):
            Wq, k = W_of_q(V, q)
            if Wq > d_gap(q):
                kill = (q, k, Wq)
                break
        if kill:
            print(f"V={V}: FOLD-KILL q={kill[0]} (W={kill[2]}>{d_gap(kill[0])})"
                  f" -> M >= {kill[2]}/{kill[0]} = {kill[2]/kill[0]:.5f}")
        else:
            M, (q, k) = exact_M(V)
            tag = "*** IN GAP ***" if FLOOR < M < THR else ""
            print(f"V={V}: EXACT M = {M} = {float(M):.6f} wit {k}/{q} {tag}")

if __name__ == "__main__":
    main()
