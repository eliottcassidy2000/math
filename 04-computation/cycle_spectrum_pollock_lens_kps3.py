#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Pollock-completeness lens on tournament cycle-count spectra.
kind-pasteur-2026-06-13-S3.

Pollock's conjecture asks: is a fixed family (tetrahedral/octahedral numbers) an
ADDITIVE BASIS of bounded order — is EVERY integer representable?  The repo's
analogue: is a tournament invariant's value-set GAP-FREE (every value in [0,max]
realized = "complete")?  THM-462 PROVED the cubic channel c3 is gap-free [0,M(n)]
(its engine: Lagrange four-square on the score deviation f = m + (sum e_i^2)/2 +
induction = a Waring/Pollock-style arity argument).  The degree-n invariant H has
GAPS (forbidden 7, 21).  So Pollock-completeness HOLDS at degree 3 and FAILS at
degree n; WHERE does it break?  The c5 (directed 5-cycle count) channel is the next
probe — and c5 is NOT score-determined (unlike c3), so this is genuinely new.

This script (exact, exhaustive small n + validation):
 (A) recompute the c3 spectrum, confirm gap-free [0,M(n)] (validates pipeline vs THM-462).
 (B) the c5 spectrum: exhaustive over labeled tournaments n<=6, heavy sample n=7;
     report value-set, max, and ALL gaps (missing values in [0,max]).
 (C) the Pollock verdict per channel: complete (gap-free) or where the first gap is.
"""

import sys, itertools, random
from math import comb
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None


def edge(adj, i, j):
    return adj[i][j]


def c3_count(adj, n):
    """directed 3-cycles = cyclic triples."""
    c = 0
    for a, b, cc in itertools.combinations(range(n), 3):
        # cyclic iff a->b->c->a or a->c->b->a
        if adj[a][b] and adj[b][cc] and adj[cc][a]:
            c += 1
        if adj[a][cc] and adj[cc][b] and adj[b][a]:
            c += 1
    return c


def c5_count(adj, n):
    """total directed 5-cycles: over each 5-subset, count directed Ham cycles
    (cyclic sequences with all 5 consecutive arcs present)."""
    c = 0
    for S in itertools.combinations(range(n), 5):
        s0 = S[0]
        rest = S[1:]
        # directed Ham cycles fixing start s0: 4! orderings of rest, cycle s0->p0->..->p3->s0
        for perm in itertools.permutations(rest):
            ok = adj[s0][perm[0]]
            if not ok: continue
            good = True
            for x in range(4 - 1):
                if not adj[perm[x]][perm[x + 1]]:
                    good = False; break
            if good and adj[perm[3]][s0]:
                c += 1
    return c


def all_tournaments(n):
    """yield adjacency matrices (adj[i][j]=True iff i->j) over all labeled tournaments."""
    pairs = list(itertools.combinations(range(n), 2))
    for bits in range(1 << len(pairs)):
        adj = [[False] * n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1:
                adj[i][j] = True
            else:
                adj[j][i] = True
        yield adj


def random_tournament(n, rng):
    adj = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5:
                adj[i][j] = True
            else:
                adj[j][i] = True
    return adj


def M3(n):
    return (n**3 - n) // 24 if n % 2 else (n**3 - 4 * n) // 24


def gaps(values):
    vmax = max(values)
    return [v for v in range(0, vmax + 1) if v not in values]


def main():
    print("=== (A) c3 spectrum — must be gap-free [0,M(n)] (THM-462 validation) ===", flush=True)
    for n in range(3, 7):
        vals = set()
        for adj in all_tournaments(n):
            vals.add(c3_count(adj, n))
        g = gaps(vals)
        print(f"   n={n}: c3 in [0,{max(vals)}], M(n)={M3(n)}, "
              f"max-correct={max(vals)==M3(n)}, gaps={len(g)} "
              f"=> {'GAP-FREE (Pollock-complete)' if not g else 'GAPS '+str(g[:10])}", flush=True)

    print("\n=== (B) c5 spectrum — the next Pollock channel (NOT score-determined) ===", flush=True)
    for n in range(5, 7):
        vals = set()
        for adj in all_tournaments(n):
            vals.add(c5_count(adj, n))
        g = gaps(vals)
        print(f"   n={n} (exhaustive {1<<comb(n,2)} labeled): c5 in [0,{max(vals)}], "
              f"#values={len(vals)}, gaps={len(g)} "
              f"=> {'GAP-FREE' if not g else 'FIRST GAPS '+str(g[:12])}", flush=True)

    # n=7 heavy sample (exhaustive too slow: 2^21=2.1M * cost) — sample + structured
    print("   n=7: heavy random sample (exhaustive 2^21 feasible but slow; sampling for the spectrum):", flush=True)
    rng = random.Random(2026613)
    vals7 = set()
    for _ in range(200000):
        adj = random_tournament(7, rng)
        vals7.add(c5_count(adj, 7))
    g7 = gaps(vals7)
    print(f"      c5 in [0,{max(vals7)}] (sampled), #values={len(vals7)}, "
          f"missing-in-range={len(g7)} (sample lower bound on gaps): "
          f"{'no gaps seen' if not g7 else 'candidate gaps '+str(sorted(g7)[:15])}", flush=True)

    print("\n=== (C) Pollock verdict per channel ===", flush=True)
    print("   c3 (degree-3): gap-free at all n (THM-462, Lagrange-four-square engine) = COMPLETE.", flush=True)
    print("   c5 (degree-5): see (B) — first gaps locate where additive completeness weakens.", flush=True)
    print("   H  (degree-n): GAPS at 7, 21 (THM-029/079) = INCOMPLETE. The Pollock-onset is between.", flush=True)


if __name__ == "__main__":
    main()
