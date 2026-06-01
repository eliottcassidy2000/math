#!/usr/bin/env python3
"""
Tournament Analysis of the n=14 sieve witnesses (project directive).
opus-2026-06-01-S551 (remote-control).

At a sieve witness t=a/q the 13 runners sit at distinct circle positions
x_i = ((v_i a) mod q)/q, ALL inside the safe band [1/n, 1 - 1/n] (they avoid the
two danger caps of total length 2/n around the observer at 0).

The half-turn runner tournament:  i -> j  iff  (x_i - x_j) mod 1 in (0, 1/2).
S522-S525 established that the lonely runner sub-tournament is always ROUND
(points on a circle) with #SCC in {1, m}.  The SIEVE adds a sharp refinement:

  loneliness <=> there is a witness arc, i.e. a GAP of length >= 2/n straddling
  the observer.  So the realizable lonely tournaments are exactly the round
  tournaments whose largest gap (>= 2/n) is the one containing the observer's
  danger zone -- a strongly RESTRICTED iso-class set inside A000568(13).

This script quantifies that restriction: over many n=14 lonely witnesses it
records the runner tournament's score sequence (out-degrees) and #SCC, and
compares the realized variety to the unrestricted count.
"""

from fractions import Fraction
from math import gcd
import random
import importlib.util, os

_spec = importlib.util.spec_from_file_location(
    "sieve_mod",
    os.path.join(os.path.dirname(__file__), "lrc_n14_multiprime_sieve_s551.py"))
S = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(S)


def witness(V, n, qmax=200):
    for q in range(2, qmax + 1):
        for a in range(1, q):
            if gcd(a, q) != 1:
                continue
            if all(S.safe_residue((v * a) % q, q, n) for v in V):
                return a, q
    return None, None


def half_turn_tournament(positions):
    """positions: list of Fractions in [0,1). Returns out-degree sequence + SCC."""
    m = len(positions)
    adj = [[False] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i == j:
                continue
            d = (positions[i] - positions[j]) % 1
            adj[i][j] = (0 < d < Fraction(1, 2))
    score = sorted(sum(adj[i]) for i in range(m))
    # Tarjan-free SCC via reachability (m small)
    def reach(start):
        seen = {start}
        stack = [start]
        while stack:
            u = stack.pop()
            for v in range(m):
                if adj[u][v] and v not in seen:
                    seen.add(v)
                    stack.append(v)
        return seen
    R = [reach(i) for i in range(m)]
    # i,j in same SCC iff j in R[i] and i in R[j]
    comp = [-1] * m
    c = 0
    for i in range(m):
        if comp[i] == -1:
            members = [j for j in range(m) if (j in R[i] and i in R[j])]
            for j in members:
                comp[j] = c
            c += 1
    n_scc = len(set(comp))
    return tuple(score), n_scc


def largest_gap_straddles_observer(positions, n):
    """Is the largest circular gap the one containing the observer (point 0)?"""
    pts = sorted(positions)
    m = len(pts)
    gaps = []
    for k in range(m):
        a = pts[k]
        b = pts[(k + 1) % m] + (1 if k == m - 1 else 0)
        gaps.append((b - a, a, b % 1 if k == m - 1 else b))
    biggest = max(gaps)
    # the observer-straddling gap is the one wrapping past 1 -> 0 (last gap)
    obs_gap = gaps[-1][0]
    return obs_gap == biggest[0], biggest[0], obs_gap


def run(n=14, n_sets=300, seed=2024):
    N = n - 1
    rng = random.Random(seed)
    print(f"====== Tournament Analysis of LRC sieve witnesses, n={n} ======\n")
    scores = {}
    sccs = {}
    straddle_ok = 0
    total = 0
    band_ok = 0
    seen = set()
    tries = 0
    while total < n_sets and tries < n_sets * 40:
        tries += 1
        g = rng.randint(0, 3)
        if g == 0:
            V = tuple(sorted(rng.sample(range(1, 80), N)))
        elif g == 1:                              # loaded-ish
            V = tuple(sorted(rng.sample(range(1, 200), N)))
        elif g == 2:                              # near-AP perturbation
            base = list(range(1, N + 1))
            for _ in range(rng.randint(1, 4)):
                base[rng.randrange(N)] += rng.choice([n, 2 * n, -n]) or n
            V = tuple(sorted(set(abs(x) + 1 for x in base)))
            if len(V) != N:
                continue
        else:
            V = tuple(sorted(rng.sample(range(1, 400), N)))
        g0 = 0
        for v in V:
            g0 = gcd(g0, v)
        V = tuple(sorted(v // g0 for v in V))
        if V in seen:
            continue
        seen.add(V)
        a, q = witness(V, n)
        if q is None:
            continue
        total += 1
        pos = [Fraction((v * a) % q, q) for v in V]
        # band check (tautological for a witness, but verify)
        if all(Fraction(1, n) <= p <= 1 - Fraction(1, n) for p in pos):
            band_ok += 1
        sc, nscc = half_turn_tournament(pos)
        scores[sc] = scores.get(sc, 0) + 1
        sccs[nscc] = sccs.get(nscc, 0) + 1
        ok, big, obs = largest_gap_straddles_observer(pos, n)
        if ok:
            straddle_ok += 1

    print(f"  witnessed lonely sets analysed: {total}")
    print(f"  runner points inside safe band [1/{n},1-1/{n}]: {band_ok}/{total} "
          f"(must be {total})")
    print(f"  #SCC histogram of runner half-turn tournament: "
          f"{dict(sorted(sccs.items()))}")
    print(f"     (S525 predicts #SCC in {{1, {N}}} only)")
    print(f"  largest circular gap straddles the observer: "
          f"{straddle_ok}/{total}")
    print(f"  distinct realised score sequences: {len(scores)}")
    top = sorted(scores.items(), key=lambda kv: -kv[1])[:6]
    print(f"  most common score sequences (out-degrees, sorted):")
    for sc, ct in top:
        print(f"     {ct:4d}x  {sc}")
    # the transitive score seq (0,1,...,N-1) = semicircle packing
    trans = tuple(range(N))
    print(f"  transitive score seq {trans} realised: "
          f"{scores.get(trans, 0)} times (= runners in a semicircle)")


if __name__ == "__main__":
    run(n=14)
