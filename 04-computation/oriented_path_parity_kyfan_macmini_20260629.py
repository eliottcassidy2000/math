#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Is there a 'Redei for ORIENTED paths' -- a parity per oriented type -- and is it
the Ky Fan / signed-cycle-index shadow?  (mac-mini-2026-06-29-S12)

MERGE: arXiv:2512.09332 (El Sahili-El Zein, 'Oriented Hamiltonian Paths in Tournaments:
Stability under Arc Deletion') opens with Redei (every tournament has an ODD number of
DIRECTED Ham paths) and studies arbitrary ORIENTED types under arc deletion.  klein THM-584:
the arc-flip graph is the hypercube Q_d, complement = the ANTIPODAL map.  An oriented TYPE
tau in {+,-}^{n-1} is a corner of a 'type hypercube'; reversal/complement act on it.

QUESTION: for which oriented types tau is N_tau(T) = #{Ham paths of type tau in T} a FIXED
PARITY across all tournaments T?  (Redei: tau=all-+ is always ODD.)  A type-indexed parity
law = a 'Redei for oriented paths', the Ky-Fan/Borsuk-Ulam odd-count graded by type.

We compute N_tau(T) mod 2 over ALL tournaments (n=4,5) and a sample (n=6), per type, and
report which types have constant parity, and the structure under reversal-flip rho(tau) and
complement -tau.
"""
from __future__ import annotations
import functools, itertools, random
print = functools.partial(print, flush=True)


def all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(pairs)):
        arc = [[False] * n for _ in range(n)]
        for b, (i, j) in enumerate(pairs):
            if (bits >> b) & 1:
                arc[i][j] = True
            else:
                arc[j][i] = True
        yield arc


def count_by_type(arc, n):
    """For each path (perm), its type tau in {0,1}^{n-1} (1=forward arc v_k->v_{k+1});
    returns dict type-tuple -> count.  Every perm IS a path in the tournament (complete),
    with some arcs forward (in T) and some backward -- the 'type' records which."""
    from collections import Counter
    cnt = Counter()
    for perm in itertools.permutations(range(n)):
        tau = tuple(1 if arc[perm[k]][perm[k + 1]] else 0 for k in range(n - 1))
        cnt[tau] += 1
    return cnt


def main():
    print("=" * 80)
    print("Redei for ORIENTED paths? parity of N_tau(T) per type (mac-mini-S12)")
    print("=" * 80)
    print("\nNote: a 'type tau' here = which of the n-1 consecutive arcs are forward in T.")
    print("tau=all-1 (=all-forward) is the DIRECTED Ham path (Redei: ODD). We scan all types.\n")

    for n in (4, 5, 6):
        types = list(itertools.product((0, 1), repeat=n - 1))
        parities = {t: set() for t in types}
        if n <= 5:
            Ts = all_tournaments(n)
            label = "ALL"
        else:
            rng = random.Random(7)
            def sample():
                for _ in range(3000):
                    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
                    arc = [[False] * n for _ in range(n)]
                    for (i, j) in pairs:
                        if rng.random() < 0.5:
                            arc[i][j] = True
                        else:
                            arc[j][i] = True
                    yield arc
            Ts = sample()
            label = "3000-sample"
        ntour = 0
        for arc in Ts:
            ntour += 1
            cnt = count_by_type(arc, n)
            for t in types:
                parities[t].add(cnt.get(t, 0) % 2)
        const_odd = [t for t in types if parities[t] == {1}]
        const_even = [t for t in types if parities[t] == {0}]
        variable = [t for t in types if parities[t] == {0, 1}]
        def fwd(t):
            return sum(t)
        print(f"n={n} ({label}, {ntour} tournaments): {len(types)} types")
        print(f"  CONSTANT-ODD types ({len(const_odd)}): " +
              ", ".join("".join('+' if b else '-' for b in t) + f"(fwd={fwd(t)})" for t in const_odd))
        print(f"  CONSTANT-EVEN types ({len(const_even)}): " +
              ", ".join("".join('+' if b else '-' for b in t) for t in const_even[:12]) +
              (" ..." if len(const_even) > 12 else ""))
        print(f"  VARIABLE-parity types: {len(variable)}")
        # structure: directed (all-1) and its complement (all-0) and reverse
        allf = tuple([1] * (n - 1)); allb = tuple([0] * (n - 1))
        print(f"  directed all-+ parity: {parities[allf]} (Redei=odd); all-- parity: {parities[allb]}")

    print("\n" + "=" * 80)
    print("If a NONTRIVIAL set of types is constant-ODD, that's a 'Redei for oriented paths' --")
    print("the graded/Ky-Fan odd-count by type, on the type-hypercube {+,-}^{n-1}.")
    print("=" * 80)


if __name__ == "__main__":
    main()
