#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Working the open question (S30->S31): does the metagraph reproduce the LRC cusp value 4cos^2(3pi/7)?

ANSWER (this script): the metagraph's APEX LAYER (Z_7 circulant TOURNAMENTS) carries the cyclotomic
gap-values {0.308 (generic), 2.0 (Paley/Fano = flat = OPTIMAL, |Gauss|^2=(1+7)/4)}, with the Paley the
H-MAXIMIZER (H=189 vs 175, THM-586). But the LRC BINDING value 0.198 = 4cos^2(3pi/7) is the DOUBLET gap
(size 2 = THM-578 R-tail), which is SUB-TOURNAMENT -- it sits BELOW the minimal circulant-tournament gap
(0.308). So: the metagraph rehearses the tournament hierarchy and the OPTIMAL (Paley=flat) and GENERIC
values, but the binding NUMBER is one level below the tournament floor, at the 2-element resonance.

Hierarchy of Z_7 autocorrelation gaps by structure SIZE:
  size 1 singleton           1.000
  size 2 DOUBLET (R-tail)    0.198  <- LRC BINDING (global min), SUB-tournament
  size 3 generic circulant   0.308  <- minimal TOURNAMENT gap (6 of 8 circulants)
  size 3 Paley/Fano          2.000  <- OPTIMAL/flat (2 of 8 circulants, H-max)
mac-mini-2026-06-29-S31
"""
from __future__ import annotations
import functools, math, cmath, itertools
print = functools.partial(print, flush=True)
W = cmath.exp(2j * math.pi / 7)


def gap(S):
    return min(abs(sum(W**(k * s) for s in S))**2 for k in range(1, 7))


def H_circ(S):
    """Ham-path count of the Z_7 circulant tournament: arc i->j iff (j-i)%7 in S."""
    fwd = set(S)
    return sum(1 for p in itertools.permutations(range(7))
               if all((p[i+1]-p[i]) % 7 in fwd for i in range(6)))


def H_metagraph_lowtail(n):
    """The low-H tail of the n-tournament metagraph (the cusp): transitive=1, first excited = 3-cycle."""
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    perms = list(itertools.permutations(range(n)))
    seen = {}
    for bits in range(1 << len(pairs)):
        arc = [[False]*n for _ in range(n)]
        for b, (i, j) in enumerate(pairs):
            if (bits >> b) & 1: arc[i][j] = True
            else: arc[j][i] = True
        canon = min(tuple(1 if arc[s[i]][s[j]] else 0 for i in range(n) for j in range(n) if i != j)
                    for s in perms)
        if canon in seen: continue
        seen[canon] = sum(1 for p in perms if all(arc[p[k]][p[k+1]] for k in range(n-1)))
    return sorted(seen.values())


def main():
    print("=" * 80)
    print("CUSP VALUE vs METAGRAPH APEX LAYER -- does the metagraph reproduce 4cos^2(3pi/7)?")
    print("=" * 80)

    # ---- the 8 Z_7 circulant tournaments ----
    print("\n[1] The Z_7 circulant TOURNAMENTS (the metagraph's apex-7 layer):")
    print(f"    {'S (connection set)':>20} {'gap':>8} {'H':>5} {'H/7':>5}  note")
    opts = [(1, 6), (2, 5), (3, 4)]
    for choice in itertools.product(*opts):
        S = set(choice)
        qr = S in ({1, 2, 4}, {3, 5, 6})
        print(f"    {str(sorted(S)):>20} {gap(S):>8.4f} {H_circ(S):>5} {H_circ(S)//7:>5}  "
              f"{'PALEY/Fano = FLAT/OPTIMAL, H-max (THM-586)' if qr else 'generic'}")
    print("    => circulant tournaments span gaps {0.308 (x6), 2.000 (x2 Paley)}; H in {175, 189}.")
    print(f"    Paley gap 2.0 = |Gauss sum|^2 split = (1+7)/4 = 2 (the flat/optimal core, octonion).")

    # ---- the size hierarchy: where 4cos^2(3pi/7) sits ----
    print("\n[2] The Z_7 gap hierarchy by structure SIZE (where the LRC binding value sits):")
    rows = [
        ("size 1  singleton {0}", gap({0}), "trivial"),
        ("size 2  DOUBLET {0,1}", gap({0, 1}), "<- LRC BINDING = 4cos^2(3pi/7) = THM-578 R-tail (SUB-tournament)"),
        ("size 3  generic {1,2,3}", gap({1, 2, 3}), "<- minimal TOURNAMENT gap (the metagraph apex floor)"),
        ("size 3  Paley {1,2,4}", gap({1, 2, 4}), "<- OPTIMAL/flat (the metagraph apex ceiling)"),
    ]
    for name, g, note in rows:
        print(f"    {name:>26}: gap = {g:.4f}  {note}")
    print(f"\n    THE KEY FACT: 0.198 (doublet) < 0.308 (minimal tournament). The LRC binds BELOW the")
    print(f"    metagraph's TOURNAMENT floor -- at the 2-element resonance (the R-tail), one level down.")
    print(f"    A doublet = a single DIFFERENCE d=b-a => the EVEN-GRAPH/Cayley-edge object, not a tournament.")

    # ---- the metagraph cusp (universal, H->1) vs the apex layer ----
    print("\n[3] The metagraph cusp (universal H->1) vs the apex-7 value:")
    for n in (4, 5, 6, 7):
        lt = H_metagraph_lowtail(n) if n <= 6 else None
        if lt is None:
            print(f"    n={n}: (skipped, too large for exhaustive iso-canon; Paley-7 lives here, H=189)")
        else:
            print(f"    n={n}: H low tail = {lt[:6]}...  (transitive H=1 = cusp; 3-cycle H=3 = 1st excited)")
    print("    => the UNIVERSAL cusp neighbor (3-cycle) is a Z_3 object (gap 1.0), NOT Z_7. The apex VALUE")
    print("    4cos^2(3pi/7) is Z_7-specific (from 14=2.7) and lives in the n=7 layer, not the universal cusp.")

    print("\n" + "=" * 80)
    print("CONCLUSION (honest split of the rehearsal):")
    print("  STRUCTURAL (transfers): cusp binds at the minimal/degenerate end; OPTIMAL = the regular Paley")
    print("    (flat, H-max). The metagraph EXPLAINS why the cusp binds low and the optimum is the Paley.")
    print("  NUMERICAL (apex-specific): the values {0.308, 2.0} ARE the metagraph's Z_7 circulant-tournament")
    print("    gaps; but the BINDING 0.198 is the DOUBLET (R-tail), SUB-tournament, below the apex floor.")
    print("  NEW OBJECT TO TRACK: the doublet = a single Z_7 DIFFERENCE = an EVEN-GRAPH/Cayley edge; the")
    print("    even-graph dual E_n is the natural home of the sub-tournament binding value.")
    print("=" * 80)


if __name__ == "__main__":
    main()
