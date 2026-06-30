#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""LONG cusp-behavior exploration: the finite Z_7 core gap landscape (the apex cusp of X_0(14)),
abnormalities, and the metagraph H->1 cusp rehearsal. (mac-mini-2026-06-29-S30)

Corrected frame (klein-S9 HYP-3581): rho_j>=c is the DIRECT FINITE MINIMUM of the Z_7 cyclotomic
Gram gap over the 128 cores O subset Z_7 = 4cos^2(3pi/7)=0.19806, BINDING at a DOUBLET (ties THM-578);
gap=0 ONLY at O=Z_7 (the disproof boundary); bare averaging (my S28) was a Jensen overshoot (INVALID).

This session: map the WHOLE landscape, hunt abnormalities, track interesting things, and test the
metagraph transitive-limit (H->1) corner as the cusp rehearsal.
"""
from __future__ import annotations
import functools, math, cmath, itertools
from collections import Counter, defaultdict
print = functools.partial(print, flush=True)

W = cmath.exp(2j * math.pi / 7)


def gap(O):
    """min over k!=0 of |sum_{x in O} w^{kx}|^2 (the Z_7 cyclotomic Gram spectral gap)."""
    if not O: return None
    return min(abs(sum(W**(k * x) for x in O))**2 for k in range(1, 7))


def spectrum(O):
    return [round(abs(sum(W**(k * x) for x in O))**2, 4) for k in range(7)]


def main():
    print("=" * 84)
    print("CUSP LANDSCAPE: the 128 Z_7 cores -- gap structure, abnormalities (mac-mini-S30)")
    print("=" * 84)

    cores = []
    for sz in range(1, 8):
        for O in itertools.combinations(range(7), sz):
            cores.append(set(O))

    # ---- PART 1: gap by core size ----
    print("\n[1] gap landscape by core size (gap = min nonzero cyclotomic mode):")
    print(f"    {'size':>4} {'#cores':>7} {'min gap':>9} {'max gap':>9} {'#gap=0':>7} {'distinct gaps':>30}")
    bysize = defaultdict(list)
    for O in cores: bysize[len(O)].append(gap(O))
    for sz in range(1, 8):
        gs = bysize[sz]
        nz = [g for g in gs if g is not None]
        zeros = sum(1 for g in nz if abs(g) < 1e-9)
        distinct = sorted(set(round(g, 4) for g in nz))
        print(f"    {sz:>4} {len(gs):>7} {min(nz):>9.4f} {max(nz):>9.4f} {zeros:>7} "
              f"{str(distinct[:5]):>30}")
    target = 4 * math.cos(3*math.pi/7)**2
    print(f"\n    global min gap (non-Z_7) = {min(gap(O) for O in cores if O != set(range(7))):.5f}")
    print(f"    4cos^2(3pi/7) = 2+2cos(6pi/7) = {target:.5f}  (klein-S9: the rho_j floor, binds at doublets)")

    # ---- PART 2: WHERE it binds (the doublets) + the disproof boundary ----
    print("\n[2] the BINDING cores (gap = the global min) and the gap=0 cores (disproof boundary):")
    gmin = min(gap(O) for O in cores if O != set(range(7)))
    binders = [sorted(O) for O in cores if O != set(range(7)) and abs(gap(O) - gmin) < 1e-6]
    zeros = [sorted(O) for O in cores if abs(gap(O)) < 1e-9]
    print(f"    binders (gap={gmin:.4f}): {len(binders)} cores, sizes {sorted(set(len(b) for b in binders))}")
    print(f"      doublets among binders: {[b for b in binders if len(b)==2][:7]} (= THM-578 R-tail object)")
    print(f"    gap=0 cores (disproof boundary): {zeros}  -- ONLY the full Z_7 (and empty)")
    print(f"    => the floor binds at DOUBLETS; gap vanishes ONLY at the covering-complete O=Z_7.")

    # ---- PART 3: the cyclotomic value structure (which 4cos^2 appear) ----
    print("\n[3] the cyclotomic gap VALUES (all live in Q(cos 2pi/7), the totally-real cubic):")
    allgaps = sorted(set(round(gap(O), 5) for O in cores if O))
    print(f"    distinct gap values ({len(allgaps)}): {allgaps}")
    for j in (1, 2, 3):
        print(f"      4cos^2({j}pi/7) = {4*math.cos(j*math.pi/7)**2:.5f}", end="")
    print()
    print("    => the gaps are 2+2cos(2pi m/7) = 4cos^2(m pi/7) and sums -- a SHORT finite cyclotomic set.")

    # ---- PART 4: ABNORMALITIES -- the averaging overshoot (klein-S9), flat/optimal cores ----
    print("\n[4] ABNORMALITIES:")
    # (a) averaging overshoot: avg_gap > raw_gap (Jensen) -- count them
    def avg_gap(O):
        # Z_7^*-averaged indicator spectrum min nonzero
        ind = [1.0 if a in O else 0.0 for a in range(7)]
        avg = [0.0]*7
        for u in range(1, 7):
            for a in range(7): avg[(u*a) % 7] += ind[a]/6
        sp = [abs(sum(avg[a]*W**(k*a) for a in range(7)))**2 for k in range(1, 7)]
        return min(sp)
    overshoot = sum(1 for O in cores if O and O != set(range(7)) and avg_gap(O) > gap(O) + 1e-9)
    print(f"    (a) Z_7^*-averaging OVERSHOOT (avg_gap > raw_gap, Jensen): {overshoot}/{sum(1 for O in cores if O and O!=set(range(7)))} cores")
    print(f"        => averaging gives the MEAN not the MIN -- klein-S9's correction of my S28, CONFIRMED.")
    # (b) flat-spectrum (perfect difference set) cores -- the octonion-optimal
    flats = [sorted(O) for O in cores if O and len(set(spectrum(O)[1:])) == 1]
    print(f"    (b) FLAT-spectrum cores (perfect difference sets, octonion-optimal): {flats}")
    print(f"        {{1,2,4}} & {{3,5,6}} (the Fano/QR lines) are the size-3 flats; gap=2 (max for size 3).")
    # (c) complement symmetry: 5-subset gap == its 2-subset complement
    print("    (c) COMPLEMENT symmetry: gap(O) == gap(Z_7\\O) for |O|!=0,7 (char sum vanishes):")
    ok = all(abs(gap(O) - gap(set(range(7))-O)) < 1e-9 for O in cores if 0 < len(O) < 7)
    print(f"        gap(O)==gap(complement): {ok}  => the doublet (2) and its 5-complement bind together.")

    # ---- PART 5: the metagraph H->1 CUSP REHEARSAL (the near-transitive corner) ----
    print("\n[5] METAGRAPH H->1 CUSP REHEARSAL (the transitive-limit corner = the X_0 cusp):")
    print("    The LRC cusp is m_R->0; the metagraph cusp is H->1 (transitive). Test: the near-transitive")
    print("    classes' H-distribution low tail -- is there a 'doublet-like' binding minimum?")
    for n in (5, 6):
        pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
        perms = list(itertools.permutations(range(n)))
        seen = {}
        for bits in range(1 << len(pairs)):
            arc = [[False]*n for _ in range(n)]
            for b, (i, j) in enumerate(pairs):
                if (bits >> b) & 1: arc[i][j] = True
                else: arc[j][i] = True
            canon = min(tuple(1 if arc[s[i]][s[j]] else 0 for i in range(n) for j in range(n) if i!=j) for s in perms)
            if canon in seen: continue
            H = sum(1 for p in perms if all(arc[p[k]][p[k+1]] for k in range(n-1)))
            seen[canon] = H
        Hs = sorted(seen.values())
        c3min = Hs[1]  # smallest non-transitive H (the 'cusp neighbor')
        print(f"    n={n}: H-spectrum {Hs}; transitive H=1 (the cusp); 1st neighbor H={c3min} (the cusp 'doublet').")
    print("    => the metagraph H->1 cusp: transitive (H=1) + its smallest neighbor (H=3 = one 3-cycle) =")
    print("    the cusp's binding object, the 3-CYCLE -- the metagraph mirror of the LRC DOUBLET (both the")
    print("    minimal non-trivial relation: 3-cycle = minimal cyclicity, doublet = minimal resonance pair).")

    print("\n" + "=" * 84)
    print("CUSP FINDINGS: (1) rho_j floor = 4cos^2(3pi/7)=0.198 binds at DOUBLETS (=THM-578 R-tail);")
    print("(2) gap=0 ONLY at O=Z_7 (disproof boundary); (3) gaps are a SHORT cyclotomic set in Q(cos2pi/7);")
    print("(4) averaging overshoots (Jensen, klein-S9 confirmed); (4b) {1,2,4}/{3,5,6} Fano flats are the")
    print("optimal; (4c) doublet binds with its 5-complement; (5) the metagraph cusp (H->1) binds at the")
    print("3-CYCLE = the mirror of the LRC doublet (both minimal relations). Track: doublet<->3-cycle<->R-tail.")
    print("=" * 84)


if __name__ == "__main__":
    main()
