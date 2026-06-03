#!/usr/bin/env python3
r"""
lrc_64_classes_antipodal_tieresolution_s577.py    oracle-2026-06-03-S577o

EXTENDING the 64 self-converse round classes of LRC n=14 (the worry-set, HYP-2094/2086).

THE KEY SUBTLETY (and the new idea). HYP-2091 says even n=14 is the CLEAN polygon ladder
because the RUNNER tournament size m=n-1=13 is ODD (no antipodal ties at generic open
times). BUT the TIGHT configs sit at the n-clock time t=j/n: the AP {1..13} at t=1/14 puts
the 13 runners at positions {1,...,13}/14 -- a 14-GON. Since n=14 is EVEN, these pair up
ANTIPODALLY: {k, 14-k} for k=1..6 (6 pairs, distance 7/14=1/2), plus the self-antipodal
7/14=1/2 (the apex). So the tight boundary REINTRODUCES the tie-wall, and there are exactly
6 tie-pairs -> 2^6 = 64 resolutions.

CONJECTURE (this script tests it):
  the 64 SELF-CONVERSE round classes  ==  the 2^6 TIE-RESOLUTIONS of the 6 antipodal pairs
  of the 14-gon  ==  the 2^6 orientation choices of the 6 mod-7 CRT pair-classes.
Because position k/14 and (14-k)/14 share residue k mod 7 (k and -k), the 6 antipodal pairs
(1,8),(2,9),(3,10),(4,11),(5,12),(6,13) ARE the 6 nonzero mod-7 classes; position 7 (=apex)
is residue 0. This would unify HYP-2091 (tie-mesh), HYP-2086 (self-converse), oracle-S552o
(mod-7 CRT), and opus-S568 (antipodal transversals mod 2n-1).

We build all 2^6 tie-resolved tournaments, canonicalize (reuse S574 canon), and compare the
set to the generator's 64 self-converse round classes.
"""
import importlib.util
from itertools import product
from pathlib import Path

def load_s574():
    spec = importlib.util.spec_from_file_location(
        "s574", Path("04-computation/lrc_round_count_m89_s574.py").resolve())
    M = importlib.util.module_from_spec(spec); spec.loader.exec_module(M)
    return M

MOD = load_s574()
canon = MOD.canon

N = 14
M = N - 1                      # 13 runners
PTS = list(range(1, N))        # positions 1..13 (units of 1/14) on the 14-gon

def half_turn_adj(resolution):
    """tournament on the 13 runners at 14-gon positions 1..13.
    Non-antipodal pairs: i beats j if (pos_j - pos_i) mod 14 in (0,7).
    The 6 antipodal pairs {k,14-k} (k=1..6) are TIES resolved by `resolution` (6 bits):
    bit=0 -> the smaller-position runner beats the larger; bit=1 -> reverse."""
    idx = {p: r for r, p in enumerate(PTS)}     # position -> runner index 0..12
    adj = [[0] * M for _ in range(M)]
    anti = [(k, N - k) for k in range(1, N // 2)]   # (1,13)? no: distance 7 pairs
    anti_pairs = [(k, k + N // 2) for k in range(1, N // 2)]   # (1,8),(2,9),...,(6,13)
    bit = {}
    for b, (a, c) in zip(resolution, anti_pairs):
        bit[(a, c)] = b
    for a in PTS:
        for c in PTS:
            if a >= c:
                continue
            diff = (c - a) % N
            if diff == N // 2:                  # antipodal tie -> resolve
                b = bit[(a, c)]
                if b == 0:
                    adj[idx[a]][idx[c]] = 1
                else:
                    adj[idx[c]][idx[a]] = 1
            elif 0 < diff < N // 2:
                adj[idx[a]][idx[c]] = 1
            else:
                adj[idx[c]][idx[a]] = 1
    return adj

def opp(a):
    m = len(a)
    return [[0 if i == j else a[j][i] for j in range(m)] for i in range(m)]

def selfconverse_round_classes():
    reps = {}
    for d in MOD.valid_dvectors(M):
        a = MOD.build_adj(d, M); reps.setdefault(canon(a, M), a)
    sc = {c for c, a in reps.items() if c == canon(opp(a), M)}
    return reps, sc

def main():
    print("=" * 80)
    print("LRC n=14: are the 64 SELF-CONVERSE round classes the 2^6 antipodal tie-resolutions?")
    print("=" * 80)

    print("\n(1) The 6 antipodal pairs of the 14-gon (tight time t=1/14) and their mod-7 class:")
    for k in range(1, N // 2):
        a, c = k, k + N // 2
        print(f"   pair (pos {a:2d}, pos {c:2d})  distance {(c-a)}/14 = 1/2   "
              f"both = residue {a % (N//2)} mod 7")
    print(f"   self-antipodal: pos 7 = 1/2 = residue 0 mod 7 = the APEX (multiple of 7)")
    print(f"   => 6 tie-pairs => 2^6 = {2**(N//2 - 1)} resolutions; 6 nonzero mod-7 classes.")

    print("\n(2) Build all 2^6 tie-resolved tournaments, canonicalize:")
    tie_classes = {}
    sc_among_ties = 0
    for res in product((0, 1), repeat=N // 2 - 1):
        adj = half_turn_adj(res)
        c = canon(adj, M)
        tie_classes.setdefault(c, []).append(res)
        if c == canon(opp(adj), M):
            sc_among_ties += 1
    n_tie_classes = len(tie_classes)
    print(f"   2^6 = 64 resolutions -> {n_tie_classes} distinct canonical classes"
          f"   (self-converse among the 64 resolutions: counted below)")

    print("\n(3) Compare to the generator's self-converse round classes:")
    reps, sc = selfconverse_round_classes()
    print(f"   round classes (A000016) = {len(reps)};  self-converse = {len(sc)}")
    tie_set = set(tie_classes)
    print(f"   tie-resolution classes  = {n_tie_classes}")
    print(f"   tie-resolution set == self-converse set ?  {tie_set == sc}")
    print(f"   tie ⊆ self-converse ? {tie_set <= sc}    self-converse ⊆ tie ? {sc <= tie_set}")
    if tie_set != sc:
        print(f"   in tie not sc: {len(tie_set - sc)};  in sc not tie: {len(sc - tie_set)}")
        # are all round? check tie classes are round
        roundset = set(reps)
        print(f"   tie classes that ARE round: {len(tie_set & roundset)}/{n_tie_classes}")

    print("\n(4) Are the 64 resolutions all DISTINCT classes (free 6-bit index)?")
    sizes = sorted(len(v) for v in tie_classes.values())
    from collections import Counter
    print(f"   resolution-count per class histogram: {dict(Counter(sizes))}")
    print(f"   (all 1 => the 6 bits are a free coordinate; the 64 classes <-> 2^6 mod-7 orientations)")

    print("\n" + "=" * 80)
    print("READING")
    print("=" * 80)
    print(f"""  If the tie-resolution set equals the 64 self-converse classes, then the LRC n=14
  worry-set is EXACTLY the 2^6 orientations of the 6 mod-7 antipodal pairs of the 14-gon
  at the tight time t=1/14. The 'clean even ladder' (HYP-2091) is clean only at GENERIC
  open times; the TIGHT boundary is a STRUCTURED tie-wall with 6 independent pairs = the
  6 mod-7 CRT classes (oracle-S552o), the apex being the self-antipodal residue 0. This
  unifies the self-converse count (HYP-2086), the parity ladder (HYP-2091), the mod-7
  reduction, and the antipodal transversals mod 2n-1 (opus-S568): they are one object.""")

if __name__ == "__main__":
    main()
