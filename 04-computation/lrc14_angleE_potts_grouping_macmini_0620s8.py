#!/usr/bin/env python3
"""
lrc14_angleE_potts_grouping_macmini_0620s8.py  (mac-mini-2026-06-20-S8)  ANGLE E

The Potts / coloring signed partition-function aggregate.

  measS7(E) = Sum_{M subset Z/7} (-1)^{|M|} J(M,E),
  J(M,E)    = meas{ x in [0,1): NO clock e*x (e in E) lands in any arc [j/7,(j+1)/7), j in M }.

Since e=0 in E always lands residue 0, J(M,E)=0 whenever 0 in M; the sum reduces to
M subset {1,...,6} (64 terms).  C1 (PROVED dead) says: per-M extremality FALSE -- consec is
neither per-M max nor min, ~half signed differences negative.

GOAL: find a GROUPING / telescoping of the 64 signed terms that makes consec-maximality
TRUE group-wise.  We test groupings:
  (G1) by |M|  (Bonferroni levels)
  (G2) by orbit of M under the affine/scaling stabilizer of the 7-clock
  (G3) by Mobius / cumulant restructuring (the connected/free-energy log)
  (G4) by gap-structure of the arc-union (number of contiguous runs of arcs)
  (G5) partial-sum prefixes in various canonical orders
For each grouping G partitioning {M}, define Group_g(E)=sum_{M in g}(-1)^|M| J(M,E).
A grouping is "consec-extremal" if for EVERY non-consec k-set E (over the certified box),
Group_g(consec) >= Group_g(E) for every group g  (then summing proves consec max).
"""
import sys, itertools
from fractions import Fraction as F
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

# ---------- core engine ----------
def J(E, M):
    """meas{x: no clock lands in any residue arc in M}."""
    Mset = set(M)
    if 0 in Mset:
        return F(0)  # e=0 always lands residue 0
    E = sorted(set(E)); bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for m in range(0, 7 * e + 1): bps.add(F(m, 7 * e))
    bps = sorted(bps); tot = F(0)
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0: continue
        xm = (x0 + x1) / 2; ok = True
        for e in E:
            if e == 0: continue
            if int(((e * xm) % 1) * 7) in Mset: ok = False; break
        if ok: tot += x1 - x0
    return tot

ALL_M = [tuple(M) for r in range(7) for M in itertools.combinations(range(1, 7), r)]
ALL_M = [()] + [m for m in ALL_M if m != ()]  # ensure empty set first

def Jvec(E):
    """dict M -> J(E,M) for all M subset {1..6}."""
    return {M: J(E, M) for M in ALL_M}

def measS7_from_J(Jv):
    return sum(F((-1) ** len(M)) * v for M, v in Jv.items())

# ---------- groupings ----------
def group_by_size(M):
    return len(M)

# scaling stabilizer of the 7-clock: multiplication by units of Z/7 maps arc j -> (a*j mod 7).
# But measS7 is scale-invariant in E (E and c*E give same value for c coprime to denom),
# and the residue arcs transform; the natural symmetry on M is the multiplicative group (Z/7)^*.
UNITS = [1, 2, 3, 4, 5, 6]
def mul_orbit(M):
    Ms = frozenset(M)
    orb = set()
    for a in UNITS:
        orb.add(frozenset((a * j) % 7 for j in Ms))
    return min(tuple(sorted(o)) for o in orb)  # canonical rep

# gap structure: number of contiguous arc-runs of M on the cycle Z/7 (treating 0 as covered/wall).
def num_runs_cyclic(M):
    if not M: return 0
    s = sorted(M); runs = 1
    for i in range(1, len(s)):
        if s[i] != s[i-1] + 1: runs += 1
    # cyclic wrap on Z/7: but residue 0 is a permanent "wall" (always covered), so do NOT wrap through 0.
    return runs

def main():
    print("=" * 96)
    print("ANGLE E: Potts/coloring signed-partition grouping search")
    print("=" * 96)

    # consec references and full-box adversaries per k.
    # We rely on the PROVED B2 span certificate for span<=N*; here we focus on grouping STRUCTURE
    # using a representative adversary set + exhaustive small-span sets.
    for k in [8]:
        consec = list(range(k))
        Jc = Jvec(consec)
        s7c = measS7_from_J(Jc)
        print(f"\nk={k}: measS7(consec)={s7c}={float(s7c):.6f}")

        # build adversary pool: all primitive k-subsets of {0..N} with 0 in E, up to a small N.
        N = 11
        pool = []
        for combo in itertools.combinations(range(1, N + 1), k - 1):
            E = (0,) + combo
            pool.append(E)
        print(f"  adversary pool size (0..{N}, |E|={k}): {len(pool)}")

        # ---- (G1) by |M| ----
        print("\n  (G1) group by |M| (Bonferroni levels):")
        lvl_consec = defaultdict(lambda: F(0))
        for M, v in Jc.items():
            lvl_consec[len(M)] += F((-1) ** len(M)) * v
        # check group-wise dominance over pool
        g1_fail = 0; g1_worst = {}
        for E in pool:
            Jv = Jvec(E)
            lvl = defaultdict(lambda: F(0))
            for M, v in Jv.items():
                lvl[len(M)] += F((-1) ** len(M)) * v
            for L in range(7):
                if lvl[L] > lvl_consec[L]:
                    g1_fail += 1
                    if L not in g1_worst or lvl[L] - lvl_consec[L] > g1_worst[L][1]:
                        g1_worst[L] = (E, lvl[L] - lvl_consec[L])
                    break
        print(f"    consec level contributions: " + ", ".join(f"L{L}={float(lvl_consec[L]):+.4f}" for L in range(7)))
        print(f"    G1 group-wise dominance FAILS for {g1_fail}/{len(pool)} adversaries")
        for L, (E, d) in sorted(g1_worst.items()):
            print(f"      worst at level {L}: E={E} excess={float(d):+.5f}")

        # ---- (G2) multiplicative orbit ----
        print("\n  (G2) group by (Z/7)^* multiplicative orbit of M:")
        orb_consec = defaultdict(lambda: F(0))
        for M, v in Jc.items():
            orb_consec[mul_orbit(M)] += F((-1) ** len(M)) * v
        g2_fail = 0; g2_worst = None
        for E in pool:
            Jv = Jvec(E)
            orb = defaultdict(lambda: F(0))
            for M, v in Jv.items():
                orb[mul_orbit(M)] += F((-1) ** len(M)) * v
            bad = False
            for g in orb_consec:
                if orb[g] > orb_consec[g]:
                    bad = True
                    if g2_worst is None or orb[g] - orb_consec[g] > g2_worst[2]:
                        g2_worst = (E, g, orb[g] - orb_consec[g])
            if bad: g2_fail += 1
        print(f"    #orbits = {len(orb_consec)}")
        print(f"    G2 group-wise dominance FAILS for {g2_fail}/{len(pool)} adversaries")
        if g2_worst: print(f"      worst: E={g2_worst[0]} orbit={g2_worst[1]} excess={float(g2_worst[2]):+.5f}")

        # ---- (G4) gap structure (#runs) ----
        print("\n  (G4) group by #contiguous arc-runs of M:")
        run_consec = defaultdict(lambda: F(0))
        for M, v in Jc.items():
            run_consec[num_runs_cyclic(M)] += F((-1) ** len(M)) * v
        g4_fail = 0; g4_worst = None
        for E in pool:
            Jv = Jvec(E)
            run = defaultdict(lambda: F(0))
            for M, v in Jv.items():
                run[num_runs_cyclic(M)] += F((-1) ** len(M)) * v
            bad = False
            for g in run_consec:
                if run[g] > run_consec[g]:
                    bad = True
                    if g4_worst is None or run[g] - run_consec[g] > g4_worst[2]:
                        g4_worst = (E, g, run[g] - run_consec[g])
            if bad: g4_fail += 1
        print(f"    consec run contributions: " + ", ".join(f"r{r}={float(run_consec[r]):+.4f}" for r in sorted(run_consec)))
        print(f"    G4 group-wise dominance FAILS for {g4_fail}/{len(pool)} adversaries")
        if g4_worst: print(f"      worst: E={g4_worst[0]} #runs={g4_worst[1]} excess={float(g4_worst[2]):+.5f}")

if __name__ == "__main__":
    main()
