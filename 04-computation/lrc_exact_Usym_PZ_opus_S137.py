"""
lrc_exact_Usym_PZ_opus_S137.py

REVERSAL-SYMMETRIZED U MOMENTS, EXACT — feeding mac-mini's PZ-on-U machinery and the
two hard cores (owner worklist items: "reversal-symmetrized moments", "the two hard cores
as the concrete targets for opus's avoidance-kernel and mac-mini's PZ-on-U machinery").

U(x) = sum_gaps (gap - 1/7)_+ (opus-S131's uncovered length; mu_{1/7} >= P(U>0) and
PZ: mu >= E[U]^2/E[U^2], mac-mini-S41's CV form).  KEY OBSERVATION: U(1-x) = U(rev E, x),
so the pointwise symmetrization
      U_sym(x) = (U(x) + U(1-x))/2
has E[U_sym] = E[U] (reversal-invariance of the MEAN) but
      E[U_sym^2] = (E[U^2] + E[U(x)U(1-x)])/2  <=  E[U^2]
whenever the cross-moment is below the second moment (Cauchy-Schwarz), so
      PZ_sym = E[U]^2 / E[U_sym^2]  >=  PZ = E[U]^2 / E[U^2].
The GAIN is strict unless U is a.e. reversal-symmetric.  mac-mini-S43 found the E[U]-min
family is a mirror PAIR (palindromy broken) — exactly the situation where symmetrization
buys the most.

This script computes EXACT (rational) E[U], E[U^2], E[U(x)U(1-x)], PZ, PZ_sym for the
record/board families, and the same on monad-S3's TWO HARD CORES with the G_P-restricted
versions (restrict all integrals to G_P: the intersected-ledger coupling): exact
      rho* >= E_GP[U]^2 / E_GP[U^2]-type floors on the load-bearing domain.

Machinery: the order-cell engine gives U as an explicit piecewise-affine function
(per micro-cell after refining by each gap's theta-crossing); moments = exact piecewise
polynomial integration; the cross moment uses the common refinement with the reflected
breakpoints; G_P-restriction intersects pieces with the G_P interval list.
"""
from fractions import Fraction as F
import time, sys

sys.path.insert(0, ".")
from lrc_exact_mu_ordercells_opus_S136 import order_cells, cell_gap_affines
from lrc_exact_rhostar_hardcores_opus_S137 import GP_intervals, meas

THETA = F(1, 7)

def U_pieces(E, theta=THETA):
    """Exact piecewise-affine representation of U(x) = sum (g_i(x)-theta)_+ on [0,1].
       Returns list of (lo, hi, c, b) with U = c*x + b on (lo, hi)."""
    pieces = []
    for a, b_ in order_cells(E):
        gaps = cell_gap_affines(E, a, b_)
        # refine by each gap's theta-crossing (activation boundaries)
        bp = set([a, b_])
        for (c, b0) in gaps:
            if c != 0:
                xc = (theta - b0) / c
                if a < xc < b_:
                    bp.add(xc)
        bp = sorted(bp)
        for u, v in zip(bp, bp[1:]):
            m2 = (u + v) / 2
            cs = F(0); bs = F(0)
            for (c, b0) in gaps:
                if c * m2 + b0 > theta:
                    cs += c; bs += b0 - theta
            pieces.append((u, v, cs, bs))
    return pieces

def integrate_affine(pieces, power=1, restrict=None):
    """Exact integral of U^power over [0,1] (power 1 or 2), optionally restricted to
       an interval list `restrict` (sorted disjoint)."""
    tot = F(0)
    for (u, v, c, b) in pieces:
        segs = [(u, v)]
        if restrict is not None:
            segs = []
            for (a, h) in restrict:
                lo = max(a, u); hi = min(h, v)
                if hi > lo: segs.append((lo, hi))
        for (lo, hi) in segs:
            if power == 1:
                # ∫ (c x + b) dx
                tot += c * (hi * hi - lo * lo) / 2 + b * (hi - lo)
            else:
                # ∫ (c x + b)^2 dx = c^2 x^3/3 + c b x^2 + b^2 x
                tot += (c * c * (hi**3 - lo**3) / 3 + c * b * (hi * hi - lo * lo)
                        + b * b * (hi - lo))
    return tot

def cross_moment(pieces, restrict=None):
    """Exact ∫ U(x)·U(1−x) dx (optionally restricted, restriction assumed symmetric)."""
    # build reflected pieces: U(1−x) = c*(1−x)+b = −c x + (c+b) on (1−v, 1−u)
    refl = [(1 - v, 1 - u, -c, c + b) for (u, v, c, b) in pieces]
    refl.sort()
    pieces_sorted = sorted(pieces)
    # common refinement by two-pointer over breakpoints
    bps = sorted(set([p[0] for p in pieces_sorted] + [p[1] for p in pieces_sorted]
                     + [p[0] for p in refl] + [p[1] for p in refl]))
    # index pieces by interval for lookup: walk both lists
    tot = F(0)
    i = j = 0
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        mid = (lo + hi) / 2
        while i < len(pieces_sorted) and pieces_sorted[i][1] <= mid: i += 1
        while j < len(refl) and refl[j][1] <= mid: j += 1
        if i >= len(pieces_sorted) or j >= len(refl): break
        (u1, v1, c1, b1) = pieces_sorted[i]
        (u2, v2, c2, b2) = refl[j]
        if not (u1 <= mid <= v1 and u2 <= mid <= v2): continue
        segs = [(lo, hi)]
        if restrict is not None:
            segs = []
            for (a, h) in restrict:
                l2 = max(a, lo); h2 = min(h, hi)
                if h2 > l2: segs.append((l2, h2))
        for (l2, h2) in segs:
            # ∫ (c1 x + b1)(c2 x + b2) dx
            A = c1 * c2; Bc = c1 * b2 + c2 * b1; C = b1 * b2
            tot += A * (h2**3 - l2**3) / 3 + Bc * (h2 * h2 - l2 * l2) / 2 + C * (h2 - l2)
    return tot

def report(E, name, P=None):
    t0 = time.time()
    pieces = U_pieces(E)
    restrict = None; gpm = None
    if P is not None:
        restrict = GP_intervals(P)
        gpm = meas(restrict)
    EU = integrate_affine(pieces, 1, restrict)
    EU2 = integrate_affine(pieces, 2, restrict)
    EUx = cross_moment(pieces, restrict)
    if EU2 > 0:
        PZ = EU * EU / EU2
        sym2 = (EU2 + EUx) / 2
        PZs = EU * EU / sym2 if sym2 > 0 else None
    else:
        PZ = PZs = None
    scope = f" | G_P-restricted, meas(G_P)={float(gpm):.4f}" if P is not None else ""
    print(f" {name:<24}{scope}")
    print(f"    E[U] = {str(EU)[:30]:>30} = {float(EU):.6f}")
    print(f"    E[U^2] = {float(EU2):.6f}   E[U(x)U(1-x)] = {float(EUx):.6f}"
          f"   (cross/second = {float(EUx/EU2) if EU2>0 else float('nan'):.3f})")
    if PZ is not None:
        gain = float(PZs / PZ) if PZs else float("nan")
        print(f"    PZ = {float(PZ):.6f}   PZ_sym = {float(PZs):.6f}   GAIN x{gain:.3f}"
              f"   [{time.time()-t0:.0f}s]")

def main():
    print("=" * 100)
    print("(1) EXACT U-moments + reversal-symmetrized PZ on the board (theta = 1/7)")
    print("    [PZ floors mu directly: mu >= P(U>0) >= E[U]^2/E[U_sym^2] = PZ_sym]")
    print("=" * 100)
    board = {
        "AP_13": list(range(1, 14)),
        "monad rec": sorted([2,4,6,8,10,12,14,16,18,20,22,11,13]),
        "2AP+13": sorted([2,4,6,8,10,12,14,16,18,20,22,24,13]),
        "S41 EU-min (mirror A)": [1,2,3,4,5,6,7,8,9,10,11,12,20],
        "S41 EU-min mirror B": [1,9,10,11,12,13,14,15,16,17,18,19,20],
    }
    for name, E in board.items():
        report(E, name)

    print()
    print("=" * 100)
    print("(2) THE TWO HARD CORES: G_P-restricted exact U-moments (the intersected PZ floor)")
    print("    rho*(P,E) >= P(U>0 on G_P) >= E_GP[U]^2/E_GP[U_sym^2]; bar m_P = 0.056487")
    print("=" * 100)
    CORE8 = [9,10,11,12,13]; CORE9 = [10,11,12,13]
    for name, (P, E) in {
        "CORE8 + AP co-offset": (CORE8, list(range(0, 8))),
        "CORE8 + CRT-block": (CORE8, [0, 15, 45, 62, 81, 100, 111, 118]),
        "CORE8 + two-cluster": (CORE8, [0,1,2,3,109,110,111,112]),
        "CORE9 + AP co-offset": (CORE9, list(range(0, 9))),
        "CORE9 + CRT-block": (CORE9, [0, 15, 45, 62, 81, 100, 111, 118, 125]),
    }.items():
        report(E, name, P)

if __name__ == "__main__":
    main()
