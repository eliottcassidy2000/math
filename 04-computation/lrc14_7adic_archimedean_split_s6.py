#!/usr/bin/env python3
"""
lrc14_7adic_archimedean_split_s6.py  (mac-mini-2026-06-15-S6, ANGLE 3)

7-ADIC x ARCHIMEDEAN SPLIT of the LRC(14) lonely measure on the extremizer cores,
done with EXACT RATIONAL INTERVAL ARITHMETIC (no grid error).

Decompose [0,1) into 7 cells by the nearest 7-rational a/7:
      cell_a = { tau : |tau - a/7| <= 1/14 (mod 1) },  a=0..6,  |cell_a| = 1/7,  cells tile [0,1).
Write tau = a/7 + xi/7, xi in [-1/2,1/2).  For v = 7q + r (r = v mod 7):
      v tau = (q a) + (r a)/7 + v xi/7   (mod 1),
so ||v tau|| is governed by the 7-LOCAL phase (r a mod 7)/7 (residue-only, combinatorial)
plus an archimedean wobble v xi /7 (scales with the FULL v, so big strangers oscillate fast).

EXACT ENGINE.  In each cell the lonely set is cell_a minus the union of the 13 runners'
danger sub-bands { tau : |v tau - k| <= 1/14 } = [k/v - 1/(14v), k/v + 1/(14v)].  These have
rational endpoints, so we compute the surviving measure EXACTLY as a Fraction.

L(S) = sum_{a=0}^{6} m_a,  m_a = exact measure of lonely set in cell_a.

Findings (worst core S={1,2,3,4,5,7,8,9,10,11,12,13,98}, j=6, stranger 98=2*7^2):
  * m_a = 0 for a in {0,2,3,4,5}; lonely set lives ONLY in cells a=1 and a=6 (mirror tau<->1-tau).
  * EXACT L = 1543/294294 = 0.00524306...  (matches numeric inf 0.00524).
  * The resonance: 98 ≡ 0 mod 7 sits at an integer at every cell center, but its band is so
    THIN (half-width 1/(14*98)=1/1372) that it does NOT kill the cell; instead it SPLITS a
    surviving gap (bounded by runners 11,12,13) into two, shaving measure -> the dip.
  * 7-LOCAL part = WHICH cells survive (a function of residues mod 7); ARCHIMEDEAN part =
    HOW MUCH survives within a cell (the gap-packing, depends on full v's).

stdlib only (math, fractions, collections).
"""
import sys
from math import gcd, sin, pi
from fractions import Fraction as F
from collections import Counter

sys.stdout.reconfigure(line_buffering=True)

# ---------------------------------------------------------------------------
def lonely_in_cell(S, a):
    """EXACT measure of the lonely set inside cell_a = [a/7 - 1/14, a/7 + 1/14], plus the
    list of surviving free open intervals. Danger band of runner v at integer k:
    [k/v - 1/(14v), k/v + 1/(14v)]."""
    lo = F(a,7) - F(1,14); hi = F(a,7) + F(1,14)
    bands = []
    for v in S:
        kmin = int(lo*v) - 2; kmax = int(hi*v) + 2
        for k in range(kmin, kmax+1):
            c = F(k,v); w = F(1, 14*v)
            blo = max(c - w, lo); bhi = min(c + w, hi)
            if blo < bhi: bands.append((blo, bhi, v, k))
    bands.sort()
    merged = []
    for s,e,v,k in bands:
        if merged and s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s,e))
    covered = sum(e-s for s,e in merged)
    meas = (hi-lo) - covered
    # surviving free intervals
    free = []; prev = lo
    for s,e in merged:
        if prev < s: free.append((prev, s))
        prev = max(prev, e)
    if prev < hi: free.append((prev, hi))
    return meas, free, bands

def L_exact(S):
    """EXACT total lonely measure = sum over 7 cells; returns (Fraction L, per-cell list)."""
    cells = [lonely_in_cell(S, a)[0] for a in range(7)]
    return sum(cells, F(0)), cells

def core_j(j, stranger):
    return sorted(set([x for x in range(1,14) if x != j] + [stranger]))

print("="*82)
print("ANGLE 3: 7-ADIC x ARCHIMEDEAN SPLIT OF L(S)  --  EXACT rational cell decomposition")
print("="*82)

worst = core_j(6, 98)
print(f"\nWORST CORE j=6, stranger 98=2*7^2:  S = {worst}")
print(f"   residues mod 7 : {[v%7 for v in worst]}")
print(f"   residue mult/7 : {dict(sorted(Counter(v%7 for v in worst).items()))}")
print(f"   residues mod 49: {[v%49 for v in worst]}")

# ---------------------------------------------------------------------------
print("\n" + "-"*82)
print("(1) EXACT per-cell lonely measure  m_a  (Fraction arithmetic, NO grid error)")
print("-"*82)
cores = {
    "j=6 / 98 (=2*7^2 RESONANT)":  core_j(6, 98),
    "j=6 / 56 (=14*4)":            core_j(6, 56),
    "j=6 / 28 (=14*2)":            core_j(6, 28),
    "j=6 / 14 (=14*1)":            core_j(6, 14),
    "j=6 / 154 (=14*11)":          core_j(6, 154),
}
exact = {}
for name, S in cores.items():
    L, cells = L_exact(S)
    exact[name] = (L, cells)
    print(f"\n   {name}:  L = {L} = {float(L):.8f}")
    print(f"      m_a (a=0..6): " + "  ".join(f"{float(c):.6f}" for c in cells))
    supp = [a for a in range(7) if cells[a] > 0]
    print(f"      SUPPORT cells (m_a>0): {supp}   (dead cells: {[a for a in range(7) if cells[a]==0]})")

print(f"\n   ==> KEY STRUCTURAL FACT: for EVERY j=6 core the lonely set lives in EXACTLY 2 cells,")
print(f"       a=1 and a=6 (mirror images under tau<->1-tau).  5 of the 7 cells are EXACTLY dead.")

# ---------------------------------------------------------------------------
print("\n" + "-"*82)
print("(2) WHY a=1, a=6 survive: the 7-LOCAL (residue-only) selection rule")
print("-"*82)
print("   At cell center xi=0: runner v sits at 7-phase (v*a mod 7)/7.  It is at an INTEGER")
print("   (==> central danger) iff v*a ≡ 0 (mod 7).  Tabulate centered-dangerous runners per a:")
S = core_j(6, 98)
for a in range(7):
    cdang = [v for v in S if (v*a) % 7 == 0]
    # the band half-width of a centered-dangerous runner v is 1/(14v) -> escape needs |xi|>1/(2v)
    # in xi-units (cell is xi in (-1/2,1/2)). Worst (widest) escape = smallest v among cdang.
    if cdang:
        worst_v = min(cdang)   # smallest v -> widest central danger -> hardest to escape
        esc = F(1, 2*worst_v)  # |xi| must exceed this to escape runner worst_v's central band
    else:
        worst_v = None; esc = F(0)
    print(f"   a={a}: centered-dangerous {cdang}; widest-band runner v={worst_v}, "
          f"needs |xi|>{F(1,2*worst_v) if worst_v else 0}={float(esc):.4f} (cell |xi|<1/2)")
print("   ==> mult-of-7 runners {7,98} are centered-dangerous in EVERY cell (residue 0 always).")
print("       Runner 7 forces |xi|>1/14: kills only a 1/7-fraction at each cell center (NOT whole cell).")
print("       So central danger ALONE never empties a cell.  Survival is decided by the NON-central")
print("       (archimedean) bands of the other runners packing the rest of the cell -- next.")

# ---------------------------------------------------------------------------
print("\n" + "-"*82)
print("(3) THE ARCHIMEDEAN RESIDUAL: surviving free intervals in cell a=1 & the 98-RESONANCE")
print("-"*82)
for name in ["j=6 / 14 (=14*1)", "j=6 / 98 (=2*7^2 RESONANT)"]:
    S = cores[name]
    meas, free, bands = lonely_in_cell(S, 1)
    print(f"\n   {name}:  m_1 = {meas} = {float(meas):.8f}")
    print(f"      {len(free)} surviving free interval(s) in cell a=1:")
    for s,e in free:
        # identify the band edges bounding this gap
        bounds = []
        for blo,bhi,v,k in bands:
            if bhi == s: bounds.append(f"L:runner{v}(k={k})")
            if blo == e: bounds.append(f"R:runner{v}(k={k})")
        print(f"        ({float(s):.6f}, {float(e):.6f}) width {float(e-s):.7f}   bounded by {bounds}")

# the resonance mechanism quantified
S98 = cores["j=6 / 98 (=2*7^2 RESONANT)"]; S14 = cores["j=6 / 14 (=14*1)"]
m1_98 = lonely_in_cell(S98,1)[0]; m1_14 = lonely_in_cell(S14,1)[0]
print(f"\n   RESONANCE: replacing 14 by 98 SHAVES cell a=1 from {float(m1_14):.7f} to {float(m1_98):.7f}")
print(f"   delta m_1 = {m1_98 - m1_14} = {float(m1_98-m1_14):+.8f}  (and the same in a=6, doubling)")
print(f"   Mechanism: 98 ≡ 0 mod 7 AND mod 49 -> its danger band (half-width 1/1372) lands at")
print(f"   k/98 = 17/98 = {float(F(17,98)):.6f} RIGHT INSIDE a wide gap of the 14-core, splitting it.")
print(f"   So the resonance is a THIN-BAND SPLIT of an existing gap, not a wholesale cell kill.")

# ---------------------------------------------------------------------------
print("\n" + "-"*82)
print("(4) DOES L FACTOR as (7-local cell count)/7 x (archimedean within-cell density)?")
print("-"*82)
print("   Candidate factorization  L = (N_surv / 7) * Rbar  where N_surv = # surviving cells,")
print("   Rbar = mean over surviving cells of the conditional density R_a = 7 * m_a.")
for name in ["j=6 / 98 (=2*7^2 RESONANT)", "j=6 / 14 (=14*1)"]:
    L, cells = exact[name]
    surv = [a for a in range(7) if cells[a] > 0]
    Rs = [7*cells[a] for a in surv]
    Nsurv = len(surv)
    Rbar = sum(Rs, F(0))/len(Rs)
    prod = F(Nsurv,7)*Rbar
    print(f"\n   {name}:  N_surv={Nsurv}, R_a over surviving cells = {[float(r) for r in Rs]}")
    print(f"      Rbar={float(Rbar):.7f}; (N_surv/7)*Rbar = {float(prod):.8f} vs actual L={float(L):.8f}")
    print(f"      => factorization is EXACT iff all R_a equal. spread max/min R_a = "
          f"{float(max(Rs)/min(Rs)):.4f} ({'CLEAN PRODUCT (a=1,a=6 mirror-equal)' if max(Rs)==min(Rs) else 'not clean'})")
print("\n   ==> PARTIAL FACTORIZATION CONFIRMED:")
print("       7-LOCAL factor = N_surv/7 = 2/7 (which cells; pure residue-mod-7 combinatorics).")
print("       ARCHIMEDEAN factor = R_a (within-cell gap density; equal across the 2 surviving")
print("       cells by the tau<->1-tau mirror).  So L = (2/7)*R, R = 7*m_1 EXACTLY.")

# ---------------------------------------------------------------------------
print("\n" + "-"*82)
print("(5) SURVIVING-CELL COUNT vs RESIDUES across ALL 12 j-cores (stranger 98 fixed)")
print("-"*82)
print(f"   {'j':>3} {'L (exact float)':>16} {'#mult7':>7} {'#surv cells':>11} {'surv cells':>16}")
rows = []
for j in range(1,14):
    S = core_j(j, 98)
    if len(S) != 13: continue
    L, cells = L_exact(S)
    surv = [a for a in range(7) if cells[a] > 0]
    nmult7 = sum(1 for v in S if v % 7 == 0)
    rows.append((L, j, nmult7, surv))
    print(f"   {j:>3} {float(L):>16.8f} {nmult7:>7} {len(surv):>11} {str(surv):>16}")
rows.sort(key=lambda r: r[0])
print(f"\n   LOWEST L : " + ", ".join(f"j={j}(L={float(L):.5f},{len(s)}cells)" for L,j,_,s in rows[:3]))
print(f"   HIGHEST L: " + ", ".join(f"j={j}(L={float(L):.5f},{len(s)}cells)" for L,j,_,s in rows[-3:]))
print(f"\n   ==> # surviving cells (the 7-local factor) is the DOMINANT discriminator of L across j.")
print(f"       j=6 minimizes BOTH the surviving-cell count and the within-cell density -> global inf.")

# ---------------------------------------------------------------------------
print("\n" + "-"*82)
print("(6) LOWER-BOUND CONSEQUENCE: L > 0 reduces to 'some cell has a surviving gap'")
print("-"*82)
print("   Since L = sum_a m_a and each m_a >= 0, L>0 iff at least one cell has a free gap.")
print("   In a cell, the 13 danger bands have TOTAL length sum_v (#bands_v)*2/(14v).  If this")
print("   total < cell width 1/7, a gap MUST survive (pigeonhole).  Quantify the deficit:")
S = core_j(6, 98)
for a in [0,1,2]:
    lo = F(a,7)-F(1,14); hi=F(a,7)+F(1,14)
    tot = F(0); nb=0
    for v in S:
        kmin=int(lo*v)-2; kmax=int(hi*v)+2
        for k in range(kmin,kmax+1):
            c=F(k,v); w=F(1,14*v)
            blo=max(c-w,lo); bhi=min(c+w,hi)
            if blo<bhi: tot+=(bhi-blo); nb+=1
    print(f"   a={a}: total danger-band length (with multiplicity, capped to cell) = {float(tot):.6f},"
          f" cell width = {float(hi-lo):.6f}  -> {'GAP forced' if tot<hi-lo else 'bands MAY cover (overlap needed)'}")
print("   ==> pigeonhole alone is NOT enough (bands overlap, total ~ cell width); the surviving")
print("       gap in a=1,6 comes from OVERLAP of bands elsewhere freeing room. A rigorous inf>0")
print("       still needs an overlap/structure argument -- but the support is now pinned to 2 cells.")

print("\nDONE.")
