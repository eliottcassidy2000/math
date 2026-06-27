"""
lrc14_pairwise_coemptiness_spectrum_kps.py  (kind-pasteur-2026-06-27-S31ag)

CRUX 1 support: the eigenstructure of the 6x6 pairwise sector co-emptiness
(Gram) matrix M[a][b] = meas{tau: inner sectors a AND b both empty},  a,b in {1..6},
with diagonal M[a][a] = meas{sector a empty}.  S2 = sum_{a<b} M[a][b].

mac-mini-S60 claims M is "reflection-symmetric + dominant-Perron, NOT 4I+2J,
NOT circulant" => a 3x3 reflection-half-block Perron bound is the CRUX-1 route.
The four-faces reflection's suggested route instead invokes the Clebsch/biplane
4I+2J SRG eigenstructure.  These are in tension.  This script computes M exactly
for consec_k and competitors, reports eigenvalues, tests reflection symmetry
(sector a <-> 6-a) and circulant/4I+2J structure, to settle which route is viable.
"""
import sys, itertools
from fractions import Fraction as F
import numpy as np

def sector_of(p):  # p a Fraction in [0,1); returns 0..6
    return int((p % 1) * 7)

def empty_indicator_measure_matrix(E):
    """Return exact 6x6 M[a][b]=meas{a,b both empty}, a,b in 1..6 (index 0..5)."""
    E = sorted(set(E)); b = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for m in range(0, 7 * e + 1):
            b.add(F(m, 7 * e))
    b = sorted(b)
    M = [[F(0)] * 6 for _ in range(6)]
    for i in range(len(b) - 1):
        x0, x1 = b[i], b[i + 1]
        if x1 <= x0: continue
        w = x1 - x0
        covered = set(sector_of(e * ((x0 + x1) / 2)) for e in E)
        empty = [s for s in range(1, 7) if s not in covered]  # inner empty sectors
        for a in empty:
            for c in empty:
                M[a - 1][c - 1] += w
    return M

def analyze(name, E):
    M = empty_indicator_measure_matrix(E)
    Mf = np.array([[float(x) for x in row] for row in M])
    evals = np.sort(np.linalg.eigvalsh(Mf))[::-1]
    diag = [M[i][i] for i in range(6)]
    S2 = sum(M[i][j] for i in range(6) for j in range(i + 1, 6))
    # reflection symmetry test: sector a <-> 6-a  (1<->5, 2<->4, 3<->3, 6<->0[anchor])
    # use involution on {1..6}: r(a) = 6-a for a in 1..5, r(6)=6 (6 has no inner partner)
    def refl(a):  # 0-indexed sector (0..5 == sectors 1..6)
        s = a + 1
        rs = 6 - s
        return (rs - 1) if 1 <= rs <= 6 else None
    refl_ok = True
    for a in range(6):
        ra = refl(a)
        if ra is None:
            continue
        for c in range(6):
            rc = refl(c)
            if rc is None:
                continue
            if M[a][c] != M[ra][rc]:
                refl_ok = False
    # circulant test (M[a][b] depends only on (a-b) mod 6)
    circ = all(M[a][c] == M[0][(c - a) % 6] for a in range(6) for c in range(6))
    # 4I+2J-proportional test: off-diagonal all equal AND diagonal all equal?
    offvals = set(M[a][c] for a in range(6) for c in range(6) if a != c)
    diagvals = set(diag)
    print(f"\n=== {name}  E={tuple(E)} ===")
    print("  M (6x6, inner sectors 1..6), entries x1000:")
    for a in range(6):
        print("   " + " ".join(f"{float(M[a][c])*1000:6.1f}" for c in range(6)))
    print(f"  diag(P[sector a empty]) = {[round(float(d),4) for d in diag]}")
    print(f"  S2 = sum_pairs = {float(S2):.5f}   trace = {float(sum(diag)):.5f}")
    print(f"  eigenvalues (desc): {[round(float(e),5) for e in evals]}")
    print(f"  Perron gap lambda1-lambda2 = {float(evals[0]-evals[1]):.5f}")
    print(f"  reflection-symmetric (a<->6-a): {refl_ok}")
    print(f"  circulant (depends on a-b mod 6): {circ}")
    print(f"  #distinct off-diag = {len(offvals)}, #distinct diag = {len(diagvals)} "
          f"(4I+2J would need 1 and 1)")
    return evals

if __name__ == "__main__":
    sys.stdout.reconfigure(line_buffering=True)
    print("PAIRWISE SECTOR CO-EMPTINESS MATRIX M — eigenstructure & symmetry")
    print("="*70)
    for k in (8, 9, 11):
        analyze(f"consec_{k}", tuple(range(k)))
    analyze("single-far(8+21)", tuple(range(7)) + (21,))
    analyze("even-AP_8", tuple(2*i for i in range(8)))
    analyze("wide-rand_8", (0,1,3,7,12,20,33,50))
    print("\nNOTE: J=all-ones; 4I+2J SRG eigenstructure needs M = aI+bJ form")
    print("(all off-diag equal, all diag equal). Reflection-Perron route needs only")
    print("the reflection symmetry + a dominant Perron eigenvalue.")
