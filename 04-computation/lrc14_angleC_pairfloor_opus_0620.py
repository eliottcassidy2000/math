#!/usr/bin/env python3
"""
lrc14_angleC_pairfloor_opus_0620.py   (opus-2026-06-20)

Angle C: the j-HOLE MONOTONICITY FLOOR.

RIGOROUS FACT (CA monotonicity of cover under clock removal, fully PROVED):
   If H' subset H (fewer holes), then measS7(full \ H) <= measS7(full \ H').
   In particular for ANY size-j subset J of H:
        measS7(E) = measS7(full \ H) <= measS7(full \ J).
   Hence  measS7(E) <= min_{J subset H, |J|=j} measS7(full \ J).
   The adversary (choosing H to maximise measS7(E)) is forced to leave SOME J of size
   j inside H; the safe bound is the value the adversary cannot avoid:
        bound_j(N, nholes) = max over H of [ min over J subset H, |J|=j of measS7(full \ J) ]
   Since the adversary wants this max large, it chooses H so that EVERY size-j subset J
   has large measS7(full\J).  The worst (largest) such value over all H with |H|=nholes is
        bound_j(N,nholes) = max_{H, |H|=nholes} min_{J in C(H,j)} measS7(full \ J).
   This equals: the (nholes choose j choose ...) -- compute directly.

   For j = nholes (use ALL holes as the single subset): bound = measS7(E) itself -> exact
   (trivial, == enumeration). For j < nholes it is a genuine LOWER-complexity floor:
   we only ever evaluate measS7 on full-boards-minus-j-clocks (C(N, j) of them), NOT on
   all C(N,nholes) shapes. The certificate closes if bound_j(N,nholes) <= cap_k.

KEY POINT vs single-hole (j=1): j=1 was too weak (closed only 1-2 holes).
j=2 uses PAIR criticality, which is supermodular-strong. Test how large j must be.

This is the honest CA "defect-cluster" floor: removing a j-cluster of defects already
saturates the coverage loss. If a SMALL fixed j (independent of nholes) closes all
residual spans, that is a genuinely cheaper certificate than full enumeration.

EXACT Fractions; stdlib only.
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

def measS7(E):
    E = sorted(set(E)); Enz = [e for e in E if e != 0]
    bps = set([F(0), F(1)])
    for e in Enz:
        for m in range(0, 7*e+1):
            bps.add(F(m, 7*e))
    bps = sorted(bps); total = F(0)
    for i in range(len(bps)-1):
        x0, x1 = bps[i], bps[i+1]
        if x1 <= x0: continue
        xm = (x0+x1)/2
        res = set(int(7*e*xm) % 7 for e in E)
        if len(res) == 7: total += x1 - x0
    return total

CAP = {8: F(2243,5880), 9: F(1979,4004), 10: F(55,91)}

def main():
    print("="*96)
    print("ANGLE C j-HOLE MONOTONICITY FLOOR  (opus-2026-06-20)")
    print("measS7(E) <= min_{J subset H,|J|=j} measS7(full\\J).  Adversary max over H => bound_j.")
    print("Certificate closes if bound_j(N,nholes) <= cap_k. Find smallest sufficient j.")
    print("="*96)
    regimes = {8: range(8, 14), 9: range(9, 14), 10: range(10, 14)}

    for k in [8, 9, 10]:
        ck = CAP[k]
        print(f"\n--- k={k}, cap_k={str(ck)}={float(ck):.5f} ---")
        for N in regimes[k]:
            nholes = (N+1) - k
            if nholes < 1: continue
            full = tuple(range(N+1)); mfull = measS7(full)
            hole_positions = list(range(1, N))   # 0 pinned, N anchors span
            # precompute measS7(full \ J) for all J of size up to nholes? expensive.
            # We need, for each j: bound_j = max_{H} min_{J in C(H,j)} measS7(full\J).
            # Precompute g(J)=measS7(full\J) for |J|=j (only need up to the j we test).
            results = {}
            for j in range(1, nholes+1):
                # g over all j-subsets of hole_positions
                gJ = {}
                for J in itertools.combinations(hole_positions, j):
                    gJ[J] = measS7(tuple(e for e in full if e not in set(J)))
                # bound_j = max over H (size nholes) of min over J in C(H,j) of gJ[J]
                bound = F(-1)
                for H in itertools.combinations(hole_positions, nholes):
                    mn = min(gJ[J] for J in itertools.combinations(H, j))
                    if mn > bound: bound = mn
                results[j] = bound
                if bound <= ck:
                    break
            # smallest j that closes
            jclose = next((j for j in sorted(results) if results[j] <= ck), None)
            line = "  ".join(f"j={j}:{float(results[j]):.4f}{'<=cap' if results[j]<=ck else ''}" for j in sorted(results))
            print(f"  N={N:2d} holes={nholes}: {line}")
            if jclose is not None:
                print(f"        --> smallest j closing certificate: j={jclose} (bound {float(results[jclose]):.5f} <= cap {float(ck):.5f})")
            else:
                print(f"        --> NO j<nholes closes; need full enumeration (j=nholes is exact).")
    print("\nDONE.")

if __name__ == "__main__":
    main()
