# -*- coding: utf-8 -*-
# kps-2026-07-11-S127cont56: THE SHAPE OF THE TARGET AND ITS TWO DOORWAYS.
#
# Target = the lonely box L = [1/N, 1-1/N]^{N-1} on the torus T^{N-1}; LRC = the closed geodesic in direction
# v=(v_1..v_{N-1}) hits L. The two-bucket dispatch splits families by WHICH doorway (base) they clear at:
#   - TIGHT bucket (non-DC, THM-366): extremal AP {1..N-1} clears at base q=N.
#   - COVERING bucket (DC): extremal deep well {1..N-2,(N-1)N} clears at base q=Phi6(N)=N^2-N+1.
# CLAIM (new): the two doorways have OPPOSITE PARITY, and that flips the geometric obstruction:
#   - base N: parity of N. For EVEN N the 7 (=N/2) diameter pairs (j, j+N/2) sit at distance EXACTLY 1/2 --
#     the order-2 antipodal symmetry x->x+1/2 (klein-S269's tournament lever). For ODD N: no exact antipode,
#     the AP is the rotational tournament R_N (regularity).
#   - base Phi6(N)=N^2-N+1 is ALWAYS ODD (N(N-1) even) -> NO exact antipode -> the floor is NOT antipodal;
#     it is a THREE-GAP BAND-PACKING: the residues fill [mu, q-mu] tightly with the AP endpoints at the edges.
from math import gcd
from fractions import Fraction as F

def norm(x): r = x - int(x); r = r + 1 if r < 0 else r; return min(r, 1 - r)

def phi6(N): return N * N - N + 1

def deepwell_residues(N):
    # {1..N-2,(N-1)N} at t*=N/Phi6(N); residues (positions * q)
    q = phi6(N); p = N
    v = list(range(1, N - 1)) + [(N - 1) * N]
    res = sorted((vi * p) % q for vi in v)
    mu = (q + N - 1) // N  # ceil(q/N)
    return q, mu, res, v

def antipodal_pairs(res, q):
    # TRUE antipodes: r' - r == q/2 (mod q), i.e. distance EXACTLY 1/2. Needs q even.
    if q % 2 != 0: return []
    S = set(res); pairs = []
    for r in res:
        if (r + q // 2) in S and r < r + q // 2 <= q: pairs.append((r, r + q // 2))
    return pairs

def reflection_pairs(res, q):
    # reflections across 0: r + r' == q (r' = -r mod q). Symmetric across the origin, NOT distance 1/2.
    S = set(res); pairs = []
    for r in res:
        if r != 0 and (q - r) in S and r < q - r: pairs.append((r, q - r))
    return pairs

def main():
    print("(1) THE COVERING DOORWAY Phi6(N)=N^2-N+1 IS ALWAYS ODD (so no exact antipode / order-2 symmetry):")
    print(f"    {'N':>2} | N (tight base) parity | Phi6(N) (covering base) | parity")
    for N in range(3, 21):
        print(f"    {N:>2} | {'even' if N%2==0 else 'odd ':>19} | {phi6(N):>23} | {'EVEN?!' if phi6(N)%2==0 else 'odd'}")

    print("\n(2) THE COVERING EXTREMAL is a BAND-PACKING, NOT antipodal. Deep well {1..N-2,(N-1)N} (with observer 0):")
    print(f"    {'N':>2} | base q=Phi6 | band [mu,q-mu] | #TRUE-antipodal (dist=1/2) | #reflection-across-0")
    for N in [8, 10, 12, 14]:
        q, mu, res, v = deepwell_residues(N)
        pts = sorted(set(res) | {0})  # include the observer at 0
        ap = antipodal_pairs(pts, q); rf = reflection_pairs(pts, q)
        print(f"    {N:>2} | {q:>10} | [{mu},{q-mu}] | {len(ap):>26} | {len(rf)}  {rf}")
    print("    => q ODD => ZERO true antipodal pairs. The only symmetry is reflection-across-0 at the band edge")
    print("       (mu,q-mu)=the two AP-endpoints; the interior residues have NO partner => one-sided band-packing.")

    print("\n(3) THE TIGHT EXTREMAL AP {1..N-1} at t=1/N (positions {0..N-1}/N, with observer 0):")
    print(f"    {'N':>2} | parity | #TRUE-antipodal (dist=1/2) pairs (j, j+N/2) | pairs")
    for N in [7, 8, 13, 14]:
        q = N; pts = list(range(0, N))  # all N grid points {0,..,N-1} (observer + 13 runners)
        ap = antipodal_pairs(pts, q)
        print(f"    {N:>2} | {'even' if N%2==0 else 'odd':>6} | {len(ap):>3} | {ap}")

    print("\n(4) THE PARITY DUALITY, one line per bucket at N=14:")
    print("    TIGHT   bucket: base 14 EVEN -> 7 diameter pairs tied at 1/2 -> order-2 antipodal obstruction")
    print("                    (klein-S269: tournament has odd |Aut|, cannot carry x->x+1/2; M>=1/14).")
    print("    COVERING bucket: base 183 ODD  -> no diameter pair -> obstruction is three-gap BAND-PACKING")
    print("                    (residues 14*{1..12,-1} fill [14,169]; AP endpoints 14,169 at the edges).")
    print("\n=> THE TARGET BOX L HAS TWO DOORWAYS OF OPPOSITE PARITY. Even doorway (tight) = antipodal/Borsuk-Ulam")
    print("   flavor; odd doorway (covering) = packing/three-gap flavor. A proof needs BOTH tools, one per bucket.")

if __name__ == "__main__":
    main()
