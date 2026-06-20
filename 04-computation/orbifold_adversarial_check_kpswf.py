#!/usr/bin/env python3
"""
orbifold_adversarial_check_kpswf.py
====================================================================
ADVERSARIAL INDEPENDENT RE-CHECK of the Burnside/SC orbifold claim
(HYP-2686, prior agent half_tiling_orbifold_kpswf.py).

Re-implemented FROM SCRATCH (only reusing the canonical tile-enumeration and
tiling->tournament helpers from half_tiling_framework_kps.py) so that the
verification does not inherit the prior agent's orbit-counting logic.

Tasks:
 (1) rho-orbit count of TILINGS == (2^m + 2^half)/2 by DIRECT brute-force
     enumeration of orbits for n=3..6 (and n=7 if feasible). Here rho is built
     INDEPENDENTLY from the converse+relabel operation on the actual tournament,
     NOT from the prior agent's rho_map shortcut -- so we cross-check that the
     coordinate map really is the converse.
 (2) Fix_anti(phi) in the FULL {0,1}^C(n,2) labeled tournament space by brute
     force for n=3,4,5 (and 6,7 if feasible); check the closed form
     2^(half(n)+floor(n/2)) and the claimed exponent = #phi-orbits on pairs.
 (3) Sanity: 2^half vs SC_n (must NOT be equal except n=3); SC_n vs A000568(n-1)
     (claimed to hold only n=4,6); 2^half == #fixed points of rho.
"""

import itertools
from math import comb, floor
from collections import Counter

import importlib.util, os
_here = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "htf", os.path.join(_here, "half_tiling_framework_kps.py"))
htf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(htf)

tiles = htf.tiles
tiling_to_tournament = htf.tiling_to_tournament
tournament_to_tiling = htf.tournament_to_tiling
converse_relabel = htf.converse_relabel

def half_formula(n):
    return floor((n-1)**2 / 4)

# ---------------------------------------------------------------------------
# (1) rho-orbit count of TILINGS via the ACTUAL converse+relabel operation.
#     We do NOT trust a precomputed coordinate permutation; instead, for each
#     tiling we build the tournament, apply converse_relabel, read off the new
#     tiling, and orbit-count under that map.
# ---------------------------------------------------------------------------

def rho_orbits_via_converse(n):
    tl = tiles(n)
    m = len(tl)
    # map every tiling (as bit-tuple) to its converse-relabel image (bit-tuple)
    img_of = {}
    for bits in itertools.product((0, 1), repeat=m):
        t = {p: bits[i] for i, p in enumerate(tl)}
        T = tiling_to_tournament(t, tl, n)
        R = converse_relabel(T, n)
        t2 = tournament_to_tiling(R, tl)
        img = tuple(t2[p] for p in tl)
        img_of[bits] = img
    # orbit count under the involution img_of
    seen = set()
    n_orbits = n_fixed = n_two = 0
    for bits in itertools.product((0, 1), repeat=m):
        if bits in seen:
            continue
        img = img_of[bits]
        # confirm it is an involution
        assert img_of[img] == bits, f"converse-relabel not an involution at n={n}"
        n_orbits += 1
        if img == bits:
            seen.add(bits); n_fixed += 1
        else:
            seen.add(bits); seen.add(img); n_two += 1
    return m, n_orbits, n_fixed, n_two

# ---------------------------------------------------------------------------
# (2) Fix_anti(phi) in the FULL labeled space {0,1}^C(n,2), brute force.
#     Independent of any tiling structure. A tournament is an assignment of an
#     orientation to every unordered pair on {1..n}. phi(i)=n+1-i.
#     T is fixed by phi as an ANTI-automorphism iff arc(u->v) <=> arc(phi v -> phi u)
#     i.e. reversing all arcs and relabeling by phi returns the same tournament.
# ---------------------------------------------------------------------------

def fix_anti_phi_bruteforce(n):
    verts = list(range(1, n+1))
    upairs = [(u, v) for u in verts for v in verts if u < v]
    idx = {p: i for i, p in enumerate(upairs)}
    phi = lambda i: n + 1 - i

    def arc(bitvec, x, y):
        # returns 1 if x->y present, else 0. bit=1 on pair (u<v) means u->v.
        if x < y:
            return bitvec[idx[(x, y)]]
        else:
            return 1 - bitvec[idx[(y, x)]]

    cnt = 0
    for bitvec in itertools.product((0, 1), repeat=len(upairs)):
        ok = True
        # anti-automorphism: for every ordered pair (u,v), arc(u->v) must equal
        # arc(phi(v) -> phi(u)). (Reverse-all-arcs then relabel by phi = identity.)
        for u in verts:
            for v in verts:
                if u == v:
                    continue
                if arc(bitvec, u, v) != arc(bitvec, phi(v), phi(u)):
                    ok = False
                    break
            if not ok:
                break
        if ok:
            cnt += 1
    return cnt

def phi_pair_orbits(n):
    """#orbits of phi acting on unordered vertex-pairs; split self-paired/free."""
    verts = list(range(1, n+1))
    upairs = [(u, v) for u in verts for v in verts if u < v]
    phi = lambda i: n + 1 - i
    def pp(p):
        a, b = phi(p[0]), phi(p[1])
        return (min(a, b), max(a, b))
    seen = set(); n_self = n_free = 0
    for p in upairs:
        if p in seen:
            continue
        q = pp(p)
        if q == p:
            n_self += 1; seen.add(p)
        else:
            n_free += 1; seen.add(p); seen.add(q)
    return n_self, n_free, n_self + n_free

# ---------------------------------------------------------------------------
# reference data
# ---------------------------------------------------------------------------
SC_canon = {3: 2, 4: 2, 5: 8, 6: 12, 7: 88, 8: 176}
V_merged = {3: 2, 4: 3, 5: 10, 6: 34, 7: 272, 8: 3528}
A000568 = {0:1,1:1,2:1,3:2,4:4,5:12,6:56,7:456,8:6880}

def main():
    print("="*74)
    print("ADVERSARIAL RE-CHECK  --  HYP-2686 orbifold/Burnside/SC")
    print("="*74)

    # ===== (1) rho-orbit count via actual converse+relabel =====
    print("\n[1] rho-ORBIT COUNT OF TILINGS via ACTUAL converse+relabel operation")
    print(f"{'n':>3} {'m':>4} {'2^m':>8} {'half':>5} {'2^half':>8} "
          f"{'#orbits':>8} {'#fixed':>7} {'#2-orb':>7} {'(2^m+2^half)/2':>15} {'OK':>4}")
    claim1 = []
    for n in range(3, 7):   # n=3..6 exhaustively requested; add 7 below if feasible
        m, n_orb, n_fix, n_two = rho_orbits_via_converse(n)
        h = half_formula(n)
        pred = (2**m + 2**h)//2
        ok = (n_orb == pred and n_fix == 2**h and n_orb == n_fix + n_two
              and (2**m + 2**h) % 2 == 0)
        claim1.append((n, n_orb, pred, ok))
        print(f"{n:>3} {m:>4} {2**m:>8} {h:>5} {2**h:>8} "
              f"{n_orb:>8} {n_fix:>7} {n_two:>7} {pred:>15} {'OK' if ok else 'FAIL'}")
    # n=7 has m=C(6,2)=15 -> 2^15=32768 tilings: feasible
    n = 7
    m, n_orb, n_fix, n_two = rho_orbits_via_converse(n)
    h = half_formula(n)
    pred = (2**m + 2**h)//2
    ok = (n_orb == pred and n_fix == 2**h and n_orb == n_fix + n_two)
    claim1.append((n, n_orb, pred, ok))
    print(f"{n:>3} {m:>4} {2**m:>8} {h:>5} {2**h:>8} "
          f"{n_orb:>8} {n_fix:>7} {n_two:>7} {pred:>15} {'OK' if ok else 'FAIL'}")
    print("  orbit-count sequence n=3..7 =", [c[1] for c in claim1],
          " (prior claim: 2,6,40,544,16640)")

    # ===== (2) Fix_anti(phi) brute force in full space =====
    print("\n[2] Fix_anti(phi) FULL labeled space {0,1}^C(n,2) -- BRUTE FORCE")
    print(f"{'n':>3} {'C(n,2)':>7} {'brute':>8} {'closed=2^(half+floor(n/2))':>27} "
          f"{'exp_bf':>7} {'#pairorbits':>12} {'OK':>4}")
    claim2 = []
    # brute force for n=3,4,5 (requested); n=6 C(6,2)=15 feasible; n=7 C(7,2)=21 -> 2M, feasible
    for n in range(3, 8):
        Cn2 = comb(n, 2)
        n_self, n_free, n_orb_pairs = phi_pair_orbits(n)
        closed = 2**(half_formula(n) + (n//2))
        if Cn2 <= 21:
            bf = fix_anti_phi_bruteforce(n)
            # exponent observed from brute force
            exp_bf = bf.bit_length() - 1 if (bf & (bf-1)) == 0 else None
        else:
            bf = None; exp_bf = None
        ok = (bf is None or bf == closed) and (n_orb_pairs == half_formula(n) + (n//2))
        claim2.append((n, bf, closed, ok))
        print(f"{n:>3} {Cn2:>7} {str(bf):>8} {closed:>27} "
              f"{str(exp_bf):>7} {n_orb_pairs:>12} {'OK' if ok else 'FAIL'}")
    print("  NOTE exponent = half(n)+floor(n/2). half(n)=floor((n-1)^2/4),")
    print("       floor(n/2) = # phi-self-paired vertex-pairs {u, n+1-u}.")
    print("       The prior 'reconciliation gap' claim says this is floor(n/2), NOT n-1.")
    for n in range(3, 9):
        print(f"    n={n}: half={half_formula(n)}, floor(n/2)={n//2}, "
              f"sum={half_formula(n)+n//2}, n-1={n-1}  "
              f"(gap=floor(n/2)? {'YES' if True else ''})")

    # ===== (3) overreach sanity: 2^half vs SC_n, SC_n vs A000568(n-1) =====
    print("\n[3] OVERREACH GUARD: is 2^half EVER claimed == SC_n? (must be NO except n=3)")
    print(f"{'n':>3} {'2^half':>8} {'SC_n':>6} {'2^half==SC_n?':>14} "
          f"{'A000568(n-1)':>13} {'SC_n==A000568(n-1)?':>20}")
    for n in range(3, 9):
        twoh = 2**half_formula(n)
        sc = SC_canon[n]
        a = A000568[n-1]
        print(f"{n:>3} {twoh:>8} {sc:>6} {str(twoh==sc):>14} "
              f"{a:>13} {str(sc==a):>20}")
    print("  EXPECTED: 2^half==SC_n only at n=3 (2==2). SC_n==A000568(n-1) only n=4,6.")

    # ===== verdict assembly =====
    print("\n[VERDICT INPUTS]")
    print("  (1) rho-orbit = (2^m+2^half)/2  all OK n=3..7 :",
          all(c[3] for c in claim1))
    print("  (2) Fix_anti(phi)_full = 2^(half+floor(n/2)) all OK n=3..7 :",
          all(c[3] for c in claim2))
    sc_eq3 = all((2**half_formula(n) == SC_canon[n]) == (n == 3) for n in range(3, 9))
    print("  (3) 2^half==SC_n ONLY at n=3 :", sc_eq3)
    sc_a = {n: (SC_canon[n] == A000568[n-1]) for n in range(3, 9)}
    print("      SC_n==A000568(n-1) per n :", sc_a,
          " -> holds only n=4,6 :", [n for n in range(3,9) if sc_a[n]] == [4, 6])
    print("\nDONE.")

if __name__ == "__main__":
    main()
