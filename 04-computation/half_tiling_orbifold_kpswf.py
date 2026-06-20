#!/usr/bin/env python3
"""
half_tiling_orbifold_kpswf.py
================================================================
THREAD 2 -- HYP-2686: the orbifold Q_m / rho and the Burnside/SC connection.
(kind-pasteur half-tiling work-fork; reuses half_tiling_framework_kps.py helpers)

rho = converse+relabel reflection on tilings (THM-280): a PURE coordinate
involution rho(a,b)=(n+1-b, n+1-a) on the m=C(n-1,2) tiles. It acts on the
2^m TILINGS with exactly 2^half(n) fixed points, half(n)=floor((n-1)^2/4).

This script answers:
 (a) merged tiling count = #rho-orbits = (2^m + 2^half)/2 for n=3..8;
     compare to V_merged (merged metagraph node count) = 2,3,10,34,272,3528.
 (b) Fix_anti(phi) in the FULL tournament space {0,1}^{C(n,2)}: count labeled
     tournaments T with T(phi u, phi v) = T(v,u). Determine it as 2^(something)
     and relate to 2^half(n). Reconcile fixed-HP-frame vs full-frame.
 (c) Characterize the 2:1 fold (orbifold): free tilings pair, grid-sym tilings
     are the branch locus. Verify (2^m+2^half)/2 = #orbits by direct enumeration
     for n=3..6 (and n=7 where feasible).
 (d) Does half-tiling give an explicit SC-halving bijection |SC(n)| <-> tournaments
     on (n-1)/2 vertices?  Compare 2^half vs SC_n; report bijection vs Burnside only.
"""

import itertools
from math import comb, floor, gcd
from collections import Counter, defaultdict
from functools import reduce

# reuse canonical helpers
import importlib.util, os
_here = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "htf", os.path.join(_here, "half_tiling_framework_kps.py"))
htf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(htf)

tiles = htf.tiles
rho = htf.rho
involution_structure = htf.involution_structure
half_formula = htf.half_formula

# ----------------------------------------------------------------------
# (a) merged tiling count via Burnside on the tiling cube Q_m under <rho>
# ----------------------------------------------------------------------

def merged_tiling_count(n):
    tl, m, d, fixed, pairs, half = involution_structure(n)
    fix_id = 2**m                 # identity fixes everything
    fix_rho = 2**half             # rho fixes 2^half tilings (grid-sym)
    orbits = (fix_id + fix_rho) // 2
    assert (fix_id + fix_rho) % 2 == 0
    return m, half, fix_id, fix_rho, orbits

def enumerate_rho_orbits(n):
    """Direct enumeration of rho-orbits of the 2^m tilings.
       Returns (#orbits, #fixed_pts, #2-orbits)."""
    tl = tiles(n)
    m = len(tl)
    rho_map = {p: rho(*p, n) for p in tl}
    # represent a tiling as a tuple of bits in tile order
    idx = {p: i for i, p in enumerate(tl)}
    perm = [idx[rho_map[tl[i]]] for i in range(m)]   # bit i moves to position perm[i]
    seen = set()
    n_orbits = n_fixed = n_two = 0
    for bits in itertools.product((0, 1), repeat=m):
        if bits in seen:
            continue
        img = tuple(bits[perm[j]] for j in range(m))  # apply rho: new[j]=old[perm[j]]
        n_orbits += 1
        if img == bits:
            seen.add(bits)
            n_fixed += 1
        else:
            seen.add(bits)
            seen.add(img)
            n_two += 1
    return n_orbits, n_fixed, n_two

# ----------------------------------------------------------------------
# (b) Fix_anti(phi) in the FULL labeled tournament space {0,1}^{C(n,2)}
#     phi(i) = n+1-i.  T is fixed-anti iff for every ordered pair (u,v),
#     T(phi u, phi v) = T(v,u). Equivalently the arc direction on pair
#     {u,v} is reversed when transported to pair {phi u, phi v}.
# ----------------------------------------------------------------------

def fix_anti_phi_full(n):
    """Count labeled tournaments T on vertices 1..n with phi an anti-automorphism,
       phi(i)=n+1-i.  Exhaustive for small n; closed form for all n.
       Returns (count_exhaustive_or_None, closed_form, exponent, details)."""
    verts = list(range(1, n+1))
    pairs = [(u, v) for u in verts for v in verts if u < v]   # unordered pairs
    phi = lambda i: n+1-i
    # induced action on unordered pairs: {u,v} -> {phi u, phi v}
    def ppair(p):
        u, v = p
        a, b = phi(u), phi(v)
        return (min(a, b), max(a, b))
    pairset = set(pairs)
    # orbit structure of phi on pairs (orbits of size 1 = self-paired, size 2 = free)
    seen = set()
    self_paired = []     # pairs {u,v} with {phi u, phi v} = {u,v}
    free_orbits = []     # 2-cycles of pairs
    for p in pairs:
        if p in seen:
            continue
        q = ppair(p)
        if q == p:
            self_paired.append(p); seen.add(p)
        else:
            free_orbits.append((p, q)); seen.add(p); seen.add(q)

    # The anti-automorphism constraint:
    #   On a FREE orbit {p, q}: arc on p reversed -> arc on q. The relation
    #   T(q) is DETERMINED by T(p) (and the orientation bookkeeping below).
    #   So a free orbit contributes exactly 1 free bit (choose T on p; q forced).
    #   => 2^(#free_orbits) choices from free orbits, each consistent.
    #
    #   On a SELF-PAIRED pair {u,v} with {phi u, phi v} = {u,v}:
    #     Either (phi u, phi v) = (u, v)  [phi fixes both vertices] -- impossible
    #         since phi has at most one fixed vertex (the center, odd n), and a
    #         pair needs two distinct vertices, so phi can't fix BOTH.
    #     Or (phi u, phi v) = (v, u)  [phi swaps the two vertices, i.e. v = phi u].
    #       Then the constraint T(phi u, phi v) = T(v,u) reads T(v,u) = T(v,u):
    #       AUTOMATICALLY satisfied for either orientation => this pair is FREE
    #       (1 free bit), NOT constrained.
    #   Hence EVERY orbit (self-paired or free) contributes exactly 1 free bit,
    #   and free 2-orbits also contribute 1 bit each (one determines the other).
    #   => total free bits = #orbits of phi on pairs = #self_paired + #free_orbits.
    n_self = len(self_paired)
    n_free = len(free_orbits)
    exponent = n_self + n_free          # = number of phi-orbits on pairs

    closed = 2**exponent

    # exhaustive check (only when C(n,2) small enough)
    Cn2 = comb(n, 2)
    exhaustive = None
    if Cn2 <= 21:                      # up to n=7: 2^21 = 2.1M, feasible
        pair_index = {p: i for i, p in enumerate(pairs)}
        # For each unordered pair store bit b: b=1 means u->v (u<v), b=0 means v->u.
        def is_fixed_anti(bitvec):
            # T(x,y)=1 iff arc x->y present.
            # build directed lookup
            def arc(x, y):
                if x < y:
                    return bitvec[pair_index[(x, y)]]           # 1 => x->y
                else:
                    return 1 - bitvec[pair_index[(y, x)]]       # x>y: x->y iff NOT(y->x)
            for (u, v) in [(a, b) for a in verts for b in verts if a != b]:
                # constraint: T(phi u, phi v) == T(v, u)
                if arc(phi(u), phi(v)) != arc(v, u):
                    return False
            return True
        cnt = 0
        for bitvec in itertools.product((0, 1), repeat=Cn2):
            if is_fixed_anti(bitvec):
                cnt += 1
        exhaustive = cnt

    details = dict(Cn2=Cn2, n_self_paired=n_self, n_free_orbits=n_free,
                   self_paired=self_paired, n_fixed_vertices=(1 if n % 2 == 1 else 0))
    return exhaustive, closed, exponent, details

def half_value(n):
    return half_formula(n)

# ----------------------------------------------------------------------
# (d) SC_n vs 2^half : is there a clean bijection?
# ----------------------------------------------------------------------

# Canon SC_n (self-complementary tournament ISO classes), n=3..8:
SC_canon = {3: 2, 4: 2, 5: 8, 6: 12, 7: 88, 8: 176}
# Canon V_merged (merged metagraph node count), n=3..8:
V_merged = {3: 2, 4: 3, 5: 10, 6: 34, 7: 272, 8: 3528}

# A000568 (# tournaments up to iso) and A002785/relevant references for context.

def sc_relation_table():
    rows = []
    for n in range(3, 9):
        h = half_value(n)
        twoh = 2**h
        sc = SC_canon[n]
        # candidate "tournaments on (n-1)/2 vertices": only meaningful for odd n
        # number of LABELED tournaments on k vertices = 2^C(k,2)
        if n % 2 == 1:
            k = (n-1)//2
            lab_k = 2**comb(k, 2)        # labeled tournaments on k vertices
            iso_k = None                  # iso classes A000568(k): 1,1,2,4,12,56,456...
        else:
            k = None; lab_k = None
        rows.append(dict(n=n, half=h, two_half=twoh, SC=sc, V_merged=V_merged[n],
                         k=k, labeled_tour_k=lab_k,
                         ratio_2half_over_SC=(twoh/sc)))
    return rows

# A000568: number of tournaments on k nodes (iso classes)
A000568 = {0:1,1:1,2:1,3:2,4:4,5:12,6:56,7:456,8:6880}
# A002785: self-converse tournament iso classes; SC_canon(n) = A002785(n+1).

# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------

def main():
    print("="*74)
    print("THREAD 2  HYP-2686  --  Orbifold Q_m/rho and Burnside/SC connection")
    print("="*74)

    # ============ (a) merged tiling count ============
    print("\n[a] MERGED TILING COUNT  = #rho-orbits = (2^m + 2^half)/2")
    print(f"{'n':>3} {'m':>4} {'half':>5} {'2^m':>10} {'2^half':>8} "
          f"{'(2^m+2^half)/2':>16} {'V_merged':>9}")
    merged_seq = []
    for n in range(3, 9):
        m, half, fix_id, fix_rho, orbits = merged_tiling_count(n)
        merged_seq.append(orbits)
        vm = V_merged.get(n, '-')
        print(f"{n:>3} {m:>4} {half:>5} {fix_id:>10} {fix_rho:>8} {orbits:>16} {str(vm):>9}")
    print("  merged-tiling-count(n=3..8) =", merged_seq)
    print("  V_merged (metagraph)        =", [V_merged[n] for n in range(3,9)])
    print("  INTERPRETATION: the merged-tiling count is the fixed-HP-FRAME analogue")
    print("  of V_merged. Both are Burnside orbit counts under the converse/complement")
    print("  Z_2, but on DIFFERENT spaces: tilings (labeled, base-path fixed) vs")
    print("  iso classes (S_n-quotient). They are NOT equal: tiling count >> V_merged,")
    print("  because the tiling cube has NOT been quotiented by S_n relabelings.")

    # ============ (b) Fix_anti(phi) in the full space ============
    print("\n[b] Fix_anti(phi) in FULL labeled tournament space {0,1}^C(n,2)")
    print("    phi(i)=n+1-i. Count T with T(phi u,phi v)=T(v,u).")
    print(f"{'n':>3} {'C(n,2)':>7} {'#selfpair':>9} {'#freeorb':>8} {'exp':>4} "
          f"{'2^exp':>10} {'exhaustive':>11} {'half':>5} {'2^half':>8}")
    for n in range(3, 8):
        exh, closed, expo, det = fix_anti_phi_full(n)
        h = half_value(n)
        flag = "OK" if (exh is None or exh == closed) else "MISMATCH"
        print(f"{n:>3} {det['Cn2']:>7} {det['n_self_paired']:>9} "
              f"{det['n_free_orbits']:>8} {expo:>4} {closed:>10} "
              f"{str(exh):>11} {h:>5} {2**h:>8}  {flag}")
    # n=8 closed form only
    exh, closed, expo, det = fix_anti_phi_full(8)
    print(f"{8:>3} {det['Cn2']:>7} {det['n_self_paired']:>9} "
          f"{det['n_free_orbits']:>8} {expo:>4} {closed:>10} "
          f"{'(skip)':>11} {half_value(8):>5} {2**half_value(8):>8}  closed-form")

    print("\n  EXPONENT formula for Fix_anti(phi)_full:")
    print("    exponent = (#phi-orbits on the C(n,2) unordered vertex-pairs)")
    print("    Derivation: pairs split into (i) self-paired pairs {u, n+1-u}")
    print("    [phi swaps the two endpoints -> constraint is vacuous, 1 free bit]")
    print("    and (ii) free 2-orbits {p, phi(p)} [T on one determines the other,")
    print("    1 free bit]. So EVERY phi-orbit on pairs = exactly 1 free bit.")
    # closed form exponent: exp = #self_paired + #free2orbits = floor(n/2) + (C(n,2)-floor(n/2))/2
    print(f"\n  {'n':>3} {'C(n,2)':>7} {'exp=Fix_anti':>12} {'self=floor(n/2)':>15} "
          f"{'half(n)':>8} {'exp-half':>9} {'half+floor(n/2)':>16}")
    closed_ok = True
    for n in range(3, 13):
        _, closed, expo, det = fix_anti_phi_full(n) if comb(n,2) <= 21 else (None, None, None, None)
        Cn2 = comb(n, 2)
        sp = n // 2                       # #self-paired pairs {u, n+1-u}
        free2 = (Cn2 - sp) // 2
        expo_cf = sp + free2              # closed-form exponent
        h = half_value(n)
        match = (expo_cf == h + sp)
        closed_ok &= match
        # also confirm exhaustive expo matches closed form where available
        if expo is not None:
            assert expo == expo_cf, f"exhaustive/closed mismatch n={n}"
        print(f"  {n:>3} {Cn2:>7} {expo_cf:>12} {sp:>15} {h:>8} {expo_cf-h:>9} {h+sp:>16}")
    print("\n  RECONCILIATION (fixed-HP frame vs full frame), CORRECTED:")
    print("    full-frame exponent = half(n) + floor(n/2).")
    print(f"    Identity  exp == half(n)+floor(n/2)  for n=3..12:  {closed_ok}")
    print("    => Fix_anti(phi)_full = 2^half(n) * 2^floor(n/2).")
    print("    The extra +floor(n/2) (NOT n-1) is the number of phi-SELF-PAIRED")
    print("    vertex-pairs {u, n+1-u}: each carries 1 free bit because the anti-")
    print("    automorphism constraint on a self-paired pair is VACUOUS (phi swaps")
    print("    its two endpoints, so T(v,u)=T(v,u) is automatic). The fixed-HP frame")
    print("    PINS exactly these floor(n/2)... no: it pins the n-1 base-path arcs and")
    print("    leaves m=C(n-1,2) tile bits, of which 2^half are grid-sym. The clean")
    print("    statement is the COUNT identity exp = half + floor(n/2); the half-tiling")
    print("    (grid-sym) part is 2^half, the remaining 2^floor(n/2) lives on the")
    print("    self-paired (anti-diagonal) pairs that the half-tiling frame collapses.")

    # ============ (c) the 2:1 orbifold fold ============
    print("\n[c] THE 2:1 FOLD (orbifold Q_m/rho): branch locus = grid-sym tilings")
    print(f"{'n':>3} {'2^m':>8} {'#orbits':>8} {'#fixed(branch)':>14} "
          f"{'#2-orbits(free)':>15} {'(2^m+2^half)/2':>15} {'OK':>4}")
    for n in range(3, 8):
        m, half, fix_id, fix_rho, predicted = merged_tiling_count(n)
        n_orb, n_fix, n_two = enumerate_rho_orbits(n)
        ok = (n_orb == predicted and n_fix == 2**half and
              n_orb == n_fix + n_two)
        print(f"{n:>3} {2**m:>8} {n_orb:>8} {n_fix:>14} {n_two:>15} "
              f"{predicted:>15} {'OK' if ok else 'FAIL'}")
    print("    free tilings pair up (2-orbits); grid-sym tilings are the branch")
    print("    locus (fixed points of the involution).  Euler-orbifold style:")
    print("    #orbits = #2-orbits + #fixed = (2^m - 2^half)/2 + 2^half = (2^m+2^half)/2.")

    # ============ (d) SC bijection? ============
    print("\n[d] SC-halving:  2^half  vs  SC_n  vs  V_merged")
    print(f"{'n':>3} {'half':>5} {'2^half':>8} {'SC_n':>6} {'V_merged':>9} "
          f"{'2^half/SC':>10} {'k=(n-1)/2':>9} {'A000568(k)':>11} {'2^C(k,2)':>9}")
    for r in sc_relation_table():
        n = r['n']
        k = r['k']
        a = A000568.get(k, '-') if k is not None else '-'
        lab = r['labeled_tour_k'] if r['labeled_tour_k'] is not None else '-'
        print(f"{n:>3} {r['half']:>5} {r['two_half']:>8} {r['SC']:>6} "
              f"{r['V_merged']:>9} {r['ratio_2half_over_SC']:>10.4f} "
              f"{str(k):>9} {str(a):>11} {str(lab):>9}")
    print("\n  open-Q#4 target is |SC(n)| = A000568(n-1) (iso classes of (n-1)-vertex")
    print("  tournaments), per two-models-staircase-recursion.md. TEST that identity:")
    print(f"  {'n':>3} {'SC_n':>6} {'A000568(n-1)':>13} {'equal?':>7} {'2^half':>8} {'2^half==SC':>11}")
    sc_eq_a = True
    for n in range(3, 9):
        a = A000568[n-1]
        h = 2**half_value(n)
        sc = SC_canon[n]
        eq = (sc == a)
        sc_eq_a &= eq
        print(f"  {n:>3} {sc:>6} {a:>13} {str(eq):>7} {h:>8} {str(h==sc):>11}")
    print(f"  => SC_n == A000568(n-1) for ALL n=3..8?  {sc_eq_a}")
    print("     (holds ONLY at n=4 and n=6; FAILS at n=3,5,7,8 -> the canon")
    print("      'SC(n)=A000568(n-1)' claim is a small-n coincidence, REFUTED as a")
    print("      general identity. SC_n = A002785(n+1), self-converse iso classes.)")

    print("\n  Test exact bijection candidates for ODD n (where (n-1)/2 is integer):")
    for n in [3, 5, 7]:
        k = (n-1)//2
        h = half_value(n)
        twoh = 2**h
        sc = SC_canon[n]
        labk = 2**comb(k, 2)
        isok = A000568[k]
        print(f"    n={n}: 2^half={twoh}, SC_n={sc}, "
              f"labeled-tour({k})=2^C({k},2)={labk}, iso-tour({k})=A000568({k})={isok}")
        print(f"          2^half == SC_n? {twoh==sc};  "
              f"2^half == labeled-tour({k})? {twoh==labk}")
    print("\n  CONCLUSION (d): 2^half(n) counts GRID-SYM TILINGS (labeled, base-path-")
    print("  fixed phi-self-converse tournaments), NOT SC iso classes. 2^half grows")
    print("  as 2^floor((n-1)^2/4) (super-exponential in n), while SC_n grows like")
    print("  A000568-style.  2^half != SC_n for all n>=4 (equal only at n=3: 2=2).")
    print("  The relationship is the BURNSIDE/fixed-point identity (2^half = #fixed")
    print("  points of rho), NOT a clean |SC(n)| <-> tournaments-on-(n-1)/2 bijection.")
    print("  The half-tiling does NOT by itself realize the SC-halving bijection of")
    print("  open-Q#4: that bijection lives at the ISO-CLASS (S_n-quotient) level,")
    print("  whereas 2^half is a LABELED/frame-fixed count. Half-tiling gives the")
    print("  correct fixed-point COUNT feeding Burnside, but the explicit class-level")
    print("  bijection still requires the S_n quotient (remains the open problem).")

    print("\nDONE.")

if __name__ == "__main__":
    main()
