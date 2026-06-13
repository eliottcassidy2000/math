"""Signed-LRC: does C=81 (=3^4) / C=105 (=3.5.7) have SIZE->=4 homometry classes?
(monad-explorer-S710)  -- answers S708b open handoff via the affine silent(H_d) subspaces.

silent(H_3) (H_3={C/3}, single order-3 flip) is an AFFINE F2-SUBSPACE (HYP-2291) with the
THM-413 value-pairing basis:
    v0 = {m}              (m=C/3)
    T_a = {a, m-a, m+a}   for a=2..(m-1)/2   (a=1 excluded: magnitude 1 is gauge-fixed)
so dim silent(H_3) = (m-1)/2 = (C/3-1)/2.  For C=81: m=27, dim 13 -> only 8192 cuts.

We CONSTRUCT silent(H_3) explicitly, then EXACT-test each member for H_d-silence (other single
moves) via the difference multiset.  Co-silence with an independent move => a SIZE->=4 class.
This reaches C=81/105 where 2^{n-2} brute force is impossible.

Validation: value-pairing construction == brute silent(H_3) for all composite 3|C, C<=39.
"""
import sys
from collections import Counter

def proper_divisors(C):
    return [d for d in range(2, C) if C % d == 0]

def order_d_halfsystem(C, d):
    K = set(((C // d) * j) % C for j in range(d))
    return frozenset(x for x in K if 1 <= x <= (C - 1) // 2)

def diff_multiset(signed_set, C):
    cnt = [0]*C
    L = list(signed_set)
    for a in L:
        for b in L:
            if a != b:
                cnt[(a-b)%C] += 1
    return tuple(cnt)

def cut_signed_set(C, n, neg_mags):
    """signed half-system for cut that negates magnitudes in neg_mags (magnitude 1 always +)."""
    return [(-(i) if i in neg_mags else i) % C for i in range(1, n)]

def homometric_after_flip(C, n, neg_mags, D):
    """True iff flipping move D (set of magnitudes) is silent at the cut neg_mags."""
    s1 = cut_signed_set(C, n, neg_mags)
    neg2 = neg_mags ^ D                      # symmetric difference
    s2 = cut_signed_set(C, n, neg2)
    return diff_multiset(s1, C) == diff_multiset(s2, C)

def silent_H3_basis(C):
    """value-pairing basis (list of frozenset magnitudes) for the linear part of silent(H_3)."""
    m = C // 3
    assert C % 3 == 0
    basis = [frozenset({m})]
    for a in range(2, (m-1)//2 + 1):
        trip = frozenset({a, m-a, m+a})
        assert len(trip) == 3 and all(1 <= x <= (C-1)//2 for x in trip)
        basis.append(trip)
    return basis

def silent_H3_offset(C):
    """a canonical silent cut for flipping H_3 = the interval [(m+1)/2, m-1], m=C/3
    (verified the coset representative; size (m-1)/2)."""
    m = C // 3
    return frozenset(range((m + 1) // 2, m))

def enumerate_subspace(basis, offset=frozenset()):
    """yield all (offset XOR XOR-combinations) (frozenset of magnitudes)."""
    bl = [set(b) for b in basis]
    pts = [frozenset(offset)]
    for b in bl:
        pts = pts + [frozenset(p ^ b) for p in pts]
    return pts

# ---- brute silent(H_3) for validation ----
def brute_silent(C, n, D):
    from itertools import product
    res = set()
    free = list(range(2, n))   # magnitudes 2..n-1
    for bits in range(1 << len(free)):
        neg = frozenset(free[k] for k in range(len(free)) if (bits>>k)&1)
        if homometric_after_flip(C, n, neg, D):
            res.add(neg)
    return res

def validate():
    print("=== validate value-pairing silent(H_3) construction vs brute (3|C, C<=39) ===")
    ok = True
    for n in [5, 8, 11, 14, 17, 20]:        # C=9,15,21,27,33,39 (all 3|C)
        C = 2*n-1
        if C % 3: continue
        m = C//3
        Hd = order_d_halfsystem(C, 3)
        D = frozenset(Hd)
        basis = silent_H3_basis(C)
        constructed = set(enumerate_subspace(basis, silent_H3_offset(C)))
        bruted = brute_silent(C, n, D)
        match = constructed == bruted
        ok = ok and match
        print(f"  C={C:3d} m={m:2d}: |constructed|={len(constructed):5d} |brute|={len(bruted):5d} "
              f"dim={len(basis)} {'OK' if match else '*** MISMATCH ***'}")
    print(f"  ALL MATCH: {ok}\n")
    return ok

def hunt(C, n, report=True):
    """Enumerate silent(H_3), test co-silence with every OTHER single subgroup-half move."""
    assert C % 3 == 0
    m = C//3
    basis = silent_H3_basis(C)
    pts = enumerate_subspace(basis, silent_H3_offset(C))   # all of silent(H_3)
    divs = [d for d in proper_divisors(C) if d != 3]
    other_moves = {}
    for d in divs:
        Hd = order_d_halfsystem(C, d)
        if 1 in Hd:    # never (proper d)
            continue
        other_moves[d] = frozenset(Hd)
    # for each cut in silent(H_3), which other moves are ALSO silent?
    cosilent = Counter()       # d -> #cuts in silent(H_3) where H_d also silent
    examples = {}
    for neg in pts:
        for d, D in other_moves.items():
            if homometric_after_flip(C, n, neg, D):
                cosilent[d] += 1
                if d not in examples:
                    examples[d] = sorted(neg)
    if report:
        print(f"=== C={C} (n={n}, m=C/3={m}) : silent(H_3) has {len(pts)} cuts (dim {len(basis)}) ===")
        print(f"  other single moves H_d, d|C: {sorted(other_moves)}")
        if not cosilent:
            print(f"  NO other single move is ever co-silent with H_3 on silent(H_3).")
            print(f"  ==> NO size->=4 class arising from two co-silent single subgroup-halves "
                  f"(via the H_3 axis).")
        for d in sorted(cosilent):
            print(f"  H_{d} co-silent with H_3 at {cosilent[d]} of {len(pts)} cuts "
                  f"=> SIZE->=4 classes EXIST. example neg-mags={examples[d]}")
    return cosilent, len(pts)

if __name__ == "__main__":
    okv = validate()
    print("###### deploying at the brute-infeasible moduli ######\n")
    # C=27 (sanity: known 1 size-4 class), C=81=3^4, C=63=3^2*7, C=45=3^2*5, C=105=3*5*7, C=99=3^2*11
    for C, n in [(27,14),(45,23),(63,32),(81,41),(99,50),(105,53)]:
        hunt(C, n)
        print()
