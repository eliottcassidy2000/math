"""Signed-LRC: EXACT homometry-class sizes among silent(H_3)-touching cuts, reaching C=81,99,105.
(monad-explorer-S710)  Refines signed_lrc_size4_hunt_s710c.py.

For every cut in the explicit affine subspace silent(H_3) (offset interval [(m+1)/2,m-1] XOR the
value-pairing linear part), compute its FULL silent group  G_eps = { D in V : D silent at eps }
by exact-testing each element of the order-block lattice V.  The class size = |G_eps|.  Report
the size histogram among the size->=4 (rank>=2) classes -- answering S708b: do C=81/105 have
size-4 (and size-8) classes, and how many (restricted to the H_3 axis)?

|V| is small (dim tau(C)-2), so testing all D in V per cut is cheap.  silent(H_3) enumerated
explicitly (8192 cuts at C=81), so NO 2^{n-2} brute force.
"""
import sys
from collections import Counter

def proper_divisors(C):
    return [d for d in range(2, C) if C % d == 0]

def order_d_halfsystem(C, d):
    K = set(((C // d) * j) % C for j in range(d))
    return frozenset(x for x in K if 1 <= x <= (C - 1) // 2)

def diff_multiset(L, C):
    cnt = [0]*C
    for a in L:
        for b in L:
            if a != b:
                cnt[(a-b)%C] += 1
    return tuple(cnt)

def signed_set(C, n, neg):
    return [(-(i) if i in neg else i) % C for i in range(1, n)]

def silent_H3_basis(C):
    m = C // 3
    basis = [frozenset({m})]
    for a in range(2, (m-1)//2 + 1):
        basis.append(frozenset({a, m-a, m+a}))
    return basis

def silent_H3_offset(C):
    m = C // 3
    return frozenset(range((m + 1) // 2, m))

def enumerate_silentH3(C):
    basis = silent_H3_basis(C)
    offs = silent_H3_offset(C)
    pts = [frozenset(offs)]
    for b in basis:
        bb = set(b)
        pts = pts + [frozenset(p ^ bb) for p in pts]
    return pts

def V_moves(C, n):
    """nonzero elements of V = span_F2{H_d}, as frozensets of magnitudes."""
    mags = list(range(1, n))
    idx = {x:i for i,x in enumerate(mags)}
    gens = []
    for d in proper_divisors(C):
        Hd = order_d_halfsystem(C, d)
        if 1 in Hd: continue
        msk = 0
        for x in Hd: msk |= 1 << idx[x]
        gens.append((d, msk))
    basis = []
    for _, msk in gens:
        cur = msk
        for b in basis: cur = min(cur, cur^b)
        if cur: basis.append(cur); basis.sort(reverse=True)
    elems = [0]
    for b in basis: elems = elems + [x^b for x in elems]
    elems = sorted(set(e for e in elems if e))
    def mask_to_set(msk):
        return frozenset(mags[i] for i in range(len(mags)) if (msk>>i)&1)
    return [mask_to_set(e) for e in elems], len(basis), gens

def analyze(C, n, cap=None):
    pts = enumerate_silentH3(C)
    if cap and len(pts) > cap:
        print(f"=== C={C} n={n}: silent(H_3)={len(pts)} > cap {cap}, skipping ==="); return
    moves, dimV, gens = V_moves(C, n)
    sizehist = Counter()        # |G_eps| -> #cuts
    rank_ge2_examples = {}
    for neg in pts:
        base_ms = diff_multiset(signed_set(C, n, neg), C)
        G = 0   # count silent moves (incl identity)
        silent_moves = []
        for D in moves:
            neg2 = neg ^ D
            if diff_multiset(signed_set(C, n, neg2), C) == base_ms:
                silent_moves.append(D)
        gsize = len(silent_moves) + 1     # +identity
        sizehist[gsize] += 1
        if gsize >= 4 and gsize not in rank_ge2_examples:
            rank_ge2_examples[gsize] = (sorted(neg), [sorted(s) for s in silent_moves])
    # classes among these cuts: each size-s class contributes s cuts -> #classes = cuts/s
    print(f"=== C={C}={'x'.join(str(p) for p in _factor(C))} n={n}: silent(H_3) has {len(pts)} cuts, dimV={dimV} ===")
    for s in sorted(sizehist):
        cuts = sizehist[s]
        ncl = cuts / s
        print(f"  |G_eps|={s}: {cuts} cuts  (=> {ncl:.3g} classes of size {s} touching silent(H_3))")
    for s in sorted(rank_ge2_examples):
        neg, sm = rank_ge2_examples[s]
        print(f"    size-{s} example: silent moves (besides H_3) = {sm}")

def _factor(C):
    f,m,d=[],C,2
    while d*d<=m:
        while m%d==0: f.append(d); m//=d
        d+=1
    if m>1: f.append(m)
    return f or [C]

if __name__ == "__main__":
    for C, n in [(27,14),(45,23),(63,32),(81,41),(99,50)]:
        analyze(C, n)
        print()
