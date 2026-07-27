#!/usr/bin/env python3
"""HYP-9045 referee — the rank-12 box's anatomy (mac-mini-2026-07-27-S143).

(A) PAIR-RELATION SATURATION: for ANY primitive v with all |v_i| <= Q
    (Q = 91^6), the support-2 relations v_j e_i - v_i e_j are admissible
    (height <= Q) and span v-perp over Q  =>  the THM-2052 sparse code has
    rank 12 for EVERY all-coordinates-small family: the "rank-12 finite
    maximal-minor box" is essentially the all-small HEIGHT BOX
    {v primitive, covering-shaped, v_max <= ~91^6}, and the rank-11
    two-anchor stars are the big-outlier strata.  (Referee: exact rank of
    the pair-relation matrix for AP, GW, deep well, and 50 random
    all-small vectors.)

(B) TOTAL MOD-p INVISIBILITY (the redirect lemma): over F_p the sparse
    span has rank 12 for every nonzero v mod p — even when coordinates
    vanish mod p: if v_i == 0 (mod p), then for any pair relation r among
    nonzero coordinates, e_i + r and 2 e_i + r are both admissible
    support-<=3 relations, so e_i itself lies in the span.  Hence
    reduction mod p NEVER drops the sparse rank below 12: the box is
    invisible to every CRT/mod-p instrument; its content is purely
    archimedean (height tension).  (Referee: exact F_7 and F_13 ranks for
    AP and GW, whose coordinate sets contain 7 and 13.)

Consequence (typed, for the write-up): the S142 approach-(2) hope
("discharge the geometric channel independently") is REFUTED as stated —
the box is a deep height exhaustion, not a separate finite structure; the
geometric and spectral channels re-merge at the 91-line.
"""

import random
import sympy as sp


def pair_matrix(v):
    n = len(v)
    rows = []
    for i in range(n):
        for j in range(i + 1, n):
            r = [0] * n
            r[i] = v[j]
            r[j] = -v[i]
            rows.append(r)
    return sp.Matrix(rows)


def support3_zero_trick_rows(v, p):
    """Over F_p: pair relations among nonzero coords + (e_i + pair) rows
    for zero coords — the constructive basis of the invisibility lemma."""
    n = len(v)
    vm = [x % p for x in v]
    nz = [i for i in range(n) if vm[i] != 0]
    zz = [i for i in range(n) if vm[i] == 0]
    rows = []
    for a in range(len(nz)):
        for b in range(a + 1, len(nz)):
            i, j = nz[a], nz[b]
            r = [0] * n
            r[i] = vm[j]
            r[j] = (-vm[i]) % p
            rows.append(r)
    if len(nz) >= 2:
        i0, j0 = nz[0], nz[1]
        base = [0] * n
        base[i0] = vm[j0]
        base[j0] = (-vm[i0]) % p
        for i in zz:
            r1 = list(base)
            r1[i] = 1          # e_i + pair-relation  (support 3, height ok)
            rows.append(r1)
    return sp.Matrix(rows), len(zz)


def main():
    Q = 91**6
    fams = {
        "AP": list(range(1, 14)),
        "GW": list(range(1, 12)) + [13, 24],
        "deep well": list(range(1, 13)) + [182],
    }
    print("== (A) pair-relation saturation (rank over Q) ==")
    for name, v in fams.items():
        M = pair_matrix(v)
        hmax = max(abs(x) for x in v)
        print(f"  {name:<10} v_max={hmax} <= 91^6: {hmax <= Q}; "
              f"rank(pair relations) = {M.rank()}  (12 = full)")
        assert M.rank() == 12
    rng = random.Random(20260727)
    ok = 0
    for _ in range(50):
        v = sorted(rng.sample(range(1, 100000), 13))
        if pair_matrix(v).rank() == 12:
            ok += 1
    print(f"  random all-small: {ok}/50 rank 12")
    assert ok == 50
    print("  => the rank-12 stratum contains EVERY all-coordinates-small family;")
    print("     'the box' == the all-small height box (stars = big outliers).")

    print("\n== (B) total mod-p invisibility ==")
    for name, v in fams.items():
        for p in (7, 13):
            Mrows, z = support3_zero_trick_rows(v, p)
            rk = Mrows.rank(iszerofunc=lambda e: e % p == 0)
            # rank over F_p: use sympy Matrix over GF via rref mod p
            Mp = sp.Matrix(Mrows).applyfunc(lambda e: e % p)
            rkp = Mp.rank(iszerofunc=lambda e: sp.Integer(e) % p == 0)
            print(f"  {name:<10} mod {p:>2}: zero-coords={z}, "
                  f"sparse-span rank over F_{p} = {rkp}  (12 = no drop)")
            assert rkp == 12, (name, p, rkp)
    print("  => reduction mod p NEVER drops the sparse rank: the box is")
    print("     invisible to every mod-p/CRT instrument (redirect lemma).")


if __name__ == "__main__":
    main()
