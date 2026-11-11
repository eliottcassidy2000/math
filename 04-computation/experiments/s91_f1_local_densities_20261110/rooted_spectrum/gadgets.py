# gadgets.py - closing gadgets with at most two new vertices (question Q3).
#
# (T, v) rooted, B = T - v.  A gadget vertex g is of type O (g -> every vertex of B) or I (every vertex of B -> g),
# with an arbitrary arc to v; two gadget vertices carry an arbitrary arc between them.  For every gadget the count
# H(T + G) is fitted (exact rational elimination on 8 independent rooted tournaments) as an integer combination of
#     (H, s, e, Yplus, Yminus, d, P2T, P2B)      d = H(T - v), P2 = ordered covers by two non-empty paths,
# and the fitted identity is then verified on every rooted tournament of the universe:
#     all labelled tournaments of orders 1..5 with every root (5,405 rooted tournaments),
#     all classes of orders 6 and 7 with every root (336 + 3,192).
# The fitted coefficients are compared with the closed forms derived by hand:
#     one vertex:   O, g -> v: H            O, v -> g: H - s + Yplus      I, g -> v: H - e + Yminus     I, v -> g: H
#     O, O:         v source of the triple: H - s - 2d + 3 Yplus;  3-cycle: H + 2 Yplus;  others reduce to one vertex
#     I, I:         v sink of the triple:   H - e - 2d + 3 Yminus; 3-cycle: H + 2 Yminus; others reduce to one vertex
#     O = o, I = i, al = [v -> o], be = [i -> v], ga = [i -> o]:
#        (1 + 2 ga) H - (al + be ga) s - (be + al ga) e + (al - be ga) Yplus + (be - al ga) Yminus
#        + 2 al be d + ga P2T + al be P2B
import os
import sys
from fractions import Fraction
import numpy as np
from numba import njit
import core

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = []
NAMES = ["H", "s", "e", "Y+", "Y-", "d", "P2(T)", "P2(T-v)"]


def check(name, ok, detail=""):
    print(("ok   " if ok else "FAIL ") + name + ((" : " + detail) if detail else ""), flush=True)
    if not ok:
        FAIL.append(name)


@njit(cache=True)
def gadget_count(n, code, v, k, typ, atv, internal):
    # typ[g] = 0 (O: g -> B) or 1 (I: B -> g); atv[g] = 1 iff v -> g; internal = 1 iff g0 -> g1
    N = n + k
    out = np.zeros(N, np.int64)
    core.out_masks(code, n, out)
    for g in range(k):
        gb = 1 << (n + g)
        for u in range(n):
            if u == v:
                if atv[g] == 1:
                    out[u] |= gb
                else:
                    out[n + g] |= (1 << u)
            else:
                if typ[g] == 1:
                    out[u] |= gb
                else:
                    out[n + g] |= (1 << u)
    if k == 2:
        if internal == 1:
            out[n] |= (1 << (n + 1))
        else:
            out[n + 1] |= (1 << n)
    size = 1 << N
    f = np.zeros((size, N), np.int64)
    g2 = np.zeros((size, N), np.int64)
    hS = np.zeros(size, np.int64)
    core.dp_tables(N, out, f, g2, hS)
    return hS[size - 1]


@njit(cache=True)
def all_counts(n, codes, k, typ, atv, internal, res):
    r = 0
    for c in range(codes.shape[0]):
        for v in range(n):
            res[r] = gadget_count(n, codes[c], v, k, typ, atv, internal)
            r += 1
    return r


def solve_exact(Arows, b):
    # Gaussian elimination over Q; Arows: list of 8 integer rows of length 8
    n = len(Arows)
    M = [[Fraction(int(x)) for x in row] + [Fraction(int(bb))] for row, bb in zip(Arows, b)]
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return None
        M[c], M[p] = M[p], M[c]
        piv = M[c][c]
        M[c] = [x / piv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                fac = M[r][c]
                M[r] = [x - fac * y for x, y in zip(M[r], M[c])]
    return [M[r][n] for r in range(n)]


def rank_int(rows):
    M = [[Fraction(int(x)) for x in row] for row in rows]
    rk = 0
    ncol = len(M[0])
    for c in range(ncol):
        p = next((r for r in range(rk, len(M)) if M[r][c] != 0), None)
        if p is None:
            continue
        M[rk], M[p] = M[p], M[rk]
        piv = M[rk][c]
        M[rk] = [x / piv for x in M[rk]]
        for r in range(len(M)):
            if r != rk and M[r][c] != 0:
                fac = M[r][c]
                M[r] = [x - fac * y for x, y in zip(M[r], M[rk])]
        rk += 1
    return rk


def fmt(coef):
    out = []
    for c, nm in zip(coef, NAMES):
        if c == 0:
            continue
        c = int(c)
        sgn = "+" if c > 0 else "-"
        a = abs(c)
        out.append("%s %s%s" % (sgn, "" if a == 1 else "%d " % a, nm))
    s = " ".join(out)
    return s[2:] if s.startswith("+ ") else s


def closed_form(k, typ, atv, internal):
    # hand-derived coefficients over NAMES
    if k == 1:
        if typ[0] == 0:
            return [1, 0, 0, 0, 0, 0, 0, 0] if atv[0] == 0 else [1, -1, 0, 1, 0, 0, 0, 0]
        return [1, 0, -1, 0, 1, 0, 0, 0] if atv[0] == 0 else [1, 0, 0, 0, 0, 0, 0, 0]
    if typ[0] == typ[1]:
        # tournament on {v, g0, g1}
        beats_v = [1 - atv[0], 1 - atv[1]]            # g -> v
        g0g1 = internal
        if typ[0] == 0:                               # O, O (dual of the tail law)
            if atv[0] == 1 and atv[1] == 1:           # v beats both: v source of the triple
                return [1, -1, 0, 3, 0, -2, 0, 0]
            if atv[0] == 0 and atv[1] == 0:           # both beat v: both are sources of T + G in turn
                return [1, 0, 0, 0, 0, 0, 0, 0]
            # exactly one, say a, is beaten by v; the other, b, beats v
            a = 0 if atv[0] == 1 else 1
            a_beats_b = (g0g1 == 1) if a == 0 else (g0g1 == 0)
            if a_beats_b:                             # v -> a -> b -> v: 3-cycle
                return [1, 0, 0, 2, 0, 0, 0, 0]
            return [1, -1, 0, 1, 0, 0, 0, 0]          # b is a source: one-vertex gadget a
        else:                                         # I, I (tail law)
            if atv[0] == 0 and atv[1] == 0:           # both beat v: v sink of the triple
                return [1, 0, -1, 0, 3, -2, 0, 0]
            if atv[0] == 1 and atv[1] == 1:
                return [1, 0, 0, 0, 0, 0, 0, 0]
            a = 0 if atv[0] == 0 else 1               # a beats v, b is beaten by v
            b_beats_a = (g0g1 == 0) if a == 0 else (g0g1 == 1)
            if b_beats_a:                             # a -> v -> b -> a: 3-cycle
                return [1, 0, 0, 0, 2, 0, 0, 0]
            return [1, 0, -1, 0, 1, 0, 0, 0]          # b is a sink: one-vertex gadget a
    # mixed: put o = the O vertex, i = the I vertex
    o = 0 if typ[0] == 0 else 1
    i = 1 - o
    al = atv[o]
    be = 1 - atv[i]
    i_beats_o = (internal == 0) if o == 0 else (internal == 1)
    ga = 1 if i_beats_o else 0
    return [1 + 2 * ga, -(al + be * ga), -(be + al * ga), al - be * ga, be - al * ga, 2 * al * be, ga, al * be]


def main():
    z = np.load(os.path.join(HERE, "reps_le8.npz"))
    universe = []          # (n, codes)
    for n in range(1, 6):
        universe.append((n, np.arange(1 << (n * (n - 1) // 2), dtype=np.int64)))
    universe.append((6, z["n6"]))
    universe.append((7, z["n7"]))
    inv = []
    for n, codes in universe:
        buf = np.zeros((len(codes) * n, 10), np.int64)
        r = core.rooted_rows_direct(n, codes, buf)
        inv.append(buf[:r][:, [0, 1, 2, 3, 4, 5, 8, 9]].copy())
    INV = np.concatenate(inv)
    print("universe: %d rooted tournaments (labelled orders 1-5, classes of orders 6, 7)" % len(INV))
    rk = rank_int(np.unique(INV, axis=0))
    check("the eight invariants are linearly independent on the universe (rank %d)" % rk, rk == 8)
    # pick 8 independent rows greedily
    pick = []
    for idx in range(len(INV)):
        cand = pick + [idx]
        if rank_int([INV[i] for i in cand]) == len(cand):
            pick = cand
        if len(pick) == 8:
            break
    gadgets = []
    for t0 in (0, 1):
        for a0 in (0, 1):
            gadgets.append((1, (t0, 0), (a0, 0), 0))
    for t0 in (0, 1):
        for t1 in (0, 1):
            if t0 > t1:
                continue                       # (I, O) is (O, I) with the names exchanged
            for a0 in (0, 1):
                for a1 in (0, 1):
                    for it in (0, 1):
                        gadgets.append((2, (t0, t1), (a0, a1), it))
    print("gadgets: %d labelled (4 with one vertex, 8 + 8 + 8 with two)" % len(gadgets))
    coefs = []
    tn = {0: "O", 1: "I"}
    for (k, typ, atv, it) in gadgets:
        typ_a = np.array(typ, np.int64)
        atv_a = np.array(atv, np.int64)
        vals = []
        for n, codes in universe:
            res = np.zeros(len(codes) * n, np.int64)
            all_counts(n, codes, k, typ_a, atv_a, it, res)
            vals.append(res)
        V = np.concatenate(vals)
        sol = solve_exact([INV[i] for i in pick], [V[i] for i in pick])
        ok_int = sol is not None and all(c.denominator == 1 for c in sol)
        coef = np.array([int(c) for c in sol], np.int64) if ok_int else None
        ok_all = ok_int and np.array_equal(INV @ coef, V)
        cf = closed_form(k, typ, atv, it)
        desc = "%s%s" % (tn[typ[0]], "" if k == 1 else tn[typ[1]])
        arcs = []
        for g in range(k):
            arcs.append(("v->g%d" % g) if atv[g] == 1 else ("g%d->v" % g))
        if k == 2:
            arcs.append("g0->g1" if it == 1 else "g1->g0")
        line = "%-2s %-22s H(T+G) = %s" % (desc, " ".join(arcs), fmt(coef) if ok_int else "NO INTEGER FIT")
        print(line)
        if not ok_all:
            check("identity verified on the whole universe for gadget " + desc + " " + " ".join(arcs), False)
        if not (ok_int and list(coef) == list(cf)):
            check("closed form agrees for gadget " + desc + " " + " ".join(arcs), False, str(cf))
        if ok_int:
            coefs.append(coef)
    check("all %d gadget identities fitted with integer coefficients and verified on %d rooted tournaments each"
          % (len(gadgets), len(INV)), len(coefs) == len(gadgets) and not FAIL)
    C = np.array(coefs)
    distinct = np.unique(C, axis=0)
    print("distinct identities: %d" % len(distinct))
    for row in distinct:
        print("   ", fmt(row))
    rk = rank_int(distinct)
    check("the coefficient vectors span a space of rank 8 (every one of the eight invariants is needed)", rk == 8,
          "rank %d" % rk)
    check("coefficient of H is at least 1 in every identity", int(C[:, 0].min()) >= 1, "min %d" % C[:, 0].min())
    # no identity of the form (single rooted count) + constant or 2 (rooted count) + 1: the H-coefficient never vanishes
    pure = [row for row in distinct if row[0] == 0]
    check("no gadget identity is free of H", len(pure) == 0)
    # the interior coefficient: write H = (H - s - e) + s + e
    cint = C[:, 0]
    cs = C[:, 0] + C[:, 1]
    ce = C[:, 0] + C[:, 2]
    print("coefficients of (interior paths, paths starting at v, paths ending at v): ",
          sorted(set((int(a), int(b), int(c)) for a, b, c in zip(cint, cs, ce))))
    check("coefficients of the three one-path classes are non-negative, interior >= 1",
          int(cint.min()) >= 1 and int(cs.min()) >= 0 and int(ce.min()) >= 0)
    if FAIL:
        print("FAILED:", FAIL)
        sys.exit(1)
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
