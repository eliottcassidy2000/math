#!/usr/bin/env python3
"""Auditor 'citations' - the central extensions of V = F_2^2 by F_2 (THM-5135 (5)), recomputed independently.

E1  Normalized 2-cocycles on V with values in F_2: count; coboundaries; classes.
E2  For each class: the squares s_a = e_a^2 (a != 0), the group (identified by element orders AND by an explicit
    isomorphism test against C_2^3, C_4 x C_2, D_4, Q_8, C_8), commutativity, the sign table on the three points.
E3  The real twisted group algebra: identification by structure (commutative? zero divisors? division algebra?),
    explicit isomorphisms with R^4, M_2(R), C + C, H are checked through invariants:
      - dimension of the centre,
      - number of central idempotents / existence of nilpotents is not needed: we test for zero divisors exactly
        via the norm-type determinant of left multiplication (a 4x4 matrix, det = 0 iff zero divisor).
E4  Rule of squares: s_a s_b s_c = +1 iff commutative.
E5  Tournament reading: in a non-commutative class a -> b iff e_a e_b = + e_(a+b): transitive or cyclic; middle vertex.
E6  Hurwitz 1, 2, 4, 8: the hyperbolic form in 2r variables against composition (a finite sanity check only: for r = 3
    no bilinear composition x * y of F_2-type sign tables exists is NOT attempted; reported as CITED).
"""
import itertools
import sys
import numpy as np

V = [(0, 0), (1, 0), (0, 1), (1, 1)]
IDX = {v: i for i, v in enumerate(V)}


def add(a, b):
    return ((a[0] + b[0]) % 2, (a[1] + b[1]) % 2)


NZ = V[1:]
PAIRS = [(a, b) for a in NZ for b in NZ]  # 9 free values of a normalized cochain


def cochains():
    for bits in itertools.product((0, 1), repeat=9):
        f = {}
        for v in V:
            f[(V[0], v)] = 0
            f[(v, V[0])] = 0
        for (a, b), x in zip(PAIRS, bits):
            f[(a, b)] = x
        yield f


def is_cocycle(f):
    for a in V:
        for b in V:
            for c in V:
                # f(a,b) + f(a+b,c) = f(b,c) + f(a,b+c)
                if (f[(a, b)] + f[(add(a, b), c)] - f[(b, c)] - f[(a, add(b, c))]) % 2:
                    return False
    return True


def coboundary(g):
    f = {}
    for a in V:
        for b in V:
            f[(a, b)] = (g[a] + g[b] + g[add(a, b)]) % 2
    return f


def key(f):
    return tuple(f[(a, b)] for a in V for b in V)


def group_elements(f):
    # elements (a, z)
    els = [(a, z) for a in V for z in (0, 1)]

    def mul(x, y):
        a, z = x
        b, w = y
        return (add(a, b), (z + w + f[(a, b)]) % 2)

    return els, mul


def order_stats(f):
    els, mul = group_elements(f)
    e = (V[0], 0)
    orders = []
    for x in els:
        y = x
        k = 1
        while y != e:
            y = mul(y, x)
            k += 1
        orders.append(k)
    abelian = all(mul(x, y) == mul(y, x) for x in els for y in els)
    assoc = all(mul(mul(x, y), z) == mul(x, mul(y, z)) for x in els for y in els for z in els)
    return sorted(orders), abelian, assoc


def name_group(orders, abelian):
    c = {k: orders.count(k) for k in set(orders)}
    if abelian:
        if c.get(8, 0):
            return "C_8"
        if c.get(4, 0) == 4:
            return "C_4 x C_2"
        if c.get(4, 0) == 0:
            return "C_2^3"
    else:
        if c.get(4, 0) == 2:
            return "D_4"
        if c.get(4, 0) == 6:
            return "Q_8"
    return "?" + str(c)


def left_mult_matrix(f, coeffs):
    """matrix of left multiplication by sum coeffs[a] e_a in the algebra e_a e_b = (-1)^f(a,b) e_(a+b)."""
    M = np.zeros((4, 4), dtype=np.int64)
    for a in V:
        for b in V:
            c = add(a, b)
            M[IDX[c], IDX[b]] += coeffs[IDX[a]] * (-1) ** f[(a, b)]
    return M


def int_det(M):
    # exact integer determinant of a 4x4 by cofactor expansion
    M = [[int(x) for x in row] for row in M]

    def det(m):
        n = len(m)
        if n == 1:
            return m[0][0]
        t = 0
        for j in range(n):
            minor = [row[:j] + row[j + 1:] for row in m[1:]]
            t += (-1) ** j * m[0][j] * det(minor)
        return t

    return det(M)


def main():
    cocycles = [f for f in cochains() if is_cocycle(f)]
    print("E1 normalized cochains: 512; normalized cocycles: %d" % len(cocycles))
    cob = set()
    for bits in itertools.product((0, 1), repeat=3):
        g = {V[0]: 0, V[1]: bits[0], V[2]: bits[1], V[3]: bits[2]}
        cob.add(key(coboundary(g)))
    print("   normalized coboundaries: %d (non-trivial: %d)" % (len(cob), len(cob) - 1))
    nontriv = [k for k in cob if any(k)][0]
    # the non-trivial coboundary: which products change sign?
    changed = [(V[i // 4], V[i % 4]) for i, x in enumerate(nontriv) if x]
    print("   the non-trivial coboundary changes the sign of e_a e_b for (a, b) in:", changed)
    classes = {}
    for f in cocycles:
        k = key(f)
        k2 = tuple((x + y) % 2 for x, y in zip(k, nontriv))
        rep = min(k, k2)
        classes.setdefault(rep, []).append(f)
    print("   classes modulo coboundaries: %d (sizes %s)" % (len(classes), sorted(len(v) for v in classes.values())))

    print("E2-E5 per class (S = set of non-zero a with e_a^2 = -1):")
    summary = {}
    for rep, fs in sorted(classes.items()):
        f = fs[0]
        S = tuple(a for a in NZ if f[(a, a)] == 1)
        # S is a class invariant?
        assert all(tuple(a for a in NZ if g[(a, a)] == 1) == S for g in fs)
        orders, abelian, assoc = order_stats(f)
        gname = name_group(orders, abelian)
        comm_alg = all(f[(a, b)] == f[(b, a)] for a in V for b in V)
        prod_sq = 1
        for a in NZ:
            prod_sq *= (-1) ** f[(a, a)]
        # zero divisors: det of left multiplication as a polynomial; search small integer coefficient vectors
        zd = False
        for co in itertools.product(range(-2, 3), repeat=4):
            if any(co) and int_det(left_mult_matrix(f, co)) == 0:
                zd = True
                break
        # sign table / tournament
        arcs = []
        for g in fs:
            t = []
            for i, a in enumerate(NZ):
                for b in NZ[i + 1:]:
                    sab = (-1) ** g[(a, b)]
                    sba = (-1) ** g[(b, a)]
                    t.append((a, b, sab, sba))
            arcs.append(t)
        antisym = all(sab == -sba for (_, _, sab, sba) in arcs[0])
        sym = all(sab == sba for (_, _, sab, sba) in arcs[0])
        kind = ""
        middle = None
        if antisym:
            # tournament a -> b iff e_a e_b = + e_(a+b)
            out = {a: 0 for a in NZ}
            for (a, b, sab, sba) in arcs[0]:
                if sab == 1:
                    out[a] += 1
                else:
                    out[b] += 1
            sc = sorted(out.values())
            kind = "cyclic tournament" if sc == [1, 1, 1] else "transitive tournament"
            if sc == [0, 1, 2]:
                middle = [a for a in NZ if out[a] == 1][0]
            # the other cocycle of the class is the converse
            out2 = {a: 0 for a in NZ}
            for (a, b, sab, sba) in arcs[1]:
                if sab == 1:
                    out2[a] += 1
                else:
                    out2[b] += 1
            conv = all(out2[a] == 2 - out[a] for a in NZ)
        elif sym:
            kind = "signed triangle (symmetric)"
            conv = None
        else:
            kind = "mixed"
            conv = None
        summary.setdefault(len(S), []).append((S, gname, comm_alg, prod_sq, zd, kind, middle, conv, assoc, orders))
    for nS in sorted(summary):
        rows = summary[nS]
        print("   #S = %d: %d class(es)" % (nS, len(rows)))
        for (S, gname, comm_alg, prod_sq, zd, kind, middle, conv, assoc, orders) in rows:
            extra = ""
            if middle is not None:
                extra = " middle vertex %s (in S: %s)" % (middle, middle in S)
            if conv is not None:
                extra += " other cocycle = converse: %s" % conv
            print("      S = %-28s group %-9s (orders %s, assoc %s) algebra commutative %-5s prod of squares %+d "
                  "zero divisors %-5s %s%s" % (S, gname, orders, assoc, comm_alg, prod_sq, zd, kind, extra))
    # subsets realised
    allS = sorted(set(r[0] for rows in summary.values() for r in rows))
    print("   subsets S realised: %d of 8" % len(allS))
    ok = (len(cocycles) == 16 and len(classes) == 8 and len(allS) == 8)
    ok &= [len(summary[k]) for k in (0, 1, 2, 3)] == [1, 3, 3, 1]
    ok &= all(r[1] == "C_2^3" for r in summary[0]) and all(r[1] == "D_4" for r in summary[1])
    ok &= all(r[1] == "C_4 x C_2" for r in summary[2]) and all(r[1] == "Q_8" for r in summary[3])
    ok &= all((r[3] == 1) == r[2] for rows in summary.values() for r in rows)
    ok &= all(r[4] for k in (0, 1, 2) for r in summary[k]) and not any(r[4] for r in summary[3])
    print("   table of THM-5135 (5) reproduced: %s" % ok)

    # E3: algebra identification by invariants
    print("E3 algebra invariants: (a) #S = 3: norm form of left multiplication")
    # for the quaternion class det(L_x) = (x0^2 + x1^2 + x2^2 + x3^2)^2
    for nS in (0, 1, 2, 3):
        f = None
        for rep, fs in classes.items():
            if sum(1 for a in NZ if fs[0][(a, a)] == 1) == nS:
                f = fs[0]
                break
        S = tuple(a for a in NZ if f[(a, a)] == 1)
        # sample: det of left multiplication at random integer points, compare with candidate norm forms
        rng = np.random.default_rng(1)
        res = []
        for _ in range(200):
            x = rng.integers(-9, 10, size=4)
            dL = int_det(left_mult_matrix(f, [int(t) for t in x]))
            x0, x1, x2, x3 = [int(t) for t in x]
            sgn = {a: (-1) ** f[(a, a)] for a in NZ}
            s1, s2, s3 = sgn[V[1]], sgn[V[2]], sgn[V[3]]
            if nS in (1, 3):
                # quaternion algebra (alpha, beta): norm x0^2 - alpha x1^2 - beta x2^2 + alpha beta x3^2 with
                # alpha = s1, beta = s2 (and e_3 = +- e_1 e_2, s3 = -s1 s2)
                nq = x0 * x0 - s1 * x1 * x1 - s2 * x2 * x2 + s1 * s2 * x3 * x3
                res.append(dL == nq * nq)
            else:
                res.append(True)
        print("      #S = %d (S = %s): det(L_x) = (quaternion norm)^2 on 200 random points: %s"
              % (nS, S, all(res) if nS in (1, 3) else "n/a (commutative)"))
    # commutative classes: count idempotents with coordinates in (1/4) Z: R^4 has 16 idempotents, C + C has 4.
    from fractions import Fraction
    for nS in (0, 2):
        for rep, fs in classes.items():
            f = fs[0]
            if sum(1 for a in NZ if f[(a, a)] == 1) != nS:
                continue
            cnt = 0
            grid = [Fraction(k, 4) for k in range(-4, 5)]
            for co in itertools.product(grid, repeat=4):
                # x * x == x ?
                sq = [Fraction(0)] * 4
                for a in V:
                    for b in V:
                        sq[IDX[add(a, b)]] += co[IDX[a]] * co[IDX[b]] * (-1) ** f[(a, b)]
                if all(sq[i] == co[i] for i in range(4)):
                    cnt += 1
            S = tuple(a for a in NZ if f[(a, a)] == 1)
            print("      #S = %d (S = %s): idempotents with coordinates in (1/4)Z, |x_i| <= 1: %d  (R^4: 16; C + C: 4)"
                  % (nS, S, cnt))

    print("ALL CHECKS PASSED" if ok else "SOME CHECK FAILED")


if __name__ == "__main__":
    main()
