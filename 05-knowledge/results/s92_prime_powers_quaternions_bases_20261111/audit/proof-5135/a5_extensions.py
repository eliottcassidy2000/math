#!/usr/bin/env python3
"""Auditor proof-5135, script a5: THM-5135 (5), the central extensions of F_2^2 by F_2.

Everything from scratch (written before reading the session's script):
  [Z] normalized 2-cochains f on V = F_2^2 (2^9), the cocycle condition on all triples, coboundaries of the
      normalized 1-cochains, classes; the invariant S = {a != 0 : f(a,a) = 1}.
  [G] the group on V x F_2: associativity, orders of elements, abelian or not, identification among the five
      groups of order 8 (C_8, C_4 x C_2, C_2^3, D_4, Q_8) by (abelian, order statistics); the number of elements
      of order 4 alone does NOT separate C_8 from D_4 (both 2), so exponent is recorded too.
  [A] the real algebra e_a e_b = (-1)^f(a,b) e_(a+b): associativity, commutativity, centre, trace form
      tr(L_x L_y) and its signature, explicit isomorphisms:
        #S = 0 : four orthogonal idempotents (R^4);
        #S = 2 : two orthogonal idempotents e+, e- with e+ A = C, e- A = C (C + C), an element of square -1 in each;
        #S = 1 : a faithful 2 x 2 real matrix model (M_2(R));
        #S = 3 : Hamilton relations i^2 = j^2 = k^2 = ijk = -1 (up to the sign convention), norm multiplicative.
  [T] sign tables: tournament (non-commutative) or signed triangle (commutative); cyclic / transitive, middle
      vertex; effect of the non-trivial coboundary; the rule s_a s_b s_c = eps.
  [0] zero divisors: explicit ones for #S = 0, 1, 2; for #S = 3 the determinant of left multiplication is
      (sum of squares)^2 (checked as a polynomial identity on a grid of integer points of degree-exceeding size).
  [P] the action of GL_2(F_2) = S_3 on the classes.
"""
import sys
import itertools
from fractions import Fraction

import numpy as np

V = [0, 1, 2, 3]            # F_2^2 as bit masks, addition = xor
NZ = [1, 2, 3]


def cochains():
    """all normalized cochains: dict (a,b) -> bit, zero when a or b is 0."""
    pairs = [(a, b) for a in NZ for b in NZ]
    for bits in itertools.product((0, 1), repeat=9):
        f = {(a, b): 0 for a in V for b in V}
        for (pr, bit) in zip(pairs, bits):
            f[pr] = bit
        yield f


def is_cocycle(f):
    for a in V:
        for b in V:
            for c in V:
                if (f[(a, b)] + f[(a ^ b, c)] + f[(b, c)] + f[(a, b ^ c)]) % 2:
                    return False
    return True


def key(f):
    return tuple(f[(a, b)] for a in NZ for b in NZ)


def group_data(f):
    els = [(a, z) for a in V for z in (0, 1)]

    def mul(x, y):
        return (x[0] ^ y[0], (x[1] + y[1] + f[(x[0], y[0])]) % 2)

    # associativity
    assoc = all(mul(mul(x, y), z) == mul(x, mul(y, z)) for x in els for y in els for z in els)
    e = (0, 0)
    orders = []
    for x in els:
        k = 1
        y = x
        while y != e:
            y = mul(y, x)
            k += 1
        orders.append(k)
    abel = all(mul(x, y) == mul(y, x) for x in els for y in els)
    stat = tuple(orders.count(o) for o in (1, 2, 4, 8))
    names = {(True, (1, 7, 0, 0)): "C_2^3", (True, (1, 3, 4, 0)): "C_4 x C_2", (True, (1, 1, 2, 4)): "C_8",
             (False, (1, 5, 2, 0)): "D_4", (False, (1, 1, 6, 0)): "Q_8"}
    return assoc, abel, stat, names.get((abel, stat), "?")


def algebra(f):
    """structure: left multiplication matrices L[a] of e_a on the basis (e_0, e_1, e_2, e_3)."""
    L = []
    for a in V:
        M = np.zeros((4, 4), dtype=np.int64)
        for b in V:
            M[a ^ b, b] = (-1) ** f[(a, b)]
        L.append(M)
    return L


def alg_mul(f, x, y):
    out = [Fraction(0)] * 4
    for a in V:
        if x[a] == 0:
            continue
        for b in V:
            if y[b] == 0:
                continue
            out[a ^ b] += x[a] * y[b] * (-1) ** f[(a, b)]
    return out


def signature(M):
    w = np.linalg.eigvalsh(M.astype(float))
    return int((w > 1e-9).sum()), int((w < -1e-9).sum())


def main():
    out = []
    P = out.append

    coc = [f for f in cochains() if is_cocycle(f)]
    P("[Z] normalized cochains: %d; cocycles: %d" % (2 ** 9, len(coc)))
    # coboundaries of normalized 1-cochains
    cobs = set()
    for g1, g2, g3 in itertools.product((0, 1), repeat=3):
        g = {0: 0, 1: g1, 2: g2, 3: g3}
        d = {(a, b): (g[a] + g[b] + g[a ^ b]) % 2 for a in V for b in V}
        cobs.add(key(d))
    P("    normalized 1-cochains: 8; distinct coboundaries: %d: %s" % (len(cobs), sorted(cobs)))
    nontriv = [c for c in cobs if any(c)][0]
    pairs = [(a, b) for a in NZ for b in NZ]
    P("    the non-trivial coboundary is 1 exactly on the ordered pairs %s (six ordered pairs a != b, i.e. the three unordered products, both orders)"
      % [pr for pr, bit in zip(pairs, nontriv) if bit])
    # classes
    cockeys = {key(f): f for f in coc}
    classes = []
    seen = set()
    for k in sorted(cockeys):
        if k in seen:
            continue
        other = tuple((x + y) % 2 for x, y in zip(k, nontriv))
        assert other in cockeys
        seen.add(k)
        seen.add(other)
        classes.append((k, other))
    P("    classes modulo coboundaries: %d (each of size 2)" % len(classes))
    Ssets = {}
    for (k, other) in classes:
        f = cockeys[k]
        f2 = cockeys[other]
        S = tuple(a for a in NZ if f[(a, a)])
        S2 = tuple(a for a in NZ if f2[(a, a)])
        assert S == S2
        Ssets.setdefault(S, []).append((k, other))
    P("    the set S (units of square -1) by class: %s" % sorted(Ssets))
    P("    every subset of the three points occurs exactly once: %s"
      % (sorted(Ssets) == sorted(tuple(s) for r in range(4) for s in itertools.combinations(NZ, r)) and all(len(v) == 1 for v in Ssets.values())))

    P("")
    P("[G][A][T] by class")
    summary = {}
    rule_ok = True
    for S in sorted(Ssets, key=lambda s: (len(s), s)):
        (k, other) = Ssets[S][0]
        rows = []
        for kk in (k, other):
            f = cockeys[kk]
            assoc, abel, stat, name = group_data(f)
            L = algebra(f)
            # associativity and commutativity of the algebra on basis elements
            a_assoc = all((L[a] @ L[b] == ((-1) ** f[(a, b)]) * L[a ^ b]).all() for a in V for b in V)
            a_comm = all((L[a] @ L[b] == L[b] @ L[a]).all() for a in V for b in V)
            # trace form
            TF = np.array([[int(np.trace(L[a] @ L[b])) for b in V] for a in V])
            sig = signature(TF)
            # centre dimension: x with x e_b = e_b x for all b -> linear system
            rowsC = []
            for b in V:
                # commutator [e_a, e_b] coefficients
                for c in V:
                    rowsC.append([int((L[a] @ L[b] - L[b] @ L[a])[c, 0]) for a in V])
            rankC = np.linalg.matrix_rank(np.array(rowsC, dtype=float))
            centre_dim = 4 - rankC
            # sign table
            sq = {a: (-1) ** f[(a, a)] for a in NZ}
            eps = {}
            for a in NZ:
                for b in NZ:
                    if a < b:
                        eps[(a, b)] = 1 if f[(a, b)] == f[(b, a)] else -1
            if len(set(eps.values())) != 1 or sq[1] * sq[2] * sq[3] != list(eps.values())[0]:
                rule_ok = False
            if not a_comm:
                arcs = [(a, b) for a in NZ for b in NZ if a != b and f[(a, b)] == 0]   # a -> b iff e_a e_b = + e_(a+b)
                assert len(arcs) == 3 and all((b, a) not in arcs for (a, b) in arcs)
                outdeg = {a: sum(1 for (x, y) in arcs if x == a) for a in NZ}
                if sorted(outdeg.values()) == [1, 1, 1]:
                    table = "cyclic tournament %s" % arcs
                else:
                    mid = [a for a in NZ if outdeg[a] == 1][0]
                    table = "transitive tournament %s, middle vertex %d (in S: %s)" % (arcs, mid, mid in S)
            else:
                signs = {(a, b): (-1) ** f[(a, b)] for a in NZ for b in NZ if a < b}
                table = "signed triangle %s" % signs
            rows.append((name, stat, abel, a_assoc, a_comm, sig, centre_dim, table))
        P("  S = %s (#S = %d):" % (list(S), len(S)))
        for i, r in enumerate(rows):
            P("    cocycle %d: group %s (orders 1,2,4,8: %s, abelian %s); algebra associative %s, commutative %s, trace-form signature %s, centre dimension %d"
              % (i + 1, r[0], r[1], r[2], r[3], r[4], r[5], r[6]))
            P("               sign table: %s" % r[7])
        summary.setdefault(len(S), []).append((rows[0][0], rows[0][4], rows[0][5], rows[0][6]))
    P("  rule s_a s_b s_c = eps (the same eps for the three pairs), all 16 cocycles: %s" % rule_ok)
    P("  summary by #S: %s" % {n: (len(v), sorted(set(v))) for n, v in summary.items()})
    P("  identification: commutative + trace form (4,0) -> R^4; commutative + (2,2) -> C + C; central simple with trace form")
    P("  (3,1) -> M_2(R); with (1,3) -> Hamilton (the trace form of the regular representation of M_2(R) is 2 tr(xy), signature (3,1);")
    P("  of H it is 4 Re(xy), signature (1,3); of R^r x C^s it has signature (r+s, s)).")

    # ------------------------------------------------------------- explicit models
    P("")
    P("[A'] explicit models")
    one = [Fraction(1), Fraction(0), Fraction(0), Fraction(0)]
    half = Fraction(1, 2)

    def vec(**kw):
        v = [Fraction(0)] * 4
        for kk, val in kw.items():
            v[int(kk[1:])] = Fraction(val)
        return v

    for S in sorted(Ssets, key=lambda s: (len(s), s)):
        f = cockeys[Ssets[S][0][0]]
        mul = lambda x, y: alg_mul(f, x, y)
        if len(S) == 0:
            # idempotents (1 +- e1)(1 +- e2)/4 up to signs fixed by the cocycle: search signs
            found = None
            for s1, s2, s3 in itertools.product((1, -1), repeat=3):
                e = [Fraction(1, 4), Fraction(s1, 4), Fraction(s2, 4), Fraction(s3, 4)]
                if mul(e, e) == e:
                    found = (s1, s2, s3) if found is None else found
            cnt = sum(1 for s in itertools.product((1, -1), repeat=3)
                      if mul([Fraction(1, 4)] + [Fraction(x, 4) for x in s], [Fraction(1, 4)] + [Fraction(x, 4) for x in s])
                      == [Fraction(1, 4)] + [Fraction(x, 4) for x in s])
            P("  #S = 0: primitive idempotents (1 +- e_1 +- e_2 +- e_3)/4: %d of the 8 sign choices are idempotent (R^4 needs 4)" % cnt)
        elif len(S) == 2:
            c = [a for a in NZ if a not in S][0]         # the unit of square +1
            ep = [half, 0, 0, 0]
            ep = [Fraction(x) for x in ep]
            ep[c] = half
            em = [Fraction(x) for x in [half, 0, 0, 0]]
            em[c] = -half
            ok = mul(ep, ep) == ep and mul(em, em) == em and mul(ep, em) == [0, 0, 0, 0]
            a = S[0]
            ea = [Fraction(0)] * 4
            ea[a] = Fraction(1)
            ia = mul(ep, ea)                              # an element of e+ A
            sq = mul(ia, ia)
            P("  S = %s: e+- = (1 +- e_%d)/2 orthogonal idempotents: %s; (e+ e_%d)^2 = -e+: %s  => e+ A = C, e- A = C, A = C + C"
              % (list(S), c, ok, a, sq == [-x for x in ep]))
        elif len(S) == 1:
            # 2x2 real model: e_s -> J = [[0,-1],[1,0]], another unit -> diag(1,-1), third -> +- their product
            s = S[0]
            others = [a for a in NZ if a != s]
            J = np.array([[0, -1], [1, 0]])
            Dg = np.array([[1, 0], [0, -1]])
            rep = {0: np.eye(2, dtype=int), s: J, others[0]: Dg}
            rep[others[1]] = ((-1) ** f[(s, others[0])]) * (J @ Dg)      # e_s e_o0 = (-1)^f e_o1
            ok = all((rep[a] @ rep[b] == ((-1) ** f[(a, b)]) * rep[a ^ b]).all() for a in V for b in V)
            mats = np.array([rep[a].ravel() for a in V])
            P("  S = %s: matrices e_%d -> [[0,-1],[1,0]], e_%d -> diag(1,-1): all 16 products correct: %s; the four matrices span M_2(R): rank %d"
              % (list(S), s, others[0], ok, np.linalg.matrix_rank(mats)))
        else:
            i = vec(e1=1)
            j = vec(e2=1)
            kq = mul(i, j)
            m1 = [Fraction(-1), 0, 0, 0]
            ok = mul(i, i) == m1 and mul(j, j) == m1 and mul(kq, kq) == m1 and mul(i, j) == [-x for x in mul(j, i)]
            ijk = mul(mul(i, j), kq)
            P("  #S = 3: i^2 = j^2 = (ij)^2 = -1 and ij = -ji: %s; (ij)(ij) computed %s" % (ok, [str(x) for x in ijk]))

    # ------------------------------------------------------------- zero divisors
    P("")
    P("[0] zero divisors")
    for S in sorted(Ssets, key=lambda s: (len(s), s)):
        f = cockeys[Ssets[S][0][0]]
        L = algebra(f)
        # exact: det of left multiplication by x = sum x_a e_a, evaluated on the integer grid {-3..3}^4
        zero = None
        minabs = None
        allpos = True
        normsq_ok = True
        for x in itertools.product(range(-3, 4), repeat=4):
            if not any(x):
                continue
            M = sum(int(xa) * L[a] for a, xa in zip(V, x))
            d = int(round(np.linalg.det(M.astype(float))))
            if d == 0 and zero is None:
                # find y with x y = 0
                zero = x
            if d <= 0:
                allpos = False
            if d != (sum(t * t for t in x)) ** 2:
                normsq_ok = False
        if zero is not None:
            M = sum(int(xa) * L[a] for a, xa in zip(V, zero))
            # kernel vector with small integer coordinates
            y = None
            for cand in itertools.product(range(-2, 3), repeat=4):
                if any(cand) and not (M @ np.array(cand)).any():
                    y = cand
                    break
            P("  S = %s: zero divisor x = %s, y = %s with x y = 0" % (list(S), zero, y))
        else:
            P("  S = %s: det L_x > 0 on all 7^4 - 1 nonzero grid points: %s; det L_x = (x_0^2+x_1^2+x_2^2+x_3^2)^2 on the grid (a degree-4 identity, 7 > 4 nodes per variable, hence a polynomial identity): %s"
              % (list(S), allpos, normsq_ok))

    # ------------------------------------------------------------- GL_2(F_2)
    P("")
    P("[P] GL_2(F_2) = S_3 permutes the three points; classes with the same #S form one orbit: sizes %s"
      % {n: len([s for s in Ssets if len(s) == n]) for n in range(4)})

    sys.stdout.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
