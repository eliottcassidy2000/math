#!/usr/bin/env python3
"""S93 lane "seeds", deliverable (a): brute-force verification of the seed formulas.

V1  Q, U (subset recursion of seedlib) = literal enumeration (all vertex orders, all cut points), every labelled
    tournament of order <= 5.
V2  closing law H(B => TT_t => v, v -> O) = r + 2^t M, r = H(B) - W_B(O), M = U_B(O):
    (i) every labelled B of order <= 4, every O (empty included), t = 1, 2, 3, H of the closed tournament by literal
        enumeration of all vertex orders; (ii) every labelled B of order 5 and 6000 random B of order 6, 7, every O,
        t = 1..4, by an independent subset recursion.
V3  |O| = 2, e -> f in B, B' = B with that arc reversed:  r(B, {e, f}) = H(B') - 2 W_B'(f -> e), where W_T(x -> y) is
    the number of Hamiltonian paths of T through the arc x -> y; every labelled B of order <= 6, literal counts.
    Also r = A - C with A = number of Hamiltonian paths of B in which e, f are not consecutive, C = Q_B(f, e).
V4  F1: U_B(o) = H(B) + sum over out-neighbours x of o of Q_B(x, o); sink rule (r, 2M); source rule.
V5  one-apex closed forms (spine z_1 -> ... -> z_n, apex a, a -> z_i iff y_i = 1), all words of length <= 10:
    H(B) = [y_1 = 1] + [y_n = 0] + #{j : y_j = 0, y_(j+1) = 1} + sum_{i<k, y_i=1, y_k=0} 2^(k-i-1)
    U_B(a) = 1 + sum_i y_i 2^(n-i)
    W_{a z_j} = sum_{i<j, y_i=1} 2^(j-1-i) + sum_{k>j, y_k=0} 2^(k-j-1) + [j = n or y_(j+1) = 1] + [j = 1 or y_(j-1) = 0]
V6  kappa identity: for a pair {e, f} of X with e -> f:  Pi_2(X) - 2 W_ef - 2 W3_ef = 2 Pi_2(X - ef) - Pi_2(X'),
    X' = X with the arc reversed, Pi_2(X - ef) = covers (A, B) by two possibly empty paths in which e, f are not
    consecutive; every labelled X of order <= 5, literal enumeration.
"""
import itertools
import sys

import numpy as np
from numba import njit

from seedlib import dp_tables, pair_mats, subset_sums, hpaths, closed_out


@njit(cache=True)
def next_perm(p, n):
    i = n - 2
    while i >= 0 and p[i] > p[i + 1]:
        i -= 1
    if i < 0:
        return False
    j = n - 1
    while p[j] < p[i]:
        j -= 1
    p[i], p[j] = p[j], p[i]
    a, b = i + 1, n - 1
    while a < b:
        p[a], p[b] = p[b], p[a]
        a += 1
        b -= 1
    return True


@njit(cache=True)
def literal(n, out, Q, U, arcuse, H3):
    """literal enumeration: H; Q[e, f]; U[o]; arcuse[x, y] = Hamiltonian paths through the arc x -> y;
    returns (H, Pi2, P3count) and fills H3[e, f] = triples (D1, D2, D3) nonempty, D1 ends e, D2 starts f."""
    for e in range(n):
        U[e] = 0
        for f in range(n):
            Q[e, f] = 0
            arcuse[e, f] = 0
            H3[e, f] = 0
    p = np.arange(n)
    H = 0
    pi2 = 0
    good = np.zeros(n, np.bool_)          # good[i]: p[i] -> p[i+1]
    while True:
        for i in range(n - 1):
            good[i] = ((out[p[i]] >> p[i + 1]) & 1) == 1
        # prefix-good and suffix-good
        allgood = True
        for i in range(n - 1):
            if not good[i]:
                allgood = False
        if allgood:
            H += 1
            for i in range(n - 1):
                arcuse[p[i], p[i + 1]] += 1
        # cut after position k-1 (first part p[0..k-1], second p[k..n-1]); k = 0: first part empty
        for k in range(0, n + 1):
            ok = True
            for i in range(n - 1):
                if i != k - 1 and not good[i]:
                    ok = False
                    break
            if ok:
                pi2 += 1
                if k < n:
                    U[p[k]] += 1
                    if k >= 1:
                        Q[p[k - 1], p[k]] += 1
        # two cuts k < l: D1 = p[0..k-1], D2 = p[k..l-1], D3 = p[l..n-1], all nonempty
        for k in range(1, n - 1):
            for l in range(k + 1, n):
                ok = True
                for i in range(n - 1):
                    if i != k - 1 and i != l - 1 and not good[i]:
                        ok = False
                        break
                if ok:
                    H3[p[k - 1], p[k]] += 1
        if not next_perm(p, n):
            break
    return H, pi2


@njit(cache=True)
def literal_H(n, out):
    p = np.arange(n)
    H = 0
    while True:
        ok = True
        for i in range(n - 1):
            if ((out[p[i]] >> p[i + 1]) & 1) == 0:
                ok = False
                break
        if ok:
            H += 1
        if not next_perm(p, n):
            break
    return H


@njit(cache=True)
def literal_pi2_avoid(n, out, e, f):
    """covers (A, B) by two possibly empty paths in which e and f are not consecutive (either order)."""
    p = np.arange(n)
    tot = 0
    while True:
        for k in range(0, n + 1):
            ok = True
            for i in range(n - 1):
                if i != k - 1:
                    if ((out[p[i]] >> p[i + 1]) & 1) == 0:
                        ok = False
                        break
                    if (p[i] == e and p[i + 1] == f) or (p[i] == f and p[i + 1] == e):
                        ok = False
                        break
            if ok:
                tot += 1
        if not next_perm(p, n):
            break
    return tot


def labelled(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for code in range(1 << len(pairs)):
        out = np.zeros(n, np.int64)
        for b, (i, j) in enumerate(pairs):
            if (code >> b) & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
        yield out


def state_tables(n, out):
    fe = np.zeros((1 << n, n), np.int64)
    gs = np.zeros((1 << n, n), np.int64)
    hS = np.zeros(1 << n, np.int64)
    Q = np.zeros((n, n), np.int64)
    U = np.zeros(n, np.int64)
    dp_tables(n, out, fe, gs, hS)
    pair_mats(n, fe, gs, hS, Q, U)
    Wsub = np.zeros(1 << n, np.int64)
    Msub = np.zeros(1 << n, np.int64)
    subset_sums(n, Q, U, Wsub, Msub)
    return fe, gs, hS, Q, U, Wsub, Msub


def main():
    rng = np.random.default_rng(93)
    allok = True

    # V1
    cnt = 0
    ok = True
    for n in range(1, 6):
        for out in labelled(n):
            fe, gs, hS, Q, U, Wsub, Msub = state_tables(n, out)
            Q2 = np.zeros((n, n), np.int64)
            U2 = np.zeros(n, np.int64)
            au = np.zeros((n, n), np.int64)
            H3 = np.zeros((n, n), np.int64)
            H, pi2 = literal(n, out, Q2, U2, au, H3)
            full = (1 << n) - 1
            pi2dp = sum(int(hS[S]) * int(hS[full ^ S]) for S in range(full + 1))
            if H != hS[full] or not (Q == Q2).all() or not (U == U2).all() or pi2 != pi2dp:
                ok = False
            cnt += 1
    print(f"V1  Q, U, H, Pi_2: subset recursion = literal enumeration on {cnt} labelled tournaments of order <= 5: {ok}")
    allok &= ok

    # V2 (i)
    cnt = 0
    ok = True
    for n in range(1, 5):
        for out in labelled(n):
            fe, gs, hS, Q, U, Wsub, Msub = state_tables(n, out)
            full = (1 << n) - 1
            h = int(hS[full])
            for O in range(full + 1):
                r = h - int(Wsub[O])
                M = int(Msub[O])
                for t in (1, 2, 3):
                    res = np.zeros(n + t + 1, np.int64)
                    N = closed_out(n, out, O, t, res)
                    if literal_H(N, res) != r + (1 << t) * M:
                        ok = False
                    cnt += 1
    print(f"V2i closing law, literal enumeration of all vertex orders, labelled B of order <= 4, every O, t = 1, 2, 3: {cnt} cases: {ok}")
    allok &= ok

    # V2 (ii)
    cnt = 0
    ok = True
    todo = list(labelled(5))
    for n in (6, 7):
        for _ in range(3000):
            code = int(rng.integers(0, 1 << (n * (n - 1) // 2)))
            out = np.zeros(n, np.int64)
            b = 0
            for i in range(n):
                for j in range(i + 1, n):
                    if (code >> b) & 1:
                        out[i] |= 1 << j
                    else:
                        out[j] |= 1 << i
                    b += 1
            todo.append(out)
    for out in todo:
        n = len(out)
        fe, gs, hS, Q, U, Wsub, Msub = state_tables(n, out)
        full = (1 << n) - 1
        h = int(hS[full])
        Os = range(full + 1) if n == 5 else [int(x) for x in rng.integers(0, full + 1, 6)]
        for O in Os:
            r = h - int(Wsub[O])
            M = int(Msub[O])
            for t in (1, 2, 3, 4):
                res = np.zeros(n + t + 1, np.int64)
                N = closed_out(n, out, O, t, res)
                if hpaths(N, res) != r + (1 << t) * M:
                    ok = False
                cnt += 1
    print(f"V2ii closing law, independent subset recursion, all labelled B of order 5 (every O) and 6000 random B of order 6, 7: {cnt} cases: {ok}")
    allok &= ok

    # V3
    cnt = 0
    ok = True
    ok2 = True
    for n in range(2, 7):
        gen = labelled(n) if n <= 5 else (o for k, o in enumerate(labelled(n)) if k % 7 == 0)
        for out in gen:
            fe, gs, hS, Q, U, Wsub, Msub = state_tables(n, out)
            full = (1 << n) - 1
            h = int(hS[full])
            Q2 = np.zeros((n, n), np.int64)
            U2 = np.zeros(n, np.int64)
            au = np.zeros((n, n), np.int64)
            H3 = np.zeros((n, n), np.int64)
            literal(n, out, Q2, U2, au, H3)
            for e in range(n):
                for f in range(n):
                    if e != f and (out[e] >> f) & 1:
                        r = h - int(Q[e, f]) - int(Q[f, e])
                        outp = out.copy()
                        outp[e] &= ~(1 << f)
                        outp[f] |= 1 << e
                        Qp = np.zeros((n, n), np.int64)
                        Up = np.zeros(n, np.int64)
                        aup = np.zeros((n, n), np.int64)
                        H3p = np.zeros((n, n), np.int64)
                        Hp, _ = literal(n, outp, Qp, Up, aup, H3p)
                        if r != Hp - 2 * int(aup[f, e]):
                            ok = False
                        # r = A - C
                        A = h - int(au[e, f])
                        if r != A - int(Q[f, e]) or int(Q[e, f]) != int(au[e, f]) or int(Q[f, e]) != int(aup[f, e]):
                            ok2 = False
                        cnt += 1
    print(f"V3  r(B, {{e, f}}) = H(B') - 2 W_B'(f -> e) (literal), labelled B of order <= 5 and every 7th of order 6, every arc: {cnt} cases: {ok};  r = A - Q(f, e), Q(e, f) = W_B(e -> f): {ok2}")
    allok &= ok and ok2

    # V4
    cnt = 0
    ok = True
    for n in range(1, 6):
        for out in labelled(n):
            fe, gs, hS, Q, U, Wsub, Msub = state_tables(n, out)
            full = (1 << n) - 1
            h = int(hS[full])
            for o in range(n):
                s = h + sum(int(Q[x, o]) for x in range(n) if (out[o] >> x) & 1)
                if s != int(U[o]):
                    ok = False
            # sink and source rules
            outs = np.zeros(n + 1, np.int64)
            outq = np.zeros(n + 1, np.int64)
            for v in range(n):
                outs[v] = out[v] | (1 << n)       # new vertex n is a sink
                outq[v] = out[v]
            outq[n] = full                         # new vertex n is a source
            a = state_tables(n + 1, outs)
            b = state_tables(n + 1, outq)
            for O in range(1, full + 1):
                r, M = h - int(Wsub[O]), int(Msub[O])
                if (int(a[2][2 * full + 1]) - int(a[5][O]), int(a[6][O])) != (r, 2 * M):
                    ok = False
                if (int(b[2][2 * full + 1]) - int(b[5][O]), int(b[6][O])) != (r, M):
                    ok = False
                cnt += 1
    print(f"V4  F1, sink rule (r, 2M), source rule (r, M): labelled B of order <= 5, every nonempty O: {cnt} cases: {ok}")
    allok &= ok

    # V5
    cnt = 0
    ok = True
    for n in range(1, 11):
        for word in range(1 << n):
            y = [(word >> (n - 1 - i)) & 1 for i in range(n)]      # y[0] = y_1
            out = np.zeros(n + 1, np.int64)
            for i in range(n):
                for j in range(i + 1, n):
                    out[i] |= 1 << j
                if y[i]:
                    out[n] |= 1 << i
                else:
                    out[i] |= 1 << n
            fe, gs, hS, Q, U, Wsub, Msub = state_tables(n + 1, out)
            full = (1 << (n + 1)) - 1
            Hf = (1 if y[0] == 1 else 0) + (1 if y[n - 1] == 0 else 0)
            Hf += sum(1 for j in range(n - 1) if y[j] == 0 and y[j + 1] == 1)
            Hf += sum(1 << (k - i - 1) for i in range(n) for k in range(i + 1, n) if y[i] == 1 and y[k] == 0)
            if Hf != int(hS[full]):
                ok = False
            if int(U[n]) != 1 + word:
                ok = False
            for j in range(n):
                w = sum(1 << (j - 1 - i) for i in range(j) if y[i] == 1)
                w += sum(1 << (k - j - 1) for k in range(j + 1, n) if y[k] == 0)
                w += 1 if (j == n - 1 or y[j + 1] == 1) else 0
                w += 1 if (j == 0 or y[j - 1] == 0) else 0
                if w != int(Q[n, j]) + int(Q[j, n]):
                    ok = False
            cnt += 1
    print(f"V5  one-apex closed forms for H, U(apex) = 1 + val(y), W(apex, z_j): all {cnt} words of length <= 10: {ok}")
    allok &= ok

    # V6
    cnt = 0
    ok = True
    for n in range(2, 6):
        for out in labelled(n):
            Q2 = np.zeros((n, n), np.int64)
            U2 = np.zeros(n, np.int64)
            au = np.zeros((n, n), np.int64)
            H3 = np.zeros((n, n), np.int64)
            H, pi2 = literal(n, out, Q2, U2, au, H3)
            for e in range(n):
                for f in range(n):
                    if e != f and (out[e] >> f) & 1:
                        kap = pi2 - 2 * int(Q2[e, f] + Q2[f, e]) - 2 * int(H3[e, f] + H3[f, e])
                        outp = out.copy()
                        outp[e] &= ~(1 << f)
                        outp[f] |= 1 << e
                        Qp = np.zeros((n, n), np.int64)
                        Up = np.zeros(n, np.int64)
                        aup = np.zeros((n, n), np.int64)
                        H3p = np.zeros((n, n), np.int64)
                        Hp, pi2p = literal(n, outp, Qp, Up, aup, H3p)
                        if kap != 2 * literal_pi2_avoid(n, out, e, f) - pi2p:
                            ok = False
                        cnt += 1
    print(f"V6  kappa_X({{e, f}}) = 2 Pi_2(X - ef) - Pi_2(X'): labelled X of order <= 5, every arc: {cnt} cases: {ok}")
    allok &= ok
    print("ALL CHECKS PASS" if allok else "SOME CHECK FAILED")
    if not allok:
        sys.exit(1)


if __name__ == "__main__":
    main()
