#!/usr/bin/env python3
"""Audit of THM-4665, item 3: independent recount of faithful EI17 realizations, by direct subgraph search.

Method B (main): search for injective homomorphisms EI17 -> Cay(A, H) with the gauge v17 = 0, v11 = -1, vertex by
vertex along a Henneberg-type order (two 'free' placements, every other vertex adjacent to >= 2 placed ones), checking
every edge as soon as both ends are placed and injectivity at every step.  Nothing from THM-4665 is used: no
parallelogram forcing, no direction classes, no relations, no quadratic formula.
    split model:   A = F_p x F_p, H = {(x, 1/x)}, N(z, w) = z w;
    twisted model: A = F_p[t]/(t^2 - t - b) (a different model of F_(p^2) from the S39 script's F_p(sqrt n)),
                   H = {norm one} = mu_(p+1), N(a0 + a1 t) = a0^2 + a0 a1 - b a1^2.
Apex step (vertex with placed neighbours a, b): v = a + u with u in H and v - b in H.  Two ways:
    'table': all pairs (u, w) in H x H are tabulated by c = u + w (the table asserts at most two u per c != 0,
             an independent check of the codegree lemma in both models);
    'brute': scan all of H (no table, no lemma), for small p.
Method C (cross-check of the direction parametrization in a different elimination order): free U2, U3; U6 from
|U10| = 1, U1 from |U12| = 1, check |U7| = 1, U0 from |U9| = 1, U4 from |U8| = 1, check |U11| = 1; degenerate
(c = 0) steps are enumerated in full, not skipped; faithful = 17 distinct VERT points.
Usage: python -B a2_embedding_search.py MODEL METHOD PRIMES     MODEL in {split, twisted}, METHOD in {table, brute, C}
       PRIMES comma-separated, or a range a-b (all primes a <= p < b).
"""
import os
import sys
import time

import numpy as np
from numba import njit, prange
from sympy import isprime, primerange

PE = [(1, 4), (1, 7), (1, 8), (2, 3), (2, 6), (2, 10), (3, 5), (3, 7), (3, 12), (4, 6), (4, 13), (4, 14), (5, 8),
      (5, 9), (5, 14), (6, 11), (6, 16), (7, 9), (7, 15), (8, 13), (8, 15), (9, 10), (9, 11), (10, 12), (10, 16),
      (11, 12), (11, 17), (13, 17), (14, 15), (14, 17), (16, 17)]
ADJ = {v: set() for v in range(1, 18)}
for i, j in PE:
    ADJ[i].add(j); ADJ[j].add(i)


# ----------------------------------------------------------------------------- orders
def closure(placed, steps):
    placed = list(placed)
    while True:
        best = None
        for v in range(1, 18):
            if v in placed:
                continue
            nb = [w for w in placed if w in ADJ[v]]
            if len(nb) >= 2 and (best is None or len(nb) > best[1]):
                best = (v, len(nb), nb)
        if best is None:
            return placed, steps
        v, _, nb = best
        steps.append((v, 1, nb[0], nb[1], nb[2:]))
        placed.append(v)


def orders_with_two_free():
    out = []
    for v1 in range(1, 18):
        for a1 in (11, 17):
            if v1 in (11, 17) or a1 not in ADJ[v1]:
                continue
            pl, st = closure([11, 17, v1], [(v1, 0, a1, -1, [])])
            for v2 in range(1, 18):
                if v2 in pl:
                    continue
                for a2 in sorted(ADJ[v2]):
                    if a2 not in pl:
                        continue
                    if len([w for w in pl if w in ADJ[v2]]) != 1:
                        continue
                    pl2, st2 = closure(pl + [v2], st + [(v2, 0, a2, -1, [])])
                    if len(pl2) == 17:
                        out.append(st2)
    return out


def pick_order():
    # S24 used the free placements 6 (from 11) and 2 (from 6); take an order with a different free pair
    cands = orders_with_two_free()
    for st in cands:
        free = [(s[0], s[2]) for s in st if s[1] == 0]
        if free[0][0] not in (6, 2) and free[1][0] not in (6, 2):
            return st, len(cands)
    return cands[0], len(cands)


def order_arrays(st):
    n = len(st)
    sv = np.array([s[0] for s in st], dtype=np.int64)
    sk = np.array([s[1] for s in st], dtype=np.int64)
    s1 = np.array([s[2] for s in st], dtype=np.int64)
    s2 = np.array([s[3] for s in st], dtype=np.int64)
    sc = -np.ones((n, 4), dtype=np.int64); nc = np.zeros(n, dtype=np.int64)
    for k, s in enumerate(st):
        for t, w in enumerate(s[4]):
            sc[k, t] = w
        nc[k] = len(s[4])
    # every edge is used exactly once: base edge (11,17) + parents + checks
    used = {(11, 17)}
    for s in st:
        for w in [s[2]] + ([s[3]] if s[1] == 1 else []) + list(s[4]):
            used.add((min(s[0], w), max(s[0], w)))
    assert used == set(PE), 'order does not cover the edge set'
    return sv, sk, s1, s2, sc, nc


# ----------------------------------------------------------------------------- models
@njit(cache=True)
def nrm(model, p, b, x0, x1):
    x0 %= p; x1 %= p
    if model == 0:
        return x0 * x1 % p
    return (x0 * x0 + x0 * x1 - b * (x1 * x1 % p)) % p


@njit(cache=True)
def units(model, p, b):
    if model == 0:
        u0 = np.zeros(p - 1, np.int64); u1 = np.zeros(p - 1, np.int64)
        for x in range(1, p):
            u0[x - 1] = x
            # inverse by exhaustive search (no Fermat), independent of the scripts' pow
            for y in range(1, p):
                if x * y % p == 1:
                    u1[x - 1] = y
                    break
        return u0, u1
    cnt = 0
    for a0 in range(p):
        for a1 in range(p):
            if nrm(model, p, b, a0, a1) == 1:
                cnt += 1
    u0 = np.zeros(cnt, np.int64); u1 = np.zeros(cnt, np.int64)
    k = 0
    for a0 in range(p):
        for a1 in range(p):
            if nrm(model, p, b, a0, a1) == 1:
                u0[k] = a0; u1[k] = a1; k += 1
    return u0, u1


@njit(cache=True)
def build_table(p, u0, u1):
    nu = u0.shape[0]
    T = -np.ones((p * p, 2), np.int16)
    overflow = 0
    for i in range(nu):
        for j in range(nu):
            c0 = (u0[i] + u0[j]) % p; c1 = (u1[i] + u1[j]) % p
            if c0 == 0 and c1 == 0:
                continue
            idx = c0 * p + c1
            if T[idx, 0] == -1:
                T[idx, 0] = i
            elif T[idx, 1] == -1:
                T[idx, 1] = i
            else:
                overflow += 1
    return T, overflow


# ----------------------------------------------------------------------------- method B
@njit(cache=True)
def search_one(first, model, p, b, u0, u1, T, brute, sv, sk, s1, s2, sc, nc, rec, maxrec):
    nst = sv.shape[0]; nu = u0.shape[0]
    pos0 = np.zeros(18, np.int64); pos1 = np.zeros(18, np.int64)
    pos0[11] = p - 1
    pos1[11] = (p - 1) if model == 0 else 0
    width = nu if brute else 2
    cand = np.zeros((nst, width), np.int64); ncand = np.zeros(nst, np.int64); ptr = np.zeros(nst, np.int64)
    count = 0
    d = 0
    v = sv[0]; a = s1[0]
    pos0[v] = (pos0[a] + u0[first]) % p; pos1[v] = (pos1[a] + u1[first]) % p
    if (pos0[v] == pos0[11] and pos1[v] == pos1[11]) or (pos0[v] == 0 and pos1[v] == 0):
        return 0
    d = 1
    ptr[1] = 0
    # prepare depth 1
    prep = True
    while d >= 1:
        if d == nst:
            if count < maxrec:
                for w in range(1, 18):
                    rec[count, w - 1, 0] = pos0[w]; rec[count, w - 1, 1] = pos1[w]
            count += 1
            d -= 1
            prep = False
            continue
        if prep:
            ptr[d] = 0
            if sk[d] == 1:
                a = s1[d]; bb = s2[d]
                c0 = (pos0[bb] - pos0[a]) % p; c1 = (pos1[bb] - pos1[a]) % p
                n = 0
                if brute:
                    for k in range(nu):
                        if nrm(model, p, b, pos0[a] + u0[k] - pos0[bb], pos1[a] + u1[k] - pos1[bb]) == 1:
                            cand[d, n] = k; n += 1
                else:
                    idx = c0 * p + c1
                    for s in range(2):
                        if T[idx, s] >= 0:
                            cand[d, n] = T[idx, s]; n += 1
                ncand[d] = n
            else:
                ncand[d] = nu
            prep = False
        if ptr[d] >= ncand[d]:
            d -= 1
            continue
        k = ptr[d]; ptr[d] += 1
        uidx = k if sk[d] == 0 else cand[d, k]
        v = sv[d]; a = s1[d]
        x0 = (pos0[a] + u0[uidx]) % p; x1 = (pos1[a] + u1[uidx]) % p
        ok = True
        if (x0 == pos0[11] and x1 == pos1[11]) or (x0 == 0 and x1 == 0):
            ok = False
        if ok:
            for e in range(d):
                w = sv[e]
                if x0 == pos0[w] and x1 == pos1[w]:
                    ok = False
                    break
        if ok and sk[d] == 1 and brute == 0:
            # the table guarantees the two parent edges; re-check them anyway
            if nrm(model, p, b, x0 - pos0[s2[d]], x1 - pos1[s2[d]]) != 1 or nrm(model, p, b, x0 - pos0[a], x1 - pos1[a]) != 1:
                ok = False
        if ok:
            for t in range(nc[d]):
                w = sc[d, t]
                if nrm(model, p, b, x0 - pos0[w], x1 - pos1[w]) != 1:
                    ok = False
                    break
        if not ok:
            continue
        pos0[v] = x0; pos1[v] = x1
        d += 1
        prep = True
    return count


@njit(parallel=True, cache=True)
def search_all(model, p, b, u0, u1, T, brute, sv, sk, s1, s2, sc, nc):
    nu = u0.shape[0]
    counts = np.zeros(nu, np.int64)
    for first in prange(nu):
        rec = np.zeros((1, 17, 2), np.int64)
        counts[first] = search_one(first, model, p, b, u0, u1, T, brute, sv, sk, s1, s2, sc, nc, rec, 0)
    return counts


# ----------------------------------------------------------------------------- method C
VERT = {1: (-1, 0, -1, 0, -1, 0, -1), 2: (0, 0, -1, -1, -1, 0, -1), 3: (0, 0, -1, -1, 0, 0, -1), 4: (-1, 1, -1, 0, -1, 0, -1),
        5: (0, 0, -1, -1, 0, 0, 0), 6: (0, 0, -1, -1, -1, 1, -1), 7: (-1, 0, 0, 0, -1, 0, -1), 8: (-1, 0, -1, 1, -1, 0, -1),
        9: (-1, 0, 0, 0, -1, 0, 0), 10: (1, 0, -1, -1, -1, 0, -1), 11: (-1, 0, 0, 0, 0, 0, 0), 12: (1, 0, -1, -1, 0, 0, -1),
        13: (-1, 1, -1, 1, -1, 0, -1), 14: (0, 0, 0, -1, 0, 0, 0), 15: (-1, 0, 0, 1, -1, 0, -1), 16: (1, 0, -1, -1, -1, 1, -1),
        17: (0, 0, 0, 0, 0, 0, 0)}
VM = np.array([VERT[v] for v in range(1, 18)], dtype=np.int64)   # verified against a1_structure.py's derivation


@njit(cache=True)
def sols_for(model, p, b, u0, u1, T, c0, c1, out):
    """unit indices k with c - u_k in H; all units if c = 0."""
    c0 %= p; c1 %= p
    nu = u0.shape[0]
    if c0 == 0 and c1 == 0:
        for k in range(nu):
            out[k] = k
        return nu
    idx = c0 * p + c1
    n = 0
    for s in range(2):
        if T[idx, s] >= 0:
            out[n] = T[idx, s]; n += 1
    return n


@njit(parallel=True, cache=True)
def method_c(model, p, b, u0, u1, T, VM):
    nu = u0.shape[0]
    res = np.zeros((nu, 2), np.int64)      # faithful, all
    gm0 = p - 1; gm1 = (p - 1) if model == 0 else 0   # the element -1
    for i2 in prange(nu):
        o6 = np.zeros(nu, np.int64); o1 = np.zeros(nu, np.int64); o0 = np.zeros(nu, np.int64); o4 = np.zeros(nu, np.int64)
        z0 = np.zeros(17, np.int64); z1 = np.zeros(17, np.int64)
        fa = 0; al = 0
        for i3 in range(nu):
            # U10 = (2U2 - U3 - 1) - U6
            n6 = sols_for(model, p, b, u0, u1, T, 2 * u0[i2] - u0[i3] + gm0, 2 * u1[i2] - u1[i3] + gm1, o6)
            for a6 in range(n6):
                i6 = o6[a6]
                # U12 = (2 - U2 - U6) - U1      (2 = -2 * (-1))
                n1 = sols_for(model, p, b, u0, u1, T, -2 * gm0 - u0[i2] - u0[i6], -2 * gm1 - u1[i2] - u1[i6], o1)
                for a1 in range(n1):
                    i1 = o1[a1]
                    # check U7 = U1 + U2 - U3 - 1
                    if nrm(model, p, b, u0[i1] + u0[i2] - u0[i3] + gm0, u1[i1] + u1[i2] - u1[i3] + gm1) != 1:
                        continue
                    # U9 = (U1 - U2 + U3 + U6 + 1) - U0
                    n0 = sols_for(model, p, b, u0, u1, T, u0[i1] - u0[i2] + u0[i3] + u0[i6] - gm0,
                                  u1[i1] - u1[i2] + u1[i3] + u1[i6] - gm1, o0)
                    for a0 in range(n0):
                        i0 = o0[a0]
                        # U8 = U4 - (U0 + U2 - 1): U4 - e in H iff e - U4 in H (H = -H)
                        n4 = sols_for(model, p, b, u0, u1, T, u0[i0] + u0[i2] + gm0, u1[i0] + u1[i2] + gm1, o4)
                        for a4 in range(n4):
                            i4 = o4[a4]
                            # check U11 = U1 + U2 + U3 - U4 + U6 - 1
                            if nrm(model, p, b, u0[i1] + u0[i2] + u0[i3] - u0[i4] + u0[i6] + gm0,
                                   u1[i1] + u1[i2] + u1[i3] - u1[i4] + u1[i6] + gm1) != 1:
                                continue
                            al += 1
                            # vertices: c0 * 1 + c1 U0 + c2 U1 + c3 U2 + c4 U3 + c5 U4 + c6 U6 ; 1 = -(-1)
                            for v in range(17):
                                z0[v] = (-VM[v, 0] * gm0 + VM[v, 1] * u0[i0] + VM[v, 2] * u0[i1] + VM[v, 3] * u0[i2]
                                         + VM[v, 4] * u0[i3] + VM[v, 5] * u0[i4] + VM[v, 6] * u0[i6]) % p
                                z1[v] = (-VM[v, 0] * gm1 + VM[v, 1] * u1[i0] + VM[v, 2] * u1[i1] + VM[v, 3] * u1[i2]
                                         + VM[v, 4] * u1[i3] + VM[v, 5] * u1[i4] + VM[v, 6] * u1[i6]) % p
                            ok = True
                            for v in range(17):
                                for w in range(v + 1, 17):
                                    if z0[v] == z0[w] and z1[v] == z1[w]:
                                        ok = False
                            if ok:
                                fa += 1
        res[i2, 0] = fa; res[i2, 1] = al
    return res


# ----------------------------------------------------------------------------- driver
def setup(model, p):
    b = 0
    if model == 1:
        b = next(x for x in range(1, p) if pow((1 + 4 * x) % p, (p - 1) // 2, p) == p - 1)
    u0, u1 = units(model, p, b)
    assert u0.shape[0] == (p - 1 if model == 0 else p + 1), ('unit count', p, u0.shape[0])
    return b, u0, u1


def run(model, method, p, st_arrays):
    t0 = time.time()
    b, u0, u1 = setup(model, p)
    if method in ('table', 'C'):
        T, ov = build_table(p, u0, u1)
        assert ov == 0, ('codegree > 2 somewhere', p, ov)
    else:
        T = -np.ones((1, 2), np.int16)
    if method == 'C':
        r = method_c(model, p, b, u0, u1, T, VM)
        f, al = int(r[:, 0].sum()), int(r[:, 1].sum())
        return f, al, time.time() - t0
    sv, sk, s1, s2, sc, nc = st_arrays
    cnt = search_all(model, p, b, u0, u1, T, 1 if method == 'brute' else 0, sv, sk, s1, s2, sc, nc)
    return int(cnt.sum()), -1, time.time() - t0


def solutions(model, p, st_arrays, maxrec=64):
    """all faithful embeddings (positions) for a prime, by method B with the table."""
    b, u0, u1 = setup(model, p)
    T, ov = build_table(p, u0, u1)
    assert ov == 0
    sv, sk, s1, s2, sc, nc = st_arrays
    cnt = search_all(model, p, b, u0, u1, T, 0, sv, sk, s1, s2, sc, nc)
    out = []
    for first in np.nonzero(cnt)[0]:
        rec = np.zeros((maxrec, 17, 2), np.int64)
        n = search_one(int(first), model, p, b, u0, u1, T, 0, sv, sk, s1, s2, sc, nc, rec, maxrec)
        for k in range(min(n, maxrec)):
            out.append(rec[k].copy())
    return b, out


def parse_primes(s):
    out = []
    for tok in s.split(','):
        if '-' in tok:
            a, c = tok.split('-')
            out += [p for p in primerange(int(a), int(c))]
        else:
            out.append(int(tok))
    return out


if __name__ == '__main__':
    model = {'split': 0, 'twisted': 1}[sys.argv[1]]
    method = sys.argv[2]
    plist = parse_primes(sys.argv[3])
    st, ncand = pick_order()
    print('# %d Henneberg orders with two free placements after the gauge edge (11,17); using:' % ncand)
    for s in st:
        print('#   %s' % ('free %2d from %2d' % (s[0], s[2]) if s[1] == 0 else 'apex %2d from (%2d,%2d) checks %s' % (s[0], s[2], s[3], s[4])))
    arrs = order_arrays(st)
    print('# model %s, method %s' % (sys.argv[1], method))
    print('    p  p%%55  faithful  %s  seconds' % ('all(C)' if method == 'C' else ''))
    for p in plist:
        assert isprime(p) and p > 2
        f, al, dt = run(model, method, p, arrs)
        print('%6d  %3d  %6d  %s  %.1f' % (p, p % 55, f, ('%6d' % al) if method == 'C' else '', dt), flush=True)
