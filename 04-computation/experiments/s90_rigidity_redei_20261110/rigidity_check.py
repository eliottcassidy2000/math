#!/usr/bin/env python3
"""THM-5124 (additive rigidity of the Hamiltonian-path count on presence cubes): exhaustive and sampled checks.

Notation: T a tournament on V, Hs[mask] = number of Hamiltonian paths of T[mask] (Hs[0] = 1).

(P)  Lemma P.  For a pair {a, b} and R = V - {a, b}:
        [Hs[R+a] > Hs[R], Hs[R+b] > Hs[R], Hs[V] - Hs[R+a] - Hs[R+b] + Hs[R] = 0]
     iff  T = Y => Q => X with Q = {x0, y0, a, b} the strong 4-tournament whose two one-triangle vertices are a, b
          (x0 -> y0, y0 -> a, y0 -> b, a -> x0, b -> x0).           Universe: every labelled tournament, n <= 7.
(R)  Theorem R.  For every ordered partition V = B + W1 + W2 (W1, W2 nonempty) with both modules active over B and
     zero interaction: values (c, 3c, 3c, 5c), exactly one active vertex in each module, they form the Lemma P
     structure over B, and every other module vertex is count-inert in T.   Universe: n <= 6.
(K)  No additive cube of dimension 3: no partition V = B + W1 + W2 + W3 with f affine on the 8 corners and all
     three increments positive.  Universe: n <= 6 (modules), n <= 7 (single vertices, every base).
(S)  Sampled random tournaments of orders 8-10: Lemma P equivalence on every pair; Theorem R on random partitions.

Exit status 0 iff every check passes.  Run:  python3 rigidity_check.py > rigidity_check.out
"""
import sys
import time

import numpy as np
from numba import njit


@njit(cache=True)
def hs_table(n, beats, cnt, Hs):
    size = 1 << n
    for m in range(size):
        for v in range(n):
            cnt[m, v] = 0
    for v in range(n):
        cnt[1 << v, v] = 1
    Hs[0] = 1
    for mask in range(1, size):
        tot = 0
        for v in range(n):
            c = cnt[mask, v]
            if c:
                tot += c
                nxt = beats[v] & ~mask
                idx = 0
                while nxt:
                    if nxt & 1:
                        cnt[mask | (1 << idx), idx] += c
                    nxt >>= 1
                    idx += 1
        Hs[mask] = tot


@njit(cache=True)
def structure_P(n, beats, a, b, Rmask):
    """T[Rmask + a + b] = Y => Q => X, Q = {x0, y0, a, b}, x0->y0, y0->a, y0->b, a->x0, b->x0."""
    full = Rmask | (1 << a) | (1 << b)
    for x0 in range(n):
        if not (Rmask >> x0) & 1:
            continue
        for y0 in range(n):
            if y0 == x0 or not (Rmask >> y0) & 1:
                continue
            if not ((beats[x0] >> y0) & 1 and (beats[y0] >> a) & 1 and (beats[y0] >> b) & 1
                    and (beats[a] >> x0) & 1 and (beats[b] >> x0) & 1):
                continue
            Q = (1 << x0) | (1 << y0) | (1 << a) | (1 << b)
            ok = True
            Y = 0
            X = 0
            for r in range(n):
                if not (full >> r) & 1 or (Q >> r) & 1:
                    continue
                if (beats[r] & Q) == Q:
                    Y |= 1 << r
                elif (beats[r] & Q) == 0:
                    X |= 1 << r
                else:
                    ok = False
                    break
            if not ok:
                continue
            for r in range(n):
                if (Y >> r) & 1 and (beats[r] & X) != X:
                    ok = False
                    break
            if ok:
                return True
    return False


@njit(cache=True)
def set_tournament(n, bits, beats):
    for i in range(n):
        beats[i] = 0
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (bits >> k) & 1:
                beats[i] |= 1 << j
            else:
                beats[j] |= 1 << i
            k += 1


@njit(cache=True)
def check_pairs(n, beats, Hs, res):
    """res[0] += #pairs satisfying the hypotheses; res[1] += #pairs with the structure; res[2] += mismatches."""
    full = (1 << n) - 1
    for a in range(n):
        for b in range(a + 1, n):
            R = full ^ (1 << a) ^ (1 << b)
            hyp = (Hs[R | (1 << a)] > Hs[R]) and (Hs[R | (1 << b)] > Hs[R]) and \
                  (Hs[full] - Hs[R | (1 << a)] - Hs[R | (1 << b)] + Hs[R] == 0)
            st = structure_P(n, beats, a, b, R)
            if hyp:
                res[0] += 1
                if not (Hs[R | (1 << a)] == 3 * Hs[R] and Hs[R | (1 << b)] == 3 * Hs[R] and Hs[full] == 5 * Hs[R]):
                    res[2] += 1
            if st:
                res[1] += 1
            if hyp != st:
                res[2] += 1


@njit(cache=True)
def exhaustive_pairs(n):
    npairs = n * (n - 1) // 2
    size = 1 << n
    cnt = np.zeros((size, n), dtype=np.int64)
    Hs = np.zeros(size, dtype=np.int64)
    beats = np.zeros(n, dtype=np.int64)
    res = np.zeros(3, dtype=np.int64)
    for bits in range(1 << npairs):
        set_tournament(n, bits, beats)
        hs_table(n, beats, cnt, Hs)
        check_pairs(n, beats, Hs, res)
    return res


@njit(cache=True)
def check_partition2(n, beats, Hs, B, W1, W2, res):
    """Theorem R on one ordered partition; res[0] hypotheses hold, res[1] failures."""
    c = Hs[B]
    f1 = Hs[B | W1]
    f2 = Hs[B | W2]
    f12 = Hs[B | W1 | W2]
    if not (f1 > c and f2 > c and f12 + c == f1 + f2):
        return
    res[0] += 1
    bad = False
    if not (f1 == 3 * c and f2 == 3 * c and f12 == 5 * c):
        bad = True
    na = 0
    nb = 0
    a = -1
    b = -1
    full = B | W1 | W2
    for v in range(n):
        if (W1 >> v) & 1 and Hs[B | (1 << v)] > c:
            na += 1
            a = v
        if (W2 >> v) & 1 and Hs[B | (1 << v)] > c:
            nb += 1
            b = v
    if na != 1 or nb != 1:
        bad = True
    else:
        if not structure_P(n, beats, a, b, B):
            bad = True
        for v in range(n):
            if ((W1 | W2) >> v) & 1 and v != a and v != b:
                if Hs[full] != Hs[full ^ (1 << v)]:
                    bad = True
    if bad:
        res[1] += 1


@njit(cache=True)
def exhaustive_modules(n):
    npairs = n * (n - 1) // 2
    size = 1 << n
    cnt = np.zeros((size, n), dtype=np.int64)
    Hs = np.zeros(size, dtype=np.int64)
    beats = np.zeros(n, dtype=np.int64)
    res2 = np.zeros(2, dtype=np.int64)
    res3 = np.zeros(2, dtype=np.int64)      # [additive active 3-cubes found, module triples examined]
    n3 = 1
    for i in range(n):
        n3 *= 3
    n4 = 1 << (2 * n)
    for bits in range(1 << npairs):
        set_tournament(n, bits, beats)
        hs_table(n, beats, cnt, Hs)
        for code in range(n3):
            B = 0
            W1 = 0
            W2 = 0
            t = code
            for v in range(n):
                d = t % 3
                t //= 3
                if d == 0:
                    B |= 1 << v
                elif d == 1:
                    W1 |= 1 << v
                else:
                    W2 |= 1 << v
            if W1 == 0 or W2 == 0:
                continue
            check_partition2(n, beats, Hs, B, W1, W2, res2)
        for code in range(n4):
            B = 0
            W1 = 0
            W2 = 0
            W3 = 0
            t = code
            for v in range(n):
                d = t & 3
                t >>= 2
                if d == 0:
                    B |= 1 << v
                elif d == 1:
                    W1 |= 1 << v
                elif d == 2:
                    W2 |= 1 << v
                else:
                    W3 |= 1 << v
            if W1 == 0 or W2 == 0 or W3 == 0:
                continue
            if not ((W1 & -W1) < (W2 & -W2) and (W2 & -W2) < (W3 & -W3)):
                continue
            c = Hs[B]
            d1 = Hs[B | W1] - c
            d2 = Hs[B | W2] - c
            d3 = Hs[B | W3] - c
            res3[1] += 1
            if d1 > 0 and d2 > 0 and d3 > 0:
                if (Hs[B | W1 | W2] == c + d1 + d2 and Hs[B | W1 | W3] == c + d1 + d3 and
                        Hs[B | W2 | W3] == c + d2 + d3 and Hs[B | W1 | W2 | W3] == c + d1 + d2 + d3):
                    res3[0] += 1
    return res2, res3


@njit(cache=True)
def exhaustive_triples_single(n):
    npairs = n * (n - 1) // 2
    size = 1 << n
    full = size - 1
    cnt = np.zeros((size, n), dtype=np.int64)
    Hs = np.zeros(size, dtype=np.int64)
    beats = np.zeros(n, dtype=np.int64)
    found = 0
    for bits in range(1 << npairs):
        set_tournament(n, bits, beats)
        hs_table(n, beats, cnt, Hs)
        for a in range(n):
            for b in range(a + 1, n):
                for d in range(b + 1, n):
                    rest = full ^ (1 << a) ^ (1 << b) ^ (1 << d)
                    B = rest
                    while True:
                        c = Hs[B]
                        da = Hs[B | (1 << a)] - c
                        db = Hs[B | (1 << b)] - c
                        dd = Hs[B | (1 << d)] - c
                        if da > 0 and db > 0 and dd > 0:
                            if (Hs[B | (1 << a) | (1 << b)] == c + da + db
                                    and Hs[B | (1 << a) | (1 << d)] == c + da + dd
                                    and Hs[B | (1 << b) | (1 << d)] == c + db + dd
                                    and Hs[B | (1 << a) | (1 << b) | (1 << d)] == c + da + db + dd):
                                found += 1
                        if B == 0:
                            break
                        B = (B - 1) & rest
    return found


@njit(cache=True)
def sampled(n, ntrials, seed, nparts):
    np.random.seed(seed)
    size = 1 << n
    cnt = np.zeros((size, n), dtype=np.int64)
    Hs = np.zeros(size, dtype=np.int64)
    beats = np.zeros(n, dtype=np.int64)
    res = np.zeros(3, dtype=np.int64)
    res2 = np.zeros(2, dtype=np.int64)
    for t in range(ntrials):
        for i in range(n):
            beats[i] = 0
        # half of the samples are sparse-feedback tournaments (few reversed arcs), where the structure is likelier
        sparse = (t % 2 == 0)
        for i in range(n):
            for j in range(i + 1, n):
                r = np.random.random()
                fwd = (r < 0.5) if not sparse else (r < 0.85)
                if fwd:
                    beats[i] |= 1 << j
                else:
                    beats[j] |= 1 << i
        hs_table(n, beats, cnt, Hs)
        check_pairs(n, beats, Hs, res)
        for q in range(nparts):
            B = 0
            W1 = 0
            W2 = 0
            for v in range(n):
                d = np.random.randint(0, 4)
                if d <= 1:
                    B |= 1 << v
                elif d == 2:
                    W1 |= 1 << v
                else:
                    W2 |= 1 << v
            if W1 == 0 or W2 == 0:
                continue
            check_partition2(n, beats, Hs, B, W1, W2, res2)
    return res, res2


def main():
    ok = True
    t0 = time.time()
    # number of labelled (tournament, pair) with the Lemma P structure:
    #   C(n,4) * 24 * sum_k C(n-4,k) 2^C(k,2) 2^C(n-4-k,2)
    expected_pairs = {3: 0, 4: 24, 5: 240, 6: 2160, 7: 23520}
    for n in range(3, 8):
        res = exhaustive_pairs(n)
        good = (res[2] == 0 and res[0] == res[1] == expected_pairs[n])
        print(f"(P) n={n}: pairs with hypotheses {res[0]}, pairs with structure {res[1]}, mismatches {res[2]}",
              "OK" if good else "FAIL", flush=True)
        ok &= good
    for n in range(3, 7):
        res2, res3 = exhaustive_modules(n)
        good = (res2[1] == 0 and res3[0] == 0)
        print(f"(R) n={n}: active zero-interaction ordered module pairs {res2[0]}, violations of Theorem R {res2[1]}; "
              f"(K) module triples examined {res3[1]}, additive active 3-cubes {res3[0]}",
              "OK" if good else "FAIL", flush=True)
        ok &= good
    for n in range(3, 8):
        found = exhaustive_triples_single(n)
        good = (found == 0)
        print(f"(K) n={n}: additive active triples of single vertices over every base: {found}",
              "OK" if good else "FAIL", flush=True)
        ok &= good
    for n, trials, parts in ((8, 20000, 60), (9, 6000, 60), (10, 1500, 60)):
        res, res2 = sampled(n, trials, 20261110 + n, parts)
        good = (res[2] == 0 and res2[1] == 0)
        print(f"(S) n={n}: {trials} random tournaments: Lemma P pairs with hypotheses {res[0]}, with structure {res[1]}, "
              f"mismatches {res[2]}; Theorem R instances {res2[0]}, violations {res2[1]}",
              "OK" if good else "FAIL", flush=True)
        ok &= good
    print("ALL CHECKS PASSED" if ok else "SOME CHECK FAILED", f"[{time.time()-t0:.0f}s]")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
