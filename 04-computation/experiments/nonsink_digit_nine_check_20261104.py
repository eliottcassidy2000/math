# nonsink_digit_nine_check_20261104.py  (opus-2026-11-04-S78, after the audit)
# Author-side check of the audit's upgrade U1 and its exact residual formula U4.
# U1: in T1 = C2(B, a), every old vertex w != a with d_B(w) >= 1 has sigma_w(T1) >= 9, with equality iff
#     B is transitive, a is its sink and w is its penultimate vertex.
# U4: for B = TT_m + w with w -> u_j only (1 <= j < m) and a = u_m,
#     sigma_w(C2(B, a)) = 4 + 2[j = m - 1] + 3 * 2^(m - j).
# Method: every LABELLED tournament B of order 2..6 and every mark a (no isomorphism reduction), direct subset DP.
import itertools


def tables(n, out):
    end = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        end[1 << v][v] = 1
    for S in range(1, 1 << n):
        row = end[S]
        for v in range(n):
            c = row[v]
            if c:
                nb = out[v] & ~S
                while nb:
                    w = (nb & -nb).bit_length() - 1
                    nb &= nb - 1
                    end[S | (1 << w)][w] += c
    return end


def digit(n, out, t):
    """sigma_t = R_t - h with R_t = sum_{X subset V - t} h(T[X]) s_t(T[V - X])."""
    full = (1 << n) - 1
    end = tables(n, out)
    hs = [sum(end[S]) if S else 1 for S in range(1 << n)]
    inn = [0] * n
    for p in range(n):
        for q in range(n):
            if out[p] >> q & 1:
                inn[q] |= 1 << p
    st = tables(n, inn)          # paths of the converse ending at t = paths starting at t
    R = 0
    for X in range(1 << n):
        if X >> t & 1:
            continue
        R += hs[X] * st[full & ~X][t]
    return R - hs[full]


def c2(n, outB, a):
    """C2(B, a): x = n beaten by B - a, x -> a; y = n + 1 beaten by all of B, y -> x."""
    out = list(outB) + [0, 0]
    x, y = n, n + 1
    out[x] |= 1 << a
    for b in range(n):
        if b != a:
            out[b] |= 1 << x
        out[b] |= 1 << y
    out[y] |= 1 << x
    return out


def transitive_order(n, out):
    sc = sorted(range(n), key=lambda v: -bin(out[v]).count('1'))
    if [bin(out[v]).count('1') for v in sc] != list(range(n - 1, -1, -1)):
        return None
    return sc                    # sc[0] the source, sc[-1] the sink


pairs = list(itertools.combinations(range(6), 2))
for n in range(2, 7):
    arcs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    mn = None
    eq_hits = 0
    eq_bad = 0
    cnt = 0
    for bits in range(1 << len(arcs)):
        out = [0] * n
        for k, (i, j) in enumerate(arcs):
            if bits >> k & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
        tor = transitive_order(n, out)
        for a in range(n):
            T1 = c2(n, out, a)
            for w in range(n):
                if w == a or out[w] == 0:
                    continue
                cnt += 1
                s = digit(n + 2, T1, w)
                mn = s if mn is None else min(mn, s)
                predicted_eq = tor is not None and tor[-1] == a and tor[-2] == w
                if s == 9:
                    eq_hits += 1
                if (s == 9) != predicted_eq or s < 9:
                    eq_bad += 1
    print('order', n, ': labelled (B, a, w) with w != a non-sink:', cnt, '| min digit', mn,
          '| digit-9 cases', eq_hits, '| mismatches with the equality characterisation', eq_bad)

bad4 = 0
cnt4 = 0
for m in range(2, 9):
    for j in range(1, m):
        n = m + 1                                # vertices u_1..u_m = 0..m-1, w = m
        out = [0] * n
        for p in range(m):
            for q in range(p + 1, m):
                out[p] |= 1 << q
        w = m
        for i in range(m):
            if i == j - 1:
                out[w] |= 1 << i
            else:
                out[i] |= 1 << w
        T1 = c2(n, out, m - 1)                   # a = u_m
        s = digit(n + 2, T1, w)
        cnt4 += 1
        if s != 4 + 2 * (j == m - 1) + 3 * 2 ** (m - j):
            bad4 += 1
print('U4 residual formula: cases (m <= 8)', cnt4, 'failures', bad4)
