# one_apex_triple_check_20261104.py  (opus-2026-11-04-S78)
# Brute-force checks:
#   (1) the one-apex apex triple: T(w) = chain 1 -> ... -> n plus apex x with x -> i iff bit i of w is 1 (most significant
#       first). Then H(T(w)) = g(v), e_x = e(v), R_x = 1 + v, where v = value of w,
#       g(0) = 1, g(2v + 1) = g(v), g(2v) = g(v) + v + (v mod 2), e(v) = [v = 0 or v even]. All words of length <= 7.
#   (2) the one-apex x one-apex one-arc join T(w1) => T(w2)^op with x2 -> x1 reversed has
#       H = (1 + u)(1 + v) + g(u) g(v) - e(u) e(v), v = val(w1), u = val(w2). All 961 pairs of words of length <= 4.
import itertools


def path_tables(n, out):
    full = (1 << n) - 1
    end = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        end[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            c = end[S][v]
            if not c:
                continue
            nb = out[v] & ~S
            while nb:
                w = (nb & -nb).bit_length() - 1
                nb &= nb - 1
                end[S | (1 << w)][w] += c
    return end


def apex_triple(n, out, t):
    full = (1 << n) - 1
    end = path_tables(n, out)
    hs = [sum(end[S]) if S else 1 for S in range(1 << n)]
    inn = [0] * n
    for a in range(n):
        for b in range(n):
            if out[a] >> b & 1:
                inn[b] |= 1 << a
    st = path_tables(n, inn)  # paths in the converse ending at v = paths starting at v
    R = sum(hs[S] * st[full & ~S][t] for S in range(1 << n) if not S >> t & 1)
    return hs[full], end[full][t], R


def hcount(n, out):
    end = path_tables(n, out)
    return sum(end[(1 << n) - 1])


def g(v):
    if v == 0:
        return 1
    if v & 1:
        return g(v >> 1)
    u = v >> 1
    return g(u) + u + (u & 1)


def e(v):
    return 1 if (v == 0 or v % 2 == 0) else 0


def one_apex(w):
    n = len(w)
    out = [0] * (n + 1)
    for i in range(n):
        for j in range(i + 1, n):
            out[i] |= 1 << j
        if w[i]:
            out[n] |= 1 << i
        else:
            out[i] |= 1 << n
    return out


bad1 = cnt1 = 0
for n in range(0, 8):
    for w in itertools.product([0, 1], repeat=n):
        v = int(''.join(map(str, w)), 2) if n else 0
        h, ex, R = apex_triple(n + 1, one_apex(w), n)
        cnt1 += 1
        if (h, ex, R) != (g(v), e(v), 1 + v):
            bad1 += 1
print('(1) one-apex triple: words checked', cnt1, 'failures', bad1)

bad2 = cnt2 = 0
for n in range(0, 5):
    for m in range(0, 5):
        for w in itertools.product([0, 1], repeat=n):
            for w2 in itertools.product([0, 1], repeat=m):
                A = one_apex(w)
                Bt = one_apex(w2)
                N = n + 1 + m + 1
                out = [0] * N
                for a in range(n + 1):
                    out[a] = A[a]
                for i in range(m + 1):          # converse of T(w2) on vertices n+1 .. n+1+m
                    for j in range(m + 1):
                        if Bt[i] >> j & 1:
                            out[n + 1 + j] |= 1 << (n + 1 + i)
                X1, X2 = n, n + 1 + m
                for a in range(n + 1):
                    for b in range(n + 1, N):
                        if a == X1 and b == X2:
                            out[X2] |= 1 << X1
                        else:
                            out[a] |= 1 << b
                v = int(''.join(map(str, w)), 2) if n else 0
                u = int(''.join(map(str, w2)), 2) if m else 0
                cnt2 += 1
                if hcount(N, out) != (1 + u) * (1 + v) + g(u) * g(v) - e(u) * e(v):
                    bad2 += 1
print('(2) one-apex x one-apex joins: pairs checked', cnt2, 'failures', bad2)
