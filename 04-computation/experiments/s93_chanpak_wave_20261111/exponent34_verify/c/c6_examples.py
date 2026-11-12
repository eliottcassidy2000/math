#!/usr/bin/env python3
"""Verifier C: explicit examples (exact big-integer arithmetic, written out from the definitions)."""
import itertools


def q(h):
    t = 0
    while h:
        if h & 1:
            t += ((h >> 1) + 1) // 2
        h >>= 1
    return t


def xv(h):
    return 2 * h + 1 - 2 * q(h)


def g(v):
    if v == 0:
        return 1
    u = v >> 1
    return g(u) if v & 1 else g(u) + u + (u & 1)


def show(a, b):
    n = len(a)
    assert len(b) == n and a[0] == b[0] == "1" and a != b
    h, hp = int(a, 2), int(b, 2)
    tp = 0
    while a[tp] == b[tp]:
        tp += 1
    s = 0
    while a[n - 1 - s] == b[n - 1 - s]:
        s += 1
    w = n - tp - s
    P = int(a[:tp], 2)
    xs, xps, Zs = a[tp:tp + w], b[tp:tp + w], a[tp + w:]
    x, xp = int(xs, 2), int(xps, 2)
    c = 2 ** (s + 1) - (int(Zs[::-1], 2) if s else 0)
    X = x - xp
    Xi = int(xs[::-1], 2) - int(xps[::-1], 2)
    gam = 2 * (q(x) - q(xp)) + (P % 2) * (int(xs[0]) - int(xps[0])) + (c % 2) * (int(xs[-1]) - int(xps[-1]))
    m = max(max(len(list(gr)) for _, gr in itertools.groupby(a)), max(len(list(gr)) for _, gr in itertools.groupby(b)))
    mn = min(abs(X), abs(Xi))
    rho = w - (mn.bit_length() - 1)
    print("h  = %s  x(h)  = %d = g(2h)  = %d" % (a, xv(h), g(2 * h)))
    print("h' = %s  x(h') = %d = g(2h') = %d   collision: %s" % (b, xv(hp), g(2 * hp), xv(h) == xv(hp)))
    print("   n = %d, t' = %d, s = %d, w = %d, P = %d, c = %d, x = %s, x' = %s, X = %d, Xi = %d" % (n, tp, s, w, P, c, xs, xps, X, Xi))
    print("   gamma = %d, c X - P Xi = %d, 2^w = %d, rho = %d (lemma: %d <= s - t' = %d <= %d), longest run m = %d, |s - t'| - m = %d"
          % (gam, c * X - P * Xi, 2 ** w, rho, -(rho + 2), s - tp, rho - 1, m, abs(s - tp) - m))
    print()


print("pair with |s - t'| = R + 2 for R = m = 4 (n = 22):")
show("1000101111011000100001", "1000100001010101011110")
print("pair attaining s = t' + rho - 1 (the second carries bound is sharp), n = 22:")
show("1000100100000000011111", "1000110100111111111111")
print("pair with t' = s + rho + 1 (the first carries bound t' <= s + rho + 2 is never attained for n <= 22), n = 22:")
show("1000000000000110100111", "1000001111111110110111")
print("smallest collision:")
show("10", "11")
print("a pair with |gamma| = 2^w (w = 1):")
show("100000001110111111", "100000001010111111")
