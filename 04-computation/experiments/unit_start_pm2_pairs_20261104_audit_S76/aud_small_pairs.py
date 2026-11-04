"""aud_small_pairs.py - do Session189-wrapper +-2 pairs (equal F, h' = h + 2, both strong) exist below order 7 when the
unit-start condition s = 1 is dropped?

Exhaustive over all labelled tournaments on {0..n-1} with root 0, n = 4, 5, 6.  Since h differs inside a pair, two
labelled tournaments with profiles (h, F) and (h + 2, F) are automatically non-isomorphic, so pair EXISTENCE only needs
the set of realised (h, s, F) profiles.  Prints every realised pair of profiles (h, s, F) / (h + 2, s', F).
"""
import itertools


def profile(n, out):
    full = (1 << n) - 1
    end = [[0] * n for _ in range(1 << n)]
    beg0 = [[0] * n for _ in range(1 << n)]
    for x in range(n):
        end[1 << x][x] = 1
    beg0[1][0] = 1
    for S in range(1, full + 1):
        for x in range(n):
            a, b = end[S][x], beg0[S][x]
            if not a and not b:
                continue
            nb = out[x] & ~S
            y = 0
            while nb:
                if nb & 1:
                    end[S | (1 << y)][y] += a
                    beg0[S | (1 << y)][y] += b
                nb >>= 1
                y += 1
    HS = [1] + [sum(end[S]) for S in range(1, full + 1)]
    h = HS[full]
    s = sum(beg0[full])
    F = sum(end[S][0] * HS[full ^ S] for S in range(1, full + 1, 2))
    return h, s, F


def strong(n, out):
    full = (1 << n) - 1
    f = 1
    while True:
        g = f
        for x in range(n):
            if (f >> x) & 1:
                g |= out[x]
        if g == f:
            break
        f = g
    if f != full:
        return False
    b = 1
    while True:
        g = b
        for y in range(n):
            if not (b >> y) & 1 and out[y] & b:
                g |= 1 << y
        if g == b:
            break
        b = g
    return b == full


for n in (4, 5, 6):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    profs = set()
    for code in range(1 << len(pairs)):
        out = [0] * n
        for k, (i, j) in enumerate(pairs):
            if (code >> k) & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
        if not strong(n, out):
            continue
        profs.add(profile(n, out))
    byF = {}
    for (h, s, F) in profs:
        byF.setdefault(F, set()).add((h, s))
    found = []
    for F, hs in byF.items():
        for (h, s) in hs:
            for (h2, s2) in hs:
                if h2 == h + 2:
                    found.append(((h, s, F), (h2, s2, F)))
    found.sort()
    print("n=%d: %d realised strong rooted profiles (h,s,F); +-2 profile pairs with equal F: %d" % (n, len(profs), len(found)))
    for a, b in found:
        print("   (h,s,F) = %s / %s%s" % (a, b, "   [s = s' = 1]" if a[1] == 1 and b[1] == 1 else ""))
