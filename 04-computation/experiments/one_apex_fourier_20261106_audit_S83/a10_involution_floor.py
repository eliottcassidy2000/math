# Audit S83, step 10: the forced energy floor.  The converse of T(w) is T(w*), w* = reverse complement, so H(w) = H(w*)
# (Session 16 / Session 191).  On X_n the map iota(x' 1^t) = (x')* 1^t (x' strong), iota(1^n) = 1^n, is an H-preserving
# involution; hence E_n >= 2 |X_n| - #Fix(iota), and #Fix <= 2^(n/2) + 1 (x' self-reverse-complementary has even length).
import numpy as np
from oa_engine import g_auto, g_block, energy

def iota(v, n):
    # v in [2^(n-1), 2^n): word x of length n; split x = x' 1^t with x' ending in 0 (or x = 1^n)
    if v == (1 << n) - 1:
        return v
    t = 0
    while (v >> t) & 1:
        t += 1
    m = n - t; xp = v >> t                          # x' has length m, ends in 0
    rc = 0
    for j in range(m):                              # reverse complement of x'
        rc = (rc << 1) | (1 - ((xp >> j) & 1))
    return (rc << t) | ((1 << t) - 1)

for n in (8, 12, 16, 18):
    lo = 1 << (n - 1)
    ok = True; fix = 0
    for v in range(lo, 1 << n):
        w = iota(v, n)
        ok &= (lo <= w < (1 << n)) and iota(w, n) == v and g_auto(w) == g_auto(v)
        fix += (w == v)
    vals = g_block(n)
    E, D, mx = energy(vals, n)
    X = 1 << (n - 1)
    print(f'n={n}: iota maps X_n to X_n, is an involution and preserves H: {ok}; #Fix = {fix} (bound 2^(n/2)+1 = {2 ** (n // 2) + 1}); '
          f'E_n/|X_n| = {E / X:.4f} >= 2 - #Fix/|X_n| = {2 - fix / X:.4f}')
