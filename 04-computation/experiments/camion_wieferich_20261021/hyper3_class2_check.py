# Is the forced class-2 root of connected 3-uniform hypergraphs the point n = 2 (b_3(2) = 0)?  Fit the mod-p^2 cubic of the class
# n = 2 (mod p-1) from l = 0..3 at base m = p + 1 and evaluate at l = -1 (n = 2) and at the other small points n = 2 - (p-1)k.
# Also test the classes of n = 1 and n = 3, 4 (b_3(1) = 1, b_3(3) = 1, b_3(4) = 4+... ) for extrapolation down to small n.
# mac-mini-2026-10-21-fanopsp.
from math import comb
def b3_mod(M, q):
    T = [pow(2, comb(n, 3), q) for n in range(M + 1)]; b = [0] * (M + 1); row = [1]
    for n in range(1, M + 1):
        if n >= 2: row = [1] + [(row[i - 1] + row[i]) % q for i in range(1, n - 1)] + [1]
        b[n] = (T[n] - sum(row[k - 1] * b[k] * T[n - k] for k in range(1, n))) % q
    return b
for p in (23, 31, 47, 71):
    q = p * p; M = p + 1 + 4 * (p - 1) + 5
    b = b3_mod(M, q)
    out = []
    for r in (1, 2, 3, 4, 5):
        m = r + (p - 1) * (2 if r < 6 else 1)      # base well above the threshold, in the class of r
        A = [b[m + l * (p - 1)] for l in range(4)]
        d1 = (A[1] - A[0]) % q; d2 = (A[2] - 2 * A[1] + A[0]) % q; d3 = (A[3] - 3 * A[2] + 3 * A[1] - A[0]) % q
        def F(l):   # Newton form, l may be negative
            return (A[0] + d1 * l + d2 * (l * (l - 1) // 2) + d3 * (l * (l - 1) * (l - 2) // 6)) % q
        lsmall = -2                                  # n = m - 2(p-1) = r
        out.append((r, b[r] % q, F(lsmall), F(lsmall) == b[r] % q))
    print(f"p={p}: extrapolating each class cubic down to n = r: " + "; ".join(f"n={r}: b={br}, cubic={fv} {'=' if ok else '!='}" for r, br, fv, ok in out))
