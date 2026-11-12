# Verifier (lane architecture): spot check of two S93 two-apex statements that route R1 cites as proved / exact.
# Family G(2,0): transitive spine z_1 -> ... -> z_n, apexes u -> v, every apex-spine arc free (4 letters per spine
# vertex).  (a) first moment: sum of H over the 4^n words = 4^n (3^n + 3 2^n + 3)/8 (S93 TA-10);
# (b) first misses (S93 TA-06): smallest odd N > 21 with no word of length <= L is 61 (L = 5), 297 (L = 6).
import sys

def ham_count(n, out):
    full = 1 << n
    dp = [[0] * n for _ in range(full)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, full):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if c:
                nxt = out[v] & ~mask
                while nxt:
                    b = nxt & -nxt
                    dp[mask | b][b.bit_length() - 1] += c
                    nxt ^= b
    return sum(dp[full - 1])

vals = set()
ok = True
for n in range(1, 7):
    N = n + 2
    u, v = n, n + 1
    tot = 0
    for word in range(4 ** n):
        out = [0] * N
        for i in range(n):
            for j in range(i + 1, n):
                out[i] |= 1 << j
        out[u] |= 1 << v
        w = word
        for i in range(n):
            l = w & 3; w >>= 2
            if l & 1: out[u] |= 1 << i
            else: out[i] |= 1 << u
            if l & 2: out[v] |= 1 << i
            else: out[i] |= 1 << v
        h = ham_count(N, out)
        tot += h
        vals.add(h)
    pred = 4 ** n * (3 ** n + 3 * 2 ** n + 3) // 8
    first = next(x for x in range(23, 10 ** 6, 2) if x not in vals)
    print("n = %d: sum H = %d, formula = %d, %s ; smallest odd N > 21 not attained with length <= %d: %d ; 7 or 21 attained: %s"
          % (n, tot, pred, "OK" if tot == pred else "MISMATCH", n, first, (7 in vals) or (21 in vals)))
    if tot != pred: ok = False
print("CHECK first moment", "PASS" if ok else "FAIL")
