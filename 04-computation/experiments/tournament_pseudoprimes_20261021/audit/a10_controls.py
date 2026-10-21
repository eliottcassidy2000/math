# AUDIT a10 (independent auditor, 2026-10-21): THM-4746 statement 7 and its controls.
# (A) semicomplete digraphs: brute force n <= 5 (3^C(n,2) digraphs), count strong ones and check that the counts
#     satisfy 3^C(m,2) = sum_k C(m,k) s3(k) 3^C(m-k,2) (linearly ordered strong components);
# (B) s3(m+p-1) = (3/p) s3(m) mod p for 5 <= p <= 61, 1 <= m <= 4p; and p = 3 separately;
# (C) connected labelled graphs mod p: the session's "no (anti)period L < 2p" and the auditor's finer claim:
#     cg(m+p-1) = eps cg(m) (mod p) whenever m != 0, 1 (mod p) and m >= 2, eps = (2/p)
#     (it can only fail when m or m+p-1 is a multiple of p).
import itertools
from math import comb

def strong_bits(n, out):
    full = (1 << n) - 1
    def reach(a):
        seen = 1; fr = 1
        while fr:
            nx = 0
            for v in range(n):
                if fr >> v & 1: nx |= a[v]
            fr = nx & ~seen; seen |= nx
        return seen == full
    inn = [0] * n
    for v in range(n):
        for w in range(n):
            if out[v] >> w & 1: inn[w] |= 1 << v
    return reach(out) and reach(inn)

bf = {}
for n in range(1, 6):
    pairs = list(itertools.combinations(range(n), 2)); cnt = 0
    for code in itertools.product(range(3), repeat=len(pairs)):
        out = [0] * n
        for (i, j), cde in zip(pairs, code):
            if cde in (0, 2): out[i] |= 1 << j
            if cde in (1, 2): out[j] |= 1 << i
        cnt += strong_bits(n, out)
    bf[n] = cnt
def strong_rec(base, M, mod=None):
    T = [base ** (j * (j - 1) // 2) for j in range(M + 1)]
    if mod: T = [t % mod for t in T]
    s = [0] * (M + 1)
    for m in range(1, M + 1):
        v = T[m] - sum(comb(m, k) * s[k] * T[m - k] for k in range(1, m))
        s[m] = v % mod if mod else v
    return s
s3 = strong_rec(3, 5)
print("(A) strong semicomplete digraphs, brute force n=1..5:", [bf[n] for n in range(1, 6)], " recurrence:", s3[1:6],
      " equal:", [bf[n] for n in range(1, 6)] == s3[1:6])
def leg(a, p):
    return 0 if a % p == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)
primes = [q for q in range(3, 62) if all(q % d for d in range(2, int(q ** .5) + 1))]
res = []
for p in primes:
    s = strong_rec(3, 5 * p + 2, p); e = leg(3, p)
    ok = all((s[m + p - 1] - e * s[m]) % p == 0 for m in range(1, 4 * p + 1))
    res.append((p, e, ok))
print("(B) s3(m+p-1) = (3/p) s3(m) mod p, 1<=m<=4p:", res)
def conn_mod(M, p):
    g = [pow(2, j * (j - 1) // 2, p) for j in range(M + 1)]; c = [0] * (M + 1)
    for m in range(1, M + 1):
        c[m] = (g[m] - sum(comb(m - 1, k - 1) * c[k] * g[m - k] for k in range(1, m))) % p
    return c
cg_exact = conn_mod(12, 10 ** 30)
print("(C) connected labelled graphs 1..8 (A001187):", cg_exact[1:9])
for p in [q for q in primes if q >= 5][:9]:
    M = 6 * p; cg = conn_mod(M + p, p); e = leg(2, p)
    per = [L for L in range(1, 2 * p) if all((cg[m + L] - cg[m]) % p == 0 for m in range(1, M - L + 1)) or
           all((cg[m + L] + cg[m]) % p == 0 for m in range(1, M - L + 1))]
    good = [m for m in range(2, M + 1) if m % p not in (0, 1)]
    viol = [m for m in good if (cg[m + p - 1] - e * cg[m]) % p]
    fails_at = [m for m in range(1, M + 1) if (cg[m + p - 1] - e * cg[m]) % p]
    print(f"p={p}: session test (anti)periods L<2p: {per};  auditor: quasi-period p-1 with eps={e:+d} on m != 0,1 mod p "
          f"(m<= {M}): violations {len(viol)} of {len(good)};  all m where it fails: {fails_at[:12]}")
