#!/usr/bin/env python3
"""Owner's odd-square bracket observations (2026-09-23), tested exactly (vectorized segmented sieve).
Bracket B_n = ((2n-1)^2, (2n+1)^2], n >= 1 (B_1 = (1,9]).
(1) Theorem: the primes p with 2p in the same bracket as p are exactly {2,3,11} (proof in the note; exhaustive check).
(2) Prime counts per bracket ('lengths' in the owner's encoding = count + 1 with the square appended).
(3) d_hi(n) = (2n+1)^2 - max prime in B_n;  d_lo(n) = min prime in B_n - (2n-1)^2.
(4) 'Something every 5th occurrence': autocorrelation at lags 1..10 and means by n mod 5."""
import math, sys
import numpy as np
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 30000
LIM = (2*NMAX+1)**2
r = math.isqrt(LIM)+1
base = np.ones(r+1, dtype=bool); base[:2] = False
for i in range(2, math.isqrt(r)+1):
    if base[i]: base[i*i::i] = False
small = np.flatnonzero(base)
cnt = np.zeros(NMAX+2, dtype=np.int64); mn = np.full(NMAX+2, -1, dtype=np.int64); mx = np.full(NMAX+2, -1, dtype=np.int64)
seg = 1 << 24; start = 2
while start <= LIM:
    end = min(start+seg, LIM+1)
    a = np.ones(end-start, dtype=bool)
    for p in small:
        p = int(p)
        if p*p >= end: break
        f = max(p*p, ((start+p-1)//p)*p)
        a[f-start::p] = False
    q = np.flatnonzero(a).astype(np.int64) + start
    s = np.floor(np.sqrt(q.astype(np.float64))).astype(np.int64)
    s -= (s*s > q); s += ((s+1)*(s+1) <= q)          # exact isqrt
    n = (s + 1)//2                                      # (2n-1)^2 < q <= (2n+1)^2  <=>  n = ceil((sqrt q - 1)/2); q non-square
    n = np.where((2*n+1)**2 < q, n+1, n)
    keep = n <= NMAX; q = q[keep]; n = n[keep]
    if q.size:
        cnt += np.bincount(n, minlength=NMAX+2)
        u, first = np.unique(n, return_index=True); last = np.r_[first[1:], len(n)] - 1
        newmin = mn[u] < 0; mn[u[newmin]] = q[first[newmin]]; mx[u] = q[last]
    start = end
N = np.arange(1, NMAX+1)
counts = cnt[1:NMAX+1]; dhi = (2*N+1)**2 - mx[1:NMAX+1]; dlo = mn[1:NMAX+1] - (2*N-1)**2
print(f"brackets n = 1..{NMAX} (odd squares up to {LIM}); every bracket nonempty: {bool((counts>0).all())}; min count {counts.min()} (n>=2: {counts[1:].min()})")
esc = []
for p in np.flatnonzero(base[:10**6]):
    p = int(p); nn = max(1, (math.isqrt(p-1)+1 + 1)//2)
    while (2*nn+1)**2 < p: nn += 1
    while nn > 1 and (2*nn-1)**2 >= p: nn -= 1
    if 2*p <= (2*nn+1)**2: esc.append(p)
print("(1) primes p < 10^6 with 2p inside p's own bracket:", esc)
print("(2) counts n=1..24:", counts[:24].tolist())
fails = [int(i) for i in N[:40] if counts[i-1] != i+3]
print("    counts == n+3 exactly for n =", [int(i) for i in N[:40] if counts[i-1] == i+3], "; first failure n =", fails[0])
print("(3) d_hi n=1..24:", dhi[:24].tolist())
print("    d_lo n=1..24:", dlo[:24].tolist(), " (with 2,3,11 removed the first two become", [5-1, 13-9], ")")
v, f = np.unique(dhi, return_counts=True); top = sorted(zip(v.tolist(), (f/len(dhi)).round(4).tolist()), key=lambda t: -t[1])[:6]
print("    d_hi value shares (top 6):", top)
def acf(x, L):
    x = x - x.mean(); return float((x[:-L]*x[L:]).mean()/(x*x).mean())
expc = 4*N/np.log(2*N+1)   # ~ (8n)/(2 log(2n+1)) primes in B_n
for name, x in [("d_hi", dhi.astype(float)), ("d_lo", dlo.astype(float)), ("count/(4n/log(2n+1))", counts/expc)]:
    ac = [round(acf(x, L), 4) for L in range(1, 11)]
    means = [round(float(x[N % 5 == r].mean()), 4) for r in range(5)]
    se = float(x.std()/math.sqrt(len(x)/5))
    print(f"(4) {name}: acf lags 1..10 = {ac}")
    print(f"    mean by n mod 5 (r=0..4) = {means}  (std error per class ~ {se:.4f})")

# (5) the mechanism: local densities (singular series). Group d_lo by (2n-1)^2 mod 5 and d_hi by (2n+1)^2 mod 5,
#     and by the first even offset j at which 5 divides (2n-1)^2 + j (resp. (2n+1)^2 - j).
for name, x, sq in [("d_lo", dlo, (2*N-1)**2), ("d_hi", dhi, (2*N+1)**2)]:
    for p in (3, 5, 7):
        cls = sq % p
        out = []
        for c in sorted(set(cls.tolist())):
            sel = cls == c
            out.append((c, round(float(x[sel].mean()), 3), int(sel.sum())))
        print(f"(5) {name} by square mod {p}: (residue, mean, count) = {out}")
