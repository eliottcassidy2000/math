"""Gilbreath depth: for the first N primes, the least row k such that row k = (1, then only 0s and 2s) on the first N-k
entries (the absorbing sea; afterwards the leading 1 persists, THM-4511 / S9 note).  Also the frontier F(r) = position of
the first entry >= 4 in row r (last defect), and max entry per row.  Primes below LIMIT by a numpy sieve (odd-only)."""
import numpy as np, sys, time, math
LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 200_000_000
t0 = time.time()
n = LIMIT // 2 + 1                       # odd numbers 2i+1, i = 0..n-1
s = np.ones(n, dtype=bool); s[0] = False # 1 is not prime
for i in range(1, int(math.isqrt(LIMIT))//2 + 1):
    if s[i]:
        p = 2*i + 1
        s[(p*p)//2::p] = False
primes = np.concatenate(([2], 2*np.nonzero(s)[0] + 1)).astype(np.int64)
NP = len(primes); print(f"{NP} primes below {LIMIT} (sieve {time.time()-t0:.1f}s)")
# depth for several N (prefix lengths), computed in one pass on the longest row: row r has length NP - r; the prefix of
# length N - r of row r depends only on the first N primes, so one pass gives all depths.
Ns = [10**k for k in range(3, 8) if 10**k <= NP] + [NP]
row = primes.copy(); depth = {N: None for N in Ns}; frontier = []
r = 0
maxrows = 3000
while r < maxrows and any(v is None for v in depth.values()):
    r += 1
    row = np.abs(np.diff(row))           # row r
    big = np.nonzero(row[1:] > 2)[0]     # entries after the leading one that are not 0/2
    first_big = (big[0] + 2) if len(big) else None    # 1-indexed position in the row
    last_big = (big[-1] + 2) if len(big) else None
    if r <= 20 or r % 10 == 0:
        frontier.append((r, int(row[0]), first_big, last_big, int(row[:min(len(row), 10**7)].max())))
    for N in Ns:
        if depth[N] is None:
            L = N - r                    # prefix length of row r determined by the first N primes
            if L < 2: depth[N] = r; continue
            if row[0] == 1 and np.all(row[1:L] <= 2):
                depth[N] = r
print("rows computed:", r, f"({time.time()-t0:.1f}s)")
print("Gilbreath depth k(N) (first row that is 1 followed by 0/2 on the first N-k entries):")
for N in Ns: print(f"  N = {N:9d} primes (p_N = {primes[N-1]}): depth {depth[N]}; log p_N = {math.log(primes[N-1]):.2f}; (log p_N)^2 = {math.log(primes[N-1])**2:.0f}; max gap below p_N = {int(np.max(np.diff(primes[:N])))}")
print("frontier profile (row, leading entry, first position with entry >= 4, last such position, max entry in first 1e7):")
for f in frontier[:60]: print("  ", f)
