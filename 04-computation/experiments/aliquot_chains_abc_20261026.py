"""aliquot_chains_abc_20261026.py -- opus-2026-10-26-S60

Exact small-number structures behind the owner's merge of the 1-2-4 decomposition of 7, the untouchable primes
{2, 5}, the forbidden Hamiltonian-path counts {7, 21}, the powers/doublings chain structure and the aliquot
preimage function s(n) = sigma(n) - n.

1. Aliquot preimage table s^{-1}(v), v = 0..60 (exact: a composite n has s(n) > sqrt(n), so n < v^2), the untouchable
   values, the odd values with a single source below 200, the values with >= 2 sources and their shapes.
2. Untouchable values <= 1000 certified exactly (n < 10^6), which are prime, which are odd.
3. Closed forms on chains: s(2^k) = 2^k - 1, s(p^k) = (p^k - 1)/(p - 1), s(2^j p) = (2^j - 1) p + 2^(j+1) - 1,
   s(p^2) = p + 1, s(2p) = p + 3, s(pq) = p + q + 1; the twin-prime collision s(p^2) = s(2(p-2)) = p + 1 for
   primes p - 2, p; the fixed-point condition s(2q) = 2q iff q = 3 (so the collision is a fixed point only at p = 5).
4. Exact preimages of 7, 21, 28, 496 and of the small Mersenne values.
5. The chain structure (each n: column n, n^2, n^3, ... and row n, 2n, 4n, ...): which n and which chain values
   already occur among the rays of 1..n-1; NEW numbers = odd non-perfect-powers; the square law n = 2^a r (r odd):
   n^2 already occurs iff r^2 < n or n is a perfect power, and n^2 is a LATER doubling (of r^2 > n) iff a >= 1 and
   r > 2^a; first instance 6; the same for cubes; the first odd number whose row meets an earlier column (9: 36 = 6^2).
6. abc hits (a + b = c coprime, rad(abc) < c) with c <= 10000, each term typed new/old; the hits with c a chain value.
7. Goldbach routes for odd values: v = s(pq) iff v - 1 is a sum of two distinct primes; failures below 2000.
8. Hamiltonian path counts of the Paley tournaments P_3, P_7, P_11 and of the small transitive stacks, for the
   {7, 21} comparison (THM-1370: 7 and 21 are never path counts).
Usage: python3 aliquot_chains_abc_20261026.py   (numpy)
"""
import numpy as np
from math import isqrt, gcd, log
from collections import defaultdict
from itertools import permutations

LIMIT = 1000000

def sigma_table(N):
    s = np.zeros(N + 1, dtype=np.int64)
    for d in range(1, N + 1): s[d::d] += d
    return s

def primes_upto(N):
    sieve = np.ones(N + 1, dtype=bool); sieve[:2] = False
    for i in range(2, isqrt(N) + 1):
        if sieve[i]: sieve[i*i::i] = False
    return np.flatnonzero(sieve), sieve

def rad_table(N, primes):
    r = np.ones(N + 1, dtype=np.int64)
    for p in primes:
        if p > N: break
        r[p::p] *= p
    return r

def is_perfect_power(n):
    if n < 4: return False
    for k in range(2, n.bit_length() + 1):
        c = round(n ** (1.0 / k))
        for cc in (c - 1, c, c + 1):
            if cc >= 2 and cc ** k == n: return True
    return False

def ham_paths(n, adj):
    """number of Hamiltonian paths of the tournament adj (adj[i][j] True iff i -> j), DP over subsets."""
    full = 1 << n
    dp = [[0] * n for _ in range(full)]
    for v in range(n): dp[1 << v][v] = 1
    for S in range(full):
        row = dp[S]
        for v in range(n):
            c = row[v]
            if not c: continue
            for w in range(n):
                if not (S >> w) & 1 and adj[v][w]: dp[S | (1 << w)][w] += c
    return sum(dp[full - 1])

def paley(q):
    qr = {(x * x) % q for x in range(1, q)}
    return [[((j - i) % q) in qr for j in range(q)] for i in range(q)]

def main():
    out = []
    sig = sigma_table(LIMIT)
    s = sig - np.arange(LIMIT + 1); s[0] = 0
    primes, isprime = primes_upto(LIMIT)
    # 1. preimage table
    pre = defaultdict(list)
    for n in range(1, 3601):
        v = int(s[n])
        if v <= 60: pre[v].append(n)
    out.append("1. aliquot preimages s^{-1}(v), v = 0..60 (exact; n < v^2 for composite n):")
    for v in range(0, 61):
        L = pre[v]
        txt = "primes" if v == 1 else (str(L) if L else "NONE (untouchable)")
        out.append(f"   {v:2d} <- {txt}")
    single_odd = [v for v in range(3, 200, 2) if int(np.count_nonzero(s[1:40000] == v)) == 1]
    out.append(f"   odd values below 200 with exactly one source (n < 200^2): {single_odd}")
    multi = {v: pre[v] for v in range(2, 61) if len(pre[v]) >= 2}
    out.append(f"   values <= 60 with >= 2 sources: {multi}")
    cnt = [len(pre[v]) for v in range(2, 61)]
    out.append(f"   source counts for v = 2..60: {cnt}")
    # 2. untouchables to 1000 (exact since LIMIT = 10^6 >= 1000^2)
    V = 1000
    hit = np.zeros(V + 1, dtype=bool)
    sv = s[1:LIMIT + 1]; sv = sv[sv <= V]; hit[sv] = True
    unt = [v for v in range(2, V + 1) if not hit[v]]
    out.append(f"2. untouchable values <= {V} (exact, all n < 10^6 scanned): count {len(unt)}; list {unt}; primes among them {[v for v in unt if isprime[v]]}; odd ones {[v for v in unt if v % 2]}")
    # 3. closed forms and the twin collision
    ok = True
    for k in range(1, 20): ok &= (int(s[2**k]) == 2**k - 1)
    for p in primes[:40]:
        p = int(p)
        for k in range(2, 7):
            if p**k <= LIMIT: ok &= (int(s[p**k]) == (p**k - 1) // (p - 1))
        if p > 2:
            for j in range(1, 10):
                if 2**j * p <= LIMIT: ok &= (int(s[2**j * p]) == (2**j - 1) * p + 2**(j+1) - 1)
    for p in primes[:80]:
        for q in primes[:80]:
            p = int(p); q = int(q)
            if p < q and p * q <= LIMIT: ok &= (int(s[p*q]) == p + q + 1)
    twins = [(int(p) - 2, int(p)) for p in primes if p > 3 and isprime[p - 2] and p * p <= LIMIT]
    twinok = all(int(s[p*p]) == int(s[2*(p-2)]) == p + 1 for (q, p) in twins)
    out.append(f"3. closed forms (Mersenne, repunit, affine rows, pq): {'all verified' if ok else 'FAILED'}; twin collision s(p^2) = s(2(p-2)) = p+1: {'verified' if twinok else 'FAILED'} on {len(twins)} twin pairs with p^2 <= 10^6; values p+1 for the first pairs: {[p + 1 for (q, p) in twins[:10]]}")
    fixed = [n for n in range(1, 10000) if int(s[n]) == n]
    out.append(f"   fixed points of s below 10^4 (perfect numbers): {fixed}; s(2q) = 2q for an odd prime q iff q + 3 = 2q iff q = 3; so the twin collision value p+1 equals the doubling source 2(p-2) iff p = 5 (6 = s(6) = s(25))")
    # 4. exact preimages of 7, 21, 28, 496 and small Mersenne values
    def preim(v):
        B = min(LIMIT, v * v + 1)
        return [int(n) for n in np.flatnonzero(s[:B] == v) if n >= 1]
    for v in (3, 7, 15, 21, 28, 31, 63, 127, 496):
        out.append(f"4. s^-1({v}) = {preim(v)}" + ("  (exact, n < v^2)" if v * v < LIMIT else "  (n < 10^6)"))
    out.append("   7 = 1 + 2 + 4 = s(8) (the proper divisors of 2^3); 21 = s(18) = s(51) = s(91); 21 = 1 + 4 + 16 = (4^3 - 1)/3 is the square-divisor part of s(64) = 63 (the non-square part is 2 + 8 + 32 = 42)")
    # 5. chain structure bookkeeping
    N = 40
    rays = set()  # all ray values of 1..n-1 up to a bound
    BOUND = 10**6
    rows = []
    def add_rays(m):
        if m == 1: rays.add(1)
        else:
            x = m
            while x <= BOUND: rays.add(x); x *= m
        x = m
        while x <= BOUND: rays.add(x); x *= 2
    for n in range(1, N + 1):
        a = 0; r = n
        while r % 2 == 0: r //= 2; a += 1
        self_old = n in rays
        sq_old = (n * n) in rays
        cube_old = (n ** 3) in rays
        later_sq = (a >= 1 and r > 2 ** a)
        later_cube = (a >= 1 and r ** 2 > 2 ** a)
        # odd n: does its doubling row hit an earlier column? (2^j n a perfect power with root < n)
        row_hits = []
        if n % 2 == 1 and n > 1:
            x = 2 * n
            while x <= BOUND:
                if is_perfect_power(x) and x in rays: row_hits.append(x)
                x *= 2
        rows.append((n, a, r, self_old, sq_old, later_sq, cube_old, later_cube, row_hits))
        add_rays(n)
    out.append(f"5. chain structure, n = 1..{N} (n = 2^a r): n already present | n^2 already present | n^2 a later doubling (a>=1, r>2^a) | n^3 present | n^3 later | odd n whose row hits an earlier column:")
    for (n, a, r, so, sq, ls, co, lc, rh) in rows:
        out.append(f"   {n:2d} = 2^{a}*{r:2d}: n {'old' if so else 'NEW'} | n^2 {'old' if sq else 'new'}{' LATER' if ls else ''} | n^3 {'old' if co else 'new'}{' LATER' if lc else ''}" + (f" | row hits {rh}" if rh else ""))
    new_list = [n for (n, a, r, so, *_rest) in rows if not so]
    out.append(f"   NEW numbers <= {N}: {new_list}; first n with n^2 a later doubling: {next(n for (n,a,r,so,sq,ls,*_r) in rows if ls)}; first n with n^3 a later doubling: {next(n for (n,a,r,so,sq,ls,co,lc,rh) in rows if lc)}; claim 'NEW iff odd non-perfect-power' holds for n <= {N}: {all((not so) == (n % 2 == 1 and not is_perfect_power(n)) for (n,a,r,so,*_rest) in rows)}")
    # 6. abc hits
    C = 10000
    rad = rad_table(C, primes)
    hits = []
    def typ(x): return "new" if (x % 2 == 1 and not is_perfect_power(x)) else "old"
    for c in range(3, C + 1):
        a = np.arange(1, c // 2 + 1); b = c - a
        m = (np.gcd(a, b) == 1) & (rad[a] * rad[b] * rad[c] < c)
        for aa, bb in zip(a[m], b[m]):
            aa = int(aa); bb = int(bb); R = int(rad[aa] * rad[bb] * rad[c])
            hits.append((c, aa, bb, round(log(c) / log(R), 3), typ(aa), typ(bb), typ(c)))
    hits.sort(key=lambda t: -t[3])
    out.append(f"6. abc hits (rad(abc) < c), c <= {C}: {len(hits)}; top 12 by quality (c, a, b, q, type a/b/c): {hits[:12]}")
    out.append(f"   hits with c old: {sum(1 for h in hits if h[6]=='old')}; with b old: {sum(1 for h in hits if h[5]=='old')}; with a old: {sum(1 for h in hits if h[4]=='old')}; with a = 1: {sum(1 for h in hits if h[1]==1)}; with all three old: {sum(1 for h in hits if h[4]==h[5]==h[6]=='old')}; with no old term: {sum(1 for h in hits if 'old' not in (h[4],h[5],h[6]))}")
    small = [h for h in hits if h[0] <= 130]
    out.append(f"   all hits with c <= 130: {small}")
    # 7. Goldbach routes
    P = [int(p) for p in primes if p < 2000]
    def two_distinct(e):
        return any(isprime[e - p] and p != e - p for p in P if p < e)
    fails = [m for m in range(3, 2000, 2) if not two_distinct(m - 1)]
    out.append(f"7. odd v < 2000 with v - 1 NOT a sum of two distinct primes (no source of the form pq): {fails}; their sources of other forms: {[(v, preim(v)) for v in fails]}")
    # 8. Hamiltonian path counts
    for q in (3, 7, 11):
        out.append(f"8. H(P_{q}) = {ham_paths(q, paley(q))}")
    # transitive stacks: H multiplicative; the strong 4-tournament
    adj4 = [[False]*4 for _ in range(4)]
    for (i, j) in [(0,1),(1,2),(2,3),(3,0),(0,2),(1,3)]: adj4[i][j] = True   # one strong 4-tournament
    out.append(f"   H of the strong 4-tournament 0>1>2>3>0, 0>2, 1>3: {ham_paths(4, adj4)}; H(C3) = 3; stack products 3*3 = 9, 3*5 = 15, 5*5 = 25; 7 and 21 are not products of attainable strong values (THM-1370)")
    # 9. perfect numbers as values (exact by the parity lemma: even sources m <= 2N-2, odd sources are odd squares q^2, q <= N-1),
    #    twin-collision values and their full source sets, Mersenne abc triples (1, 2^k - 1, 2^k), refined abc tallies
    import sympy
    def sources_even_value(N):
        """all n with s(n) = N for even N: even n <= 2N-2 (n/2 is a proper divisor) and odd squares q^2 with q <= N-1."""
        if 2 * N - 1 <= LIMIT + 1:
            res = [int(n) for n in np.flatnonzero(s[:2 * N - 1] == N) if n >= 1 and n % 2 == 0]
        else:
            res = [m for m in range(2, 2 * N - 1, 2) if int(sympy.divisor_sigma(m)) - m == N]
        for q in range(3, N, 2):
            if int(sympy.divisor_sigma(q * q)) - q * q == N: res.append(q * q)
        return sorted(res)
    perfect = [6, 28, 496, 8128]
    for N in perfect:
        q = (N - 7) // 3 if (N - 7) % 3 == 0 and N > 7 else None
        out.append(f"9. perfect N = {N}: s^-1({N}) = {sources_even_value(N)}; (N-7)/3 = {q} prime: {sympy.isprime(q) if q else '-'}; N - 3 prime: {sympy.isprime(N - 3)}")
    for pexp in (13, 17, 19, 31):
        N = 2 ** (pexp - 1) * (2 ** pexp - 1); q = (N - 7) // 3
        out.append(f"   perfect 2^{pexp-1}(2^{pexp}-1) = {N}: (N-7)/3 = {q} prime: {sympy.isprime(q)} (so 4q is a source iff prime); N - 3 prime: {sympy.isprime(N-3)}")
    out.append("   row sources 2^a q of an even perfect number N = 2^e M: (2^a - 1) q = N + 2 - 2^(a+1); a = 2 gives q = (N - 7)/3, an integer for every even perfect N > 6 since N = 1 mod 3; a = e gives q = M, the number itself")
    twinvals = [int(p) + 1 for p in primes if 5 <= p <= 9999 and isprime[p - 2]]
    counts = []
    for v in twinvals:
        src = sources_even_value(v); counts.append((v, len(src), src if len(src) <= 4 else src[:4] + ['...']))
    exactly_two = [v for (v, c, src) in counts if c == 2]
    out.append(f"   twin-collision values p+1 <= 10^4: {len(twinvals)}; with exactly the two twin sources: {len(exactly_two)} = {exactly_two}; source counts of the first 14: {counts[:14]}")
    big = [c for (v, c, _) in counts if v > 1000]
    out.append(f"   source counts for twin values in (10^3, 10^4]: min {min(big)}, max {max(big)}, mean {sum(big)/len(big):.2f} over {len(big)} values")
    mh = []
    for k in range(2, 65):
        M = 2 ** k - 1; f = sympy.factorint(M); r = 1
        for pr in f: r *= pr
        if 2 * r < 2 ** k: mh.append((k, round(log(2 ** k) / log(2 * r), 3), {int(a): int(b) for a, b in f.items() if b > 1}))
    out.append(f"   Mersenne abc hits (1, 2^k - 1, 2^k) with k <= 64 (hit iff rad(2^k - 1) < 2^(k-1)): {mh}")
    out.append("   (1, s(n), n) is a coprime additive triple iff sigma(n) = 2n - 1 (almost perfect; the powers of 2 are the only known ones): the row of 1 is where the aliquot function itself produces abc triples")
    odd_pp = sum(1 for h in hits if any(x % 2 == 1 and is_perfect_power(x) for x in (h[1], h[2], h[0])))
    even_pow2 = sum(1 for h in hits if any(x % 2 == 0 and (x & (x - 1)) == 0 for x in (h[1], h[2], h[0])))
    out.append(f"   abc refined tallies (c <= {C}): exactly one term is even in every hit (parity); hits with an odd perfect power among the terms: {odd_pp} of {len(hits)}; hits whose even term is a power of 2: {even_pow2}; hits with 1 as a term: {sum(1 for h in hits if h[1] == 1)}")
    print("\n".join(out))

if __name__ == "__main__":
    main()
