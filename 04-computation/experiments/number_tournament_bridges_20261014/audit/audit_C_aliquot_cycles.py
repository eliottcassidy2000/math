#!/usr/bin/env python3
"""Independent audit of the aliquot census (referee script).

sigma(n) for n <= 1e8 by a LINEAR (multiplicative) sieve; cross-checked against a SEGMENTED divisor-pair sieve that
also extends to M = max s(n) (n <= 1e8), which decides the amicable pairs whose larger member exceeds 1e8.
Cycle detection by start-stamps (no step cap).  Then: classification, parity, squares / twice-squares, abundancy
patterns, the cycle identities, the Dickson / frustration census (G, Q, regularity, odd excess, both readings of
"odd prime divides every sigma but not every member"), sigma-divisibility percentages, and the published longer
cycles recomputed from their seeds by an independent trial-division sigma.
"""
import sys, time, math, json
from collections import Counter
import numpy as np
import numba as nb

N = 10 ** 8


@nb.njit(cache=True)
def linear_sigma(N):
    sig = np.zeros(N + 1, dtype=np.int64)
    pk = np.zeros(N + 1, dtype=np.int32)          # largest power of the smallest prime factor dividing n
    comp = np.zeros(N + 1, dtype=np.uint8)
    primes = np.zeros(6000000, dtype=np.int32)
    npr = 0
    sig[1] = 1
    pk[1] = 1
    for i in range(2, N + 1):
        if comp[i] == 0:
            primes[npr] = i
            npr += 1
            sig[i] = i + 1
            pk[i] = i
        for t in range(npr):
            p = primes[t]
            if np.int64(p) * i > N:
                break
            m = p * i
            comp[m] = 1
            if i % p == 0:
                pk[m] = pk[i] * p
                q = pk[i]
                sig[m] = sig[i // q] * (sig[q] * p + 1)
                break
            else:
                pk[m] = p
                sig[m] = sig[i] * (p + 1)
    return sig


@nb.njit(parallel=True, cache=True)
def segmented_s_u32(M, block):
    """s(m) = sigma(m) - m for 1 <= m <= M, by pairing divisors d <= sqrt(m) with m/d, block by block"""
    out = np.zeros(M + 1, dtype=np.uint32)
    nblocks = (M + block) // block
    for bi in nb.prange(nblocks):
        L = bi * block
        R = min(L + block, M + 1)          # [L, R)
        if L == 0:
            L1 = 1
        else:
            L1 = L
        acc = np.zeros(R - L, dtype=np.int64)
        dmax = int(math.sqrt(R)) + 2
        for d in range(1, dmax + 1):
            dd = np.int64(d) * d
            if dd >= R:
                break
            # multiples m of d with m >= max(L1, d*d)
            start = max(L1, dd)
            r = start % d
            if r != 0:
                start += d - r
            for m in range(start, R, d):
                q = m // d
                if q == d:
                    acc[m - L] += d
                else:
                    acc[m - L] += d + q
        for m in range(L1, R):
            out[m] = np.uint32(acc[m - L] - m)
    return out


@nb.njit(cache=True)
def to_s_inplace(sig):
    for i in range(sig.shape[0]):
        sig[i] -= i
    sig[0] = 0


@nb.njit(cache=True)
def agree_check(s, s2, N):
    bad = 0
    for i in range(1, N + 1):
        if np.int64(s2[i]) != s[i]:
            bad += 1
    return bad


@nb.njit(cache=True)
def amicable_scan(s, s2, N, Mmax):
    cnt = 0
    res = []
    for n in range(1, N):
        m = s[n]
        if m > n and m <= Mmax:
            if np.int64(s2[m]) == n:
                res.append((n, m))
    return res


@nb.njit(cache=True)
def find_cycles(N, s):
    stamp = np.zeros(N + 1, dtype=np.int32)
    starts = []
    for n0 in range(1, N + 1):
        if stamp[n0] != 0:
            continue
        x = n0
        while x >= 1 and x <= N and stamp[x] == 0:
            stamp[x] = n0
            x = s[x]
        if x >= 1 and x <= N and stamp[x] == n0:
            starts.append(x)
    return starts


def odd_part(n):
    while n % 2 == 0:
        n //= 2
    return n


def in_sq_or_2sq(v):
    r = math.isqrt(v)
    if r * r == v:
        return True
    if v % 2 == 0:
        r = math.isqrt(v // 2)
        return r * r == v // 2
    return False


def trial_sigma(n, primes):
    tot = 1
    m = n
    for p in primes:
        if p * p > m:
            break
        if m % p == 0:
            e = 1
            pe = p
            while m % p == 0:
                m //= p
                pe *= p
            tot *= (pe - 1) // (p - 1)
    if m > 1:
        if m >= primes[-1] ** 2:
            raise ValueError('cofactor too large for the trial-division bound')
        tot *= m + 1
    return tot


def small_primes(B):
    sieve = bytearray([1]) * (B + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(B ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    return [i for i in range(B + 1) if sieve[i]]


def prime_factors(n, primes):
    fs = []
    m = n
    for p in primes:
        if p * p > m:
            break
        if m % p == 0:
            fs.append(p)
            while m % p == 0:
                m //= p
    if m > 1:
        fs.append(m)
    return fs


def main():
    out = []
    t0 = time.time()
    sig = linear_sigma(N)
    out.append(f'linear sigma sieve to {N:.0e}: {time.time() - t0:.0f} s')
    s = sig.copy()
    to_s_inplace(s)
    Mmax = int(s[1:].max())
    out.append(f'max s(n) for n <= 1e8: {Mmax} (at n = {int(np.argmax(s))})')
    t0 = time.time()
    s2 = segmented_s_u32(Mmax, 1 << 20)
    out.append(f'segmented divisor-pair sieve to {Mmax}: {time.time() - t0:.0f} s')
    nbad = agree_check(s, s2, N)
    out.append(f'two sieves agree on 1..1e8: {nbad == 0} ({nbad} disagreements)')
    pairs = [(int(a), int(b)) for a, b in amicable_scan(s, s2, N, Mmax)]
    big = [(a, b) for a, b in pairs if b > N]
    out.append(f'amicable pairs with smaller member < 1e8: {len(pairs)}; of these, larger member > 1e8: {len(big)} {big}')
    oddpairs = [(a, b) for a, b in pairs if a % 2 == 1 and b <= N]
    out.append(f'odd-odd amicable pairs with both members <= 1e8: {len(oddpairs)}')
    del s2

    # ---- cycles
    t0 = time.time()
    starts = find_cycles(N, s)
    cycles = []
    for x in starts:
        c = [int(x)]
        y = int(s[x])
        while y != x:
            c.append(y)
            y = int(s[y])
        cycles.append(c)
    out.append(f'cycles with all members <= 1e8: {len(cycles)} ({time.time() - t0:.0f} s); by length {sorted(Counter(len(c) for c in cycles).items())}')
    mixed = [c for c in cycles if len({v % 2 for v in c}) > 1]
    sq = [c for c in cycles if any(in_sq_or_2sq(v) for v in c)]
    odd_c = [c for c in cycles if all(v % 2 for v in c)]
    out.append(f'mixed parity: {len(mixed)}; with a member in squares or twice squares: {len(sq)}; all-odd cycles: {len(odd_c)}'
               f' (lengths {sorted(Counter(len(c) for c in odd_c).items())})')
    two = [c for c in cycles if len(c) == 2]
    out.append(f'amicable pairs with both members <= 1e8 from the cycle census: {len(two)}; consistent with pair list: '
               f'{len(two) == len([1 for a, b in pairs if b <= N])}')
    for c in cycles:
        if len(c) >= 3:
            k = len(c)
            j = c.index(min(c))
            c2 = c[j:] + c[:j]
            pat = ''.join('A' if s[v] > v else ('D' if s[v] < v else 'P') for v in c2)
            out.append(f'   length {k}: min {min(c)}, parity {"odd" if c[0] % 2 else "even"}, pattern from min: {pat}, members {c2 if k <= 9 else str(c2[:5]) + "..."}')
    # identities
    ok1 = all(sum(int(sig[v]) for v in c) == 2 * sum(c) for c in cycles)
    ok2 = all(math.prod(int(sig[v]) for v in c) >= 2 ** len(c) * math.prod(c) for c in cycles)
    eq2 = [c for c in cycles if math.prod(int(sig[v]) for v in c) == 2 ** len(c) * math.prod(c)]
    ad = all((any(s[v] > v for v in c) and any(s[v] < v for v in c)) for c in cycles if len(c) > 1)
    out.append(f'sum sigma = 2 sum n on every cycle: {ok1}; prod sigma/n >= 2^k: {ok2}; equality only for k = 1: {all(len(c) == 1 for c in eq2)};'
               f' every nontrivial cycle has an abundant and a deficient member: {ad}')
    # abundancy patterns up to rotation for 4-cycles
    def canon_rot(pat):
        return min(pat[i:] + pat[:i] for i in range(len(pat)))
    pats4 = Counter(canon_rot(''.join('A' if s[v] > v else 'D' for v in c)) for c in cycles if len(c) == 4)
    out.append(f'4-cycle abundancy patterns up to rotation: {dict(pats4)}')

    # ---- frustration census
    primes = small_primes(3_000_000)
    seeds = (21548919483, 90632826380, 1771417411016, 4773123705616, 1095447416, 1276254780, 805984760)
    extra = {}
    for n0 in seeds:
        c = [n0]
        x = trial_sigma(n0, primes) - n0
        while x != n0 and len(c) < 40:
            c.append(x)
            x = trial_sigma(x, primes) - x
        extra[n0] = c if x == n0 else None
    out.append('published cycles from seeds: ' + str({n0: (len(c) if c else None) for n0, c in extra.items()}))
    for n0, c in extra.items():
        if c:
            out.append(f'   seed {n0}: min member {min(c)}, max {max(c)}, parity {"odd" if c[0] % 2 else "even"}, all <= 1e8: {max(c) <= N}')
    allc = cycles + [c for c in extra.values() if c and max(c) > N]
    out.append(f'cycles in the frustration census: {len(allc)}')

    rows = []
    for c in allc:
        k = len(c)
        G = 0
        for v in c:
            G = math.gcd(G, v)
        sg = [c[i] + c[(i + 1) % k] for i in range(k)]
        # verify sigma(n_i) = n_i + n_(i+1) with an independent sigma for the extra cycles
        Q = 0
        for v in sg:
            Q = math.gcd(Q, v)
        regular = all(math.gcd(G, v // G) == 1 for v in c)
        Qo = odd_part(Q)
        excess = Qo // math.gcd(Qo, G)
        oddprimes_Q = [p for p in prime_factors(Qo, primes)] if Qo > 1 else []
        prime_not_in_G = [p for p in oddprimes_Q if G % p != 0]
        if k % 2 == 1:
            lemma = (G % Qo == 0) and ((2 * G) % Q == 0)
            pollack = (G % Q == 0) or (Q % 4 == 2 and G % (Q // 2) == 0)
        else:
            lemma = all((c[i] - (-1) ** i * c[0]) % Q == 0 for i in range(k))
            pollack = None
        rows.append(dict(k=k, G=G, Q=Q, regular=regular, excess=excess, prime_not_in_G=bool(prime_not_in_G), lemma=lemma,
                         pollack=pollack, c=c))
    out.append(f'lemma holds on every cycle: {all(r["lemma"] for r in rows)}; Pollack-form (Q | G, or 2||Q and Q/2 | G) on odd cycles: '
               f'{[r["pollack"] for r in rows if r["k"] % 2 == 1]}')
    out.append('   k : cycles, regular, odd(Q) not dividing G (author\'s "excess > 1"), an odd PRIME dividing Q but not G')
    byk = {}
    for r in rows:
        byk.setdefault(r['k'], []).append(r)
    for k in sorted(byk):
        rs = byk[k]
        out.append(f'   {k:2d}: {len(rs):3d}, {sum(r["regular"] for r in rs):3d}, {sum(r["excess"] > 1 for r in rs):3d}, '
                   f'{sum(r["prime_not_in_G"] for r in rs):3d}')
        if k in (5, 9):
            for r in rs:
                out.append(f'        G = {r["G"]}, Q = {r["Q"]}')
    am = byk.get(2, [])
    ex = [r['excess'] for r in am]
    out.append(f'amicable: median excess {np.median(ex)}, share excess > 1 {np.mean([e > 1 for e in ex]):.4f}, share 3 | excess {np.mean([e % 3 == 0 for e in ex]):.4f}')
    # examples where the readings differ
    diff = [(r['k'], r['c'][:2], r['G'], r['Q']) for r in rows if (r['excess'] > 1) != r['prime_not_in_G']]
    out.append(f'cycles where "odd(Q) does not divide G" and "an odd prime divides Q but not G" differ: {len(diff)}; first few {diff[:5]}')

    # ---- sigma divisibility
    sig_all = sig[2:]
    mem24 = sorted({v for c in cycles if len(c) in (2, 4) for v in c})
    memall = sorted({v for c in cycles if len(c) > 1 for v in c})
    for q in (3, 9):
        fa = float(np.mean(sig_all % q == 0))
        f24 = float(np.mean([int(sig[v]) % q == 0 for v in mem24]))
        fall = float(np.mean([int(sig[v]) % q == 0 for v in memall]))
        out.append(f'q = {q}: all 2..1e8 {fa:.4f}; members of 2- and 4-cycles ({len(mem24)}) {f24:.4f}; members of all nontrivial cycles ({len(memall)}) {fall:.4f}')
    json.dump({'cycles': cycles, 'extra': {str(k): v for k, v in extra.items()}},
              open('audit_C_cycles.json', 'w'))
    print('\n'.join(out))


if __name__ == '__main__':
    main()
