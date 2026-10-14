#!/usr/bin/env python3
"""Aliquot parity flips only at squares and twice-squares; every aliquot cycle with members <= N (exact census),
its parity and abundancy pattern; and the 'double flips': odd a <= A with s(a^2) a square or twice a square, which are
the only possible odd members of opposite-parity amicable pairs.

Session opus-2026-10-14-S33.  Note: 05-knowledge/results/number_tournament_bridges_20261014.md
Reproduce: python3 04-computation/experiments/aliquot_parity_squares_20261014.py [N=1e8] [A=1e9]   (about 2 minutes, 5 GB)
"""
import sys, time, math
import numpy as np
import numba as nb

N = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10 ** 8
A = int(float(sys.argv[2])) if len(sys.argv) > 2 else 10 ** 9


@nb.njit(cache=True)
def sigma_sieve(N):
    s = np.zeros(N + 1, dtype=np.int64)
    for d in range(1, N // 2 + 1):
        for m in range(2 * d, N + 1, d):
            s[m] += d                 # s[m] = sum of proper divisors = sigma(m) - m
    return s


@nb.njit(cache=True)
def find_cycles(N, s, out, lens):
    # walks are followed without a step cap (audit I14: a capped walk marked done could hide a long-entry cycle);
    # the path buffer has N + 1 slots, which no walk inside [1, N] can exceed
    state = np.zeros(N + 1, dtype=np.int8)      # 0 unvisited, 1 on current path, 2 done
    path = np.zeros(N + 1, dtype=np.int64)
    nc = 0
    pos = 0
    for n0 in range(1, N + 1):
        if state[n0] != 0:
            continue
        k = 0
        x = n0
        while x >= 1 and x <= N and state[x] == 0:
            state[x] = 1
            path[k] = x
            k += 1
            x = s[x]
        if x >= 1 and x <= N and state[x] == 1:
            # cycle found: from the position of x in path
            j = 0
            while path[j] != x:
                j += 1
            L = k - j
            if nc < lens.shape[0]:
                lens[nc] = L
                for t in range(L):
                    if pos + t < out.shape[0]:
                        out[pos + t] = path[j + t]
                pos += L
                nc += 1
        for t in range(k):
            state[path[t]] = 2
    return nc, pos


@nb.njit(cache=True)
def is_square(m):
    if m < 0:
        return False
    r = np.int64(math.sqrt(m))
    while r * r > m:
        r -= 1
    while (r + 1) * (r + 1) <= m:
        r += 1
    return r * r == m


@nb.njit(cache=True)
def spf_sieve(A):
    spf = np.zeros(A + 1, dtype=np.int32)
    for i in range(2, A + 1):
        if spf[i] == 0:
            spf[i] = i
            if i <= A // i:
                for m in range(i * i, A + 1, i):
                    if spf[m] == 0:
                        spf[m] = i
    return spf


@nb.njit(cache=True)
def double_flips(A, spf, out_a, out_m):
    """odd a <= A: m = sigma(a^2) - a^2; record a if m is a square or twice a square"""
    cnt = 0
    for a in range(3, A + 1, 2):
        x = a
        sig = np.int64(1)
        while x > 1:
            p = np.int64(spf[x])
            e = 0
            while x % p == 0:
                x //= p
                e += 1
            # sigma(p^(2e)) = 1 + p + ... + p^(2e)
            t = np.int64(1)
            pw = np.int64(1)
            for _ in range(2 * e):
                pw *= p
                t += pw
            sig *= t
        a2 = np.int64(a) * np.int64(a)
        m = sig - a2
        if is_square(m) or (m % 2 == 0 and is_square(m // 2)):
            if cnt < out_a.shape[0]:
                out_a[cnt] = a
                out_m[cnt] = m
            cnt += 1
    return cnt


def sigma_py(n):
    import sympy
    tot = 1
    for p, e in sympy.factorint(n).items():
        tot *= (p ** (e + 1) - 1) // (p - 1)
    return tot


def main():
    t0 = time.time()
    s = sigma_sieve(N)
    print(f'aliquot sieve to {N:.0e} in {time.time() - t0:.1f} s')
    out = np.zeros(200000, dtype=np.int64)
    lens = np.zeros(20000, dtype=np.int64)
    nc, pos = find_cycles(N, s, out, lens)
    cycles = []
    p = 0
    for i in range(nc):
        L = int(lens[i])
        cyc = [int(v) for v in out[p:p + L]]
        p += L
        cycles.append(cyc)
    from collections import Counter
    print(f'cycles with all members <= N: {nc};  by length: {sorted(Counter(len(c) for c in cycles).items())}')
    def sq_or_twice(v):
        r = math.isqrt(v)
        if r * r == v:
            return True
        if v % 2 == 0:
            r = math.isqrt(v // 2)
            return r * r == v // 2
        return False
    mixed = [c for c in cycles if len({v % 2 for v in c}) > 1]
    sqm = [c for c in cycles if any(sq_or_twice(v) for v in c)]
    odd_cycles = [c for c in cycles if all(v % 2 for v in c)]
    print(f'mixed-parity cycles: {len(mixed)};  cycles with a member in squares or twice squares: {len(sqm)};  all-odd cycles: {len(odd_cycles)}')
    for c in cycles:
        if len(c) >= 3:
            print(f'   length {len(c)}: min member {min(c)}, parity {"odd" if c[0] % 2 else "even"}, members {c if len(c) <= 9 else str(c[:6]) + " ..."}')
    # collective perfection: sum sigma = 2 sum n
    ok = all(sum(int(s[v]) + v for v in c) == 2 * sum(c) for c in cycles)
    print(f'every cycle satisfies sum sigma(n_i) = 2 sum n_i: {ok}')
    # abundant/deficient pattern
    for c in cycles:
        if len(c) >= 4 and len(c) <= 9:
            pat = ''.join('A' if s[v] > v else ('D' if s[v] < v else 'P') for v in c)
            print(f'   pattern of the {len(c)}-cycle from {min(c)}: {pat}')
    del s
    print(f'\ndouble flips: odd a <= {A:.0e} with s(a^2) a square or twice a square')
    t1 = time.time()
    spf = spf_sieve(A)
    print(f'   spf sieve in {time.time() - t1:.1f} s')
    out_a = np.zeros(100000, dtype=np.int64)
    out_m = np.zeros(100000, dtype=np.int64)
    cnt = double_flips(A, spf, out_a, out_m)
    print(f'   {cnt} odd a found in {time.time() - t1:.1f} s')
    amicable = []
    rows = []
    for a, m in zip(out_a[:cnt], out_m[:cnt]):
        a = int(a)
        m = int(m)
        n = a * a
        sn = sigma_py(n)
        sm = sigma_py(m) if m > 0 else 0
        kind = 'square' if math.isqrt(m) ** 2 == m else 'twice a square'
        rows.append((a, n, m, kind, sm == sn))
        if sm == sn:
            amicable.append((n, m))
    import sympy
    prime_rows = [r for r in rows if sympy.isprime(r[0])]
    pp_rows = [r for r in rows if not sympy.isprime(r[0]) and len(sympy.factorint(r[0])) == 1]
    comp_rows = [r for r in rows if len(sympy.factorint(r[0])) > 1]
    print(f'   by type of a: prime {len(prime_rows)} (then s(a^2) = a + 1, far below a^2), prime power {len(pp_rows)}, '
          f'two or more primes {len(comp_rows)}')
    print(f'   ratio s(a^2)/a^2 over all double flips: max {max((r[2] / r[1] for r in rows), default=0):.4f} '
          f'(amicability needs sigma(m)/m = 1 + a^2/m, so m must be comparable to a^2)')
    for a, n, m, kind, ami in pp_rows + comp_rows:
        print(f'   a = {a:>10d} = {sympy.factorint(a)}: s(a^2) = {m} ({kind}), s/a^2 = {m / n:.4f};  sigma(m) = sigma(a^2)? {ami}')
    print(f'   opposite-parity amicable pairs with odd member <= {A}^2: {amicable}')
    import json
    json.dump({'cycles_by_length': sorted(Counter(len(c) for c in cycles).items()), 'double_flips': [(r[0], r[2], r[3]) for r in rows]},
              open('05-knowledge/results/aliquot_parity_squares_20261014.json', 'w'))


if __name__ == '__main__':
    main()
