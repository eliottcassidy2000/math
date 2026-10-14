#!/usr/bin/env python3
"""Odd aliquot cycles are frustrated: Dickson's 1913 irregularity theorem read as signed-cycle balance.

If q divides sigma(n_i) for every member of an aliquot cycle n_0 -> n_1 -> ... -> n_(k-1) -> n_0, then
n_(i+1) = sigma(n_i) - n_i = -n_i (mod q): along the cycle the aliquot map is a sign flip mod q, and a closed walk of
odd length k gives n_i = -n_i, i.e. q | 2 n_i (Dickson's argument, as reproved by Kobayashi-Pollack-Pomerance 2009,
Theorem E).  Even cycles carry no such constraint (they alternate n_i = (-1)^i n_0 mod q).  Over Q the cycle equations
n_i + n_(i+1) = sigma(n_i) are the unoriented incidence system of the cycle graph C_k: uniquely solvable (members =
half alternating sums of the sigma values) iff k is odd, i.e. iff C_k is not bipartite.

This script measures it on every aliquot cycle with all members <= N (exact sieve) plus the published longer cycles
of lengths 6, 8 and 9 (their least members from the standard tables, members recomputed here):
  * G = gcd of the members, regular (Cohen) iff G is a unitary divisor of every member;
  * Q = gcd of the sigma(n_i) and its odd part; for odd k the theorem forces odd(Q) | G and Q | 2G (checked);
    for even k the alternation n_i = (-1)^i n_0 mod Q (checked) and how much of odd(Q) is NOT in G
    (that excess is exactly what an odd cycle with the same sigma-divisibility could not afford);
  * the decade counts of the double flips (odd a with s(a^2) a square or twice a square) from the census JSON,
    against the Iannucci-Luca bound x / (log x)^(3/2 + o(1)).

Session opus-2026-10-14-S33.  Note: 05-knowledge/results/number_tournament_bridges_20261014.md
Reproduce: python3 04-computation/experiments/aliquot_dickson_frustration_20261014.py [N=1e8]   (about 1 minute, 2 GB)
"""
import sys, time, math, json
from collections import Counter
import numpy as np
import numba as nb
import sympy

N = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10 ** 8
FAIL = []


def check(label, ok):
    if not ok:
        FAIL.append(label)
    return ok


@nb.njit(cache=True)
def sigma_sieve(N):
    s = np.zeros(N + 1, dtype=np.int64)
    for d in range(1, N // 2 + 1):
        for m in range(2 * d, N + 1, d):
            s[m] += d                 # s[m] = sum of proper divisors = sigma(m) - m
    return s


@nb.njit(cache=True)
def find_cycles(N, s, out, lens):
    state = np.zeros(N + 1, dtype=np.int8)      # 0 unvisited, 1 on current path, 2 done
    path = np.zeros(1000, dtype=np.int64)
    nc = 0
    pos = 0
    for n0 in range(1, N + 1):
        if state[n0] != 0:
            continue
        k = 0
        x = n0
        while x >= 1 and x <= N and state[x] == 0 and k < 1000:
            state[x] = 1
            path[k] = x
            k += 1
            x = s[x]
        if x >= 1 and x <= N and state[x] == 1:
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


def sigma_py(n):
    tot = 1
    for p, e in sympy.factorint(n).items():
        tot *= (p ** (e + 1) - 1) // (p - 1)
    return tot


def cycle_from(n0, maxlen=40):
    cyc = [n0]
    x = sigma_py(n0) - n0
    while x != n0 and len(cyc) < maxlen:
        cyc.append(x)
        x = sigma_py(x) - x
    return cyc if x == n0 else None


def odd_part(n):
    while n % 2 == 0:
        n //= 2
    return n


def analyse(cyc):
    k = len(cyc)
    G = 0
    for v in cyc:
        G = math.gcd(G, v)
    regular = all(math.gcd(G, v // G) == 1 for v in cyc)
    sig = [cyc[i] + cyc[(i + 1) % k] for i in range(k)]          # sigma(n_i) = n_i + n_(i+1)
    Q = 0
    for v in sig:
        Q = math.gcd(Q, v)
    Qo = odd_part(Q)
    if k % 2 == 1:
        ok = (Qo != 0 and G % Qo == 0) and ((2 * G) % Q == 0)
    else:
        ok = all((cyc[i] - (-1) ** i * cyc[0]) % Q == 0 for i in range(k))
    excess = Qo // math.gcd(Qo, G)               # part of odd(Q) not dividing the members
    return dict(k=k, G=G, regular=regular, Q=Q, Qodd=Qo, excess=excess, ok=ok,
                parity='even' if cyc[0] % 2 == 0 else 'odd')


def main():
    t0 = time.time()
    s = sigma_sieve(N)
    out = np.zeros(200000, dtype=np.int64)
    lens = np.zeros(20000, dtype=np.int64)
    nc, pos = find_cycles(N, s, out, lens)
    cycles = []
    p = 0
    for i in range(nc):
        L = int(lens[i])
        cycles.append([int(v) for v in out[p:p + L]])
        p += L
    print(f'aliquot cycles with all members <= {N:.0e}: {len(cycles)} ({time.time() - t0:.1f} s); '
          f'by length {sorted(Counter(len(c) for c in cycles).items())}')
    check('census count', len(cycles) == 248 if N == 10 ** 8 else True)

    # some published longer cycles (starting members from the standard lists of known sociable cycles; the cycle is
    # recomputed and closed here, so a wrong starting value would fail the check)
    extra = {}
    for n0 in (21548919483, 90632826380, 1771417411016, 4773123705616,     # 6-cycles
               1095447416, 1276254780,                                       # 8-cycles
               805984760):                                                   # the 9-cycle
        c = cycle_from(n0, maxlen=12)
        if c is None:
            print(f'   {n0}: no cycle of length <= 12 found')
            check(f'cycle {n0}', False)
            continue
        extra[n0] = c
    print('   published longer cycles recomputed:', {n0: len(c) for n0, c in extra.items()})
    allc = cycles + [c for c in extra.values() if max(c) > N]

    rows = [analyse(c) for c in allc]
    print('\n== Dickson / signed-cycle balance on every cycle')
    print('   theorem checks (odd k: odd(Q) | G and Q | 2G;  even k: n_i = (-1)^i n_0 mod Q):',
          check('dickson', all(r['ok'] for r in rows)))
    by = {}
    for r in rows:
        by.setdefault(r['k'], []).append(r)
    for k in sorted(by):
        rs = by[k]
        nreg = sum(r['regular'] for r in rs)
        if k == 1:
            print(f'   k = 1 (perfect): {len(rs)} cycles, Q = 2n, G = n trivially')
            continue
        exc = Counter(r['excess'] for r in rs)
        big = sum(1 for r in rs if r['excess'] > 1)
        print(f'   k = {k:2d}: {len(rs):3d} cycles; regular (gcd unitary in every member) {nreg}, irregular {len(rs) - nreg};'
              f'  odd part of gcd(sigma) not dividing the members (impossible for odd k): {big} of {len(rs)} cycles'
              + (f'; values {sorted(exc.items())[:8]}' if k > 2 else ''))
        if k % 2 == 1:
            check(f'odd k={k} irregular', nreg == 0)
            for r in rs:
                print(f'      G = {r["G"]} = {sympy.factorint(r["G"])}, Q = gcd sigma = {r["Q"]} = {sympy.factorint(r["Q"])}, '
                      f'Q | 2G: {(2 * r["G"]) % r["Q"] == 0}, regular: {r["regular"]}')
    am = by.get(2, [])
    if am:
        qo = [r['excess'] for r in am]
        print(f'   amicable pairs: median of odd(Q)/gcd(odd(Q), G) = {int(np.median(qo))}, '
              f'share with an excess > 1: {sum(1 for v in qo if v > 1) / len(qo):.3f}; '
              f'share with 3 | excess: {sum(1 for v in qo if v % 3 == 0) / len(qo):.3f}')
    # members' sigma divisibility: how often q | sigma(n) for cycle members vs all n <= N (typicality)
    print('\n== how often q | sigma(n): cycle members vs all n <= N (sigma(n) = s(n) + n)')
    idx = np.arange(N + 1, dtype=np.int64)
    sig_all = s + idx
    mem = sorted({v for c in cycles for v in c if c and len(c) > 1})
    for q in (3, 4, 5, 7, 8, 9, 12, 24):
        frac_all = float(np.mean(sig_all[2:] % q == 0))
        frac_mem = float(np.mean([(s[v] + v) % q == 0 for v in mem]))
        print(f'   q = {q:2d}: all n {frac_all:.4f}, members of 2- and 4-cycles {frac_mem:.4f}')

    # double flips by decade (from the census JSON)
    print('\n== double flips by decade (odd a with s(a^2) in squares or twice squares), census JSON')
    try:
        df = json.load(open('05-knowledge/results/aliquot_parity_squares_20261014.json'))['double_flips']
        cum = []
        for e in range(1, 10):
            X = 10 ** e
            allc_ = sum(1 for a, m, t in df if a <= X)
            comp = sum(1 for a, m, t in df if a <= X and not sympy.isprime(a))
            bound = X / math.log(X) ** 1.5
            cum.append((e, allc_, comp))
            print(f'   a <= 10^{e}: {allc_:5d} double flips, {comp:2d} composite;  x/(log x)^1.5 = {bound:.3g}')
        check('double flips 1e9', cum[-1][1:] == (4476, 18))
    except FileNotFoundError:
        print('   census JSON not found (run aliquot_parity_squares_20261014.py first)')

    print('\nFAILED CHECKS:', FAIL if FAIL else 'none')
    print('ALL CHECKS PASSED' if not FAIL else 'SOME CHECKS FAILED')


if __name__ == '__main__':
    main()
