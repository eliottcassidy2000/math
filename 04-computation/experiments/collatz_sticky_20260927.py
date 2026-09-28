#!/usr/bin/env python3
"""STICKY as a barrier: size-coupled persistence of the aliquot map's 2-adic driver
against the size-free memorylessness of Collatz; random-model controls; re-check of
the nine-iteration typology of the parallel session (collatz_two_carries_typology).

Session: opus, collatz-poset-dag-20260927 (S16), 2026-09-27.
Exact integer arithmetic for all arithmetic facts; numpy sieves for divisor sums;
sympy.factorint for sampled large ranges.

Run: python 04-computation/experiments/collatz_sticky_20260927.py
"""
from __future__ import annotations

import math
import random
import sys
from collections import Counter, defaultdict

import numpy as np
from sympy import factorint, divisor_sigma


def v2(x: int) -> int:
    return (x & -x).bit_length() - 1 if x else 10 ** 9


def sigma_sieve(N: int) -> np.ndarray:
    s = np.zeros(N + 1, dtype=np.int64)
    for d in range(1, N + 1):
        s[d::d] += d
    return s


def v2_sigma_odd(m: int) -> int:
    """v_2(sigma(m)) for odd m from its factorization: for p^e || m, e odd,
    v_2(sigma(p^e)) = v_2(p+1) + v_2(e+1) - 1; e even contributes 0."""
    tot = 0
    for p, e in factorint(m).items():
        if e % 2 == 1:
            tot += v2(p + 1) + v2(e + 1) - 1
    return tot


# ----------------------------------------------------------------------------
# P1: aliquot driver persistence by size
# ----------------------------------------------------------------------------

def part1():
    print("== P1: aliquot driver persistence P(v2 s(n) = a | v2 n = a) and loss P(v2 s(n) < a | v2 n = a) by size ==")
    print("   exact ranges [N, 2N) by sieve; sampled ranges by factorization")
    ranges_exact = [10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6]
    sig = sigma_sieve(2 * 10 ** 6)
    # exact characterization checks on [2, 2*10^6]: v2(s(n)) = min(a, v2 sigma(m)) if differ, >= a+1 if equal;
    # loss (v2 s < a) iff v2 sigma(m) <= a-1; for a = 1 loss iff m is an odd square
    checked = 0
    for n in range(2, 2 * 10 ** 6 + 1, 2):
        a = v2(n)
        m = n >> a
        sm = int(sig[m])
        vs = v2(sm)
        s = int(sig[n]) - n
        if s == 0:
            continue
        v = v2(s)
        if vs != a:
            assert v == min(a, vs), (n, a, vs, v)
        else:
            assert v >= a + 1, (n, a, vs, v)
        assert (v < a) == (vs <= a - 1)
        if a == 1:
            r = math.isqrt(m)
            assert (v < 1) == (r * r == m)
        checked += 1
    print(f"   exact identities checked on {checked} even n <= 2*10^6 (v2 s = min(a, v2 sigma m) / >= a+1; loss iff v2 sigma(m) <= a-1; a=1 loss iff odd square)")
    rows = []
    for N in ranges_exact:
        cnt = Counter()
        keep = Counter()
        loss = Counter()
        for n in range(N, 2 * N, 2):
            a = v2(n)
            if a > 4:
                continue
            s = int(sig[n]) - n
            v = v2(s) if s else 99
            cnt[a] += 1
            if v == a:
                keep[a] += 1
            if v < a:
                loss[a] += 1
        rows.append((N, cnt, keep, loss))
    for N in (10 ** 7, 10 ** 8, 10 ** 9):
        rng = random.Random(20260927 + N)
        cnt = Counter()
        keep = Counter()
        loss = Counter()
        for _ in range(20000):
            a = rng.choice([1, 1, 2, 2, 3, 4])  # oversample driver classes
            m = rng.randrange(N // 2 ** a, 2 * N // 2 ** a) | 1
            vs = v2_sigma_odd(m)
            # v2(s) from the identity (exact, no need to compute s)
            if vs != a:
                v = min(a, vs)
            else:
                v = a + 1  # >= a+1; only the comparison with a matters
            cnt[a] += 1
            if v == a:
                keep[a] += 1
            if v < a:
                loss[a] += 1
        rows.append((N, cnt, keep, loss))
    print("   N        | a=1 keep  loss (N^-1/2 pred) | a=2 keep  loss  loss*lnN | a=3 keep  loss  loss*lnN/lnlnN | a=4 keep loss")
    for N, cnt, keep, loss in rows:
        lnN = math.log(N)
        l1 = loss[1] / cnt[1]
        pred1 = math.sqrt(2) / math.sqrt(N) * (math.sqrt(2) - 1) / (2 * (math.sqrt(2) - 1))  # placeholder, replaced below
        # a=1 loss = P(m odd square) for m in [N/2, N): #odd squares in [N/2, N) / #odd m in [N/2, N)
        lo, hi = N // 2, N
        nsq = sum(1 for r in range(math.isqrt(lo - 1) + 1, math.isqrt(hi - 1) + 1) if r % 2 == 1)
        nodd = (hi - lo) // 2
        pred1 = nsq / nodd
        l2 = loss[2] / cnt[2]
        l3 = loss[3] / cnt[3]
        l4 = loss[4] / cnt[4] if cnt[4] else float('nan')
        print(f"   {N:>9} | {keep[1]/cnt[1]:.4f}  {l1:.5f} ({pred1:.5f})  | {keep[2]/cnt[2]:.4f}  {l2:.4f}  {l2*lnN:.3f} | "
              f"{keep[3]/cnt[3]:.4f}  {l3:.4f}  {l3*lnN/math.log(lnN):.3f} | {keep[4]/cnt[4] if cnt[4] else float('nan'):.4f} {l4:.4f}")
    return sig


# ----------------------------------------------------------------------------
# P2: Collatz is size-free and exactly memoryless under counting measure
# ----------------------------------------------------------------------------

def U(m: int) -> tuple[int, int]:
    x = 3 * m + 1
    v = v2(x)
    return x >> v, v


def part2():
    print("\n== P2: Collatz: P(next v = 1 | v = 1) is exactly 1/2 at every size (n = 7 mod 8 among n = 3 mod 4) ==")
    for N in (10 ** 3, 10 ** 5, 10 ** 7, 10 ** 9):
        c3 = c7 = 0
        # exact count over odd n in [N, 2N)
        lo = N | 1
        for r in (3, 7):
            pass
        n3 = len(range((N + (3 - N) % 4), 2 * N, 4))
        n7 = len(range((N + (7 - N) % 8), 2 * N, 8))
        print(f"   [N, 2N) with N = {N:>10}: #n=3 mod 4: {n3}, #n=7 mod 8: {n7}, ratio {n7/n3:.6f}")
    # along orbits: persistence 0.52 comes from small values
    tot = Counter()
    keep = Counter()
    for n in range(1, 2 * 10 ** 5, 2):
        m = n
        prev = None
        while m != 1:
            m2, v = U(m)
            if prev is not None:
                band = 'small(<10^4)' if m < 10 ** 4 else 'large(>=10^4)'
                if prev == 1:
                    tot[band] += 1
                    if v == 1:
                        keep[band] += 1
            prev = v
            m = m2
    for band in sorted(tot):
        print(f"   along orbits of odd n < 2*10^5, steps with current value {band}: P(next v=1 | v=1) = {keep[band]/tot[band]:.4f} ({tot[band]} steps)")


# ----------------------------------------------------------------------------
# P3: random-model controls
# ----------------------------------------------------------------------------

def part3():
    print("\n== P3: random models: stationary memory does not change the fate; size-coupled persistence can ==")
    rng = random.Random(1)
    L3 = math.log2(3)

    def geom():
        v = 1
        while rng.random() < 0.5:
            v += 1
        return v

    for p in (0.0, 0.5, 0.9, 0.99):
        # valuation chain: with prob p keep the previous valuation, else fresh geometric(1/2);
        # stationary law is geometric(1/2) for every p
        slopes = []
        maxh = []
        for _ in range(200):
            v = geom()
            h = 0.0
            mx = 0.0
            for step in range(4000):
                if rng.random() >= p:
                    v = geom()
                h += L3 - v
                mx = max(mx, h)
            slopes.append(h / 4000)
            maxh.append(mx)
        print(f"   persistence p={p:.2f}: mean slope {sum(slopes)/len(slopes):+.4f} bits/step (i.i.d. value log2 3 - 2 = {L3-2:+.4f}); "
              f"mean max height {sum(maxh)/len(maxh):.1f}; paths ending above start: {sum(1 for s in slopes if s > 0)}/200")
    # two-state size-coupled model: growth state multiplies by 2^0.3 per step, escape probability p(N)
    # escape probabilities written in log2 N to avoid overflow: c/sqrt(N) = c 2^(-logN/2), c/ln N = c/(logN ln 2)
    for law, f in (("c/sqrt(N), c=3", lambda logN: 3.0 * 2.0 ** (-logN / 2)), ("c/ln N, c=0.6", lambda logN: 0.6 / (logN * math.log(2)))):
        survive = 0
        trials = 2000
        for _ in range(trials):
            logN = 20.0  # log2 N0 = 20
            ok = True
            for step in range(20000):
                if rng.random() < f(logN):
                    ok = False
                    break
                logN += 0.3
            survive += ok
        print(f"   growth state with escape {law}: fraction never escaping in 20000 steps from N0 = 2^20: {survive/trials:.3f}")


# ----------------------------------------------------------------------------
# P4: the nine-iteration typology, re-checked
# ----------------------------------------------------------------------------

def juggler(n: int) -> int:
    return math.isqrt(n) if n % 2 == 0 else math.isqrt(n ** 3)


def part4(sig):
    print("\n== P4: typology numbers re-checked ==")
    # Juggler
    tot = keep = 0
    ll = []
    maxsteps = maxdig = 0
    for n in range(2, 2001):
        m = n
        steps = 0
        prev_par = None
        while m != 1:
            m2 = juggler(m)
            if m2 > 1 and m > 1:
                ll.append(math.log2(math.log(m2) / math.log(m)))
            par = m % 2
            if prev_par is not None:
                tot += 1
                keep += (par == prev_par)
            prev_par = par
            maxdig = max(maxdig, len(str(m2)))
            m = m2
            steps += 1
            assert steps < 10 ** 4
        maxsteps = max(maxsteps, steps)
    print(f"   Juggler n <= 2000: all reach 1; mean log2(log m'/log m) per step {sum(ll)/len(ll):+.3f} (fair coin {0.5*math.log2(1.5)+0.5*math.log2(0.5):+.3f}); "
          f"P(parity persists) {keep/tot:.3f}; longest {maxsteps} steps; largest {maxdig} digits")
    # reverse-and-add
    def ra(n):
        return n + int(str(n)[::-1])

    cand = []
    growth = []
    for n in range(1, 10 ** 5):
        m = n
        hit = False
        for k in range(300):
            s = str(m)
            if s == s[::-1] and k > 0:
                hit = True
                break
            m2 = ra(m)
            if k < 50:
                growth.append(math.log10(m2 / m))
            m = m2
        if not hit:
            cand.append(n)
    print(f"   reverse-and-add n < 10^5, 300 steps: not palindromic: {len(cand)} (first {cand[:12]}); below 10^4: {sum(1 for c in cand if c < 10**4)}; "
          f"mean digit growth per step (first 50 steps) {sum(growth)/len(growth):+.3f}")
    # look-and-say
    s = "1"
    lens = []
    for _ in range(60):
        out = []
        i = 0
        while i < len(s):
            j = i
            while j < len(s) and s[j] == s[i]:
                j += 1
            out.append(str(j - i) + s[i])
            i = j
        s = "".join(out)
        lens.append(len(s))
    print(f"   look-and-say: length ratio at 60 steps {lens[-1]/lens[-2]:.5f} (Conway 1.30358)")
    # Ducci (length 8): all vectors with entries < 30 reach zero
    from itertools import product
    worst = 0
    for vec in product(range(4), repeat=8):
        v = list(vec)
        k = 0
        while any(v):
            v = [abs(v[i] - v[(i + 1) % 8]) for i in range(8)]
            k += 1
            assert k < 100
        worst = max(worst, k)
    print(f"   Ducci length 8, entries < 4 (4^8 vectors): all reach 0, worst {worst} steps")
    # Kaprekar 4 digits
    worst = 0
    for n in range(0, 10000):
        d = f"{n:04d}"
        if len(set(d)) == 1:
            continue
        x = n
        k = 0
        while x != 6174:
            d = f"{x:04d}"
            x = int("".join(sorted(d, reverse=True))) - int("".join(sorted(d)))
            k += 1
            assert k < 20
        worst = max(worst, k)
    print(f"   Kaprekar 4 digits: every non-repdigit reaches 6174, worst {worst} steps")
    # aliquot class drifts and abundant density
    drift = defaultdict(list)
    for n in range(2, 2 * 10 ** 6 + 1, 2):
        a = v2(n)
        s = int(sig[n]) - n
        if s > 0 and a <= 6:
            drift[a].append(math.log2(s / n))
    print("   aliquot mean log2(s(n)/n) by class a (even n <= 2*10^6): " + ", ".join(f"a={a}: {sum(v)/len(v):+.3f}" for a, v in sorted(drift.items())))
    ab = sum(1 for n in range(1, 2 * 10 ** 6 + 1) if int(sig[n]) > 2 * n)
    print(f"   abundant density to 2*10^6: {ab/(2*10**6):.4f} (Deleglise: 0.2476..0.2480)")


# ----------------------------------------------------------------------------
# P5: Erdos persistence for abundant n (fixed horizon), with the sieve
# ----------------------------------------------------------------------------

def part5():
    print("\n== P5: Erdos persistence: P(k consecutive increases | abundant), abundant n <= 10^6, sieve to 2*10^7 ==")
    M = 2 * 10 ** 7
    sig = sigma_sieve(M)
    rng = random.Random(7)
    abund = [n for n in rng.sample(range(2, 10 ** 6), 60000) if int(sig[n]) > 2 * n]
    abund = abund[:5000]
    K = 8
    alive = [0] * (K + 1)
    inside = [0] * (K + 1)
    for n in abund:
        x = n
        ok = True
        for k in range(1, K + 1):
            if x > M:
                break
            s = int(sig[x]) - x
            inside[k] += 1
            if s > x:
                alive[k] += 1
                x = s
            else:
                break
    print("   k: " + ", ".join(f"{k}: {alive[k]/len(abund):.3f} (chains still inside the sieve: {inside[k]})" for k in range(1, K + 1)))
    print(f"   ({len(abund)} abundant starts; a chain leaving the sieve range counts as not increasing)")


if __name__ == "__main__":
    sig = part1()
    part2()
    part3()
    part4(sig)
    part5()
    print("\nDONE")
