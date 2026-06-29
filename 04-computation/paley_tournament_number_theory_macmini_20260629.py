#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tournaments <-> primes: the number theory of PALEY tournament invariants.
(mac-mini-2026-06-29-S9)

Paley tournament T_p (p prime, p=3 mod4): vertices Z_p, i->j iff (j-i) is a nonzero
QR mod p.  This is THE bridge between tournaments and primes.  We compute, as
functions of the prime p:
  H(T_p)   = # directed Hamiltonian PATHS        (via subset DP)
  C(T_p)   = # directed Hamiltonian CYCLES       (via subset DP)
  f(T_p)   = # palindromic Hamiltonian paths     (THM-583 half-system DP; cheap)
and look for number-theoretic structure: parity (Redei: H odd), factorizations,
2-adic valuations, OEIS, and relations among H, C, f, and arithmetic of p.
"""
from __future__ import annotations
import functools
print = functools.partial(print, flush=True)
from sympy import factorint, isprime


def paley(p):
    qr = set((x * x) % p for x in range(1, p))
    return [[(i != j and ((j - i) % p) in qr) for j in range(p)] for i in range(p)]


def ham_paths_count(arc, n):
    """# directed Hamiltonian paths via DP over (visited bitmask, last vertex)."""
    # dp[mask][v] = # paths covering 'mask' ending at v
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        row = dp[mask]
        if not any(row):
            continue
        for v in range(n):
            c = row[v]
            if not c:
                continue
            for w in range(n):
                if not (mask >> w) & 1 and arc[v][w]:
                    dp[mask | (1 << w)][w] += c
    return sum(dp[full])


def ham_cycles_count(arc, n):
    """# directed Hamiltonian cycles (fix start=0 to avoid n-fold rotation overcount,
    count directed cycles through all vertices returning to 0)."""
    full = (1 << n) - 1
    dp = {}
    dp[(1, 0)] = 1
    # iterate
    table = [[0] * n for _ in range(1 << n)]
    table[1][0] = 1
    for mask in range(1 << n):
        if not (mask & 1):
            continue
        row = table[mask]
        for v in range(n):
            c = row[v]
            if not c:
                continue
            for w in range(n):
                if not (mask >> w) & 1 and arc[v][w]:
                    table[mask | (1 << w)][w] += c
    # close cycle: last -> 0
    return sum(table[full][v] for v in range(n) if arc[v][0])


def palindromic_count(arc, p):
    """THM-583: f = Ham-path count on the (p-1)/2 half-system, transfer DP."""
    m = (p - 1) // 2
    pairs = [(x, (-x) % p) for x in range(1, m + 1)]
    repidx = {}
    for idx, (a, b) in enumerate(pairs):
        repidx[a] = idx; repidx[b] = idx
    starts = [v for pr in pairs for v in pr]
    from functools import lru_cache
    table = {}
    def dp(last, used):
        key = (last, used)
        if key in table:
            return table[key]
        if used == (1 << m) - 1:
            r = 1 if arc[last][0] else 0
            table[key] = r; return r
        tot = 0
        for v in starts:
            idx = repidx[v]
            if used & (1 << idx):
                continue
            if arc[last][v]:
                tot += dp(v, used | (1 << idx))
        table[key] = tot; return tot
    return sum(dp(v, 1 << repidx[v]) for v in starts)


def main():
    print("=" * 78)
    print("Paley tournament invariants as functions of the prime p (mac-mini-S9)")
    print("=" * 78)
    primes3mod4 = [p for p in range(3, 32) if isprime(p) and p % 4 == 3]
    print(f"\nprimes p=3 mod4: {primes3mod4}")
    Hs, Cs, fs = [], [], []
    print(f"\n{'p':>3} {'H(T_p)':>12} {'factor H':>22} {'v2(H-1)':>7} {'f(palin)':>10} {'C(cyc)':>12}")
    for p in primes3mod4:
        arc = paley(p)
        f = palindromic_count(arc, p)
        H = ham_paths_count(arc, p) if p <= 21 else None
        C = ham_cycles_count(arc, p) if p <= 21 else None
        fs.append(f); Hs.append(H); Cs.append(C)
        if H is not None:
            v2 = 0; t = H - 1
            while t and t % 2 == 0:
                v2 += 1; t //= 2
            fac = factorint(H)
            facstr = "*".join(f"{q}^{e}" if e > 1 else f"{q}" for q, e in sorted(fac.items()))
            print(f"{p:>3} {H:>12} {facstr:>22} {v2:>7} {f:>10} {C if C else '-':>12}")
        else:
            print(f"{p:>3} {'(skip)':>12} {'':>22} {'':>7} {f:>10} {'(skip)':>12}")

    print(f"\nf-sequence (palindromic, all p=3mod4): {fs}")
    print(f"H-sequence: {[h for h in Hs if h is not None]}")
    print(f"C-sequence: {[c for c in Cs if c is not None]}")
    # number-theoretic probes
    print("\n[probes]")
    for i, p in enumerate(primes3mod4):
        if Hs[i] is None:
            continue
        H = Hs[i]
        print(f"  p={p}: H mod p = {H % p}, H mod (p-1) = {H % (p-1)}, "
              f"f mod p = {fs[i] % p}, f factor = {dict(factorint(fs[i]))}")
    print("\nWatch for: H odd (Redei), H/f patterns, p | (H - something), Gauss-sum (sqrt p)")
    print("=" * 78)


if __name__ == "__main__":
    main()
