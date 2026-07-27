#!/usr/bin/env python3
"""
ham_excess_spectra_kps_S134.py
(kind-pasteur-2026-07-26-S134; companion to HYP-9028)

Hamiltonian excess spectra of the two circulant tournament families.

H(T) vs the random-tournament mean n!/2^(n-1):
  excess(T) = H(T) / (n!/2^(n-1)).
Families: rotational R_n = circulant {1..(n-1)/2} (odd n),
Paley QR_p (p = 3 mod 4).  Also tc = H/|Aut| (|Aut(R_n)| >= n from
rotations; Paley |Aut| = p(p-1)/2), Korselt/Carmichael check on tc,
and the free-rotation lemma n | H(R_n), p | H(QR_p) verified
(the rotation group acts freely on directed Hamiltonian paths).

Data (this run): excess climbs 2.0 -> 2.522 (R_n, n<=17) and
2.4 -> 2.527 (QR_p, p<=19), consistent with excess -> e with O(1/n)
correction.  tc(QR_11) = 1729 is Carmichael (Chernick k=1) but
tc(QR_19) is NOT: the Carmichael hit is a coincidence, recorded to
prevent numerology.  OEIS lookups rate-limited at run time: novelty
of H(R_n) = 3, 15, 175, 3267, 93027, 3711175, 198464295, 13689269499
UNVERIFIED.

Reproduction: python 04-computation/ham_excess_spectra_kps_S134.py
"""
import numpy as np
from math import factorial

def ham_bitdp(n, adj):
    # adj[i] = bitmask of out-neighbors
    full = (1 << n) - 1
    dp = np.zeros((1 << n, n), dtype=np.int64)
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, 1 << n):
        row = dp[mask]
        nz = np.nonzero(row)[0]
        for last in nz:
            c = row[last]
            nxts = adj[last] & ~mask
            while nxts:
                v = (nxts & -nxts).bit_length() - 1
                nxts &= nxts - 1
                dp[mask | (1 << v)][v] += c
    return int(dp[full].sum())

def circulant(n, conn):
    return [sum(1 << ((i + d) % n) for d in conn) for i in range(n)]

def korselt(m):
    if m < 3: return False
    x = m; ps = []
    d = 2
    while d * d <= x:
        while x % d == 0:
            ps.append(d); x //= d
        d += 1
    if x > 1: ps.append(x)
    if len(set(ps)) != len(ps): return False  # not squarefree
    if len(ps) < 3: return False
    return all((m - 1) % (p - 1) == 0 for p in set(ps))

print("family n  H          excess      tc         korselt(tc)")
for n in range(3, 20, 2):
    conn = list(range(1, (n - 1) // 2 + 1))
    adj = circulant(n, conn)
    if n > 17:
        break
    H = ham_bitdp(n, adj)
    exc = H / (factorial(n) / 2 ** (n - 1))
    assert H % n == 0, "rotation lemma fails"
    # |Aut| for rotational: at least n; report tc_n = H/n (rotations only)
    print(f"R_{n:2d}     {H:<10d} {exc:<11.5f} H/n={H//n}")

for p in (7, 11, 19):
    qr = sorted({(i * i) % p for i in range(1, p)})
    adj = circulant(p, qr)
    H = ham_bitdp(p, adj)
    exc = H / (factorial(p) / 2 ** (p - 1))
    assert H % p == 0
    aut = p * (p - 1) // 2
    tc = H // aut if H % aut == 0 else None
    print(f"QR_{p:2d}    {H:<12d} {exc:<11.5f} tc={tc} korselt={korselt(tc) if tc else '-'}")
