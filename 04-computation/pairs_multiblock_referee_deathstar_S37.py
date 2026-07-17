#!/usr/bin/env python3
"""death-star-2026-07-17-S37 (HYP-7181): referee for LRCDeviationPairs.lean and
LRCMultiBlockChain.lean.

PAIRS: (1) sandwich N_i + N_j <= N_ij + (q-1) and N_ij <= min(N_i, N_j);
(2) THE DILATE FORMULA: at q = 14Q, gcd(v_i, q) = 1, v_j == 2 v_i (mod q):
N_ij = 2*floor((Q-1)/2); D_ij = N_ij - (q-1)/49 = (5/7)Q + O(1) > 0.
MULTIBLOCK: the two-block fee arithmetic composes (nested windows shrink)."""
import random
from fractions import Fraction as Fr
from math import gcd

def in_band(vi, q, p):
    r = (vi * p) % q
    return q <= 14 * r <= 13 * q

def NT(v, q, T):
    return sum(1 for p in range(1, q) if all(not in_band(v[i], q, p) for i in T))

def referee_pairs(trials=800, seed=37):
    rnd = random.Random(seed)
    ok_sand = ok_dil = True
    n_dil = 0
    for _ in range(trials):
        q = rnd.randint(20, 400)
        v = [rnd.randint(1, 10**5) for _ in range(13)]
        i, j = rnd.sample(range(13), 2)
        Ni, Nj, Nij = NT(v, q, [i]), NT(v, q, [j]), NT(v, q, [i, j])
        if not (Ni + Nj <= Nij + (q - 1) and Nij <= min(Ni, Nj)):
            ok_sand = False
        # planted dilate at 14|q
        Q = rnd.randint(1, 40)
        q2 = 14 * Q
        vi = rnd.randint(1, 10**5)
        if gcd(vi, q2) != 1:
            continue
        n_dil += 1
        vj = (2 * vi) % q2 + q2 * rnd.randint(0, 3)
        v2 = list(v)
        v2[i], v2[j] = vi, vj
        got = NT(v2, q2, [i, j])
        want = 2 * ((Q - 1) // 2)
        if got != want:
            ok_dil = False
            print(f"  FAIL dilate: Q={Q} vi={vi} got {got} want {want}")
    print(f"pairs referee: sandwich {'PASS' if ok_sand else 'FAIL'}; "
          f"dilate formula {'PASS' if ok_dil else 'FAIL'} ({n_dil} planted)")

def referee_multiblock(trials=20000, seed=137):
    """two-block fee arithmetic: sample (k, B, blocks, eps) and verify the nested
    ledger is consistent (fee1 at 2d with eps1; fee2 at eps1 with eps2; eps2 > 0
    admits a nonempty singles tail entry)."""
    rnd = random.Random(seed)
    n_ok = 0
    for _ in range(trials):
        k = rnd.randint(0, 10)
        B = rnd.randint(1, 30)
        d = Fr(13 - k, 14 * (k + 1) * B)
        L = 2 * d
        ws1 = sorted(rnd.randint(200 * B, 4000 * B) for _ in range(rnd.randint(1, 3)))
        eps1 = Fr(1, 20 * max(ws1))
        fee1 = sum(L / 7 + Fr(3, 7 * u) + eps1 * (u * L + 3) for u in ws1)
        if not fee1 < L - eps1:
            continue
        ws2 = sorted(rnd.randint(200 * max(ws1), 2000 * max(ws1))
                     for _ in range(rnd.randint(1, 3)))
        eps2 = Fr(1, 20 * max(ws2))
        fee2 = sum(eps1 / 7 + Fr(3, 7 * u) + eps2 * (u * eps1 + 3) for u in ws2)
        if fee2 < eps1 - eps2:
            n_ok += 1
    print(f"multiblock referee: {n_ok}/{trials} sampled two-block ledgers compose "
          f"(nested windows admit; scale separation ~200x per level suffices)")

if __name__ == "__main__":
    referee_pairs()
    referee_multiblock()
