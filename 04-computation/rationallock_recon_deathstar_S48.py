#!/usr/bin/env python3
"""Recon for HYP-7240 (death-star-S48): the rational-ratio lock + branches.

Speeds a = g*i', b = g*j' (reduced ratio i':j').  Both fail at p/q with
witnesses w_a, w_b (14|v*p - w*q| < q).
(a) LOCK: i'+j' <= 13  =>  j'*w_a = i'*w_b exactly.
(b) BRANCHES: i'+j' <= 27  =>  |j'*w_a - i'*w_b| <= 1.
(c) LOCKED COUNT: gcd(i',j')=1, i'+j' <= 13, gcd(g,q)=1 =>
    N_ab = 2*floor((q-1)/(14*max(i',j'))).
(d) BRIDGE: N_ab/(q-1) -> boxeph LEM-044 mu for consecutive pairs.
"""
from math import gcd
import random
random.seed(48)

def failing_witness(v, q, p):
    r = (v * p) % q
    if 14 * r < q: return (v * p) // q
    if 14 * r > 13 * q: return (v * p) // q + 1
    return None

# (a)+(b) over random reduced pairs
lock_ok = lock_bad = br_ok = br_bad = 0
for _ in range(400000):
    ip = random.randint(1, 13); jp = random.randint(1, 13)
    if gcd(ip, jp) != 1 or ip == jp: continue
    g = random.randint(1, 60)
    q = random.randint(2, 900); p = random.randint(1, q - 1)
    wa = failing_witness(g * ip, q, p); wb = failing_witness(g * jp, q, p)
    if wa is None or wb is None: continue
    k = jp * wa - ip * wb
    if ip + jp <= 13:
        lock_ok += (k == 0); lock_bad += (k != 0)
        if k != 0 and lock_bad <= 3: print("LOCK FAIL", ip, jp, g, q, p, wa, wb)
    if ip + jp <= 27:
        br_ok += (abs(k) <= 1); br_bad += (abs(k) > 1)
        if abs(k) > 1 and br_bad <= 3: print("BRANCH FAIL", ip, jp, g, q, p, k)
print(f"(a) rational lock (i'+j'<=13): {lock_ok} ok, {lock_bad} fail",
      "PASS" if lock_bad == 0 else "FAIL")
print(f"(b) branch bound (i'+j'<=27): {br_ok} ok, {br_bad} fail",
      "PASS" if br_bad == 0 else "FAIL")

# (c) the locked count over all locked pairs from (1..13)
cnt_ok = cnt_bad = 0
pairs_locked = [(a, b) for a in range(1, 14) for b in range(a + 1, 14)
                if (a // gcd(a, b)) + (b // gcd(a, b)) <= 13]
pairs_sparse = [(a, b) for a in range(1, 14) for b in range(a + 1, 14)
                if (a // gcd(a, b)) + (b // gcd(a, b)) > 13]
for (a, b) in pairs_locked:
    g = gcd(a, b); ip, jp = a // g, b // g
    for q in random.sample(range(200, 1400), 12):
        if gcd(g, q) != 1: continue
        N = sum(1 for p in range(1, q)
                if failing_witness(a, q, p) is not None
                and failing_witness(b, q, p) is not None)
        pred = 2 * ((q - 1) // (14 * max(ip, jp)))
        cnt_ok += (N == pred); cnt_bad += (N != pred)
        if N != pred and cnt_bad <= 3: print("COUNT FAIL", a, b, q, N, pred)
print(f"(c) locked count on {len(pairs_locked)} locked pairs of (1..13): "
      f"{cnt_ok} ok, {cnt_bad} fail", "PASS" if cnt_bad == 0 else "FAIL")
print(f"    locked pairs: {len(pairs_locked)}/78; sparse: {len(pairs_sparse)} {pairs_sparse}")

# (d) bridge to boxeph LEM-044 on consecutive pairs k=2..6 (locked) and k=7,8 (sparse)
from fractions import Fraction as F
q = 98000 + 1  # large, coprime to small numbers likely
print("(d) bridge: consecutive pairs (k,k+1), N/(q-1) vs LEM-044 mu:")
for k in range(2, 9):
    a, b = k, k + 1
    N = sum(1 for p in range(1, q)
            if failing_witness(a, q, p) is not None
            and failing_witness(b, q, p) is not None)
    r = k % 7
    mu = F(1, 49) + F(r * (6 - r), 49 * k * (k + 1))
    print(f"    k={k}: N/(q-1) = {N}/{q-1} = {N/(q-1):.6f}   mu = {float(mu):.6f}"
          f"   {'locked' if 2*k+1 <= 13 else 'SPARSE'}")
