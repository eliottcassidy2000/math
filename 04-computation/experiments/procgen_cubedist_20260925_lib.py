#!/usr/bin/env python3
"""
procgen_cubedist_20260925_lib.py -- shared exact helpers for the cube-distance lane
(session collatz-procgen-20260922, 2026-09-25; HYP-9138).

Setting (THM-4474): a level-k sign strategy sigma : odd residues mod 2^k -> {+1,-1};
    T_sigma(n) = n/2 (n even),  (m n + sigma(n mod 2^k))/2 (n odd),   m = 3 (Collatz sheet) or 5 (DRIFT control).
Flip set R = residues where sigma = -1 (sigma = + elsewhere); Collatz is R = {} (m = 3).
Mask encoding (as in the cube lane): bit i <-> residue 2i+1, bit 1 = minus.
Parity graph G_sigma on Z/2^k: s -> the two lifts of T(s) mod 2^(k-1); a cycle (a odd nodes, length p)
is expanding iff m^a > 2^p.  Class (i) <=> every cycle contracting (THM-4474 Theorem A).

Everything that is used as a proof step is exact (integers / Fractions); the C engine
(procgen_cubedist_20260925_engine.c) is used for speed and every certificate it returns is re-checked here.
"""
import os
import subprocess
from fractions import Fraction
from functools import lru_cache
from math import comb, gcd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
SCR = os.path.join(ROOT, 'scratch', 'procgen_cubedist')
ENGINE_SRC = os.path.join(HERE, 'procgen_cubedist_20260925_engine.c')
ENGINE = os.path.join(SCR, 'engine')


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)


def ensure_engine():
    os.makedirs(SCR, exist_ok=True)
    if (not os.path.exists(ENGINE)) or os.path.getmtime(ENGINE) < os.path.getmtime(ENGINE_SRC):
        subprocess.run(['gcc', '-O2', '-o', ENGINE, ENGINE_SRC], check=True)
    return ENGINE


# ----------------------------------------------------------------------------- basics
def mask_of(flips):
    m = 0
    for r in flips:
        m |= 1 << ((r - 1) // 2)
    return m


def flips_of(mask, k):
    return [2 * i + 1 for i in range(1 << (k - 1)) if (mask >> i) & 1]


def step(n, mul, sig):
    """one step of T_sigma on an integer n (sig = +1/-1 for odd n)"""
    return n // 2 if n % 2 == 0 else (mul * n + sig) // 2


def parity_word(r, k, mul=3):
    """all-plus (m n + 1) parity word of length k of the residue r mod 2^k (Terras: determined by r mod 2^k)"""
    w = []
    n = r
    for _ in range(k):
        w.append(n & 1)
        n = step(n, mul, 1)
    return w


def ballot(word, mul=3):
    """all prefixes j = 1..len have mul^(a_j) > 2^j"""
    a = 0
    for j, b in enumerate(word, 1):
        a += b
        if not mul ** a > 2 ** j:
            return False
    return True


def bad_set(k, mul=3):
    """Bad_k: residues r mod 2^k on which the all-plus map has no descent within k steps
    (every prefix multiplier mul^(a_j)/2^j > 1, j = 1..k)."""
    return [r for r in range(1, 1 << k, 2) if ballot(parity_word(r, k, mul), mul)]


def first_descent(r, k, mul=3):
    """least j <= k with mul^(a_j) < 2^j for the all-plus word of r, or None"""
    a = 0
    n = r
    for j in range(1, k + 1):
        a += n & 1
        n = step(n, mul, 1)
        if mul ** a < 2 ** j:
            return j
    return None


@lru_cache(maxsize=None)
def ballot_count(L, mul=3):
    """|Bad_L| = number of ballot words of length L (exact DP over a)."""
    # a_j must satisfy mul^a > 2^j; track counts by a
    cur = {0: 1}
    for j in range(1, L + 1):
        nxt = {}
        for a, c in cur.items():
            for b in (0, 1):
                na = a + b
                if mul ** na > 2 ** j:
                    nxt[na] = nxt.get(na, 0) + c
        cur = nxt
    return sum(cur.values())


def critical_ones(L, mul=3):
    """least a with mul^a > 2^L"""
    a = 0
    while not mul ** a > 2 ** L:
        a += 1
    return a


def necklaces_with_ones(n, j):
    """number of binary necklaces of length n with exactly j ones (Burnside)"""
    tot = 0
    # rotations by t: fixed words have period g = gcd(n, t); need (n/g) | j
    for t in range(n):
        g = gcd(n, t) if t else n
        if j % (n // g) == 0:
            tot += comb(g, j // (n // g))
    check(tot % n == 0, "Burnside count not divisible")
    return tot // n


def necklace_lower_bound(k, mul=3):
    """N_k = number of binary necklaces of length k whose density exceeds log_mul 2 (mul^a > 2^k)."""
    a0 = critical_ones(k, mul)
    return sum(necklaces_with_ones(k, j) for j in range(a0, k + 1))


def best_lower_approx(D, mul=3):
    """largest fraction a/p with p <= D and mul^a < 2^p (exact)."""
    best = Fraction(0, 1)
    a = 0
    for p in range(1, D + 1):
        while mul ** (a + 1) < 2 ** p:
            a += 1
        # a = max with mul^a < 2^p
        f = Fraction(a, p)
        if f > best:
            best = f
    return best


# ----------------------------------------------------------------------------- graph
def graph(k, flipset, mul=3):
    K = 1 << k
    H = K >> 1
    succ = []
    for s in range(K):
        if s % 2 == 0:
            t = s // 2
        else:
            t = (mul * s + (-1 if s in flipset else 1)) // 2
        t0 = t % H
        succ.append((t0, t0 + H))
    return succ


# ----------------------------------------------------------------------------- engine wrappers
def _run(header, k, mask):
    ensure_engine()
    txt = header + '\n' + (format(mask, 'x') if mask else '0') + '\n'
    out = subprocess.run([ENGINE], input=txt, capture_output=True, text=True, check=True).stdout
    return out


def karp_density(k, mask, mul=3):
    out = _run(f"karp {k} {mul}", k, mask)
    a, p = out.split()[1:3]
    return Fraction(int(a), int(p))


def certificate(k, mask, q, r, mul=3):
    """('OK', psi list) or ('CYCLE', cycle nodes)."""
    out = _run(f"cert {k} {mul} {q} {r}", k, mask).split('\n')
    if out[0].startswith('OK'):
        psi = [int(v) for v in out[1].split()]
        return 'OK', psi
    cyc = [int(v) for v in out[1].split(',')]
    return 'CYCLE', cyc


def disjoint_cycles(k, mask, q, r, M, mul=3):
    out = _run(f"cycles {k} {mul} {q} {r} {M}", k, mask)
    cycles = []
    for line in out.splitlines():
        if line.startswith('C '):
            cycles.append([int(v) for v in line[2:].split(',')])
    return cycles


def short_cycles(k, mask, q, r, P, cap, needflip=True, mul=3):
    """every simple cycle of G_sigma of length <= P with a/p > q/r (each once), at most cap; needflip: only
    cycles through a residue with sigma = -"""
    out = _run(f"short {k} {mul} {q} {r} {P} {cap} {1 if needflip else 0}", k, mask)
    return [[int(v) for v in line[2:].split(',')] for line in out.splitlines() if line.startswith('C ')]


def verify_certificate(k, flipset, psi, q, r, mul=3):
    """exact check: psi(t) <= psi(s) - w(s) on every edge, w(odd) = r - q, w(even) = -q.
    This proves every cycle has a/p <= q/r."""
    succ = graph(k, flipset, mul)
    for s in range(1 << k):
        w = (r - q) if s & 1 else -q
        for t in succ[s]:
            if not psi[t] <= psi[s] - w:
                return False
    return True


def verify_cycle(k, flipset, cyc, mul=3):
    """cyc is a closed walk of G_sigma; return (a, p, expanding?)"""
    succ = graph(k, flipset, mul)
    p = len(cyc)
    for i in range(p):
        check(cyc[(i + 1) % p] in succ[cyc[i]], "not a closed walk of G_sigma")
    a = sum(1 for v in cyc if v & 1)
    return a, p, mul ** a > 2 ** p


def expanding_threshold(k, mul=3):
    """q/r = best lower approximation of log_mul 2 with denominator <= 2^k (every simple cycle has p <= 2^k),
    so 'a/p > q/r' <=> 'expanding' for simple cycles."""
    return best_lower_approx(1 << k, mul)


def is_class_i(k, flipset, mul=3):
    """exact: max cycle density (Karp) < log_mul 2  (compare mul^a < 2^p)."""
    d = karp_density(k, mask_of(flipset), mul)
    return mul ** d.numerator < 2 ** d.denominator, d
