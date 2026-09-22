#!/usr/bin/env python3
"""collatz_mod6_20260917_row_braid_typing_audit.py

Adversarial verify-and-fix audit (2026-09-21) of lane row_braid_typing, session collatz-mod6-20260917.
Independent code: nothing is imported from the lane script.  Every check raises (survives python -O).
Exact integer/Fraction arithmetic throughout; RAM < 100 MB; runtime well under a minute (plus an optional
PARI/GP call, a fraction of a second).

Sections
  A1  orders of 2, power-of-two target classes, Wieferich single classes (with exact residues)
  A2  the three image rows over Z, 21 target classes mod 63
  A3  R(n)=4n+1: R^3 by exact affine composition, the exact universal modulus 42 by gcd, mod-7 orbits, mod-18 period
  A4  odd-multiplier tower, the signed diagonal solved as a unit equation, boundary p=2 / p=1,
      inverse braids R_p with the FULL orbit decomposition (not just the orbit of 1), Wieferich braid at p=1093
  A5  cycle census with the lane's caps, independently coded
  A6  rows j<=60 and j in [-40,-1], inherited JSON (path relative to this file), first appearances, exact count law
  A7  row law over Z (both signs), AP layers, index recurrence, base table
  A8  the three minus cycles and their (row,j,h) placements
  A9  the -7/4 question: closed form, coprimality lemma, reduction of BOTH finite searches to ONE cubic Thue
      equation a^3+2a^2b+ab^2+b^3 = +-1 (discriminant -23, plastic-number field), unit exponents, PARI/GP `thue`
      (certified) if gp is on PATH, convergents of rho^2
  A10 rational PCF quadratics, hostiles x^2-6 / x^2-12, THM-4146 cycle, x^2-2 step on rows, Bang exception n=6

Reproduction:
  python3 04-computation/experiments/collatz_mod6_20260917_row_braid_typing_audit.py \
      > 05-knowledge/results/collatz_mod6_20260917_row_braid_typing_audit.out
"""
import os
import json
import shutil
import subprocess
import hashlib
from fractions import Fraction
from math import gcd

CHECKS = 0


def check(cond, msg):
    global CHECKS
    CHECKS += 1
    if not cond:
        raise RuntimeError("AUDIT CHECK FAILED: " + msg)


def order(a, m):
    if gcd(a, m) != 1:
        raise RuntimeError("non-unit")
    t, x = 1, a % m
    while x != 1:
        x = x * a % m
        t += 1
    return t


def val(n, p):
    if n == 0:
        raise RuntimeError("val(0)")
    n = abs(n)
    k = 0
    while n % p == 0:
        n //= p
        k += 1
    return k


def ch(m):
    """(core, height) with sign on the core."""
    if m == 0:
        raise RuntimeError("ch(0)")
    h = val(m, 2)
    return m // 2 ** h, h


def F(n):
    if n % 2 == 0:
        raise RuntimeError("F on even")
    return (3 * n + 1) // 2


def hdr(s):
    print()
    print("-" * 78)
    print(s)
    print("-" * 78)


B = {1: 2, 3: 5, 5: 8}
RHO = {2: 1, 5: 3, 8: 5}

# ---------------------------------------------------------------------------
hdr("A1  orders of 2, power-of-two target classes, Wieferich single classes")
o = {m: order(2, m) for m in (3, 7, 9, 21, 63)}
print("ord_m(2):", o)
check(o == {3: 2, 7: 3, 9: 6, 21: 6, 63: 6}, "orders")
check(min(k for k in range(1, 100) if (2 ** k - 1) % 9 == 0) == 6, "first Mersenne divisible by 9 is 2^6-1")
check(sorted(pow(2, k, 9) for k in range(6)) == [1, 2, 4, 5, 7, 8], "2 primitive root mod 9")
check(sorted(set(pow(2, k, 63) for k in range(6))) == [1, 2, 4, 8, 16, 32], "<2> mod 63")
# CRT: 2 -> (2 mod 9, 2 mod 7); <2> mod 63 -> (Z/9)^* is a bijection, -> {1,2,4} mod 7 is 2:1
img9 = [pow(2, k, 9) for k in range(6)]
img7 = [pow(2, k, 7) for k in range(6)]
check(len(set(img9)) == 6 and sorted(set(img7)) == [1, 2, 4] and all(img7.count(v) == 2 for v in (1, 2, 4)), "CRT projections")


def target_classes(p):
    """image classes mod p^2 of odd n under (pn+1)/2: exactly m = 2^{-1} mod p (p classes mod p^2)."""
    return sorted(set(((p * n + 1) // 2) % (p * p) for n in range(1, 2 * p * p, 2)))


def pow2_classes(p):
    o1, o2 = order(2, p), order(2, p * p)
    inv2 = pow(2, -1, p)
    ks = [k for k in range(0, o2) if pow(2, k, p) == inv2]
    return o1, o2, sorted(set(pow(2, k, p * p) for k in ks))


print("p | ord_p | ord_p^2 | #target classes hit by powers of two | p | first class")
for p in (3, 5, 7, 11, 13, 17, 19, 1093, 3511):
    o1, o2, hit = pow2_classes(p)
    tc = target_classes(p) if p < 100 else None
    if tc is not None:
        check(len(tc) == p and all(x % p == pow(2, -1, p) for x in tc), "target classes = 2^{-1} mod p, p of them")
        check(set(hit) <= set(tc), "powers of two land in target classes")
    check(len(hit) == o2 // o1, "count = ord_{p^2}/ord_p")
    check(o2 in (o1, p * o1), "ord_{p^2} is ord_p or p ord_p")
    wief = pow(2, p - 1, p * p) == 1
    check((len(hit) == p) == (not wief), "count p iff non-Wieferich")
    print("%5d | %5d | %8d | %3d | %s" % (p, o1, o2, len(hit), hit[:7] if p < 100 else hit))
    if p == 3:
        check(hit == [2, 5, 8], "p=3 classes 2,5,8")
    if p == 1093:
        check(hit == [597325] and 1093 ** 2 == 1194649 and pow(2, 363, 1093 ** 2) == 597325, "1093 single class 597325 = 2^363 mod 1093^2")
        check(val(2 ** 364 - 1, 1093) == 2, "v_1093(2^364-1)=2")
    if p == 3511:
        check(hit == [6163561] and 3511 ** 2 == 12327121 and pow(2, 1754, 3511 ** 2) == 6163561, "3511 single class 6163561 = 2^1754 mod 3511^2")
        check(val(2 ** 1755 - 1, 3511) == 2, "v_3511(2^1755-1)=2")
# exponent table for p=3
tab = {}
for k in range(1, 121):
    if k % 2 == 0:
        check(pow(2, k, 3) == 1, "even k not an image")
        continue
    src = (2 ** (k + 1) - 1) // 3
    check(src % 2 == 1 and F(src) == 2 ** k, "source of 2^k")
    tab.setdefault(k % 6, set()).add((pow(2, k, 9), src % 6, pow(2, k, 7)))
check(tab == {1: {(2, 1, 2)}, 3: {(8, 5, 1)}, 5: {(5, 3, 4)}}, "k mod 6 -> (b, r, 2^k mod 7)")
print("CONFIRMED: k=1,3,5 mod 6 -> image row b=2,8,5 (source row 1,5,3), 2^k mod 7 = 2,1,4; on odd k, k mod 3 <-> k mod 6.")

# ---------------------------------------------------------------------------
hdr("A2  the three image rows over Z, 21 target classes mod 63")
for n in range(-200001, 200002, 2):
    j, r = divmod(n, 6)
    check(F(n) == 9 * j + B[r] and F(n) % 3 == 2, "F(6j+r)=9j+b_r over Z")
tc63 = sorted(set(F(n) % 63 for n in range(1, 127, 2)))
check(len(tc63) == 21 and tc63 == [x for x in range(63) if x % 3 == 2], "21 target classes mod 63")
p2 = sorted(set(pow(2, k, 63) for k in range(60)) & set(tc63))
check(p2 == [2, 8, 32] and [x % 9 for x in p2] == [2, 8, 5] and [x % 7 for x in p2] == [2, 1, 4], "powers of two in target classes mod 63")
print("CONFIRMED (odd |n|<=200001): rows 9j+2, 9j+5, 9j+8; target classes mod 63 with powers of two: [2, 8, 32].")
# boundary n=1, 9
check(ch(F(1)) == (1, 1) and ch(F(9)) == (7, 1) and 9 % 6 == 3 and ch(F(3)) == (5, 0), "n=1 -> (1,1) in row 1; n=9 -> (7,1) in row 3; n=3 -> (5,0)")
for bad in (2, 4):
    try:
        F(bad)
        raise RuntimeError("F accepted even n")
    except RuntimeError as e:
        check("even" in str(e), "F rejects even n=%d" % bad)

# ---------------------------------------------------------------------------
hdr("A3  R(n)=4n+1: exact affine composition, universal modulus, mod-7 orbits, mod-18 period")


def compose(f, g):  # affine maps as (m, k): x -> m x + k ; f o g
    return (f[0] * g[0], f[0] * g[1] + f[1])


Rm = (4, 1)
R3 = compose(Rm, compose(Rm, Rm))
check(R3 == (64, 21), "R^3 = 64n+21 exactly (integer coefficients)")
check(64 - 1 == 63 and 21 * 3 == 63 and 21 == 3 * 7, "R^3(n)-n = 63n+21 = 21(3n+1)")
# universal modulus: gcd over odd n of 21(3n+1)
g = 0
for n in range(-999, 1000, 2):
    g = gcd(g, 21 * (3 * n + 1))
check(g == 42, "gcd_n 21(3n+1) over odd n = 42 (exact universal modulus)")
check(gcd(21 * 4, 21 * 10) == 42, "already forced by n=1, n=3")
Rt = (1, 0)
for t in range(1, 13):
    Rt = compose(Rm, Rt)
    check(Rt == (4 ** t, (4 ** t - 1) // 3), "R^t = 4^t n + (4^t-1)/3")
    check(val(4 ** t - 1, 3) == 1 + val(t, 3), "v_3(4^t-1) = 1 + v_3(t)")
# periods on residues
def perm_order(m, k, M, odd_only=True):
    xs = [x for x in range(M) if (x % 2 == 1 or not odd_only)]
    best = 0
    lens = set()
    for x0 in xs:
        x, t = x0, 0
        while True:
            x = (m * x + k) % M
            t += 1
            if x == x0:
                break
        lens.add(t)
    return lens


check(perm_order(4, 1, 6) == {3}, "R period 3 on rows mod 6")
check(perm_order(4, 1, 18) == {9}, "R period 9 on odd classes mod 18")
check(perm_order(4, 1, 42) == {3}, "R period 3 on odd classes mod 42")
check(perm_order(4, 1, 7, odd_only=False) == {1, 3}, "R mod 7: fixed point and 3-cycles")
check([x for x in range(7) if (4 * x + 1) % 7 == x] == [2] and (3 * 2 + 1) % 7 == 0, "fixed class 2 mod 7 = 7 | 3n+1")
check(order(4, 9) == 3 and order(4, 7) == 3 and order(4, 63) == 3, "4 has order 3 mod 9, 7, 63")
check(order(4, 9) == order(2, 9) // gcd(2, order(2, 9)) == 3, "ord_9(4) = ord_9(2)/gcd(2,ord_9(2)) = 3")
for n in range(-2001, 2002, 2):
    x = 64 * n + 21
    check((x - n) % 126 != 0 and ((x - n) % 84 == 0) == (n % 4 == 1), "never 0 mod 126; 0 mod 84 iff n=1 mod 4")
    check(F(x) == 64 * F(n), "F(R^3 n) = 64 F(n)")
print("CONFIRMED: R^3=(64,21); universal source modulus exactly 42; R mod 7 orbits {2},{0,1,5},{3,4,6}; period 9 mod 18.")

# ---------------------------------------------------------------------------
hdr("A4  tower F_p=(pn+1)/2, signed diagonal, boundary p=1,2, inverse braids with full orbit decomposition")


def Fp(p, n):
    if (p * n + 1) % 2:
        raise RuntimeError("parity")
    return (p * n + 1) // 2


for p in range(-9, 40, 2):
    for n in range(-501, 502, 2):
        check(Fp(p + 2, n) == Fp(p, n) + n and Fp(p, n) == Fp(1, n) + (p - 1) // 2 * n, "tower identities (signed p, n)")
        check(Fp(p, n) - n == Fp(p - 2, n), "companion is F_{p-2}")
# p=2 is excluded by parity: 2n+1 is odd for every n
check(all((2 * n + 1) % 2 == 1 for n in range(-50, 50)), "p=2: pn+1 odd, F_2 undefined on every n (boundary)")
# diagonal: (p-4) n = -1 over odd integers <=> p-4 = +-1, n = -+1
diag = sorted((p, n) for p in range(-99, 100, 2) for n in range(-999, 1000, 2) if (p - 4) * n == -1)
check(diag == [(3, 1), (5, -1)], "signed diagonal = {(3,1),(5,-1)}")
check(Fp(3, -1) == -1 and Fp(1, 1) == 1, "F_3(-1)=-1, F_1(1)=1")
check(all(Fp(-1, n) <= 0 for n in range(1, 100, 2)), "p=1 companion (1-n)/2 <= 0 for n>=1")
print("CONFIRMED: diagonal (p-4)n=-1 over odd p,n has exactly (3,1),(5,-1); p=2 excluded by parity.")


def least_q(p):
    d = 1
    while (2 ** d - 1) % p:
        d += 1
    return d


print("p | d | c | r | orbit lengths of R_p on odd classes mod 2p, 2p^2, 2p^3 (full decomposition)")
for p in (3, 5, 7, 11, 13):
    d = least_q(p)
    q = 2 ** d
    c = (q - 1) // p
    r = val(q - 1, p)
    check(r == 1, "r=1")
    lens = [perm_order(q, c, 2 * p ** s) for s in (1, 2, 3)]
    check(lens == [{p}, {p * p}, {p ** 3}], "single cycle of length p^s on odd classes mod 2p^s (all classes)")
    check(all(Fp(p, q * n + c) == q * Fp(p, n) for n in range(-999, 1000, 2)), "F_p(R_p n) = q F_p(n) over Z")
    print("%2d | %2d | %3d | %d | %s" % (p, d, c, r, [sorted(l) for l in lens]))
    check(d == order(2, p), "d = ord_p(2)")
# p=1
check(perm_order(2, 1, 6) == {1, 2} and (2 * 5 + 1) % 6 == 5, "R_1 mod 6: orbits (1,3),(5)")
for u in range(1, 100, 2):
    for h in range(0, 7):
        n = 2 ** (h + 1) * u - 1
        check(ch(Fp(1, n)) == (u, h), "F_1 fibre n=2^{h+1}u-1 for every h")
# Wieferich braid p=1093: r=2, so period p^{max(0,s-1)}: 1 mod 2p, p mod 2p^2
p = 1093
d = order(2, p)
q = 2 ** d
c = (q - 1) // p
check(val(q - 1, p) == 2, "r=2 at 1093")
check(all((q * x + c - x) % (2 * p) == 0 for x in range(1, 2 * p, 2)), "R_1093 fixes every odd class mod 2p")
M = 2 * p * p
qm, cm = q % M, c % M
for x0 in (1, 3, 5, 1093 + 2, 2 * 1093 + 1, 12345):
    x, t = x0, 0
    while True:
        x = (qm * x + cm) % M
        t += 1
        if x == x0:
            break
    check(t == p, "R_1093 has period exactly p on odd classes mod 2p^2 (sample)")
print("CONFIRMED: r=1 braids are single cycles on ALL odd classes mod 2p^s; p=1093 (r=2) fixes odd classes mod 2p, period 1093 mod 2p^2.")

# ---------------------------------------------------------------------------
hdr("A5  cycle census of T_p (odd starts <= 20000, <= 3000 steps, cap 10^40), independently coded")


def Tp(p, n):
    m = p * n + 1
    return m >> val(m, 2)


census = {}
for p in (1, 3, 5, 7):
    cyc, unresolved = set(), 0
    for n0 in range(1, 20001, 2):
        path = {}
        x = n0
        t = 0
        found = False
        while t < 3000 and x < 10 ** 40:
            if x in path:
                # cycle from x
                y, c_ = x, []
                while True:
                    c_.append(y)
                    y = Tp(p, y)
                    if y == x:
                        break
                i = c_.index(min(c_))
                cyc.add(tuple(c_[i:] + c_[:i]))
                found = True
                break
            path[x] = t
            x = Tp(p, x)
            t += 1
        if not found:
            unresolved += 1
    census[p] = (sorted(cyc), unresolved)
    print("  p=%d: %s unresolved=%d" % (p, sorted(cyc), unresolved))
check(census[1] == ([(1,)], 0) and census[3] == ([(1,)], 0), "p=1,3")
check(census[5] == ([(1, 3), (13, 33, 83), (17, 43, 27)], 9605), "p=5 census incl. 9605 unresolved")
check(census[7] == ([(1,)], 9982), "p=7 census incl. 9982 unresolved")
check(all(Tp(1, n) < n for n in range(3, 100001, 2)) and Tp(1, 1) == 1, "T_1 descent")

# ---------------------------------------------------------------------------
hdr("A6  rows j<=60 / j in [-40,-1], inherited JSON, first appearances, exact count law")
rows = {r: [ch(F(6 * j + r)) for j in range(61)] for r in (1, 3, 5)}
nrows = {r: {j: ch(F(6 * j + r)) for j in range(-1, -41, -1)} for r in (1, 3, 5)}
here = os.path.dirname(os.path.abspath(__file__))
jpath = os.path.join(here, "..", "..", "05-knowledge", "results", "arithmetic_braids_20260917_collatz.json")
with open(jpath) as fh:
    inh = json.load(fh)
for r in (1, 3, 5):
    ir = [tuple(x) for x in inh["rows"][str(r)]]
    check(len(ir) == 35 and rows[r][:35] == ir, "inherited JSON agrees j<=34")
check(rows[1][:5] == [(1, 1), (11, 0), (5, 2), (29, 0), (19, 1)], "row 1 head")
check(rows[3][:5] == [(5, 0), (7, 1), (23, 0), (1, 5), (41, 0)], "row 3 head")
check(rows[5][:5] == [(1, 3), (17, 0), (13, 1), (35, 0), (11, 2)], "row 5 head")
check(rows[1][14] == (1, 7) and rows[5][56] == (1, 9) and rows[3][35] == (5, 6) and rows[1][46] == (13, 5), "spot values")
first = {}
for r in (1, 3, 5):
    for j, (u, h) in enumerate(rows[r]):
        first.setdefault(u, {}).setdefault(r, (j, h))
want_first = {1: {1: (0, 1), 3: (3, 5), 5: (0, 3)}, 5: {1: (2, 2), 3: (0, 0), 5: (8, 4)}, 7: {1: (6, 3), 3: (1, 1), 5: (24, 5)},
              11: {1: (1, 0), 3: (19, 4), 5: (4, 2)}, 13: {1: (46, 5), 3: (11, 3), 5: (2, 1)}, 17: {1: (30, 4), 3: (7, 2), 5: (1, 0)},
              19: {1: (4, 1), 5: (16, 3)}, 23: {1: (10, 2), 3: (2, 0), 5: (40, 4)}, 25: {1: (22, 3), 3: (5, 1)}}
for u, w in want_first.items():
    check(first[u] == w, "first appearance of %d" % u)
# u=19 in row 3 and u=25 in row 5 are beyond j=60: exact positions
check(ch(F(6 * 67 + 3)) == (19, 5) and ch(F(6 * 88 + 5)) == (25, 5), "19 enters row 3 at j=67 (h=5), 25 enters row 5 at j=88 (h=5)")
cnt = {}
for r in (1, 3, 5):
    for (u, h) in rows[r]:
        cnt[u] = cnt.get(u, 0) + 1
check(cnt[1] == 5 and cnt[5] == 4, "core 1 five times, core 5 four times (j<=60)")
three = sorted(u for u, k in cnt.items() if k == 3)
check(three == [7, 11, 13, 17, 23, 29], "cores appearing exactly three times for j<=60")
check(sorted(u for u, k in cnt.items() if k >= 3) == [1, 5, 7, 11, 13, 17, 23, 29], "all cores with >=3 appearances")
print("  cores with exactly 3 appearances (j<=60):", three)
# exact count law: #appearances of u in j<=J = #{admissible h : 2^h u <= 9J+8}
J = 60
for u in range(1, 400, 2):
    if u % 3 == 0:
        continue
    pred = sum(1 for h in range(0, 40) if (h % 2 == 1) == (u % 3 == 1) and 2 ** h * u <= 9 * J + 8)
    check(pred == cnt.get(u, 0), "exact count law for u=%d" % u)
pos1 = sorted((r, j, h) for r in (1, 3, 5) for j, (u, h) in enumerate(rows[r]) if u == 1)
check(pos1 == [(1, 0, 1), (1, 14, 7), (3, 3, 5), (5, 0, 3), (5, 56, 9)], "positions of core 1")
# negative rows
check([nrows[5][j] for j in (-1, -2, -3, -4)] == [(-1, 0), (-5, 1), (-19, 0), (-7, 2)], "neg row 5 head")
check([nrows[1][j] for j in (-1, -2, -3, -4)] == [(-7, 0), (-1, 4), (-25, 0), (-17, 1)], "neg row 1 head")
check([nrows[3][j] for j in (-1, -2, -3, -4)] == [(-1, 2), (-13, 0), (-11, 1), (-31, 0)], "neg row 3 head")
check([6 * j + 1 for j in (-1, -2, -3)] == [-5, -11, -17] and [6 * j + 5 for j in (-1, -2, -3)] == [-1, -7, -13], "negative sources")
print("CONFIRMED: rows, JSON agreement, first appearances, count law, core-1 positions, negative row heads.")

# ---------------------------------------------------------------------------
hdr("A7  row law over Z, AP layers, index recurrence, base table")
n_pairs = 0
for u in list(range(1, 3002, 2)) + list(range(-1, -3002, -2)):
    if u % 3 == 0:
        continue
    for h in range(0, 13):
        num = 2 ** (h + 1) * u - 1
        adm = (num % 3 == 0)
        check(adm == ((h % 2 == 1) == (u % 3 == 1)), "parity gate")
        if not adm:
            continue
        n = num // 3
        check(n % 2 == 1 and ch(F(n)) == (u, h), "fibre round trip")
        b = (2 ** h * u) % 9
        check(b in RHO and RHO[b] == n % 6, "row law")
        check(n == 6 * ((2 ** h * u - b) // 9) + n % 6, "index law")
        if h + 2 <= 12:
            m2 = (2 ** (h + 3) * u - 1)
            check(m2 % 3 == 0 and (m2 // 3) % 6 == (n % 6 + 4) % 6, "row(u,h+2) = row(u,h)+4 mod 6")
        check(F(4 * n + 1) == 4 * F(n), "R braid over Z")
        n_pairs += 1
print("  (u,h) pairs checked (both signs):", n_pairs)
# |u|<=3001 odd, 3 !| u: positive 501 (u=1 mod 3, 6 odd h) + 500 (u=2 mod 3, 7 even h) = 6506; negative 501*7+500*6 = 6507
check(n_pairs == 6506 + 6507 == 13013, "13013 admissible (u,h) pairs for |u|<=3001, h<=12 (the lane's 6500 is u<3001, positive only)")
base = {}
for um in (1, 5, 7, 11, 13, 17):
    h0 = 1 if um % 3 == 1 else 0
    base[um] = (h0, ((2 ** (h0 + 1) * um - 1) // 3) % 6)
check(base == {1: (1, 1), 7: (1, 3), 13: (1, 5), 5: (0, 3), 11: (0, 1), 17: (0, 5)}, "base table")
# AP layers: at height h in row r cores are the class b_r 2^{-h} mod 9 (odd -> one class mod 18), index step 2^{h+1}
for r in (1, 3, 5):
    for h in range(0, 7):
        cls = B[r] * pow(2, -h, 9) % 9
        js = [j for j in range(2 ** 12) if ch(9 * j + B[r])[1] == h]
        us = [ch(9 * j + B[r])[0] for j in js]
        check(all(u % 9 == cls and u % 2 == 1 for u in us), "class")
        check(all(js[i + 1] - js[i] == 2 ** (h + 1) for i in range(len(js) - 1)), "index step")
        check(all(us[i + 1] - us[i] == 18 for i in range(len(us) - 1)), "core step 18")
layer = {(1, 0): 11, (1, 1): 1, (1, 2): 5, (1, 3): 7, (1, 4): 17, (3, 0): 5, (5, 0): 17}
for (r, h), u18 in layer.items():
    check(B[r] * pow(2, -h, 9) % 9 == u18 % 9, "layer (%d,%d) -> %d mod 18" % (r, h, u18))
for r in (1, 3, 5):
    for u in range(1, 400, 2):
        if u % 3 == 0:
            continue
        h = 0
        while not ((2 ** (h + 1) * u - 1) % 3 == 0 and ((2 ** (h + 1) * u - 1) // 3) % 6 == r):
            h += 1
        j0 = (2 ** h * u - B[r]) // 9
        j1 = (2 ** (h + 6) * u - B[r]) // 9
        check(j1 == 64 * j0 + 7 * B[r], "index recurrence")
check([7 * B[r] for r in (1, 3, 5)] == [14, 35, 56], "shifts 14,35,56")
print("CONFIRMED: row law and R-braid over Z for |u|<=3001, h<=12; AP layers; recurrence j'=64j+7b_r.")

# ---------------------------------------------------------------------------
hdr("A8  the three minus cycles and their placements")


def T(n):
    return Tp(3, n)


cycles = [(-1,), (-5, -7), (-17, -25, -37, -55, -41, -61, -91)]
for cyc in cycles:
    for i, x in enumerate(cyc):
        check(T(x) == cyc[(i + 1) % len(cyc)], "cycle closes in order")
check([len(c) for c in cycles] == [1, 2, 7], "lengths")
check([sorted(set(x % 6 for x in c)) for c in cycles] == [[5], [1, 5], [1, 5]], "member rows")
check(all(x % 3 != 0 for c in cycles for x in c), "no member is 0 mod 3")
pl = {}
for r in (1, 3, 5):
    for j, (u, h) in nrows[r].items():
        pl.setdefault(u, []).append((r, j, h))
want = {-1: [(1, -2, 4), (3, -29, 8), (3, -1, 2), (5, -8, 6), (5, -1, 0)], -5: [(1, -18, 5), (3, -5, 3), (5, -2, 1)],
        -7: [(1, -1, 0), (3, -13, 4), (5, -4, 2)], -17: [(1, -4, 1), (5, -16, 3)], -25: [(1, -3, 0), (5, -12, 2)],
        -37: [(3, -17, 2), (5, -5, 0)], -55: [(3, -25, 2), (5, -7, 0)], -41: [(3, -37, 3), (5, -10, 1)],
        -61: [(1, -7, 0), (5, -28, 2)], -91: [(5, -11, 0)]}
for u, w in want.items():
    check(sorted(pl.get(u, [])) == sorted(w), "placement of %d" % u)
for j in range(-3000, 0):
    for r in (1, 3, 5):
        n = 6 * j + r
        check((-n) % 6 == {1: 5, 3: 3, 5: 1}[r] and F(n) == -((3 * (-n) - 1) // 2), "negation law")
print("CONFIRMED: cycle orders, rows {5},{1,5},{1,5}, all (row,j,h) placements for |j|<=40, negation law j>=-3000.")

# ---------------------------------------------------------------------------
hdr("A9  the -7/4 question: closed form, coprimality lemma, ONE Thue equation, units of the plastic field, PARI")


def orb(c, N):
    x, out = Fraction(0), []
    for _ in range(N):
        x = x * x + c
        out.append(x)
    return out


def no_new_prime(num, earlier):
    g = abs(num)
    if g <= 1:
        return True
    N = 1
    for e in earlier:
        N *= abs(e) if e else 1
    while True:
        d = gcd(g, N)
        if d == 1:
            return g == 1
        g //= d


c = Fraction(-7, 4)
check(c * (c + 1) ** 2 == Fraction(-63, 64) and orb(c, 3) == [c, Fraction(21, 16), Fraction(-7, 256)], "-7/4 orbit")
check(orb(c, 4)[3] == Fraction(-114639, 65536) and 114639 == 3 * 7 * 53 * 103, "term 4")
check(orb(Fraction(-29, 16), 3)[2] == Fraction(23345, 65536) and 23345 == 5 * 7 * 23 * 29, "-29/16 term 3")
check(orb(Fraction(-29, 16), 2)[1] == Fraction(377, 256) and 377 == 13 * 29, "-29/16 term 2")


def N3(a, b):  # a^3+2a^2b+ab^2+b^3 = a(a+b)^2 + b^3
    return a * (a + b) ** 2 + b ** 3


# closed form and coprimality lemma: gcd(a(a+b)^2+b^3, a(a+b)) = 1 whenever gcd(a,b)=1
for b in range(1, 121):
    for a in range(-300, 301):
        if gcd(abs(a), b) != 1:
            continue
        cc = Fraction(a, b)
        f3 = orb(cc, 3)[2]
        check(f3 == Fraction(a * N3(a, b), b ** 4), "closed form")
        check(f3.numerator == a * N3(a, b) and f3.denominator == b ** 4, "lowest terms")
        check(gcd(abs(N3(a, b)), abs(a * (a + b))) == 1, "coprimality lemma")
# symbolic proof of the lemma via resultants (sympy)
try:
    import sympy as sp
    A, Bb = sp.symbols('a b')
    N3s = A * (A + Bb) ** 2 + Bb ** 3
    check(sp.expand(N3s - (A ** 3 + 2 * A ** 2 * Bb + A * Bb ** 2 + Bb ** 3)) == 0, "N3 expanded")
    r1, r2 = sp.expand(sp.resultant(N3s, A, A)), sp.expand(sp.resultant(N3s, A + Bb, A))
    print("  resultants (sympy sign convention): Res_a(N3,a) =", r1, " Res_a(N3,a+b) =", r2)
    check(r1 in (Bb ** 3, -Bb ** 3), "Res_a(N3, a) = +-b^3")
    check(r2 in (Bb ** 3, -Bb ** 3), "Res_a(N3, a+b) = +-b^3")
    check(sp.discriminant(sp.Poly(A ** 3 + 2 * A ** 2 + A + 1, A)) == -23, "discriminant -23")
    # rho^2 (rho^3=rho+1) is a root of y^3-2y^2+y-1, i.e. -rho^2 is the real root of x^3+2x^2+x+1
    rr = sp.symbols('r')
    m2 = sp.Poly(sp.rem((rr ** 2) ** 3 - 2 * (rr ** 2) ** 2 + rr ** 2 - 1, rr ** 3 - rr - 1, rr), rr)
    check(m2.is_zero, "rho^2 satisfies y^3-2y^2+y-1 modulo rho^3=rho+1")
    print("  sympy: both resultants are +-b^3 (so gcd(a,b)=1 forces gcd(N3, a(a+b))=1); disc = -23; -rho^2 is the real root.")
except ImportError:
    print("  sympy unavailable: resultant identities not re-derived here (numeric lemma check above stands).")
print("PROVED: for c=a/b in lowest terms, f^3(0)=a*N3(a,b)/b^4 in lowest terms with gcd(N3(a,b), a(a+b))=1, so")
print("        'no new numerator prime at term 3' <=> N3(a,b) = a^3+2a^2b+ab^2+b^3 = +-1 (a cubic Thue equation).")
print("        For b=2^k, k>=1: c(c+1)^2+1 = N3/b^3 with N3 odd, so the 'strong form' +-2^{-m} also forces N3 = +-1.")
# units: rho^n = x + y rho + z rho^2 ; N3(a,b) = Norm(a + b rho^2); solutions <=> y_n = 0 with (a,b) = +-(x_n, z_n)
def mul_rho(v):
    x, y, z = v
    return (z, x + z, y)


def mul_rho_inv(v):  # rho^{-1} = rho^2 - 1
    x, y, z = v
    return (y - x, z, x)


# norm identity N3(a,b) = N_{K/Q}(a + b rho^2): determinant of multiplication by a + b rho^2 in the basis {1,rho,rho^2}
def norm_a_b_rho2(a, b):
    cols = []
    for e in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
        w = mul_rho(mul_rho(e))
        cols.append(tuple(a * e[i] + b * w[i] for i in range(3)))
    (m11, m21, m31), (m12, m22, m32), (m13, m23, m33) = cols
    return (m11 * (m22 * m33 - m23 * m32) - m12 * (m21 * m33 - m23 * m31) + m13 * (m21 * m32 - m22 * m31))


for a in range(-30, 31):
    for b in range(-30, 31):
        check(norm_a_b_rho2(a, b) == N3(a, b), "N3(a,b) = Norm(a + b rho^2)")
print("  norm identity N3(a,b) = N_{K/Q}(a + b rho^2) checked on |a|,|b| <= 30 (3x3 determinant in the basis 1, rho, rho^2).")
sols_unit = []
v = (1, 0, 0)
for n_ in range(0, 3001):
    if v[1] == 0:
        sols_unit.append((n_, v[0], v[2]))
    v = mul_rho(v)
v = (1, 0, 0)
for n_ in range(0, 3001):
    if n_ > 0 and v[1] == 0:
        sols_unit.append((-n_, v[0], v[2]))
    v = mul_rho_inv(v)
sols_unit.sort()
print("  unit exponents n with rho^n = x + 0*rho + z*rho^2, |n|<=3000: ", sols_unit)
check(sols_unit == [(-14, -7, 4), (-5, 2, -1), (-1, -1, 1), (0, 1, 0), (2, 0, 1)], "y_n = 0 exactly at n=-14,-5,-1,0,2")
for n_, a, b in sols_unit:
    check(N3(a, b) == 1, "norm of rho^n is 1")
check(N3(-7, 4) == 1 and N3(-2, 1) == -1 and N3(-1, 1) == 1 and N3(0, 1) == 1 and N3(1, 0) == 1, "known solutions")
# large direct window near the real root -rho^2 (all solutions with b>=2 must have a/b within 1 of the real root)
big = []
for b in range(1, 200001):
    a0 = -(7 * b) // 4  # crude centre; N3(a,b)=b^3 g(a/b) with g increasing near its real root, scan a window
    for a in range(a0 - 3, a0 + 4):
        if gcd(abs(a), b) == 1 and abs(N3(a, b)) == 1:
            big.append((a, b))
    if b == 1:
        for a in range(-10, 11):
            if abs(N3(a, b)) == 1 and (a, b) not in big:
                big.append((a, b))
print("  window search b<=200000, a within 3 of -7b/4, |N3(a,b)|=1, gcd(a,b)=1:", sorted(big))
check(sorted(big) == [(-7, 4), (-2, 1), (-1, 1), (0, 1)], "window search b<=200000 near the real root")
# convergents of rho^2 = 1.7548776662...: 7/4 is the third convergent
cf = []
num, den = [], []
# exact CF of the real root of y^3-2y^2+y-1 via sign changes (unique real root)
poly = [1, -2, 1, -1]  # y^3-2y^2+y-1 coefficients high->low


def peval(p_, t):
    r_ = Fraction(0)
    for co in p_:
        r_ = r_ * t + co
    return r_


def cf_step(p_):
    # floor of the unique real root (assumed > 0): find k with p(k) p(k+1) <= 0 ... unique sign change
    k = 0
    while peval(p_, k) * peval(p_, k + 1) > 0:
        k += 1
    # new polynomial for t = 1/(x-k): q(t) = t^3 p(k + 1/t)
    # p(k+1/t) = sum c_i (k+1/t)^i ; multiply by t^3
    import itertools
    deg = len(p_) - 1
    q = [0] * (deg + 1)
    for i, co in enumerate(p_):
        e = deg - i  # power of (k + 1/t)
        # (k + 1/t)^e t^3 = sum_m C(e,m) k^{e-m} t^{3-m}
        from math import comb
        for m in range(e + 1):
            q[m] += co * comb(e, m) * k ** (e - m)  # coefficient of t^{3-m}
    return k, q


pp = poly
for _ in range(8):
    k, pp = cf_step(pp)
    cf.append(k)
conv = []
h0, h1, k0, k1 = 0, 1, 1, 0
for a_ in cf:
    h0, h1 = h1, a_ * h1 + h0
    k0, k1 = k1, a_ * k1 + k0
    conv.append(Fraction(h1, k1))
print("  continued fraction of rho^2:", cf, " convergents:", [str(x) for x in conv])
check(cf[:4] == [1, 1, 3, 12] and conv[2] == Fraction(7, 4), "rho^2 = [1;1,3,12,...], third convergent 7/4")
# real root of y^3-2y^2+y-1 (= rho^2) by exact bisection, 12 decimals
lo, hi = Fraction(1), Fraction(2)
for _ in range(60):
    mid = (lo + hi) / 2
    if peval([1, -2, 1, -1], mid) < 0:
        lo = mid
    else:
        hi = mid
print("  rho^2 = %.12f (exact bisection; -rho^2 is the real root of x^3+2x^2+x+1; 7/4 = 1.75)" % float(lo))
check(Fraction(1754877, 1000000) < lo < Fraction(1754878, 1000000), "rho^2 = 1.754877...")
# PARI/GP certified Thue solve
gp = shutil.which("gp")
if gp:
    # one statement per line: gp swallows a one-line "...; quit" without flushing its output
    code = "\n".join(["default(parisize,64000000);", "P=x^3+2*x^2+x+1;", "T=thueinit(P,1);", "K=bnfinit(P,1);",
                       "print(thue(T,1));", "print(thue(T,-1));", "print(bnfcertify(K));", "print(K.no);",
                       "print(nfdisc(P));", "print(version());", "quit"]) + "\n"
    out = subprocess.run([gp, "-q", "-f"], input=code, capture_output=True, text=True, timeout=300).stdout.strip().splitlines()
    out = [l for l in out if not l.startswith("  ***")]
    print("  PARI/GP thue(x^3+2x^2y+xy^2+y^3 = +1):", out[0])
    print("  PARI/GP thue(... = -1):               ", out[1])
    print("  bnfcertify, class number, nfdisc:     ", out[2:5])
    check(out[0] == "[[-7, 4], [-1, 1], [0, 1], [1, 0], [2, -1]]", "thue +1 solutions")
    check(out[1] == "[[-2, 1], [-1, 0], [0, -1], [1, -1], [7, -4]]", "thue -1 solutions")
    check(out[2:5] == ["1", "1", "-23"], "certified, h=1, disc -23")
    print("  PARI/GP version:", out[5])
    print("CITED (PARI/GP %s, thueinit(P,1) certified, bnfcertify=1): the complete solution set of N3(a,b)=+-1 is" % out[5])
    print("        +-{(1,0),(0,1),(-1,1),(-2,1),(-7,4)}; hence c in {0,-1,-2,-7/4} and the only NONDEGENERATE c are -2 and -7/4.")
else:
    print("  gp not on PATH: PARI completeness check SKIPPED here (recorded run: solutions +-{(1,0),(0,1),(-1,1),(-2,1),(-7,4)}).")
# term 4: f^4(0) = a (a N3^2 + b^7) / b^8 and a+b ALWAYS divides a N3^2 + b^7 (N3 = b^3 mod a+b), so the term-4
# condition has a different coprimality structure from term 3 (left OPEN in the note).
for b in range(1, 60):
    for a in range(-150, 151):
        if gcd(abs(a), b) != 1:
            continue
        f4 = orb(Fraction(a, b), 4)[3]
        check(f4 == Fraction(a * (a * N3(a, b) ** 2 + b ** 7), b ** 8), "closed form f^4(0)")
        check((a * N3(a, b) ** 2 + b ** 7) % (a + b) == 0 if a + b != 0 else True, "(a+b) | a N3^2 + b^7")
print("  term 4: f^4(0) = a(a N3^2 + b^7)/b^8 and (a+b) | a N3^2 + b^7 always (b<60, |a|<=150): a different structure, left OPEN.")
# strong form for integer c
check(sorted(cint for cint in range(-50, 51) if abs(cint * (cint + 1) ** 2 + 1) == 1) == [-2, -1, 0], "integer strong form")

# ---------------------------------------------------------------------------
hdr("A10  PCF quadratics, hostiles, THM-4146 cycle, x^2-2 step, Bang exception")
pcf = []
for cint in range(-200, 201):
    x, seen = 0, set()
    while x not in seen and abs(x) < 10 ** 30:
        seen.add(x)
        x = x * x + cint
    if x in seen:
        pcf.append(cint)
check(pcf == [-2, -1, 0], "integer PCF")
for cint in (-6, -12):
    fixed = [x for x in range(-20, 21) if x * x + cint == x]
    check(len(fixed) == 2, "x^2%+d has integer fixed points %s" % (cint, fixed))
    x, seen = 0, set()
    while x not in seen and abs(x) < 10 ** 30:
        seen.add(x)
        x = x * x + cint
    check(x not in seen, "0 escapes under x^2%+d" % cint)
for b in range(2, 40):
    for a in range(-60, 61):
        if gcd(abs(a), b) != 1:
            continue
        o_ = orb(Fraction(a, b), 5)
        check([t.denominator for t in o_] == [b ** (2 ** i) for i in range(5)], "denominators b^(2^i)")
G = lambda y: Fraction(y * y - 29, 4)
check(G(-7) == 5 and G(5) == -1 and G(-1) == -7, "THM-4146 (33)")
ff = lambda x: x * x - Fraction(29, 16)
check(ff(Fraction(-7, 4)) == Fraction(5, 4) and ff(Fraction(5, 4)) == Fraction(-1, 4) and ff(Fraction(-1, 4)) == Fraction(-7, 4), "x^2-29/16 cycle")
step = {}
for b in (2, 5, 8):
    t = 2 * b % 9
    t = t if t in (2, 5, 8) else 9 - t
    step[RHO[b]] = RHO[t]
check(step == {1: 3, 3: 5, 5: 1} and {r: step[step[r]] for r in step} == {1: 5, 3: 1, 5: 3}, "x^2-2 step and its square = R on rows")
earlier, nonew = [], []
for n_ in range(2, 201):
    t = 2 ** n_ - 1
    earlier.append(2 ** (n_ - 1) - 1) if n_ == 2 else None
    if no_new_prime(t, [2 ** m - 1 for m in range(1, n_)]):
        nonew.append(n_)
check(nonew == [6], "2^n-1 without primitive prime divisor, 2<=n<=200: only n=6")
print("CONFIRMED: PCF {-2,-1,0}; x^2-6, x^2-12 hostiles; THM-4146 cycle; x^2-2 step {1:3,3:5,5:1}; Bang exception n=6 (n<=200).")

print()
print("audit checks passed:", CHECKS)
print("audit script sha256:", hashlib.sha256(open(__file__, "rb").read()).hexdigest())
print("ALL AUDIT CHECKS PASSED")
