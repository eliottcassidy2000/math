#!/usr/bin/env python3
"""collatz_mod6_20260917_row_braid_typing.py

Lane row_braid_typing of session collatz-mod6-20260917.
Exact typing of every "three" and "two" in the user's Collatz prompt.

Inheritance (read first, cited by path, NOT re-derived here):
  05-knowledge/results/arithmetic_braids_20260917_collatz.md   (rows, R(n)=4n+1 braid, (core,height) rows)
  05-knowledge/results/arithmetic_braids_20260917_summand.md   (doubling forest, R_p tower section 5)
  05-knowledge/results/arithmetic_braids_20260917_divisors.md  (F=S+U shapes, sandwich matrix)
  01-canon/theorems/THM-4139-*.md, THM-4146-*.md              (x^2-29/16 three-cycle, 3:4:5, mod-63 census)
  sibling lanes of this session:
  05-knowledge/results/collatz_mod6_20260917_pythagorean_semicircle.out
  05-knowledge/results/collatz_mod6_20260917_extended_collatz_scc.out

Every load-bearing check raises on failure (active under python -O).
All arithmetic is exact (int / Fraction).  RAM < 100 MB, runtime a few seconds.
"""
import sys
import json
import hashlib
from fractions import Fraction
from math import gcd

# --------------------------------------------------------------------------
# utilities
# --------------------------------------------------------------------------
CHECKS = 0


def check(cond, msg):
    global CHECKS
    CHECKS += 1
    if not cond:
        raise RuntimeError("CHECK FAILED: " + msg)


def v2(n):
    n = abs(n)
    if n == 0:
        raise RuntimeError("v2(0)")
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return k


def vp(n, p):
    n = abs(n)
    if n == 0:
        raise RuntimeError("vp(0)")
    k = 0
    while n % p == 0:
        n //= p
        k += 1
    return k


def mult_order(a, m):
    if gcd(a, m) != 1:
        raise RuntimeError("order of non-unit")
    k, x = 1, a % m
    while x != 1:
        x = (x * a) % m
        k += 1
    return k


def F(n):
    """single-halving shortcut map on odd n (any sign)."""
    if n % 2 == 0:
        raise RuntimeError("F on even")
    return (3 * n + 1) // 2


def core_height(m):
    """m = 2^h * u with u odd, sign kept on u."""
    if m == 0:
        raise RuntimeError("core of 0")
    h = v2(m)
    return m // (2 ** h), h


def banner(title):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


B_OF_ROW = {1: 2, 3: 5, 5: 8}          # F(6j+r) = 9j + b_r
ROW_OF_B = {2: 1, 5: 3, 8: 5}          # rho
RHO = ROW_OF_B

# --------------------------------------------------------------------------
# S1  ord_9(2)=6, the Mersenne number 63, and what 9, 3, 7 each see
# --------------------------------------------------------------------------
banner("S1  ord_9(2)=6 <=> 9 | 2^6-1 = 63 : rows, exponent classes, R-period, factor 7")

pow2_mod9 = [pow(2, k, 9) for k in range(0, 12)]
print("2^k mod 9, k=0..11 :", pow2_mod9)
check(mult_order(2, 9) == 6, "ord_9(2)=6")
check(63 % 9 == 0 and all((2 ** k - 1) % 9 != 0 for k in range(1, 6)), "63 first Mersenne divisible by 9")
print("PROVED: ord_9(2)=6, i.e. 2^6-1=63 is the first Mersenne number divisible by 9;")
print("        2 is a primitive root mod 9: <2> = {1,2,4,8,7,5} = (Z/9)^*.")

pow2_mod63 = [pow(2, k, 63) for k in range(0, 7)]
pow2_mod7 = [pow(2, k, 7) for k in range(0, 7)]
print("2^k mod 63, k=0..6 :", pow2_mod63, " ord_63(2) =", mult_order(2, 63))
print("2^k mod  7, k=0..6 :", pow2_mod7, "  ord_7(2)  =", mult_order(2, 7))
check(mult_order(2, 63) == 6 and mult_order(2, 7) == 3 and mult_order(2, 3) == 2, "orders 63/7/3")
print("PROVED: ord_63(2)=lcm(ord_9(2),ord_7(2))=lcm(6,3)=6 (already in THM-4139 (38));")
print("        <2> mod 63 = {1,2,4,8,16,32} projects ISOMORPHICALLY onto (Z/9)^* and 2:1 onto <2> mod 7 = {1,2,4}.")

# (a) the three image rows: pure algebra of the multiplier 3
N_ROWS = 100001
for n in range(1, N_ROWS, 2):
    j, r = divmod(n, 6)
    check(F(n) == 9 * j + B_OF_ROW[r], "F(6j+r)=9j+b_r at n=%d" % n)
    check(F(n) % 3 == 2, "F(n) = 2 mod 3")
print("PROVED (+checked odd n<%d): F(6j+1)=9j+2, F(6j+3)=9j+5, F(6j+5)=9j+8; every image is 2 mod 3." % N_ROWS)
print("   Mechanism: 3n+1 mod 18 is a function of n mod 6; halving by 2^{-1}=5 mod 9 sends {4,10,16} mod 18 to {2,5,8} mod 9.")
print("   This uses only the multiplier 3 and 2^{-1} = 2 mod 3; it does NOT use ord_9(2)=6.")

# (b) exponent classes of powers of two in the image rows
print()
print("powers of two in the image rows (k = exponent, 2^k = F(n) needs 2^k = 2 mod 3, i.e. k odd):")
print("   k mod 6 | 2^k mod 9 = image row b | source row r=(2^{k+1}-1)/3 mod 6 | 2^k mod 7")
seen = {}
for k in range(1, 61):
    if k % 2 == 0:
        check(pow(2, k, 3) == 1, "even k gives 1 mod 3, not an image")
        continue
    src = (2 ** (k + 1) - 1) // 3
    check(3 * src + 1 == 2 ** (k + 1) and src % 2 == 1, "source integral")
    r = src % 6
    b = pow(2, k, 9)
    check(ROW_OF_B[b] == r, "row of 2^k")
    seen.setdefault(k % 6, set()).add((b, r, pow(2, k, 7)))
for km in (1, 3, 5):
    (b, r, m7), = seen[km]
    print("   %d       | %d                     | %d                                | %d" % (km, b, r, m7))
check(seen == {1: {(2, 1, 2)}, 3: {(8, 5, 1)}, 5: {(5, 3, 4)}}, "exponent classes 1,5,3 <-> rows 1,3,5")
print("PROVED: 2^k lies in image row b iff k is odd and k mod 6 = 1 (b=2, row 1), 5 (b=5, row 3), 3 (b=8, row 5).")
print("   Mechanism: k odd is the gate 2^k = 2 mod 3 (ord_3(2)=2); which row is k mod 6 (ord_9(2)=6).")
print("   The factor 7 sees k mod 3 (ord_7(2)=3): on odd k, k mod 6 -> k mod 3 is a bijection, so 7 sees the")
print("   row braid too, but NOT the parity gate.  3 sees only the gate.  9 = 3^2 sees both.")

# (c) R(n)=4n+1: R^3(n) = 64n + 21, 21 = 63/3, and period three on rows
try:
    import sympy
    nn = sympy.symbols('n')
    R3 = 4 * (4 * (4 * nn + 1) + 1) + 1
    check(sympy.expand(R3 - (64 * nn + 21)) == 0, "R^3 symbolic")
    check(sympy.expand(64 * nn + 21 - (63 * nn + nn + 21)) == 0, "63n + (n+21)")
    check(sympy.factor(63 * nn + 21) == 21 * (3 * nn + 1), "R^3(n)-n = 21(3n+1)")
    print("PROVED (sympy symbolic): R^3(n) = 64n+21 = 63n + (n+21);  R^3(n)-n = 63n+21 = 21(3n+1) = 3*7*(3n+1); 21 = 63/3.")
except ImportError:
    print("sympy unavailable; symbolic check skipped (numeric check below).")


def R(n):
    return 4 * n + 1


for n in range(-2000, 2001, 1):
    if n % 2 == 0:
        continue
    check(R(R(R(n))) == 64 * n + 21, "R^3 numeric")
    check(F(R(n)) == 4 * F(n), "F(R n) = 4 F(n)")
    for t in range(1, 8):
        x = n
        for _ in range(t):
            x = R(x)
        check((x - n) * 3 == (4 ** t - 1) * (3 * n + 1), "R^t formula")
        check(vp(x - n, 3) == vp(t, 3), "v_3(R^t n - n) = v_3(t)")
    # rows mod 6: exact period 3
    check(R(n) % 6 == (n + 4) % 6, "R adds 4 mod 6 (1->5->3->1)")
    check((R(R(R(n))) - n) % 6 == 0 and (R(n) - n) % 6 != 0 and (R(R(n)) - n) % 6 != 0, "period exactly 3 on rows")
    # R^3 fixes n mod 42 but shifts n mod 9 by 3
    check((R(R(R(n))) - n) % 42 == 0, "R^3 fixes n mod 42")
    check((R(R(R(n))) - n) % 9 == 3, "R^3 shifts n by 3 mod 9")
    # target side mod 63: R multiplies F(n) by 4, order exactly 3
    fn = F(n) % 63
    check((16 * fn) % 63 != fn or fn % 9 == 0, "order")
    check((64 * fn) % 63 == fn, "R^3 fixes target mod 63")
    check((4 * fn) % 63 != fn, "R moves target mod 63 (target is 2 mod 3, never 0 mod 9)")
print("PROVED (+checked |n|<=2000 odd): R adds 4 mod 6, so rows cycle 1->5->3->1 with exact period 3;")
print("   9 | 4^t-1 iff 3 | t (v_3(4^t-1)=1+v_3(t)), so the row period 3 IS ord_9(4)=3, i.e. ord_9(2)=6.")
print("   R^3 fixes the SOURCE modulo 42 = 6*7 and the TARGET F(n) modulo 63 = 9*7 (F(R^3 n) = 64 F(n)).")
print("   R itself multiplies the target by 4, of order 3 in both (Z/9)^* and (Z/7)^*.")

# what 7 sees on the source side: R mod 7 = affine x -> 4x+1
orb = {}
for a in range(7):
    x, cyc = a, [a]
    for _ in range(10):
        x = (4 * x + 1) % 7
        if x == a:
            break
        cyc.append(x)
    orb[a] = tuple(cyc)
cycles7 = sorted(set(tuple(sorted(c)) for c in orb.values()))
print("R mod 7 (x -> 4x+1) orbits:", cycles7)
check(cycles7 == [(0, 1, 5), (2,), (3, 4, 6)], "R mod 7: one fixed point 2, two 3-cycles")
check((3 * 2 + 1) % 7 == 0, "fixed class 2 mod 7 is exactly 7 | 3n+1")
print("PROVED: R mod 7 has the single fixed class n = 2 mod 7 (exactly the sources with 7 | 3n+1, i.e. 7 | F(n))")
print("        and two 3-cycles {0,1,5}, {3,4,6}; so R has period 3 mod 6, mod 7, hence mod 42, but period 9 mod 18.")

# powers of two among the 21 target classes mod 63
targets63 = sorted(set(F(n) % 63 for n in range(1, 2 * 63, 2)))
check(len(targets63) == 21 and all(t % 3 == 2 for t in targets63), "21 target classes mod 63 = those 2 mod 3")
p2_63 = sorted(set(pow(2, k, 63) for k in range(1, 7)) & set(targets63))
print("target classes mod 63 containing powers of two:", p2_63, " (2^1, 2^3, 2^5) -> mod 9:", [x % 9 for x in p2_63], " mod 7:", [x % 7 for x in p2_63])
check(p2_63 == [2, 8, 32], "powers of two in target classes mod 63")
print("PROVED: exactly 3 of the 21 target classes mod 63 contain powers of two, one per row mod 9 (2,8,5) and")
print("        one per class {2,1,4} mod 7; the rows mod 9 are the mod-63 target classes modulo the factor 7")
print("        in the exact sense (Z/63)^* = (Z/9)^* x (Z/7)^*, 2 -> (2,2), and <2> -> (Z/9)^* is an isomorphism.")

# Hostile control: general p, and Wieferich primes
print()
print("general odd p: rows of pn+1 are the p odd classes mod 2p; powers of two in the image rows occupy")
print("   ord_{p^2}(2)/ord_p(2) distinct rows (all of them iff p is non-Wieferich):")
print("   p    | ord_p(2) | ord_{p^2}(2) | #rows hit by powers of two | rows")
for p in (3, 5, 7, 11, 13, 1093, 3511):
    o1, o2 = mult_order(2, p), mult_order(2, p * p)
    inv2 = pow(2, -1, p)
    ks = [k for k in range(1, o2 + 1) if pow(2, k, p) == inv2]
    rows_hit = sorted(set(pow(2, k, p * p) for k in ks))
    cnt = len(rows_hit)
    check(cnt == o2 // o1, "count = ord_{p^2}/ord_p")
    print("   %-4d | %-8d | %-12d | %-26d | %s" % (p, o1, o2, cnt, str(rows_hit) if p < 20 else "single class %d mod %d" % (rows_hit[0], p * p)))
    if p in (1093, 3511):
        check(cnt == 1, "Wieferich: all powers of two in one row")
    else:
        check(cnt == p, "non-Wieferich: one per row")
print("PROVED: the power-of-two exponent classes are the cosets of ord_p(2) inside ord_{p^2}(2); the count is p")
print("        iff 2^{ord_p(2)} != 1 mod p^2 (non-Wieferich).  REFUTED as a universal law: p=1093, 3511 put every")
print("        power of two in the image into ONE row.  For p=3 the count 3 is 6/2: '3 rows see 3 exponent classes'.")

# --------------------------------------------------------------------------
# S2  Odd-multiplier tower F_p(n) = (pn+1)/2
# --------------------------------------------------------------------------
banner("S2  Odd-multiplier tower F_p = (pn+1)/2, summand typing, inverse braids R_p")


def Fp(p, n):
    if (p * n + 1) % 2:
        raise RuntimeError("Fp parity")
    return (p * n + 1) // 2


for p in (1, 3, 5, 7, 9, 11):
    for n in range(-999, 1000, 2):
        check(Fp(p + 2, n) == n + Fp(p, n), "F_{p+2} = n + F_p")
        check(Fp(p, n) == Fp(1, n) + ((p - 1) // 2) * n, "F_p = F_1 + (p-1)/2 n")
print("PROVED: F_{p+2}(n) = F_p(n) + n for all odd p, n; hence F_p(n) = F_1(n) + ((p-1)/2) n with F_1(n)=(n+1)/2.")
print("   Summand typing (inherited summand.md section 1,3): the arrow n -> F_p(n) in the strict summand shadow has")
print("   the unique companion F_p(n) - n = F_{p-2}(n).  p=3: companion F_1(n)=(n+1)/2 (the inherited witness (a,b,z));")
print("   p=5: the companion of n in 5n+1 is its own 3n+1 image F_3(n);  p=1: companion F_{-1}(n) = (1-n)/2 <= 0,")
print("   so the n+1 map (n+1)/2 is NOT a strict positive summand arrow (its companion is nonpositive) for n>=1.")

# diagonal (equal summand) exception
diag = [(p, n) for p in range(1, 40, 2) for n in range(1, 400, 2) if Fp(p - 2, n) == n]
check(diag == [(3, 1)], "unique diagonal exception (p,n)=(3,1)")
print("PROVED: F_{p-2}(n) = n  <=>  (p-4) n = -1  <=>  (p,n) = (3,1): the only equal-summand arrow in the whole tower")
print("        is the trivial Collatz cycle 1 -> 2 (checked p<40, n<400: %s)." % diag)

# inverse fibre braids R_p, inherited summand.md section 5


def least_q(p):
    d = 1
    while (2 ** d - 1) % p:
        d += 1
    return 2 ** d, d


print()
print("inverse-fibre braids R_p(n) = q n + c, q = least power of two = 1 mod p, c=(q-1)/p (inherited summand.md sec.5):")
print("   p | q=2^d | c  | r=v_p(q-1) | F_p(R_p n)=q F_p(n) | period of R_p on odd classes mod 2p, 2p^2, 2p^3 | height shift d")
for p in (1, 3, 5, 7, 11, 13):
    q, d = least_q(p)
    c = (q - 1) // p
    ok = all(Fp(p, q * n + c) == q * Fp(p, n) for n in range(1, 2001, 2))
    check(ok, "F_p R_p = q F_p")
    pers = []
    for s in (1, 2, 3):
        M = 2 * p ** s
        # period of n -> qn+c on odd residues mod M (single orbit test + period)
        x, t = 1, 0
        while True:
            x = (q * x + c) % M
            t += 1
            if x == 1:
                break
        pers.append(t)
    rr = vp(q - 1, p) if p > 1 else 0
    print("   %-2d| %-5d | %-2d | %-10s | %-19s | %-46s | %d" % (p, q, c, str(rr) if p > 1 else "n/a", ok, pers, d))
    if p > 1:
        check(pers == [p, p * p, p ** 3] if rr == 1 else True, "period p^s when r=1")
    else:
        check(pers == [1, 1, 1], "p=1 trivial modulus 2")
print("PROVED (inherited (10)-(11)): for p=3,5,7,11,13 r=1 so R_p is ONE cycle on odd classes mod 2p^s of length p^s;")
print("   the row period equals p (3 for 3n+1, 5 for 5n+1, 7 for 7n+1).  The height shift per braid step is d=ord_p(2).")
print("   p=1: q=2, c=1, R_1(n)=2n+1, F_1(R_1 n) = 2 F_1(n).  The fibre of core u under F_1 is n=2^{h+1}u-1 for EVERY h>=0")
print("   (no congruence gate, no parity of h, height shift 1): R_1 = tau o D o tau^{-1} with tau(x)=x-1, D(x)=2x,")
print("   i.e. the p=1 'braid' is the inherited doubling forest itself, shifted by one.  Its modulus 2*1^s is 2, so")
print("   the odometer statement is empty; on odd classes mod 6 R_1 has orbits (1,3) and (5): NO single row cycle.")
for u in range(1, 200, 2):
    for h in range(0, 8):
        n = 2 ** (h + 1) * u - 1
        check(core_height(Fp(1, n)) == (u, h), "F_1 fibre")
        check(2 * n + 1 == 2 ** (h + 2) * u - 1, "R_1 climbs one height")
orb6 = sorted(set(tuple(sorted({(2 * x + 1) % 6, x})) for x in (1, 3, 5)))
check(((2 * 1 + 1) % 6, (2 * 3 + 1) % 6, (2 * 5 + 1) % 6) == (3, 1, 5), "R_1 mod 6: 1->3->1, 5->5")

# bounded cycle census of the odd maps T_p (positive starts)


def Tp(p, n):
    m = p * n + 1
    return m // 2 ** v2(m)


print()
print("FINITE-EXACT positive cycles of the odd-to-odd maps T_p (starts n<=20000 odd, <=3000 steps, values capped 10^40):")
for p in (1, 3, 5, 7):
    cycles = set()
    escaped = 0
    for n0 in range(1, 20001, 2):
        seen_ = {}
        x, t = n0, 0
        while t < 3000 and x < 10 ** 40:
            if x in seen_:
                # extract cycle
                cyc = []
                y = x
                while True:
                    cyc.append(y)
                    y = Tp(p, y)
                    if y == x:
                        break
                m = min(cyc)
                i = cyc.index(m)
                cycles.add(tuple(cyc[i:] + cyc[:i]))
                break
            seen_[x] = t
            x = Tp(p, x)
            t += 1
        else:
            escaped += 1
    print("   p=%d: cycles found = %s ; starts not resolved within bounds = %d" % (p, sorted(cycles), escaped))
    if p == 1:
        check(cycles == {(1,)} and escaped == 0, "T_1 converges (trivial descent (n+1)/2 < n)")
    if p == 3:
        check(cycles == {(1,)} and escaped == 0, "T_3 on starts <= 20000")
    if p == 5:
        check({(1, 3), (13, 33, 83), (17, 43, 27)} <= cycles, "5n+1 known cycles")
    if p == 7:
        check((1,) in cycles, "7n+1 fixes 1")
print("   PROVED: T_1 converges for every n>=1 since (n+1)/2 < n for n>1.  Others: FINITE-EXACT only (unresolved starts")
print("   are divergence candidates within the cap, not proofs).")

# --------------------------------------------------------------------------
# S3  The three (core,height) rows to j<=60, the row law, and prominence
# --------------------------------------------------------------------------
banner("S3  (core,height) rows for j<=60, exact row law row(u,h) = rho(2^h u mod 9), prominence of 1,5,7")

JMAX = 60
rows = {}
for r in (1, 3, 5):
    rows[r] = [core_height(F(6 * j + r)) for j in range(0, JMAX + 1)]
    print("row %d mod 6, j=0..%d:" % (r, JMAX))
    print("   " + ",".join("(%d,%d)" % ch for ch in rows[r]))

# cross-check against inherited JSON
json_path = "/tmp/math-wt-collatz-mod6/05-knowledge/results/arithmetic_braids_20260917_collatz.json"
with open(json_path, "r") as fh:
    inh = json.load(fh)
for r in (1, 3, 5):
    inh_rows = [tuple(x) for x in inh["rows"][str(r)]]
    check(rows[r][:len(inh_rows)] == inh_rows, "row %d agrees with inherited JSON on j<=%d" % (r, len(inh_rows) - 1))
print("FINITE-EXACT: all three rows agree with arithmetic_braids_20260917_collatz.json for j<=34 (35 entries each).")
check(rows[3][35] == (5, 6) and rows[1][14] == (1, 7) and rows[5][24] == (7, 5), "spot values")

# row law
print()
print("row law.  For odd u with 3 !| u and admissible h (h odd iff u = 1 mod 3), the source n=(2^{h+1}u-1)/3 lies in row")
print("   r = rho(2^h u mod 9), rho(2)=1, rho(5)=3, rho(8)=5,   and its index is j = (2^h u - b_r)/9.")
cnt = 0
for u in range(1, 3001, 2):
    if u % 3 == 0:
        continue
    for h in range(0, 13):
        num = 2 ** (h + 1) * u - 1
        if num % 3:
            check((h % 2 == 1) == (u % 3 == 2), "parity gate: h odd iff u=1 mod 3 is admissible")
            continue
        n = num // 3
        r = n % 6
        b = (2 ** h * u) % 9
        check(ROW_OF_B[b] == r, "row law")
        check((2 ** h * u - B_OF_ROW[r]) % 9 == 0 and 6 * ((2 ** h * u - B_OF_ROW[r]) // 9) + r == n, "index law")
        check(core_height(F(n)) == (u, h), "round trip")
        cnt += 1
print("PROVED (+checked %d (u,h) pairs, u<3001, h<=12): row(u,h) = rho(2^h u mod 9).  Proof: F(n)=2^h u and F(6j+r)=9j+b_r," % cnt)
print("   so 2^h u = b_r mod 9 with b_r in {2,5,8} determined by r; rho inverts r -> b_r.  Since 4 = 2^2 has order 3 mod 9,")
print("   the three admissible heights h0, h0+2, h0+4 send u to the three DIFFERENT rows, and h -> h+6 returns:")
print("   row(u,h+2) = row(u,h) + 4 mod 6 (1->5->3->1), which is the R-braid R(n)=4n+1 read on rows.")

# base row table by u mod 18
print()
print("   base row (least admissible height h0) as a function of u mod 18:")
print("   u mod 18 | h0 | source n   | row")
base = {}
for um in (1, 5, 7, 11, 13, 17):
    h0 = 1 if um % 3 == 1 else 0
    n = (2 ** (h0 + 1) * um - 1) // 3
    base[um] = (h0, n % 6)
    print("   %-8d | %-2d | %-10s | %d" % (um, h0, "(%d*%d-1)/3=%d" % (2 ** (h0 + 1), um, n), n % 6))
check(base == {1: (1, 1), 7: (1, 3), 13: (1, 5), 5: (0, 3), 11: (0, 1), 17: (0, 5)}, "base row table")
for u in range(1, 2000, 2):
    if u % 3 == 0:
        continue
    h0, r0 = base[u % 18]
    check(((2 ** (h0 + 1) * u - 1) // 3) % 6 == r0, "base row depends on u mod 18 only")
print("PROVED: the base row depends on u mod 18 only (u mod 9 given odd), and row(u,h0+2t) = r0 + 4t mod 6.")

# mod-18 arithmetic-progression law inside a row at fixed height
print()
print("cores at fixed height h in row r form ONE odd residue class mod 18, u = b_r * 2^{-h} mod 9, listed at index step 2^{h+1}:")
for r in (1, 3, 5):
    for h in range(0, 6):
        want = (B_OF_ROW[r] * pow(2, -h, 9)) % 9
        js = [j for j in range(0, 2 ** 14) if rows_h(j, r, h) if False] if False else None
        js = []
        us = []
        for j in range(0, 2 ** 13):
            u, hh = core_height(9 * j + B_OF_ROW[r])
            if hh == h:
                js.append(j)
                us.append(u)
        check(all(u % 9 == want and u % 2 == 1 for u in us), "AP law mod 18")
        check(all(js[i + 1] - js[i] == 2 ** (h + 1) for i in range(len(js) - 1)), "index step 2^{h+1}")
        check(all(us[i + 1] - us[i] == 18 for i in range(len(us) - 1)), "core step 18")
        check(js[0] == (2 ** h * (want if want % 2 else want + 9) - B_OF_ROW[r]) // 9, "first index")
        if h <= 4:
            print("   row %d, h=%d: u = %2d mod 18, j = %2d mod %2d, cores %s..." % (r, h, want if want % 2 else want + 9, js[0], 2 ** (h + 1), us[:5]))
print("PROVED (+checked j<2^13): row r at height h = {odd u : u = b_r 2^{-h} mod 9}, an arithmetic progression of step 18")
print("   at index positions j = j0 mod 2^{h+1}; the row is the 2-adic interleaving of these six-periodic AP layers.")
print("   This is exactly the user's 'chain segments grow a pair at a time and are interleaved'.")

# index recurrence for a fixed core: j_{t+1} = 64 j_t + 7(1+3r)/2
print()
for r in (1, 3, 5):
    shift = 7 * (1 + 3 * r) // 2
    check((1 + 3 * r) % 2 == 0, "integrality")
    for u in range(1, 300, 2):
        if u % 3 == 0:
            continue
        # find first index of u in row r
        h = 0
        while True:
            num = 2 ** (h + 1) * u - 1
            if num % 3 == 0 and (num // 3) % 6 == r:
                break
            h += 1
        j0 = (2 ** h * u - B_OF_ROW[r]) // 9
        j1 = (2 ** (h + 6) * u - B_OF_ROW[r]) // 9
        check(j1 == 64 * j0 + shift, "index recurrence")
        check(core_height(F(6 * j0 + r)) == (u, h) and core_height(F(6 * j1 + r)) == (u, h + 6), "recurrence realizes u")
    print("   row %d: indices of a fixed core satisfy j_{t+1} = 64 j_t + %d  (7(1+3r)/2), heights h0+6t." % (r, shift))
print("PROVED: j_t = (2^{h0+6t} u - b_r)/9, so j_{t+1} - 64 j_t = 63 b_r/9 = 7 b_r = 14, 35, 56 for r=1,3,5. The 63 again.")

# prominence: first appearances of small cores
print()
print("first appearance (j, h) of each core u<=25 with 3 !| u in each row (j<=60 shown; '-' = beyond 60):")
first = {}
for r in (1, 3, 5):
    for j, (u, h) in enumerate(rows[r]):
        first.setdefault(u, {}).setdefault(r, (j, h))
print("   u  | row1 (j,h) | row3 (j,h) | row5 (j,h) | u mod 18")
for u in range(1, 26, 2):
    if u % 3 == 0:
        continue
    cells = [("(%d,%d)" % first.get(u, {}).get(r)) if first.get(u, {}).get(r) else "-" for r in (1, 3, 5)]
    print("   %-2d | %-10s | %-10s | %-10s | %d" % (u, cells[0], cells[1], cells[2], u % 18))
check(first[1] == {1: (0, 1), 3: (3, 5), 5: (0, 3)}, "core 1")
check(first[5] == {1: (2, 2), 3: (0, 0), 5: (8, 4)}, "core 5")
check(first[7] == {1: (6, 3), 3: (1, 1), 5: (24, 5)}, "core 7")
check(first[11] == {1: (1, 0), 3: (19, 4), 5: (4, 2)}, "core 11")
check(first[13] == {1: (46, 5), 3: (11, 3), 5: (2, 1)}, "core 13")
# what 'prominent' means exactly: the number of appearances of u in j<=J is about log_64
appear = {}
for r in (1, 3, 5):
    for j, (u, h) in enumerate(rows[r]):
        appear[u] = appear.get(u, 0) + 1
top = sorted(appear.items(), key=lambda kv: (-kv[1], kv[0]))[:8]
print("   appearance counts over the three rows, j<=60:", top)
check(top[0][0] == 1 and top[0][1] == 4, "core 1 appears 4 times (j<=60)")
print("PROVED: the index of core u at height h in its row is (2^h u - b_r)/9, so within j<=J the core u appears at exactly")
print("   the heights with 2^h u <= 9J + 8 in its parity class: about (1/2) log_2((9J+8)/u) times per row.  Small u is")
print("   prominent because 2^h u is small.  1, 5, 7, 11, 13, 17 are the least cores of the six classes mod 18;")
print("   1 and 7 are the least 1-mod-6 cores (odd heights), 5 and 11 the least 5-mod-6 cores (even heights); 7 is the")
print("   least 1-mod-6 core >1 and 5 the least 5-mod-6 core.  'Prominence' is size, not a dynamical property.")

# --------------------------------------------------------------------------
# S4  Negative rows
# --------------------------------------------------------------------------
banner("S4  Negative rows j=-40..-1 under (3n+1)/2 with signed odd cores; negation law; where the 3 minus cycles sit")

JN = 40
nrows = {}
for r in (1, 3, 5):
    nrows[r] = [(j, core_height(F(6 * j + r))) for j in range(-1, -JN - 1, -1)]
    print("row %d (n=6j+%d, j=-1..-%d):" % (r, r, JN))
    print("   " + ",".join("(%d,%d)" % ch for (_, ch) in nrows[r]))
print("   (n = 6j+1 for j<0 gives -5,-11,-17,...; 6j+5 gives -1,-7,-13,...; 6j+3 gives -3,-9,...; this is the user's")
print("    'the 1 chain contains -5 and so on, the 5 chain contains -1 and so on'.)")


def Fminus(n):
    return (3 * n - 1) // 2


for j in range(-JN, 0):
    for r in (1, 3, 5):
        n = 6 * j + r
        m = -n
        rp = m % 6
        check(rp == {1: 5, 3: 3, 5: 1}[r], "negation swaps rows 1<->5, fixes 3")
        check(F(n) == -Fminus(m), "F(-m) = -F_-(m)")
        u, h = core_height(F(n))
        um, hm = core_height(Fminus(m))
        check((u, h) == (-um, hm), "signed core = negated 3n-1 core, same height")
print("PROVED: -(6j+r) = 6(-j-1) + (6-r) so negation swaps rows 1 and 5 and fixes row 3 (inherited summand.md sec.7);")
print("   F(-m) = -(3m-1)/2, so the negative (core,height) rows are the 3n-1 rows of the positives with cores negated,")
print("   rows 1<->5 exchanged, and index j -> -j-1.  Heights are unchanged.")

# inverse fibre theorem for negative cores; the row law is one law over Z
cnt = 0
for u in range(-1, -600, -2):
    if u % 3 == 0:
        continue
    for h in range(0, 10):
        num = 2 ** (h + 1) * u - 1
        if num % 3:
            continue
        n = num // 3
        check(n < 0 and n % 2 != 0 and core_height(F(n)) == (u, h), "negative fibre")
        check(ROW_OF_B[(2 ** h * u) % 9] == n % 6, "row law over Z")
        check(F(R(n)) == 4 * F(n), "R braid over Z")
        cnt += 1
print("PROVED (+%d checks): the inverse-fibre theorem (B1),(B2) and row(u,h)=rho(2^h u mod 9) hold verbatim for negative" % cnt)
print("   odd cores (all identities are over Z); e.g. u=-1: h even, n=-1 (row 5), -3 (row 3), -11 (row 1), -43, ...")

# where the negative cycle members appear
neg_cycles = [(-1,), (-5, -7), (-17, -25, -37, -55, -41, -61, -91)]
for cyc in neg_cycles:
    for x in cyc:
        y = F(x)
        u, h = core_height(y)
        check(u in cyc, "cycle closes")
print("negative cycle members as cores in the negative rows (row of the SOURCE n=(2^{h+1}u-1)/3, all j with |j|<=40):")
for cyc in neg_cycles:
    for u in cyc:
        hits = []
        for r in (1, 3, 5):
            for (j, (uu, h)) in nrows[r]:
                if uu == u:
                    hits.append((r, j, h))
        print("   u=%-4d source row %d (member of cycle %s): appears at (row,j,h) = %s" % (u, u % 6, cyc, sorted(hits)))
        check(len(hits) >= 1, "member appears")
member_rows = sorted(set(u % 6 for cyc in neg_cycles for u in cyc))
check(member_rows == [1, 5], "no cycle member is 3 mod 6")
print("PROVED: no cycle member of any sign is 3 mod 6 (multiples of 3 are never T-images).  The three minus cycles")
print("   have lengths 1,2,7 and their members lie in rows {5}, {1,5}, {1,5}: the count 3 of cycles is NOT indexed by")
print("   the 3 rows; no map, only the number 3 (cycle completeness is OPEN; bounded census inherited).")

# --------------------------------------------------------------------------
# S5  The typing table: real maps and the 63 <-> -7/4 identity
# --------------------------------------------------------------------------
banner("S5  Typing table: every 'three' and 'two' in the prompt; the exact identity c(c+1)^2 = -63/64 at c=-7/4")


def new_prime_free(num, earlier):
    """True iff every prime of |num| divides some earlier numerator (exact, no factoring)."""
    g = abs(num)
    if g == 0:
        return True
    if g == 1:
        return True
    N = 1
    for e in earlier:
        N *= abs(e)
    if N == 0:
        N = 1
    while True:
        d = gcd(g, N)
        if d == 1:
            break
        g //= d
    return g == 1


# 2x+1 from 0: 2^n - 1; Bang/Zsigmondy exception at n=6
print("2x+1 iterated from 0 gives 2^n-1.  Terms with NO new prime (gcd-stripping, exact):")
earlier = []
nonew = []
for n in range(1, 65):
    t = 2 ** n - 1
    if n > 1 and new_prime_free(t, earlier):
        nonew.append(n)
    earlier.append(t)
print("   n in 2..64 with no new prime:", nonew, "  (63 = 3^2 * 7, 3 | 2^2-1, 7 | 2^3-1)")
check(nonew == [6], "only n=6 up to 64")
print("   CITED: Bang (1886)/Zsigmondy (1892): 2^n-1 has a primitive prime divisor for every n except n=1 and n=6.")
print("   FINITE-EXACT here for n<=64.")

# x^2+c orbits of 0


def orbit0(c, N):
    x = Fraction(0)
    out = []
    for _ in range(N):
        x = x * x + c
        out.append(x)
    return out


print()
print("x^2+c from 0, first terms; which terms introduce no new prime in the NUMERATOR:")
for c in (Fraction(-7, 4), Fraction(-29, 16), Fraction(-1, 4), Fraction(-1), Fraction(-2), Fraction(0)):
    orb = orbit0(c, 6)
    earlier = []
    flags = []
    for i, x in enumerate(orb):
        flags.append(new_prime_free(x.numerator, earlier) if i > 0 else None)
        earlier.append(x.numerator)
    shown = ", ".join(str(x) for x in orb[:4])
    print("   c=%-7s: %s ... ; no-new-prime at terms %s (of 2..6)" % (c, shown, [i + 1 for i, f in enumerate(flags) if f]))
    if c == Fraction(-7, 4):
        check(orb[2] == Fraction(-7, 256) and flags[2] is True and flags[1] is False, "x^2-7/4 third term")
    if c == Fraction(-29, 16):
        check(flags[2] is False, "x^2-29/16 third term has new primes 5,7,23")
c = Fraction(-7, 4)
check(c * (c + 1) ** 2 == Fraction(-63, 64), "c(c+1)^2 = -63/64")
check(c * (c * (c + 1) ** 2 + 1) == Fraction(-7, 256), "f^3(0) = c/64")
print("PROVED: for f=x^2+c, f(0)=c, f^2(0)=c(c+1), f^3(0)=c(c(c+1)^2+1).  At c=-7/4: c(c+1)^2 = (-7/4)(9/16) = -63/64,")
print("   so f^3(0) = c/2^6 = -7/256: the third term repeats the prime 7 of c EXACTLY because 7 * 3^2 = 63 = 2^6 - 1.")
print("   This is the same integer factorization 63 = 3^2*7 that gives ord_9(2)=6 (the row braid) and the n=6 Bang")
print("   exception of 2x+1.  Typed map: source = {2x+1 orbit of 0, term 6}; target = {x^2-7/4 orbit of 0, term 3};")
print("   map = 63/64 = -c(c+1)^2, preserved predicate = 'numerator has no primitive prime divisor', mechanism = the")
print("   single equation 2^6 - 1 = 7 * 3^2 (7 from c, 3^2 from (c+1)^2).  Lost: the dynamics (2x+1 is affine).")

# is -7/4 the unique c=-a/2^k with f^3(0) = c * 2^{-m} (strong form)?
sols = []
for k in range(1, 31):
    b = 2 ** k
    for a in range(-b * 2, 3 * b):
        if a == 0 or gcd(abs(a), b) != 1:
            continue
        cc = Fraction(-a, b)
        val = cc * (cc + 1) ** 2 + 1
        if abs(val.numerator) == 1 and val.denominator & (val.denominator - 1) == 0:
            sols.append((cc, val))
print("   c=-a/2^k, k<=30, |a|<3*2^k, with c(c+1)^2+1 = +-2^{-m}:", sols)
check(sols == [(Fraction(-7, 4), Fraction(1, 64))] or set(sols) >= {(Fraction(-7, 4), Fraction(1, 64))}, "-7/4 found")
# weaker: any c=a/b (|a|<=120, b<=64) whose THIRD term has no new numerator prime
weak = []
for b in range(1, 65):
    for a in range(-120, 121):
        if a == 0 or gcd(abs(a), b) != 1:
            continue
        cc = Fraction(a, b)
        orb = orbit0(cc, 3)
        if orb[1].numerator == 0 or orb[0].numerator == 0:
            continue  # degenerate (0 has no primes): c=-1 only
        if new_prime_free(orb[2].numerator, [orb[0].numerator, orb[1].numerator]):
            weak.append(cc)
print("   FINITE-EXACT: c=a/b, |a|<=120, b<=64, nondegenerate, third term without new numerator prime:", weak)
check(Fraction(-7, 4) in weak and Fraction(-2) in weak, "-7/4 and -2 present")
print("   (c=-2 is the Chebyshev/PCF case 0->-2->2->2; c=-1 is degenerate since f^2(0)=0.)")

# rational PCF quadratics x^2+c: exactly c in {0,-1,-2}
print()
pcf = []
for cint in range(-60, 61):
    x, seen_ = 0, set()
    per = False
    for _ in range(200):
        if x in seen_:
            per = True
            break
        seen_.add(x)
        x = x * x + cint
        if abs(x) > 10 ** 12:
            break
    if per:
        pcf.append(cint)
check(pcf == [-2, -1, 0], "integer PCF parameters")
print("PROVED: x^2+c with c in Q has a preperiodic critical point 0 iff c in {0,-1,-2}.  (c non-integral: the orbit of 0")
print("   has denominators b, b^2, b^4,... in lowest terms, never periodic; c>=1: increasing; c<=-3: f^2(0)=c^2+c>=|c|+3")
print("   and x>=|c|+1 implies x^2+c > x, escape.)  Checked |c|<=60.  The user's three systems x^2-{0,1,2} are exactly")
print("   the rational postcritically finite quadratics: critical orbits 0 fixed; 0<->-1; 0->-2->2 fixed.")
# integer preperiodic points, finite lists
for cint, want in ((0, [-1, 0, 1]), (-1, [-1, 0, 1]), (-2, [-2, -1, 0, 1, 2])):
    pp = []
    for x0 in range(-50, 51):
        x, seen_ = x0, set()
        ok = False
        for _ in range(100):
            if x in seen_:
                ok = True
                break
            seen_.add(x)
            x = x * x + cint
            if abs(x) > 10 ** 9:
                break
        if ok:
            pp.append(x0)
    check(pp == want, "integer preperiodic points of x^2%+d" % cint)
    print("   integer preperiodic points of x^2%+d: %s" % (cint, want))
print("   (sibling lane pythagorean_semicircle.out S5 has the edge lists and the sign audit x^2+{1,2} -> x^2-{1,2}.)")

# THM-4146 check: G_D(y) = (y^2-29)/4 cycle -7 -> 5 -> -1


def G(y):
    return Fraction(y * y - 29, 4)


check(G(-7) == 5 and G(5) == -1 and G(-1) == -7, "THM-4146 (33) cycle")
ff = lambda x: x * x - Fraction(29, 16)
check(ff(Fraction(-7, 4)) == Fraction(5, 4) and ff(Fraction(5, 4)) == Fraction(-1, 4) and ff(Fraction(-1, 4)) == Fraction(-7, 4), "x^2-29/16 cycle")
print("   THM-4146 (29)-(33) re-checked: (y^2-29)/4 has the cycle -7 -> 5 -> -1 (y=4x), forced by (3,4,5): -(3+4), 5, -(4-3).")

# ord_9(2)=6 <-> 3-cycle of x^2-2 on 2cos(2 pi k/9): the row 3-cycle
print()
print("rows <-> the 3-cycle of x^2-2 (sibling lane S4 has the fold): b -> 2b mod 9 on {2,5,8} modulo sign:")
step = {}
for b in (2, 5, 8):
    tb = (2 * b) % 9
    tb = tb if tb in (2, 5, 8) else 9 - tb
    step[ROW_OF_B[b]] = ROW_OF_B[tb]
print("   one squaring step on rows:", step, " ; R (x4) on rows:", {r: (r + 4) % 6 for r in (1, 3, 5)})
check(step == {5: 1, 1: 3, 3: 5}, "squaring 3-cycle 5->1->3->5")
check({r: step[step[r]] for r in (1, 3, 5)} == {r: (r + 4) % 6 for r in (1, 3, 5)}, "R on rows = squaring step squared")
print("PROVED: on the three rows, the R-braid (multiply target by 4) is the SQUARE of the angle-doubling step of x^2-2 on")
print("   {2cos(2 pi b/9)} (multiply by 2 mod +-9); both are the group (Z/9)^*/{+-1} = Z/3.  Preserved predicate:")
print("   multiplication by 2 on exponents mod 9.  Lost: everything about n except F(n) mod 9.  This is a genuine but")
print("   low-content map (one cyclic group of order three in two guises).")

# --------------------------------------------------------------------------
# the typing table proper
# --------------------------------------------------------------------------
print()
print("TYPING TABLE.  Each 'three'/'two' of the prompt; for each pair a MAP (with preserved predicate) or 'number only'.")
T3 = [
    ("T1", "3 rows mod 6 (1,3,5)", "the p=3 odd classes mod 2p; count = multiplier p", "PROVED S1/S2"),
    ("T2", "3 image residues 2,5,8 mod 9", "j-preserving bijection with T1: F(6j+r)=9j+b_r; = classes 2 mod 3 in Z/9", "PROVED S1"),
    ("T3", "3 exponent classes 1,5,3 mod 6 of 2^k", "bijection with T2 via 2^k mod 9; needs ord_9(2)=6", "PROVED S1"),
    ("T4", "period 3 of R=4n+1 on rows", "ord_9(4)=3 = ord_9(2)/ord_3(2); R^3=64n+21, 63=3*21", "PROVED S1"),
    ("T5", "3 = 63 mod 7 sees k mod 3", "(Z/63)^*=(Z/9)^*x(Z/7)^*; on odd k, k mod 6 <-> k mod 3", "PROVED S1"),
    ("T6", "3 negative (3n-1) cycles", "lengths 1,2,7; members in rows {5},{1,5},{1,5}; number only", "FINITE-EXACT census inherited; completeness OPEN"),
    ("T7", "3 solutions of F=S+U (p, p^3, p^2qr)", "exponent-box balance prod(a_i+1)=2^r+r+1 (divisors.md DB1); number only", "PROVED inherited"),
    ("T8", "3 almost-prime classes A,B,C", "Omega=1,2,3 labels; p^2qr is a 4-almost-prime (divisors.md DB3); number only", "REFUTED as C^2B, inherited"),
    ("T9", "3 systems x^2-{0,1,2}", "exactly the rational PCF quadratics; x^2 = squaring forest (summand.md sec.2)", "PROVED S5"),
    ("T10", "3 special values 63, -7/4, -29/16", "63<->-7/4: c(c+1)^2=-63/64 (S5, new); 63<->rows: ord_9(2)=6; -29/16<->3:4:5: THM-4146", "PROVED S5"),
    ("T11", "3-cycle -7/4->5/4->-1/4", "unique AP-supported rational quadratic 3-cycle (THM-4139); 3:4:5 forced (THM-4146)", "CITED canon"),
    ("T12", "3:4:5", "signed template -(a+b),h,-(b-a) forces (3k,4k,5k), D=29k^2 (THM-4146 4.2)", "CITED canon"),
    ("T13", "3 parameters e/d, l, theta", "one-parameter chart in theta (sibling pythagorean_semicircle.out S3): no braid", "PROVED sibling"),
    ("T14", "3 of 3n+1", "the multiplier; ord_3(2)=2 gives h parity; 9=3^2 gives the row braid", "PROVED S1"),
    ("T15", "3-cycle of x^2-2 at 2cos(2 pi k/9)", "rows: R on rows = square of the angle-doubling 3-cycle (S5)", "PROVED S5 (low content)"),
]
T2 = [
    ("D1", "doubling u->2u", "complement of strict summand shadow = doubling forest = reversed halving (summand.md sec.1)", "PROVED inherited"),
    ("D2", "2-way splitting +-x / negation", "3n+1 <-> 3n-1 conjugacy swaps rows 1,5 (S4); sandwich k->-k (divisors SW2): different objects, no map between them", "PROVED both; no cross map"),
    ("D3", "2 mixed sandwich cells N_12, N_21", "differ already at K=4 (divisors SW2); number only w.r.t. Collatz", "PROVED inherited"),
    ("D4", "2-cycles: {1,2} shortcut, {-5,-7}, {0,-1} of x^2-1", "no map found; number only", "no map"),
    ("D5", "period-2 height shift h->h+2", "= ord_3(2)=2: 4 is the least power of two = 1 mod 3; general p: shift ord_p(2) (S2)", "PROVED S2"),
    ("D6", "(3n+1)/2 single halving", "2 | 3n+1 for odd n: the parity gate; companion (n+1)/2 = F_1(n) (S2)", "PROVED S2"),
    ("D7", "2 = ord_3(2) inside 6 = ord_9(2)", "6 = 2*3: gate x braid; 7 sees the 3, 3 sees the 2, 9 sees both (S1)", "PROVED S1"),
]
for row_ in T3 + T2:
    print("   %-4s %-42s | %-88s | %s" % row_)

print()
print("PAIRWISE VERDICTS among the threes (only pairs with an actual map are listed; all other pairs: number only):")
pairs = [
    ("T1<->T2<->T3<->T4<->T5<->T14", "one mechanism: multiplier 3 (count) and 63=2^6-1=3^2*7 (which row / period / what 7 sees)"),
    ("T10(63)<->T10(-7/4)", "c(c+1)^2 = -63/64: same factorization 63 = 7*3^2, predicate 'no primitive prime divisor'"),
    ("T10(-29/16)<->T11<->T12", "THM-4139/4146: AP three-cycle and 3:4:5 Pythagorean forcing"),
    ("T1<->T15", "(Z/9)^*/{+-1}: R on rows is the square of the x^2-2 angle-doubling 3-cycle"),
    ("T9(x^2)<->D1", "squaring forest <-> doubling forest via exponents (summand.md sec.2)"),
    ("T6, T7, T8, T13 vs everything", "NO MAP FOUND: only the number 3 (T7 is an exponent-box identity, T13 is one-dimensional)"),
]
for a, b in pairs:
    print("   %-36s : %s" % (a, b))

# --------------------------------------------------------------------------
print()
print("checks passed:", CHECKS)
src = open(__file__, "rb").read()
print("script sha256:", hashlib.sha256(src).hexdigest())
print("ALL CHECKS PASSED")
