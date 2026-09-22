#!/usr/bin/env python3
"""collatz_mod6_20260917_row_braid_typing_audit_proof-audit.py

Independent proof-audit of lane row_braid_typing (session collatz-mod6-20260917).
Does NOT import the explorer's script.  Focus: quantifier/boundary cases the
explorer's proofs skip, the exact criterion behind the -7/4 'no new prime'
predicate, and a certified Thue-equation solve (PARI/gp) closing open question 1.
All arithmetic exact.  RAM << 1 GB, runtime < 1 minute.
"""
import subprocess
import sys
from fractions import Fraction
from math import gcd

FAILS = []


def rep(tag, ok, msg):
    print(("OK   " if ok else "FAIL ") + tag + "  " + msg)
    if not ok:
        FAILS.append(tag)


def order(a, m):
    k, x = 1, a % m
    while x != 1:
        x = x * a % m
        k += 1
    return k


def v(n, p):
    n = abs(n)
    k = 0
    while n and n % p == 0:
        n //= p
        k += 1
    return k


def F(n):          # (3n+1)/2, odd n, any sign
    assert n % 2
    return (3 * n + 1) // 2


def core(m):       # m = 2^h u, u odd (sign on u)
    h = v(m, 2)
    return m // 2 ** h, h


B = {1: 2, 3: 5, 5: 8}
RHO = {2: 1, 5: 3, 8: 5}

print("== RB1/RB3: orders and what 3, 7, 9 see ==")
rep("RB1", (order(2, 3), order(2, 9), order(2, 7), order(2, 63)) == (2, 6, 3, 6), "ord_3,9,7,63(2) = 2,6,3,6")
rep("RB1", min(k for k in range(1, 40) if (2 ** k - 1) % 9 == 0) == 6, "first Mersenne divisible by 9 is 2^6-1")
img = sorted(set(pow(2, k, 63) for k in range(6)))
rep("RB1", sorted(x % 9 for x in img) == [1, 2, 4, 5, 7, 8] and sorted(x % 7 for x in img) == [1, 1, 2, 2, 4, 4], "<2> mod 63 -> (Z/9)^* iso, -> {1,2,4} mod 7 two-to-one")
# 2^k is an F-image of an odd integer iff (2^{k+1}-1)/3 is an odd integer  (own derivation, k up to 300)
tab = {}
for k in range(0, 301):
    num = 2 ** (k + 1) - 1
    is_img = num % 3 == 0 and (num // 3) % 2 == 1
    rep("RB3", is_img == (k % 2 == 1), "k=%d image iff k odd" % k) if k < 3 or not (is_img == (k % 2 == 1)) else None
    if is_img:
        tab.setdefault(k % 6, set()).add(((num // 3) % 6, pow(2, k, 9), pow(2, k, 7)))
rep("RB3", tab == {1: {(1, 2, 2)}, 5: {(3, 5, 4)}, 3: {(5, 8, 1)}}, "k mod 6 = 1,5,3 -> source row 1,3,5, image 2,5,8 mod 9, 2,4,1 mod 7 (k<=300)")

print("== RB4: R=4n+1 bookkeeping, boundary n even ==")
# polynomial composition by coefficients: R(n)=4n+1 -> R^3 coefficients
a, b = 1, 0
for _ in range(3):
    a, b = 4 * a, 4 * b + 1
rep("RB4", (a, b) == (64, 21), "R^3(n) = 64n+21 (coefficient composition)")
rep("RB4", all((64 * n + 21 - n) == 21 * (3 * n + 1) for n in range(-50, 50)), "R^3(n)-n = 21(3n+1) identically")
odd_ok = all((21 * (3 * n + 1)) % 42 == 0 for n in range(-999, 1000, 2))
even_bad = all((21 * (3 * n + 1)) % 42 != 0 for n in range(-998, 1000, 2))
rep("RB4", odd_ok and even_bad, "R^3 fixes n mod 42 for ODD n only (even n: 21(3n+1) is odd) -> quantifier 'odd n' needed")
rep("RB4", all((64 * F(n)) % 63 == F(n) % 63 and (4 * F(n)) % 63 != F(n) % 63 for n in range(-999, 1000, 2)), "target: x64 fixes mod 63, x4 moves it")
orb = set()
for x0 in range(7):
    x, c = x0, []
    while x not in c:
        c.append(x)
        x = (4 * x + 1) % 7
    orb.add(tuple(sorted(c)))
rep("RB4", orb == {(2,), (0, 1, 5), (3, 4, 6)} and all(((3 * n + 1) % 7 == 0) == (n % 7 == 2) for n in range(1, 1000, 2)), "R mod 7 orbits; 7 | F(n) iff n = 2 mod 7")
rep("RB4", sorted(set(pow(2, k, 63) for k in range(1, 7)) & set(F(n) % 63 for n in range(1, 126, 2))) == [2, 8, 32], "powers of two among the 21 target classes mod 63")
# period-3-on-rows mechanism: (4^t-1)(3n+1)/3 = 0 mod 6 iff 9 | 4^t-1 iff 3 | t
rep("RB4", all(((4 ** t - 1) % 9 == 0) == (t % 3 == 0) for t in range(1, 40)), "9 | 4^t-1 iff 3 | t (so row period = ord_9(4) = 3)")

print("== RB5: Wieferich boundary ==")


def ord_from_divisors(a, m, n):  # order of a mod m given multiple n of the order
    ds = sorted(d for d in range(1, n + 1) if n % d == 0) if n < 10 ** 6 else None
    if ds is None:
        # factor n crudely (n = p(p-1), small)
        fs = {}
        t = n
        d = 2
        while d * d <= t:
            while t % d == 0:
                fs[d] = fs.get(d, 0) + 1
                t //= d
            d += 1
        if t > 1:
            fs[t] = fs.get(t, 0) + 1
        o = n
        for q in fs:
            while o % q == 0 and pow(a, o // q, m) == 1:
                o //= q
        return o
    for d in ds:
        if pow(a, d, m) == 1:
            return d


for p in (3, 5, 7, 11, 13, 17, 19, 1093, 3511):
    o1 = ord_from_divisors(2, p, p - 1)
    o2 = ord_from_divisors(2, p * p, p * (p - 1))
    wief = pow(2, p - 1, p * p) == 1
    inv = pow(2, -1, p)
    rows = len(set(pow(2, k, p * p) for k in range(o1 - 1, o2, o1)))
    rep("RB5", rows == o2 // o1 and (rows == p) == (not wief) and rows in (1, p), "p=%d ord_p=%d ord_p2=%d rows=%d wieferich=%s" % (p, o1, o2, rows, wief))
rep("RB5", pow(2, -1, 1093 ** 2) == 597325 and pow(2, -1, 3511 ** 2) == 6163561, "single classes are 2^{-1} mod p^2")
rep("RB5", [p for p in range(3, 5000, 2) if all(p % q for q in range(2, int(p ** .5) + 1)) and pow(2, p - 1, p * p) == 1] == [1093, 3511], "Wieferich primes below 5000 are exactly 1093, 3511")

print("== RB6: tower, diagonal over Z ==")


def Fp(p, n):
    return (p * n + 1) // 2


rep("RB6", all(Fp(p + 2, n) == Fp(p, n) + n and Fp(p, n) - n == Fp(p - 2, n) for p in range(-9, 40, 2) for n in range(-999, 1000, 2)), "F_{p+2}=F_p+n; companion F_p(n)-n = F_{p-2}(n) = ((p-2)n+1)/2 (this is summand.md (12)'s b-coordinate)")
diag_all = sorted((p, n) for p in range(-99, 100, 2) for n in range(-999, 1000, 2) if Fp(p - 2, n) == n)
diag_pos = [(p, n) for (p, n) in diag_all if n > 0]
rep("RB6", diag_pos == [(3, 1)], "n>=1: diagonal = [(3,1)]")
rep("RB6", diag_all == [(3, 1), (5, -1)], "over Z: (p-4)n=-1 has TWO solutions %s: (5,-1) is F_3(-1)=-1, the negative trivial cycle -> RB6 quantifier must say n>=1" % diag_all)

print("== RB7: R_p, p=1 conjugacy ==")
for p, q, c in ((1, 2, 1), (3, 4, 1), (5, 16, 3), (7, 8, 1), (11, 1024, 93), (13, 4096, 315)):
    d = q.bit_length() - 1
    lq = min(dd for dd in range(1, 400) if (2 ** dd - 1) % p == 0)
    rep("RB7", 2 ** lq == q and (q - 1) // p == c and all(Fp(p, q * n + c) == q * Fp(p, n) for n in range(-499, 500, 2)) and (p == 1 or v(q - 1, p) == 1), "p=%d q=2^%d c=%d r=%s, F_p(R_p n)=q F_p(n)" % (p, d, c, "n/a" if p == 1 else v(q - 1, p)))
    for s in (1, 2, 3):
        M = 2 * p ** s
        x, t = 1, 0
        while True:
            x = (q * x + c) % M
            t += 1
            if x == 1:
                break
        rep("RB7", t == p ** s, "p=%d s=%d orbit of 1 has length %d = p^s = #odd classes -> single cycle" % (p, s, t))
rep("RB7", all(2 * (x + 1) - 1 == 2 * x + 1 for x in range(-50, 50)) and all(core(Fp(1, 2 ** (h + 1) * u - 1)) == (u, h) for u in range(-99, 100, 2) for h in range(8)), "R_1 = tau D tau^-1 (tau(x)=x-1); F_1 fibre n=2^{h+1}u-1 all h, both signs")
rep("RB7", {r: (2 * r + 1) % 6 for r in (1, 3, 5)} == {1: 3, 3: 1, 5: 5}, "R_1 mod 6 orbits (1,3),(5)")

print("== RB8: T_1 descent, bounded census recheck (p=5,7) ==")
rep("RB8", all((n + 1) // 2 ** v(n + 1, 2) < n for n in range(3, 20001, 2)) and (1 + 1) // 2 == 1, "T_1(n) = oddpart(n+1) <= (n+1)/2 < n for n>1; T_1(1)=1")


def census(p, N=20000, steps=3000, cap=10 ** 40):
    cyc, esc = set(), 0
    for n0 in range(1, N + 1, 2):
        seen = set()
        x, t = n0, 0
        while t < steps and x < cap:
            if x in seen:
                c, y = [], x
                while True:
                    c.append(y)
                    m = p * y + 1
                    y = m // 2 ** v(m, 2)
                    if y == x:
                        break
                cyc.add(tuple(sorted(c)))
                break
            seen.add(x)
            m = p * x + 1
            x = m // 2 ** v(m, 2)
            t += 1
        else:
            esc += 1
    return cyc, esc


c5, e5 = census(5)
c7, e7 = census(7)
rep("RB8", c5 == {(1, 3), (13, 33, 83), (17, 27, 43)} and e5 == 9605, "p=5 cycles %s escaped %d" % (sorted(c5), e5))
rep("RB8", c7 == {(1,)} and e7 == 9982, "p=7 cycles %s escaped %d" % (sorted(c7), e7))

print("== RB9/RB10: row law, layers, index recurrence (both signs) ==")
cnt = 0
for u in list(range(1, 2001, 2)) + list(range(-1, -2001, -2)):
    if u % 3 == 0:
        continue
    for h in range(0, 12):
        num = 2 ** (h + 1) * u - 1
        if num % 3:
            continue
        n = num // 3
        r = n % 6
        ok = RHO[(2 ** h * u) % 9] == r and 6 * ((2 ** h * u - B[r]) // 9) + r == n and core(F(n)) == (u, h)
        ok &= ((h % 2 == 1) == (u % 3 == 1))
        if h + 2 < 12:
            n2 = (2 ** (h + 3) * u - 1) // 3
            ok &= n2 % 6 == (r + 4) % 6 and n2 == 4 * n + 1
        if not ok:
            rep("RB9", False, "u=%d h=%d" % (u, h))
        cnt += 1
rep("RB9", True, "row law, index law, admissibility parity, row(u,h+2)=row+4 (=R): %d pairs, both signs" % cnt)
base = {}
for um in (1, 5, 7, 11, 13, 17):
    h0 = 1 if um % 3 == 1 else 0
    base[um] = (h0, ((2 ** (h0 + 1) * um - 1) // 3) % 6)
rep("RB9", base == {1: (1, 1), 7: (1, 3), 13: (1, 5), 5: (0, 3), 11: (0, 1), 17: (0, 5)}, "base row table by u mod 18")
for r in (1, 3, 5):
    for h in range(0, 7):
        want = B[r] * pow(2, -h, 9) % 9
        js = [j for j in range(2 ** 12) if v(9 * j + B[r], 2) == h]
        us = [(9 * j + B[r]) >> h for j in js]
        ok = all(u % 18 == (want if want % 2 else want + 9) for u in us) and all(js[i + 1] - js[i] == 2 ** (h + 1) for i in range(len(js) - 1)) and all(us[i + 1] - us[i] == 18 for i in range(len(us) - 1))
        rep("RB10", ok, "row %d h=%d: one class mod 18, index step 2^{h+1}, core step 18" % (r, h)) if not ok or h == 0 else None
rep("RB10", all((2 ** (h + 6) * u - B[r]) // 9 == 64 * ((2 ** h * u - B[r]) // 9) + 7 * B[r] for r in (1, 3, 5) for h in range(0, 8) for u in range(1, 300, 2) if (2 ** h * u - B[r]) % 9 == 0), "j_{t+1} = 64 j_t + 7 b_r (14,35,56)")

print("== RB11: rows j<=60, first appearances, counts ==")
rows = {r: [core(F(6 * j + r)) for j in range(61)] for r in (1, 3, 5)}
import json
with open("/tmp/math-wt-collatz-mod6/05-knowledge/results/arithmetic_braids_20260917_collatz.json") as fh:
    J = json.load(fh)
rep("RB11", all([tuple(x) for x in J["rows"][str(r)]] == rows[r][:len(J["rows"][str(r)])] for r in (1, 3, 5)), "agree with inherited JSON (j<=34)")
first = {}
for r in (1, 3, 5):
    for j, (u, h) in enumerate(rows[r]):
        first.setdefault(u, {}).setdefault(r, (j, h))
rep("RB11", first[1] == {1: (0, 1), 3: (3, 5), 5: (0, 3)} and first[5] == {1: (2, 2), 3: (0, 0), 5: (8, 4)} and first[7] == {1: (6, 3), 3: (1, 1), 5: (24, 5)} and first[11] == {1: (1, 0), 3: (19, 4), 5: (4, 2)} and first[13] == {1: (46, 5), 3: (11, 3), 5: (2, 1)}, "first appearances of 1,5,7,11,13")
cnt = {}
for r in (1, 3, 5):
    for (u, h) in rows[r]:
        cnt[u] = cnt.get(u, 0) + 1
rep("RB11", cnt[1] == 5 and cnt[5] == 4 and sorted(h for r in (1, 3, 5) for (u, h) in rows[r] if u == 1) == [1, 3, 5, 7, 9], "core 1 five times (h=1,3,5,7,9), core 5 four times")

print("== RB12: negative rows ==")
rep("RB12", all((-(6 * j + r)) % 6 == {1: 5, 3: 3, 5: 1}[r] and F(6 * j + r) == -((3 * (-(6 * j + r)) - 1) // 2) for j in range(-500, 500) for r in (1, 3, 5)), "negation swaps rows 1,5 fixes 3; F(-m) = -(3m-1)/2")
nrows = {r: {j: core(F(6 * j + r)) for j in range(-40, 0)} for r in (1, 3, 5)}
cycles = [(-1,), (-5, -7), (-17, -25, -37, -55, -41, -61, -91)]
for cyc in cycles:
    rep("RB12", all(core(F(cyc[i]))[0] == cyc[(i + 1) % len(cyc)] for i in range(len(cyc))), "cycle %s closes in order" % (cyc,))
pos = {u: sorted((r, j, h) for r in (1, 3, 5) for j, (uu, h) in nrows[r].items() if uu == u) for cyc in cycles for u in cyc}
rep("RB12", pos[-1] == [(1, -2, 4), (3, -29, 8), (3, -1, 2), (5, -8, 6), (5, -1, 0)] and pos[-5] == [(1, -18, 5), (3, -5, 3), (5, -2, 1)] and pos[-7] == [(1, -1, 0), (3, -13, 4), (5, -4, 2)] and pos[-91] == [(5, -11, 0)], "cycle-member placements |j|<=40")
rep("RB12", all(u % 3 != 0 for cyc in cycles for u in cyc) and all(core(F(n))[0] % 3 != 0 for n in range(-9999, 10000, 2)), "no T-image core (either sign) is divisible by 3")

print("== RB13: -7/4 identity, exact 'no new prime' criterion, Thue closure ==")


def orbit0(c, N):
    x, out = Fraction(0), []
    for _ in range(N):
        x = x * x + c
        out.append(x)
    return out


def new_prime_free(num, earlier):
    g = abs(num)
    if g <= 1:
        return True
    N = 1
    for e in earlier:
        N *= abs(e)
    while True:
        d = gcd(g, N)
        if d == 1:
            break
        g //= d
    return g == 1


c = Fraction(-7, 4)
rep("RB13", c * (c + 1) ** 2 == Fraction(-63, 64) and orbit0(c, 3)[2] == c / 64 == Fraction(-7, 256), "c(c+1)^2=-63/64, f^3(0)=c/2^6")
rep("RB13", orbit0(c, 4)[3] == Fraction(-114639, 65536) and 114639 == 3 * 7 * 53 * 103, "term 4 = -3*7*53*103/2^16 (link stops at term 3)")
# exact criterion: c=-a/b lowest terms, b>=1, c not in {0,-1}: third numerator new-prime-free iff E := b^3 - a(b-a)^2 = +-1
mism, crit_sols = [], []
for bb in range(1, 151):
    for aa in range(-300, 301):
        if gcd(abs(aa), bb) != 1 or aa == 0 or (bb == 1 and aa == 1):
            continue
        cc = Fraction(-aa, bb)
        o = orbit0(cc, 3)
        E = bb ** 3 - aa * (bb - aa) ** 2
        pred = new_prime_free(o[2].numerator, [o[0].numerator, o[1].numerator])
        if pred != (abs(E) == 1):
            mism.append(cc)
        if abs(E) == 1:
            crit_sols.append(cc)
        # sanity: third numerator is -a*E, coprime to a(b-a)
        if o[2] != Fraction(-aa * E, bb ** 4) or gcd(abs(E), abs(aa * (bb - aa))) != 1:
            mism.append(("formula", cc))
rep("RB13", mism == [] and sorted(crit_sols) == [Fraction(-2), Fraction(-7, 4)], "PROVED criterion: third numerator = -a*E, E=b^3-a(b-a)^2 coprime to a(b-a); new-prime-free iff E=+-1; box |a|<=300,b<=150 gives %s" % sorted(crit_sols))
# brute-force Thue F(b,d)=b^3-b d^2+d^3=+-1 (d=b-a), |b|,|d|<=3000
bf = sorted((bb, dd) for bb in range(-3000, 3001) for dd in range(-3000, 3001) if abs(bb ** 3 - bb * dd * dd + dd ** 3) == 1)
rep("RB13", bf == [(-4, 3), (-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1), (4, -3)], "brute force |b|,|d|<=3000: %s" % bf)
# certified solve with PARI/gp (thueinit flag 1 = unconditional)
gp_ok, gp_txt = False, ""
try:
    gp_txt = subprocess.run(["gp", "-q", "-f"], input="default(parisize,64000000);T=thueinit(x^3-x+1,1);print(thue(T,1));print(thue(T,-1));quit;\n", capture_output=True, text=True, timeout=120).stdout
    lines = [ln.strip() for ln in gp_txt.splitlines() if ln.strip().startswith("[")]
    gp_ok = lines == ["[[-1, 1], [0, 1], [1, 0], [1, 1], [4, -3]]", "[[-4, 3], [-1, -1], [-1, 0], [0, -1], [1, -1]]"]
except Exception as ex:  # gp missing
    gp_txt = "gp unavailable: %r" % ex
rep("RB13", gp_ok, "PARI thue (certified): all solutions of b^3-bd^2+d^3=+-1 are +-{(1,0),(0,1),(1,1),(-1,1),(4,-3)}; hence c in {0,-1 (degenerate), -2, -7/4}: OPEN QUESTION 1 CLOSED (field Q(plastic), disc -23)")
print("     gp raw:", " | ".join(gp_txt.split("\n")).strip())
# strong form including k=0
strong = sorted(set(Fraction(-aa, 2 ** k) for k in range(0, 13) for aa in range(-2 ** (k + 1), 3 * 2 ** k) if aa and gcd(abs(aa), 2 ** k) == 1 and abs((Fraction(-aa, 2 ** k) * (Fraction(-aa, 2 ** k) + 1) ** 2 + 1).numerator) == 1 and ((Fraction(-aa, 2 ** k) * (Fraction(-aa, 2 ** k) + 1) ** 2 + 1).denominator & ((Fraction(-aa, 2 ** k) * (Fraction(-aa, 2 ** k) + 1) ** 2 + 1).denominator - 1)) == 0))
rep("RB13", strong == [Fraction(-2), Fraction(-7, 4), Fraction(-1)], "strong form with k=0 admitted: %s (explorer's 'only -7/4' needs k>=1 / non-integral c)" % strong)

print("== RB14/RB15/RB16 ==")
earlier, nonew = [], []
for n in range(1, 201):
    t = 2 ** n - 1
    if n > 1 and new_prime_free(t, earlier):
        nonew.append(n)
    earlier.append(t)
rep("RB14", nonew == [6], "2^n-1 without primitive prime divisor, 2<=n<=200: %s (Bang 1886 / Zsigmondy 1892 exception (2,1,6) correctly cited)" % nonew)
ok = True
for bb in range(2, 40):
    for aa in range(-60, 61):
        if gcd(abs(aa), bb) != 1:
            continue
        o = orbit0(Fraction(aa, bb), 5)
        ok &= all(o[i].denominator == bb ** (2 ** i) for i in range(5))
rep("RB15", ok, "non-integral c=a/b: denominator of f^n(0) is exactly b^{2^{n-1}} (strictly increasing => never preperiodic)")
rep("RB15", all(orbit0(Fraction(cc), 3)[2] > orbit0(Fraction(cc), 3)[1] > cc for cc in range(1, 60)) and all(orbit0(Fraction(cc), 2)[1] >= -cc + 3 for cc in range(-60, -2)), "c>=1 increasing; c<=-3: f^2(0)>=|c|+3 then escape")
rep("RB15", all(x * x + cc > x for cc in range(-60, -2) for x in range(-cc + 1, -cc + 100)), "x>=|c|+1 => x^2+c > x")
pp = {}
for cc in (0, -1, -2):
    pp[cc] = [x0 for x0 in range(-50, 51) if (lambda x0: (lambda s, x: [s.add(x) or (x := x * x + cc) for _ in range(60)] and x in s)(set(), x0))(x0)]
rep("RB15", pp == {0: [-1, 0, 1], -1: [-1, 0, 1], -2: [-2, -1, 0, 1, 2]}, "integer preperiodic sets")
step = {}
for bb in (2, 5, 8):
    t = 2 * bb % 9
    step[RHO[bb]] = RHO[t if t in RHO else 9 - t]
rep("RB16", step == {1: 3, 3: 5, 5: 1} and {r: step[step[r]] for r in step} == {1: 5, 3: 1, 5: 3}, "angle-doubling on rows is a 3-cycle; its square is R (+4 mod 6)")
rep("T11", Fraction(49 - 29, 4) == 5 and Fraction(25 - 29, 4) == -1 and Fraction(1 - 29, 4) == -7, "(y^2-29)/4: -7 -> 5 -> -1 -> -7 (THM-4146 (33))")

print()
print("FAILURES:", FAILS if FAILS else "none")
