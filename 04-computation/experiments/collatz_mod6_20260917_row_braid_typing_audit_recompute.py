#!/usr/bin/env python3
"""Independent recomputation for the row_braid_typing lane (audit lens: recompute).

Does NOT import the explorer's script.  Every claim RB1..RB17 is recomputed
from scratch with its own code paths, plus hostile / boundary probes the
explorer did not run (negative n in the diagonal search, integer c in the
strong -7/4 search, larger weak search, Wieferich braid r, R^3 mod 84/126).
Exact integer / Fraction arithmetic only.  RAM tiny, runtime < 1 min.
"""
import json
from fractions import Fraction
from math import gcd

FAIL = []


def rep(tag, ok, msg=""):
    print("%s %-6s %s" % ("OK  " if ok else "FAIL", tag, msg))
    if not ok:
        FAIL.append((tag, msg))


def order(a, m):
    assert gcd(a, m) == 1
    x, k = a % m, 1
    while x != 1:
        x = x * a % m
        k += 1
    return k


def vp(n, p):
    n = abs(n)
    k = 0
    while n and n % p == 0:
        n //= p
        k += 1
    return k


def ch(m):  # core, height with sign on core
    h = vp(m, 2)
    return m // 2 ** h, h


def F(n):
    assert n % 2
    return (3 * n + 1) // 2


B = {1: 2, 3: 5, 5: 8}
RHO = {2: 1, 5: 3, 8: 5}

# ---------------- RB1 ----------------
print("== RB1 orders ==")
rep("RB1", order(2, 9) == 6 and order(2, 3) == 2 and order(2, 7) == 3 and order(2, 63) == 6,
    "ord_9=6 ord_3=2 ord_7=3 ord_63=6")
first9 = min(n for n in range(1, 100) if (2 ** n - 1) % 9 == 0)
rep("RB1", first9 == 6, "first Mersenne divisible by 9 is n=%d" % first9)
sub63 = sorted(set(pow(2, k, 63) for k in range(6)))
rep("RB1", sub63 == [1, 2, 4, 8, 16, 32], "<2> mod 63 = %s" % sub63)
rep("RB1", sorted(x % 9 for x in sub63) == [1, 2, 4, 5, 7, 8], "projects onto (Z/9)^* bijectively")
rep("RB1", sorted(x % 7 for x in sub63) == [1, 1, 2, 2, 4, 4], "2:1 onto {1,2,4} mod 7")
rep("RB1", len(set(pow(2, k, 9) for k in range(6))) == 6, "2 primitive root mod 9")

# ---------------- RB2 ----------------
print("== RB2 rows ==")
ok = True
for n in range(-200001, 200002, 2):
    j, r = divmod(n, 6)
    if F(n) != 9 * j + B[r] or F(n) % 3 != 2:
        ok = False
        break
rep("RB2", ok, "F(6j+r)=9j+b_r and F=2 mod 3 for odd |n|<=200001 (both signs)")
# 'uses only the multiplier': 3n+1 mod 18 is a function of n mod 6 -- trivially true
# general p: image of odd n under (pn+1)/2 is exactly {m : 2m = 1 mod p}; sources mod 2p <-> targets mod p^2
for p in (3, 5, 7, 11, 13):
    img = set(((p * n + 1) // 2) % (p * p) for n in range(1, 2 * p * p, 2))
    inv2 = pow(2, -1, p)
    rep("RB2", len(img) == p and all(m % p == inv2 for m in img),
        "p=%d: p target classes mod p^2, all = 2^-1 mod p" % p)

# ---------------- RB3 ----------------
print("== RB3 powers of two ==")
seen = {}
for k in range(1, 301):
    m = 2 ** k
    if m % 3 != 2:
        rep("RB3", k % 2 == 0, "2^%d not an image, k even" % k) if k % 2 else None
        continue
    src = (2 * m - 1) // 3
    if src % 2 == 0:
        rep("RB3", False, "source even?!")
    seen.setdefault(k % 6, set()).add((src % 6, m % 9, m % 7))
rep("RB3", seen == {1: {(1, 2, 2)}, 5: {(3, 5, 4)}, 3: {(5, 8, 1)}},
    "k=1,5,3 mod 6 -> rows 1,3,5; mod 7 residues 2,4,1 bijective on odd k: %s" % seen)

# ---------------- RB4 ----------------
print("== RB4 R braid ==")
R = lambda n: 4 * n + 1
ok = True
for n in range(-5001, 5002, 2):
    n3 = R(R(R(n)))
    if n3 != 64 * n + 21 or n3 - n != 21 * (3 * n + 1):
        ok = False
    if (n3 - n) % 42 != 0:
        ok = False
    if (F(n3) - 64 * F(n)) != 0 or (F(n3) - F(n)) % 63 != 0:
        ok = False
    if (R(n) - n) % 6 == 0 or (R(R(n)) - n) % 6 == 0:
        ok = False
    if (R(n) - n - 4) % 6 != 0:
        ok = False
    # period on rows mod 18 is 9, not 3
    x, t = n, 0
    while True:
        x = R(x); t += 1
        if (x - n) % 18 == 0:
            break
    if t != 9:
        ok = False
    for tt in range(1, 10):
        x = n
        for _ in range(tt):
            x = R(x)
        if vp(4 ** tt - 1, 3) != 1 + vp(tt, 3) or vp(x - n, 3) != vp(tt, 3):
            ok = False
rep("RB4", ok, "R^3=64n+21, R^3-n=21(3n+1), mod 42 fixed, target x64 (mod 63 fixed), row period 3, mod-18 period 9, v_3 law; odd |n|<=5001")
# hostile: is 42 the exact modulus?  R^3(n)-n = 21(3n+1) is divisible by 84 iff n = 1 mod 4, never by 126
n84 = sum(1 for n in range(1, 1000, 2) if (R(R(R(n))) - n) % 84 == 0)
n126 = sum(1 for n in range(1, 1000, 2) if (R(R(R(n))) - n) % 126 == 0)
rep("RB4", n84 == 250 and n126 == 0, "R^3 fixes n mod 84 for exactly n=1 mod 4 (%d/500), never mod 126 (%d): 42 is the exact universal modulus" % (n84, n126))
orbs = set()
for a in range(7):
    cyc, x = [a], (4 * a + 1) % 7
    while x != a:
        cyc.append(x); x = (4 * x + 1) % 7
    orbs.add(tuple(sorted(cyc)))
rep("RB4", orbs == {(2,), (0, 1, 5), (3, 4, 6)}, "R mod 7 orbits %s" % sorted(orbs))
rep("RB4", all((3 * n + 1) % 7 == 0 for n in range(2, 700, 7)) and sum(1 for n in range(1, 700, 2) if (3 * n + 1) % 7 == 0) == 50,
    "7 | F(n) iff n = 2 mod 7")
tg = sorted(set(F(n) % 63 for n in range(1, 127, 2)))
p2 = sorted(set(pow(2, k, 63) for k in range(1, 7)) & set(tg))
rep("RB4", len(tg) == 21 and p2 == [2, 8, 32] and sorted(x % 9 for x in p2) == [2, 5, 8] and sorted(x % 7 for x in p2) == [1, 2, 4],
    "21 target classes mod 63; powers of two in %s" % p2)
rep("RB4", order(4, 9) == 3 and order(4, 7) == 3, "ord_9(4)=ord_7(4)=3")

# ---------------- RB5 ----------------
print("== RB5 Wieferich ==")
for p in (3, 5, 7, 11, 13, 17, 1093, 3511):
    o1, o2 = order(2, p), order(2, p * p)
    inv2 = pow(2, -1, p)
    # rows hit: classes 2^k mod p^2 over k with 2^k = 2^{-1} mod p
    rows = set()
    k0 = next(k for k in range(1, o1 + 1) if pow(2, k, p) == inv2)
    for k in range(k0, k0 + o2, o1):
        rows.add(pow(2, k, p * p))
    wief = pow(2, p - 1, p * p) == 1
    rep("RB5", len(rows) == o2 // o1 and (len(rows) == (1 if wief else p)),
        "p=%d ord_p=%d ord_p2=%d rows=%d wieferich=%s" % (p, o1, o2, len(rows), wief))
    if wief:
        rep("RB5", rows == {pow(2, k0, p * p)}, "p=%d single class %d mod p^2" % (p, pow(2, k0, p * p)))

# ---------------- RB6 ----------------
print("== RB6 tower ==")
Fp = lambda p, n: (p * n + 1) // 2
ok = all(Fp(p + 2, n) == Fp(p, n) + n and Fp(p, n) == Fp(1, n) + (p - 1) // 2 * n
         for p in range(-9, 40, 2) for n in range(-999, 1000, 2))
rep("RB6", ok, "F_{p+2}=F_p+n, F_p=F_1+((p-1)/2)n, odd p in [-9,39], odd |n|<1000")
diag_pos = [(p, n) for p in range(1, 60, 2) for n in range(1, 600, 2) if Fp(p - 2, n) == n]
diag_all = [(p, n) for p in range(-59, 60, 2) for n in range(-599, 600, 2) if Fp(p - 2, n) == n]
rep("RB6", diag_pos == [(3, 1)], "positive n: diagonal %s" % diag_pos)
rep("RB6", sorted(diag_all) == [(3, 1), (5, -1)], "ALL odd p,n (both signs): diagonal %s  <- (5,-1) is a second solution of (p-4)n=-1" % sorted(diag_all))

# ---------------- RB7 ----------------
print("== RB7 inverse braids ==")


def least_q(p):
    d = 1
    while (2 ** d - 1) % p:
        d += 1
    return 2 ** d, d


for p in (1, 3, 5, 7, 11, 13):
    q, d = least_q(p)
    c = (q - 1) // p
    ok = all(Fp(p, q * n + c) == q * Fp(p, n) for n in range(-999, 1000, 2))
    pers = []
    for s in (1, 2, 3):
        M = 2 * p ** s
        # count orbits of n->qn+c on odd classes mod M
        seen_ = set(); orbits = 0; lens = set()
        for a in range(1, M, 2):
            if a in seen_:
                continue
            orbits += 1; x = a; L = 0
            while x not in seen_:
                seen_.add(x); x = (q * x + c) % M; L += 1
            lens.add(L)
        pers.append((orbits, tuple(sorted(lens))))
    r = vp(q - 1, p) if p > 1 else None
    exp = [(1, (p ** s,)) for s in (1, 2, 3)] if p > 1 else [(1, (1,))] * 3
    rep("RB7", ok and pers == exp and (p == 1 or r == 1), "p=%d q=%d c=%d r=%s d=%d orbits/lengths mod 2p^s: %s" % (p, q, c, r, d, pers))
    if p == 1:
        rep("RB7", all(ch(Fp(1, 2 ** (h + 1) * u - 1)) == (u, h) for u in range(1, 400, 2) for h in range(0, 10)),
            "F_1 fibre n=2^{h+1}u-1 every h; R_1 mod 6 orbits " + str({a: (2 * a + 1) % 6 for a in (1, 3, 5)}))
# hostile: Wieferich p: r >= 2, period drops
for p in (1093,):
    q, d = least_q(p)
    r = vp(q - 1, p)
    rep("RB7", r == 2, "p=1093: least q=2^%d has r=v_p(q-1)=%d (>1, fullness fails per inherited (11))" % (d, r))

# ---------------- RB8 ----------------
print("== RB8 cycle census ==")


def Tp(p, n):
    m = p * n + 1
    return m // 2 ** vp(m, 2)


for p in (1, 3, 5, 7):
    cycles, esc = set(), 0
    for n0 in range(1, 20001, 2):
        x, t, seen_ = n0, 0, set()
        while t < 3000 and x < 10 ** 40:
            if x in seen_:
                cyc, y = [], x
                while True:
                    cyc.append(y); y = Tp(p, y)
                    if y == x:
                        break
                cycles.add(tuple(sorted(cyc)))
                break
            seen_.add(x); x = Tp(p, x); t += 1
        else:
            esc += 1
    print("   p=%d cycles(sorted members)=%s escaped=%d" % (p, sorted(cycles), esc))
    want = {1: ({(1,)}, 0), 3: ({(1,)}, 0), 5: ({(1, 3), (13, 33, 83), (17, 27, 43)}, 9605), 7: ({(1,)}, 9982)}[p]
    rep("RB8", cycles == want[0] and esc == want[1], "p=%d matches explorer (esc %d)" % (p, esc))
rep("RB8", all(Tp(1, n) < n for n in range(3, 100001, 2)) and Tp(1, 1) == 1, "T_1 descent")

# ---------------- RB9 / RB10 ----------------
print("== RB9 row law, RB10 layers ==")
ok = True; cnt = 0
for u in list(range(-3001, 3002, 2)):
    if u % 3 == 0:
        continue
    for h in range(0, 13):
        num = 2 ** (h + 1) * u - 1
        if num % 3:
            if (h % 2 == 1) != (u % 3 == 1):
                ok = False
            continue
        n = num // 3
        r = n % 6
        if RHO[(2 ** h * u) % 9] != r or ch(F(n)) != (u, h) or (2 ** h * u - B[r]) % 9 or 6 * ((2 ** h * u - B[r]) // 9) + r != n:
            ok = False
        cnt += 1
rep("RB9", ok, "row law + index law + round trip, both signs, %d (u,h) pairs" % cnt)
rep("RB9", all(RHO[(4 * b) % 9] == (RHO[b] + 4) % 6 for b in (2, 5, 8)), "row(u,h+2)=row(u,h)+4 mod 6")
base = {}
for um in (1, 5, 7, 11, 13, 17):
    h0 = 1 if um % 3 == 1 else 0
    base[um] = (h0, ((2 ** (h0 + 1) * um - 1) // 3) % 6)
rep("RB9", base == {1: (1, 1), 7: (1, 3), 13: (1, 5), 5: (0, 3), 11: (0, 1), 17: (0, 5)}, "base table %s" % base)
ok = True
for r in (1, 3, 5):
    for h in range(0, 7):
        want = B[r] * pow(2, -h, 9) % 9
        js = [j for j in range(2 ** 12) if vp(9 * j + B[r], 2) == h]
        us = [(9 * j + B[r]) >> h for j in js]
        if not all(u % 9 == want and u % 2 for u in us):
            ok = False
        if any(js[i + 1] - js[i] != 2 ** (h + 1) for i in range(len(js) - 1)):
            ok = False
        if any(us[i + 1] - us[i] != 18 for i in range(len(us) - 1)):
            ok = False
        # first index formula (explorer): j0 = (2^h * u_min - b_r)/9 with u_min least odd rep of class
        umin = want if want % 2 else want + 9
        if js[0] != (2 ** h * umin - B[r]) // 9:
            ok = False
rep("RB10", ok, "AP layers mod 18, index step 2^{h+1}, first index, h<=6, j<2^12")
ok = True
for r in (1, 3, 5):
    for u in range(1, 400, 2):
        if u % 3 == 0:
            continue
        h = 0
        while not ((2 ** (h + 1) * u - 1) % 3 == 0 and ((2 ** (h + 1) * u - 1) // 3) % 6 == r):
            h += 1
        j0 = (2 ** h * u - B[r]) // 9
        j1 = (2 ** (h + 6) * u - B[r]) // 9
        if j1 != 64 * j0 + 7 * B[r] or ch(F(6 * j1 + r)) != (u, h + 6):
            ok = False
rep("RB10", ok and [7 * B[r] for r in (1, 3, 5)] == [14, 35, 56], "j_{t+1}=64 j_t + 7 b_r (14,35,56)")

# ---------------- RB11 ----------------
print("== RB11 rows vs JSON, prominence ==")
rows = {r: [ch(F(6 * j + r)) for j in range(61)] for r in (1, 3, 5)}
inh = json.load(open("/tmp/math-wt-collatz-mod6/05-knowledge/results/arithmetic_braids_20260917_collatz.json"))["rows"]
rep("RB11", all(rows[r][:35] == [tuple(x) for x in inh[str(r)]] and len(inh[str(r)]) == 35 for r in (1, 3, 5)), "rows agree with JSON on j<=34")
first = {}
for r in (1, 3, 5):
    for j, (u, h) in enumerate(rows[r]):
        first.setdefault(u, {}).setdefault(r, (j, h))
rep("RB11", first[1] == {1: (0, 1), 3: (3, 5), 5: (0, 3)} and first[5] == {1: (2, 2), 3: (0, 0), 5: (8, 4)}
    and first[7] == {1: (6, 3), 3: (1, 1), 5: (24, 5)} and first[11] == {1: (1, 0), 3: (19, 4), 5: (4, 2)}
    and first[13] == {1: (46, 5), 3: (11, 3), 5: (2, 1)}, "first appearances of 1,5,7,11,13")
cnt = {}
for r in (1, 3, 5):
    for (u, h) in rows[r]:
        cnt[u] = cnt.get(u, 0) + 1
top = sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0]))[:4]
h1 = sorted(h for r in (1, 3, 5) for (u, h) in rows[r] if u == 1)
rep("RB11", top[0] == (1, 5) and top[1] == (5, 4) and h1 == [1, 3, 5, 7, 9], "top counts %s; core-1 heights %s" % (top, h1))
# prediction formula: appearances of u within j<=J = #{admissible h : 2^h u <= 9J + b_{row(u,h)}}
J = 60
ok = True
for u in range(1, 200, 2):
    if u % 3 == 0:
        continue
    pred = 0
    for h in range(0, 30):
        if (2 ** (h + 1) * u - 1) % 3:
            continue
        r = RHO[(2 ** h * u) % 9]
        if 2 ** h * u <= 9 * J + B[r]:
            pred += 1
    if pred != cnt.get(u, 0):
        ok = False
rep("RB11", ok, "appearance-count formula exact for u<200, J=60")
least = {c: min(u for u in range(1, 200, 2) if u % 18 == c) for c in (1, 5, 7, 11, 13, 17)}
rep("RB11", sorted(least.values()) == [1, 5, 7, 11, 13, 17], "least cores of the six classes mod 18")

# ---------------- RB12 ----------------
print("== RB12 negative rows ==")
Fm = lambda m: (3 * m - 1) // 2
ok = True
for j in range(-2000, 0):
    for r in (1, 3, 5):
        n = 6 * j + r; m = -n
        if m % 6 != {1: 5, 3: 3, 5: 1}[r] or F(n) != -Fm(m):
            ok = False
        u, h = ch(F(n)); um, hm = ch(Fm(m))
        if (u, h) != (-um, hm):
            ok = False
        if (-n) != 6 * (-j - 1) + (6 - r):
            ok = False
rep("RB12", ok, "negation law, F(-m)=-F_-(m), signed cores, j->-j-1; j in [-2000,-1]")
nrows = {r: {j: ch(F(6 * j + r)) for j in range(-40, 0)} for r in (1, 3, 5)}
rep("RB12", [nrows[5][j] for j in (-1, -2, -3, -4)] == [(-1, 0), (-5, 1), (-19, 0), (-7, 2)] and
    [nrows[1][j] for j in (-1, -2, -3, -4)] == [(-7, 0), (-1, 4), (-25, 0), (-17, 1)], "printed negative row heads")
neg = [(-1,), (-5, -7), (-17, -25, -37, -55, -41, -61, -91)]
for cyc in neg:
    L = len(cyc)
    rep("RB12", all(ch(F(cyc[i]))[0] == cyc[(i + 1) % L] for i in range(L)), "cycle %s closes under T (ordered)" % (cyc,))
want = {-1: [(1, -2, 4), (3, -29, 8), (3, -1, 2), (5, -8, 6), (5, -1, 0)], -5: [(1, -18, 5), (3, -5, 3), (5, -2, 1)], -7: [(1, -1, 0), (3, -13, 4), (5, -4, 2)],
        -17: [(1, -4, 1), (5, -16, 3)], -25: [(1, -3, 0), (5, -12, 2)], -37: [(3, -17, 2), (5, -5, 0)], -55: [(3, -25, 2), (5, -7, 0)],
        -41: [(3, -37, 3), (5, -10, 1)], -61: [(1, -7, 0), (5, -28, 2)], -91: [(5, -11, 0)]}
ok = True
for u, w in want.items():
    hits = sorted((r, j, h) for r in (1, 3, 5) for j, (uu, h) in nrows[r].items() if uu == u)
    if hits != sorted(w):
        ok = False; print("   mismatch", u, hits, w)
rep("RB12", ok, "cycle-member placements |j|<=40 match")
rep("RB12", all(u % 6 in (1, 5) for cyc in neg for u in cyc) and 1 % 6 == 1, "no cycle member 3 mod 6")
# row law for negative u already covered in RB9 (u range includes negatives)

# ---------------- RB13 / RB14 ----------------
print("== RB13 -7/4 identity, RB14 Bang ==")
c = Fraction(-7, 4)
rep("RB13", c * (c + 1) ** 2 == Fraction(-63, 64) and c * (c * (c + 1) ** 2 + 1) == Fraction(-7, 256), "c(c+1)^2=-63/64, f^3(0)=-7/256")


def orbit(c, N):
    x, out = Fraction(0), []
    for _ in range(N):
        x = x * x + c; out.append(x)
    return out


def nonew(num, earlier):
    g = abs(num)
    if g <= 1:
        return True
    N = 1
    for e in earlier:
        N *= abs(e) or 1
    while True:
        d = gcd(g, N)
        if d == 1:
            break
        g //= d
    return g == 1


o = orbit(c, 4)
rep("RB13", o[:3] == [Fraction(-7, 4), Fraction(21, 16), Fraction(-7, 256)] and not nonew(o[1].numerator, [o[0].numerator]) and nonew(o[2].numerator, [o[0].numerator, o[1].numerator]),
    "orbit %s; term 3 no new prime, term 2 has new prime 3" % o[:4])
rep("RB13", o[3] == Fraction(-114639, 65536) and 114639 == 3 * 7 * 53 * 103, "term 4 = -3*7*53*103/2^16")
o2 = orbit(Fraction(-29, 16), 3)
rep("RB13", o2[2] == Fraction(23345, 65536) and 23345 == 5 * 7 * 23 * 29, "-29/16 third term brings 5,7,23")
# strong search: c = -a/2^k INCLUDING k=0 (integers) -- explorer started at k=1
strong = []
for k in range(0, 13):
    b = 2 ** k
    for a in range(-4 * b, 4 * b + 1):
        if a == 0 or gcd(abs(a), b) != 1:
            continue
        cc = Fraction(-a, b); v = cc * (cc + 1) ** 2 + 1
        if abs(v.numerator) == 1 and v.denominator & (v.denominator - 1) == 0:
            strong.append(cc)
rep("RB13", sorted(strong) == [Fraction(-2), Fraction(-7, 4), Fraction(-1), Fraction(0)],
    "strong form c(c+1)^2+1=+-2^-m over c=-a/2^k, k<=12 INCLUDING k=0: %s  <- integers 0,-1,-2 are solutions too" % sorted(strong))
# weak search, larger: c=a/b |a|<=300, b<=150
weak = []
for b in range(1, 151):
    for a in range(-300, 301):
        if a == 0 or gcd(abs(a), b) != 1:
            continue
        cc = Fraction(a, b)
        x1 = cc; x2 = cc * cc + cc
        if x2 == 0:
            continue
        x3 = x2 * x2 + cc
        if nonew(x3.numerator, [x1.numerator, x2.numerator]):
            weak.append(cc)
rep("RB13", sorted(weak) == [Fraction(-2), Fraction(-7, 4)], "weak (third term no new prime) |a|<=300, b<=150: %s" % sorted(weak))
# structural reason: third term numerator = a * (a(a+b)^2 + b^3) ; it has no new prime iff a(a+b)^2 + b^3 is (up to sign) a product of primes dividing a(a+b)
# hostile: c = a/b in lowest terms with f^3(0) numerator = a*(a*(a+b)^2 + b^3); check formula
ok = all((lambda cc: (orbit(cc, 3)[2]).numerator * 1 == (Fraction(cc.numerator * (cc.numerator * (cc.numerator + cc.denominator) ** 2 + cc.denominator ** 3), cc.denominator ** 4)).numerator)(Fraction(a, b))
         for a in range(-30, 31) for b in range(1, 20) if a and gcd(abs(a), b) == 1)
rep("RB13", ok, "f^3(0) = a(a(a+b)^2+b^3)/b^4 with c=a/b (numerator formula)")
earlier, nn = [], []
for n in range(1, 80):
    t = 2 ** n - 1
    if n > 1 and nonew(t, earlier):
        nn.append(n)
    earlier.append(t)
rep("RB14", nn == [6], "2^n-1 no primitive prime divisor only at n=6 for 2<=n<80")

# ---------------- RB15 ----------------
print("== RB15 PCF ==")
pcf = []
for cint in range(-200, 201):
    x, s = 0, set()
    per = False
    for _ in range(300):
        if x in s:
            per = True; break
        s.add(x); x = x * x + cint
        if abs(x) > 10 ** 15:
            break
    if per:
        pcf.append(cint)
rep("RB15", pcf == [-2, -1, 0], "integer PCF %s" % pcf)
# denominators: c=a/b, b>1 -> denominator of f^n(0) is exactly b^(2^(n-1))
ok = True
for b in range(2, 30):
    for a in range(-40, 41):
        if a == 0 or gcd(abs(a), b) != 1:
            continue
        o = orbit(Fraction(a, b), 6)
        if any(o[i].denominator != b ** (2 ** i) for i in range(6)):
            ok = False
rep("RB15", ok, "denominator of f^n(0) is b^(2^(n-1)) exactly (so never preperiodic), b<30, |a|<=40")
for cint, want in ((0, [-1, 0, 1]), (-1, [-1, 0, 1]), (-2, [-2, -1, 0, 1, 2])):
    pp = []
    for x0 in range(-200, 201):
        x, s, okk = x0, set(), False
        for _ in range(200):
            if x in s:
                okk = True; break
            s.add(x); x = x * x + cint
            if abs(x) > 10 ** 12:
                break
        if okk:
            pp.append(x0)
    rep("RB15", pp == want, "integer preperiodic of x^2%+d: %s" % (cint, pp))
# hostile: x^2-6, x^2-12: integer fixed points but 0 escapes
for cint in (-6, -12):
    fp = [x for x in range(-20, 21) if x * x + cint == x]
    x = 0
    for _ in range(6):
        x = x * x + cint
    rep("RB15", len(fp) == 2 and abs(x) > 10 ** 6, "x^2%+d fixed pts %s, 0 escapes" % (cint, fp))

# ---------------- RB16 ----------------
print("== RB16 ==")
step = {}
for b in (2, 5, 8):
    t = 2 * b % 9
    t = t if t in (2, 5, 8) else 9 - t
    step[RHO[b]] = RHO[t]
sq = {r: step[step[r]] for r in (1, 3, 5)}
rep("RB16", step == {1: 3, 3: 5, 5: 1} and sq == {1: 5, 3: 1, 5: 3} and sq == {r: (r + 4) % 6 for r in (1, 3, 5)}, "step %s, squared %s = R on rows" % (step, sq))
# angle-doubling on 2cos(2 pi b/9): y->y^2-2 sends 2cos(t)->2cos(2t), so b->2b mod +-9: {2,5,8} = {2,-4,-1} -> classes {1,2,4} mod +-9
rep("RB16", sorted((min(b, 9 - b)) for b in (2, 5, 8)) == [1, 2, 4], "rows = classes {1,2,4} mod +-9")

# ---------------- THM-4146 (33) recheck ----------------
G = lambda y: Fraction(y * y - 29, 4)
rep("T11", G(-7) == 5 and G(5) == -1 and G(-1) == -7, "(y^2-29)/4 cycle -7->5->-1")

print()
print("FAILURES:", FAIL if FAIL else "none")
