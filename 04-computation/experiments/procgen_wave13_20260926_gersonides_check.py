#!/usr/bin/env python3
"""Wave-13 orchestrator check: the known integer cycles of x -> x/2, (3x+1)/2 on Z
versus the Gersonides equation |3^a - 2^p| = 1, their densities as Stern-Brocot
approximants of log_3 2, and the identity 3 + 1 = 4 in three roles.

A parity word w of length p with a ones determines T^p(x) = (3^a x + c_w)/2^p on
its cylinder, so the periodic point of w is x_w = c_w/(2^p - 3^a) (Bohm-Sontacchi).
"""
import math
from fractions import Fraction
from itertools import combinations

C = math.log(2) / math.log(3)


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg)


def c_word(bits, q=3):
    c = 0
    for j, b in enumerate(bits):
        if b:
            c = q * c + 2 ** j
    return c


def orbit_of(x, q=3, b=1, cap=10 ** 4):
    seen = []
    while x not in seen:
        seen.append(x)
        x = (q * x + b) // 2 if x % 2 else x // 2
        if len(seen) > cap:
            return None
    i = seen.index(x)
    return seen[i:]


def canon(cyc):
    m = min(cyc)
    i = cyc.index(m)
    return tuple(cyc[i:] + cyc[:i])


print("1. Gersonides pairs |3^a - 2^p| = 1 (p <= 400)")
pairs = [(p, a) for p in range(1, 401) for a in range(0, p + 1) if abs(3 ** a - 2 ** p) == 1]
print("  ", pairs)
check(pairs == [(1, 0), (1, 1), (2, 1), (3, 2)], "only (p,a) = (1,0), (1,1), (2,1), (3,2): 2-1, 3-2, 4-3, 9-8 (Levi ben Gershon 1343)")

print("2. every word of a Gersonides shape gives an integer cycle")
free = set()
for p, a in pairs:
    D = 2 ** p - 3 ** a
    for ones in combinations(range(p), a):
        bits = [1 if j in ones else 0 for j in range(p)]
        c = c_word(bits)
        assert c % D == 0
        x = c // D
        cyc = orbit_of(x)
        # the cycle must realise the word
        y = x
        for j in range(p):
            assert (y % 2) == bits[j]
            y = (3 * y + 1) // 2 if y % 2 else y // 2
        assert y == x
        free.add(canon(cyc))
print("  ", sorted(free, key=lambda t: (len(t), t)))
check(free == {(0,), (1, 2), (-1,), (-10, -5, -7)}, "the free cycles are exactly {0}, {1,2}, {-1}, {-5,-7,-10}")

print("3. shape (p,a) = (11,7): 3^7 - 2^11 = 139")
D = 3 ** 7 - 2 ** 11
hits = []
for ones in combinations(range(11), 7):
    bits = [1 if j in ones else 0 for j in range(11)]
    c = c_word(bits)
    if c % (2 ** 11 - 3 ** 7) == 0:
        hits.append((bits, c // (2 ** 11 - 3 ** 7)))
cycs = {canon(orbit_of(x)) for _, x in hits}
print(f"   D = {D}; words with 139 | c_w: {len(hits)} of 330; cycles: {cycs}")
check(len(hits) == 11 and len(cycs) == 1 and -17 in next(iter(cycs)),
      "exactly one necklace (11 rotations) of shape (11,7) is integral: the -17 cycle (min -136), a sporadic cycle")

print("4. census of integer cycles of 3x+1 on Z for words of length p <= 20")
found = set()
for p in range(1, 21):
    for a in range(0, p + 1):
        Dp = 2 ** p - 3 ** a
        # enumerate necklace representatives by brute force over words with a ones
        # (2^20 words at p = 20)
        if Dp == 0:
            continue
        for ones in combinations(range(p), a):
            bits = [1 if j in ones else 0 for j in range(p)]
            c = c_word(bits)
            if c % Dp == 0:
                found.add(canon(orbit_of(c // Dp)))
print("   cycles:", sorted(found, key=lambda t: (len(t), t)))
check(found == free | {canon(orbit_of(-17))}, "p <= 20: exactly the five known cycles {0}, {1,2}, {-1}, {-5,-7,-10}, {-17,...}")

print("5. densities as Stern-Brocot approximants of log_3 2")
dens = {cyc: Fraction(sum(1 for v in cyc if v % 2), len(cyc)) for cyc in found}
for cyc, d in sorted(dens.items(), key=lambda t: t[1]):
    print(f"   {str(cyc[:4]):>22}...  density {d}  ({'expanding' if d > C else 'contracting'})")


def best_upper(fr):
    """fr > c and no fraction in (c, fr) has denominator <= fr.denominator."""
    if not fr > Fraction(C).limit_denominator(10 ** 12):
        return False
    for qd in range(1, fr.denominator + 1):
        for pn in range(0, qd + 1):
            f = Fraction(pn, qd)
            if C < f < fr:
                return False
    return True


def best_lower(fr):
    if not float(fr) < C:
        return False
    for qd in range(1, fr.denominator + 1):
        for pn in range(0, qd + 1):
            f = Fraction(pn, qd)
            if fr < f < C:
                return False
    return True


ups = [Fraction(x) for x in ("1/1", "2/3", "7/11")]
lows = [Fraction(0, 1), Fraction(1, 2)]
check(all(best_upper(f) for f in ups), "1/1, 2/3, 7/11 are best upper approximations of log_3 2")
check(all(best_lower(f) for f in lows), "0/1, 1/2 are best lower approximations of log_3 2")
seq_up = [f for f in (Fraction(pn, qd) for qd in range(1, 12) for pn in range(qd + 1)) if best_upper(f)]
seq_up = sorted(set(seq_up), key=lambda f: f.denominator)
check(seq_up == ups, "with denominators <= 11 the best upper approximations are exactly 1/1, 2/3, 7/11")
check(Fraction(7, 11) == Fraction(2 + 5, 3 + 8), "7/11 is the mediant of 2/3 (the -5 cycle) and 5/8 (sigma_k's critical Christoffel cycle, 319/13)")

print("6. the identity 3 + 1 = 4 in three roles")
for q in range(1, 40, 2):
    # (a) pairing ladder: T(2i-1) + T(2i) = (2i-1) + 2i for all i
    ok_pair = all((q * (2 * i - 1) + 1) // 2 + i == 4 * i - 1 for i in range(1, 50))
    # (b) second moment g_q(2) = (1 + q)/4 = 1
    ok_mom = Fraction(1 + q, 4) == 1
    # (c) trivial cycle 1 -> 2 -> 1: (q*1 + 1)/2 = 2
    ok_cyc = (q + 1) // 2 == 2 and (q + 1) % 2 == 0
    assert ok_pair == ok_mom == ok_cyc == (q == 3), q
check(True, "for odd q <= 39: pairing-sum preservation (THM-4470) <=> g_q(2) = 1 (THM-4477) <=> 1 -> 2 -> 1 is a cycle <=> q = 3")
# integer roots of g_3(s) = (1 + 3^(s-1))/2^s = 1 are s = 1, 2 <=> 2^s - 3^(s-1) = 1 (words 1^(s-1) 0)
roots = [s for s in range(1, 200) if 1 + 3 ** (s - 1) == 2 ** s]
check(roots == [1, 2], "integer roots of g_3(s) = 1 are s = 1, 2, the Gersonides pairs (2,1), (4,3), i.e. the cycles {0} and {1,2} (words 0 and 10)")

print("7. other multipliers: free cycles need |2^p - q^a| = 1")
for q in range(5, 102, 2):
    sols = [(p, a) for p in range(1, 400) for a in range(0, 60) if abs(q ** a - 2 ** p) == 1]
    assert all(a <= 1 for p, a in sols), (q, sols)
check(True, "for odd 5 <= q <= 101 every solution of |q^a - 2^p| = 1 has a <= 1 (p < 400, a < 60; Mihailescu for all): "
      "q = 3 is the only multiplier with a free cycle having two odd steps (the -5 cycle)")
sol5 = [(p, a) for p in range(1, 400) for a in range(0, 60) if abs(5 ** a - 2 ** p) == 1]
check(sol5 == [(1, 0), (2, 1)], "5x+1: free cycles only {0} and {-1,-2} (5 - 4 = 1)")
