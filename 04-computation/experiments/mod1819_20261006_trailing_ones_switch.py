#!/usr/bin/env python3
"""Uniform (run-length-independent) switches of the rewrite compiler are collisions at -1 (mac-mini, 2026-10-06).

Setting: U(x) = oddpart(3x+1) on odd x; a word w = (a_1..a_j) of exponents; f_c(x) = (3x+1)/2^c on Q, f_u the
composition along u.  A word is REDUCED if its first letter is >= 2.  Every odd n > 1 is n = 2^(r+1) t - 1 (t odd),
and its word is 1^r u with u reduced (r = number of leading exponents 1).

  1. Endpoint progression: the endpoints z = U^j(n) of the sources with word 1^r u are exactly z = f_u(-1)
     (mod 3^(r+|u|)) (positive odd z with integral inverse).  [PROVED: U_w is affine of slope 3^j/2^A and
     U_w(-1) = f_w(-1) = f_u(-1) because f_1 fixes -1.]
  2. Denominators: for reduced u, f_u(-1) = N/2^(A_u - 1) with N odd.  So colliding reduced words (f_u(-1) =
     f_u'(-1)) have equal totals A_u = A_u'.
  3. Switch: if f_u(-1) = f_u'(-1) with |u'| = |u| + D, D >= 1, then every n with word prefix 1^r u, r >= D,
     satisfies U^j(n) = U^j(m), j = r + |u|, with m = (n+1)/2^D - 1 (delete D trailing ones of n; m's word is
     1^(r-D) u').  The reset switch of checked_switch_phase19 is the collision (a) ~ (2, a-2), D = 1, m = (n-1)/2.
  4. Longer partners: a collision partner of length k > j fires only on n = -1 (mod 3^(k-j)), where n already has
     the smaller predecessor (2n-1)/3.  For reset-2 sources (word 1^r (2, c, ...)) the root collision
     (2, c, .) ~ (c+2, .) points the wrong way; shifted by 3 it fires exactly on n = 8 (mod 9).
  5. Sporadic collisions exist, e.g. (8, c) ~ (4, 1, 1, c+2) (both give 125/2^(7+c); 125 = 2^7 - 3 = 2^5 + 3*2^4 +
     9*2^3 - 27), hence (2, 6, c) ~ (4, 1, 1, c+2): reset-2 sources (1^r, 2, 6, c) switch to (n-1)/2 at depth r+3.
  6. Census: reset-2 sources n < 2*10^4: uniform trailing-ones switches for 377 of 2500; every equal-length merge
     of n with some (n+1)/2^D - 1 is a collision.  The 239 residual seeds: 16 (10 multiples of 3).  Odd Mersenne
     numbers 2^a - 1 (a reset-2 family): 37 of the 60 odd a in [3, 121] switch to 2^(a-D) - 1, always with D odd
     (an even Mersenne number, which the reset switch takes one step further).
Run: python3 mod1819_20261006_trailing_ones_switch.py   (about 1 min)
"""
import itertools, json
from fractions import Fraction
from pathlib import Path
from collections import Counter, defaultdict

OK = True


def check(cond, msg):
    global OK
    print(("  ok   " if cond else "  FAIL ") + msg, flush=True)
    OK &= bool(cond)


def U(x):
    y = 3 * x + 1
    e = (y & -y).bit_length() - 1
    return y >> e, e


def orbit(n, L):
    xs, w, x = [n], [], n
    for _ in range(L):
        if x == 1:
            break
        x, e = U(x)
        xs.append(x)
        w.append(e)
    return xs, w


def fq(u):
    x = Fraction(-1)
    for c in u:
        x = (3 * x + 1) / 2 ** c
    return x


def strip1(w):
    r = 0
    while r < len(w) and w[r] == 1:
        r += 1
    return r, tuple(w[r:])


def AB(w):
    A, B, S = sum(w), 0, 0
    for i in range(len(w)):
        B += 3 ** (len(w) - 1 - i) * 2 ** S
        S += w[i]
    return A, B


def inv(w, z):
    """I_w(z) = (2^A z - B)/3^j if it is a positive odd integer whose word is w, else None"""
    A, B = AB(w)
    num = 2 ** A * z - B
    if num <= 0 or num % 3 ** len(w):
        return None
    m = num // 3 ** len(w)
    xs, ww = orbit(m, len(w))
    return m if (ww == list(w) and xs[-1] == z) else None


# ---------------------------------------------------------------- 1. endpoint progression
print("1. the endpoints of the sources with word 1^r u are z = f_u(-1) (mod 3^(r+|u|))")
good = True
for L in range(1, 5):
    for u in itertools.product(range(1, 6), repeat=L):
        if u[0] == 1:
            continue
        for r in range(0, 4):
            w = (1,) * r + u
            A, B = AB(w)
            j = len(w)
            tau = B * pow(2, -A, 3 ** j) % 3 ** j
            q = fq(u)
            good &= tau == q.numerator * pow(q.denominator, -1, 3 ** j) % 3 ** j
check(good, "tau(1^r u) = B 2^-A = f_u(-1) mod 3^(r+|u|) for every reduced u of length <= 4 (letters <= 5), r <= 3")

# ---------------------------------------------------------------- 2. denominators and the collision census
print("2. f_u(-1) = N/2^(A_u - 1), N odd, for reduced u; collisions have equal totals")
good = True
val = defaultdict(list)
for L in range(1, 7):
    for u in itertools.product(range(1, 9), repeat=L):
        if u[0] == 1:
            continue
        q = fq(u)
        good &= q.denominator == 2 ** (sum(u) - 1) and q.numerator % 2 == 1
        val[q].append(u)
check(good, "every reduced word of length <= 6 with letters <= 8: denominator exactly 2^(A_u - 1), odd numerator")
coll = {q: v for q, v in val.items() if len(v) > 1}
check(all(len({sum(u) for u in v}) == 1 for v in coll.values()), f"{len(coll)} colliding values: all colliding words share their total")


def reset_canon(u):
    # canonical form under the root collision (2, c, rest) ~ (c+2, rest)
    return (u[1] + 2,) + u[2:] if len(u) >= 2 and u[0] == 2 else u


spor = {q: v for q, v in coll.items() if len({reset_canon(u) for u in v}) > 1}
check(fq((8, 1)) == fq((4, 1, 1, 3)) == Fraction(125, 2 ** 8) and fq((2, 6, 1)) == fq((8, 1)),
      "sporadic collision (8, c) ~ (4, 1, 1, c+2): 125/2^(7+c), 125 = 2^7 - 3 = 2^5 + 3*2^4 + 9*2^3 - 27; and (2, 6, c) ~ (8, c)")
print(f"   words of length <= 6, letters <= 8: {len(coll)} colliding values, {len(spor)} not explained by the root collision alone")

# ---------------------------------------------------------------- 3. the switch m = (n+1)/2^D - 1
print("3. a collision u ~ u' with |u'| = |u| + D gives U^j(n) = U^j((n+1)/2^D - 1) on the whole family 1^r u, r >= D")
fams = [((3,), (2, 1), 1), ((5,), (2, 3), 1), ((2, 6, 1), (4, 1, 1, 3), 1), ((2, 6, 2), (2, 2, 1, 1, 4), 2), ((8, 2), (4, 1, 1, 4), 2)]
good, tested = True, 0
for u, up, D in fams:
    assert fq(u) == fq(up) and len(up) - len(u) == D
    for r in range(D, D + 6):
        w, v = (1,) * r + u, (1,) * (r - D) + up
        A, B = AB(w)
        j = len(w)
        q = fq(u)
        z0 = q.numerator * pow(q.denominator, -1, 3 ** j) % 3 ** j
        found = 0
        for z in range(z0, z0 + 3 ** j * 400, 3 ** j):
            if z % 2 == 0 or z <= 1:
                continue
            n = inv(w, z)
            if n is None:
                continue
            m = inv(v, z)
            good &= m is not None and m == (n + 1) // 2 ** D - 1 and (n + 1) % 2 ** D == 0
            found += 1
            if found >= 25:
                break
        tested += found
        good &= found > 0
check(good, f"{tested} sources in 5 collision families x 6 run lengths: the partner is exactly m = (n+1)/2^D - 1, same endpoint, same length")
check(fq((3,)) == fq((2, 1)) and all(fq((a,)) == fq((2, a - 2)) for a in range(3, 30)),
      "the reset switch (checked_switch (3)) is the collision (a) ~ (2, a-2), D = 1, m = (n-1)/2")

# ---------------------------------------------------------------- 4. longer partners and the reset-2 root collision
print("4. a longer partner (k > j) fires only on n = -1 mod 3^(k-j); for reset 2 the shifted root collision gives n = 8 mod 9")
cls, tot = Counter(), 0
for n in range(3, 100001, 2):
    xs, w = orbit(n, 60)
    r, u = strip1(w)
    if r == 0 or not u or u[0] != 2 or len(u) < 2 or r + 2 > len(w):
        continue
    tot += 1
    m = inv((1,) * (r + 3) + (u[1] + 2,), xs[r + 2])
    if m is not None and m < n:
        cls[n % 9] += 1
check(set(cls) == {8} and abs(sum(cls.values()) / tot - 1 / 9) < 0.002,
      f"reset-2 sources n < 10^5 ({tot}): the shifted root-collision join fires on {sum(cls.values())} = {sum(cls.values()) / tot:.4f}, all = 8 mod 9 (they have (2n-1)/3 < n)")

# ---------------------------------------------------------------- 5. census on reset-2 sources, the residual seeds, Mersenne numbers
print("5. census of uniform trailing-ones switches")


def switches(n, jmax):
    xs, w = orbit(n, jmax)
    r = strip1(w)[0]
    out = []
    for D in range(1, r + 1):
        m = (n + 1) // 2 ** D - 1
        ys, v = orbit(m, jmax)
        eq = next((j for j in range(1, min(len(xs), len(ys))) if xs[j] == ys[j]), None)
        if eq is not None:
            out.append((D, eq, fq(strip1(w[:eq])[1]) == fq(strip1(v[:eq])[1]), eq - r))
    return r, w, out


pop = cov = nonuni = 0
for n in range(3, 20001, 2):
    r, w, out = switches(n, 300)
    if r == 0 or len(w) <= r or w[r] != 2:
        continue
    pop += 1
    cov += bool(out)
    nonuni += sum(1 for o in out if not o[2])
check(pop == 2500 and cov == 377 and nonuni == 0,
      f"reset-2 sources n < 2*10^4: {cov} of {pop} have a uniform trailing-ones switch; non-collision equal-length merges: {nonuni}")
root = Path(__file__).resolve().parents[2]
seeds = json.load(open(root / "05-knowledge/results/checked_switch_phase19_20261004.json"))["compiler"][-1]["seed_sources"]
sw = [n for n in seeds if switches(n, 300)[2]]
check(len(sw) == 16 and sum(n % 3 == 0 for n in sw) == 10,
      f"residual seeds with a uniform switch: {len(sw)} of 239 ({sum(n % 3 == 0 for n in sw)} multiples of 3): {sw}")
mer = {}
for a in range(3, 122, 2):
    r, w, out = switches(2 ** a - 1, 4000)
    assert r == a - 1 and w[r] == 2
    mer[a] = out
hit = {a: o for a, o in mer.items() if o}
check(len(hit) == 37 and all(o[0][0] % 2 == 1 and all(x[2] for x in o) for o in hit.values()),
      f"odd Mersenne 2^a - 1 (a reset-2 family): {len(hit)} of 60 odd a in [3, 121] switch to some 2^(a-D) - 1; the least D is always odd"
      " (an even Mersenne number, which the reset switch takes to 2^(a-D-1) - 1); all are collisions")
print("   first cases (a: least D, merge depth beyond the run):", [(a, o[0][0], o[0][3]) for a, o in sorted(hit.items())[:10]])
print("ALL CHECKS PASSED" if OK else "SOME CHECK FAILED")
