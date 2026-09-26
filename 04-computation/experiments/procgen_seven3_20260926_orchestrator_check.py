#!/usr/bin/env python3
"""Orchestrator audit of lane `seven3` (itinerary-coded 7n+-1 strategies), written from the note's statements;
the lane's scripts were not read. q = 7, max-halving (MH) sign s(x) with 7x + s = 0 mod 4.

  1. Lemma R: for 2000 random 600-bit x, one flip at x followed by MH; whenever the flipped orbit meets the MH
     orbit of x (exact integer equality) the two paths have the same numbers of odd steps and of halvings.
  2. Lemma S: for 3000 random odd z whose MH itinerary starts (-s,v_1),...,(-s,v_n),(-s,.) with v_i in {2,3},
     E = (z - s)/2 is odd, has MH valuations v_1..v_n, and MH^i(E) = (z_i - s)/2.
  3. Proposition A (lower-bound witness): for D = 3..12 the MH-periodic rational point with alternating signs and
     valuations 2^(D-1) 3 is never flipped by the rule G + A_D and has density D/(2D+1).
"""
import random
from fractions import Fraction as Fr


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)
    print("  ok:", msg, flush=True)


def v2(z):
    return (z & -z).bit_length() - 1


def mh_step(x):
    """accelerated MH step on an exact integer (an element of Z_2): returns (s, v, x'); exact arithmetic,
    so orbit points are compared exactly (no truncation bias)."""
    s = 1 if (7 * x + 1) % 4 == 0 else -1
    z = 7 * x + s
    v = v2(z)
    return s, v, (z >> v)


# ---------------------------------------------------------------- 1. Lemma R
random.seed(3)
rejoins = 0
for trial in range(2000):
    x = random.randrange(1, 1 << 600, 2)
    # MH orbit of x with cumulative (odd, halvings)
    mh = {}
    y, odd, hv = x, 0, 0
    for i in range(60):
        mh.setdefault(y, (odd, hv))
        s, v, y = mh_step(y)
        odd += 1
        hv += v
    # flipped: first step uses -s with one halving (valuation 1 for the flipped sign), then MH
    s = 1 if (7 * x + 1) % 4 == 0 else -1
    y = (7 * x - s) >> 1
    odd, hv = 1, 1
    for i in range(60):
        if y in mh and (odd, hv) != (0, 0):
            o2, h2 = mh[y]
            assert (o2, h2) == (odd, hv), (trial, (o2, h2), (odd, hv))
            rejoins += 1
            break
        s2, v, y = mh_step(y)
        odd += 1
        hv += v
check(rejoins > 0, f"Lemma R: {rejoins} of 2000 random flips rejoin the MH orbit within 60 steps, all with equal odd and halving counts")

# ---------------------------------------------------------------- 2. Lemma S
random.seed(4)
tested = 0
while tested < 3000:
    z = random.randrange(1, 1 << 400, 2)
    # read the itinerary
    it = []
    y = z
    for i in range(12):
        s, v, y2 = mh_step(y)
        it.append((s, v, y))
        y = y2
    s0 = -it[0][0]                       # the lemma's s: the itinerary starts with sign -s
    n = 0
    while n < len(it) - 1 and it[n][0] == -s0 and it[n][1] in (2, 3):
        n += 1
    if n == 0 or it[n][0] != -s0:
        continue
    E = (z - s0) >> 1
    assert E & 1
    y = E
    zi = z
    for i in range(n):
        s, v, y2 = mh_step(y)
        assert v == it[i][1], (i, v, it[i][1])
        _, _, zi = mh_step(zi)
        assert y2 == ((zi - s0) >> 1), i
        y = y2
    tested += 1
check(True, f"Lemma S: {tested} random itineraries (-s,v_1..v_n),(-s,.) with v_i in {{2,3}}: E = (z-s)/2 shadows with equal valuations and MH^i(E) = (z_i - s)/2")


# ---------------------------------------------------------------- 3. Proposition A witness
def mh_sign_frac(x):
    r = x.numerator * pow(x.denominator, -1, 1 << 64) % (1 << 64)
    return 1 if (7 * r + 1) % 4 == 0 else -1


def itinerary_frac(x, n):
    out = []
    for _ in range(n):
        s = mh_sign_frac(x)
        y = 7 * x + s
        v = 0
        while y.numerator % 2 == 0:
            y /= 2
            v += 1
        out.append((s, v))
        x = y
    return out


def flips(it, D):
    """rule G + A_D on the itinerary read from the current point."""
    (s1, v1), (s2, v2_) = it[0], it[1]
    if v1 == 2 and v2_ == 2 and s1 == s2:
        return True
    if all(v == 2 for (s, v) in it[:D]) and all(it[i][0] == -it[i + 1][0] for i in range(D - 1)):
        return True
    return False


for D in range(3, 13):
    per = D if D % 2 == 0 else 2 * D
    vals = ([2] * (D - 1) + [3]) * (per // D)
    signs = [(-1) ** i for i in range(per)]
    # fixed point of the composition of inverse branches y -> (2^v y - s)/7 (first symbol applied first)
    a, b = Fr(1), Fr(0)
    for s, v in reversed(list(zip(signs, vals))):
        a, b = a * Fr(2 ** v, 7), (b * 2 ** v - s) / 7
    x0 = b / (1 - a)
    x = x0
    odd = steps = 0
    for i in range(per):
        it = itinerary_frac(x, D + 2)
        assert it[0] == (signs[i], vals[i]), (D, i, it[0])
        assert not flips(it, D), (D, i)
        s, v = it[0]
        x = (7 * x + s) / 2 ** v
        odd += 1
        steps += v
    assert x == x0
    assert Fr(odd, steps) == Fr(D, 2 * D + 1), (D, odd, steps)
check(True, "Proposition A: for D = 3..12 the alternating MH-periodic point with valuations 2^(D-1) 3 is never flipped by G + A_D and has density D/(2D+1)")
