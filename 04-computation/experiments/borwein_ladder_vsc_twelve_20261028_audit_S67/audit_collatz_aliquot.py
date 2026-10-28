"""Independent referee re-computation (audit of opus-2026-10-28-S67) -- section 5 of the note:
the Collatz word 231 -> 233 and the aliquot / square-sigma facts about 231, 233, 420.

Run:  python audit_collatz_aliquot.py
"""
from math import isqrt

FAIL = []


def say(*a):
    print(*a, flush=True)


def check(cond, label):
    say(("  [PASS] " if cond else "  [FAIL] ") + label)
    if not cond:
        FAIL.append(label)


def syr(n):
    """Syracuse map on odd n: (3n+1)/2^v; also returns v."""
    n = 3 * n + 1
    v = 0
    while n % 2 == 0:
        n //= 2
        v += 1
    return n, v


say("=" * 90)
say("1. Collatz: the 17-step Syracuse word 231 -> 233 and its cylinder")
say("=" * 90)
w, vs = [231], []
for _ in range(17):
    nxt, v = syr(w[-1])
    w.append(nxt)
    vs.append(v)
say("  word: %s" % w)
say("  valuations: %s (sum %d)" % (vs, sum(vs)))
check(w == [231, 347, 521, 391, 587, 881, 661, 31, 47, 71, 107, 161, 121, 91, 137, 103, 155, 233],
      "the recorded word")
check(sum(vs) == 27, "27 halvings in 17 odd steps (the 27/17 clock)")
# the affine law on the 2-adic cylinder n = 231 mod 2^28: T^17(n) = (3^17 n + c)/2^27
c = 233 * 2 ** 27 - 3 ** 17 * 231
ok = True
for t in list(range(0, 2001)) + [10 ** 6 + 7, 2 ** 40 + 3]:
    x = 231 + 2 ** 28 * t
    vv = []
    for _ in range(17):
        x, v = syr(x)
        vv.append(v)
    if not (vv == vs and x == 233 + 2 * 3 ** 17 * t == (3 ** 17 * (231 + 2 ** 28 * t) + c) // 2 ** 27):
        ok = False
check(ok, "m_17(231 + 2^28 t) = 233 + 2*3^17 t with the same valuation word, t = 0..2000 and two large t")
x = 231 + 2 ** 27
vv = []
for _ in range(17):
    x, v = syr(x)
    vv.append(v)
say("  modulus check: 231 + 2^27 has valuation word %s (differs from 231's: %s)" % (vv[-3:], vv != vs))

say("")
say("=" * 90)
say("2. Aliquot / square-sigma facts")
say("=" * 90)
LIM = 420 ** 2 + 10          # covers the exhaustive ranges for 231, 233 and 420
sig = [0] * (LIM + 1)
for d in range(1, LIM + 1):
    for mlt in range(d, LIM + 1, d):
        sig[mlt] += d


def s(n):
    return sig[n] - n


def is_sq(x):
    r = isqrt(x)
    return r * r == x


check(s(345) == 231 and sig[345] == 576 == 24 ** 2, "s(345) = 231, sigma(345) = 576 = 24^2")
check(s(364) == 420 and sig[364] == 784 == 28 ** 2, "s(364) = 420, sigma(364) = 784 = 28^2")
pre233 = [n for n in range(2, 233 ** 2 + 1) if s(n) == 233]
sq233 = [n for n in pre233 if is_sq(sig[n])]
say("  all n <= 233^2 with s(n) = 233: %d of them, e.g. %s ; with square sigma: %s" % (len(pre233), pre233[:6], sq233))
check(sq233 == [], "no n with square sigma has s(n) = 233 (composite n has s(n) >= 1 + sqrt(n), so n < 233^2 is exhaustive)")
# the bound used for exhaustiveness
okb = all(s(n) >= 1 + isqrt(n) for n in range(4, 20000) if sig[n] != n + 1)
check(okb, "s(n) >= 1 + floor(sqrt n) for composite n < 20000 (sanity of the exhaustiveness bound)")
pre231 = [n for n in range(2, 231 ** 2 + 1) if s(n) == 231 and is_sq(sig[n])]
pre420 = [n for n in range(2, 420 ** 2 + 1) if s(n) == 420 and is_sq(sig[n])]
say("  complete square-sigma preimage lists: of 231 -> %s ; of 420 -> %s" % (pre231, pre420))

say("")
say("FAILURES: %s" % (FAIL if FAIL else "none"))
