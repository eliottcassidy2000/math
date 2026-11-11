"""AUDIT: section 7 (the monoid M of odd numbers other than 7 and 21), brute force below LIMIT.

 - atoms of M from the definition (no factorisation m = a b with a, b in M, a, b > 1);
 - the set of factorisation lengths of every element (dynamic programme over divisors);
 - the largest ratio max length / min length and where it is attained; the smallest element with two lengths;
 - the same for the sub-question "is 63 the cause": elements NOT divisible by 63 with two different lengths.
usage: python3 aud_elasticity.py [LIMIT]
"""
import sys
from fractions import Fraction

LIMIT = int(sys.argv[1]) if len(sys.argv) > 1 else 50000


def inM(m):
    return m % 2 == 1 and m not in (7, 21)


# divisors d of m with 1 < d <= m/d, both in M
divs = {m: [] for m in range(1, LIMIT, 2)}
for d in range(3, int(LIMIT ** 0.5) + 2, 2):
    if not inM(d):
        continue
    for e in range(d, LIMIT // d + 1, 2):
        if d * e < LIMIT and inM(e):
            divs[d * e].append(d)

atoms = [m for m in range(3, LIMIT, 2) if inM(m) and not divs[m]]
atomset = set(atoms)


def is_prime(x):
    if x < 2:
        return False
    i = 2
    while i * i <= x:
        if x % i == 0:
            return False
        i += 1
    return True


# predicted atoms: odd primes other than 7; 7p with p an odd prime other than 3; 63; 343
pred = set()
for m in range(3, LIMIT, 2):
    if is_prime(m) and m != 7:
        pred.add(m)
    elif m % 7 == 0 and is_prime(m // 7) and m // 7 not in (2, 3):
        pred.add(m)
    elif m in (63, 343):
        pred.add(m)
print("limit %d: atoms found %d; equal to the predicted list (odd primes != 7; 7p, p odd prime != 3; 63; 343): %s" % (
    LIMIT, len(atoms), atomset == pred))
print("atoms divisible by 9: %s" % [a for a in atoms if a % 9 == 0])

# length sets
lengths = {1: frozenset([0])}
for m in range(3, LIMIT, 2):
    if not inM(m):
        continue
    if m in atomset:
        lengths[m] = frozenset([1])
        continue
    s = set()
    for d in divs[m]:
        e = m // d
        for a in lengths[d]:
            for b in lengths[e]:
                s.add(a + b)
    lengths[m] = frozenset(s)

best = Fraction(1)
arg = []
first_two = None
first_two_no63 = None
count_two = 0
count_two_no63 = 0
best_no63 = Fraction(1)
arg_no63 = []
for m in sorted(lengths):
    L = lengths[m]
    if m == 1:
        continue
    r = Fraction(max(L), min(L))
    if len(L) > 1:
        count_two += 1
        if first_two is None:
            first_two = (m, sorted(L))
        if m % 63:
            count_two_no63 += 1
            if first_two_no63 is None:
                first_two_no63 = (m, sorted(L))
    if r > best:
        best, arg = r, [m]
    elif r == best:
        arg.append(m)
    if m % 63:
        if r > best_no63:
            best_no63, arg_no63 = r, [m]
        elif r == best_no63 and r > 1:
            arg_no63.append(m)
print("largest ratio of two factorisation lengths below %d: %s, attained at %s" % (LIMIT, best, arg))
print("smallest element with two different lengths: %s" % (first_two,))
print("elements with at least two lengths: %d; among them NOT divisible by 63: %d; the smallest of those: %s" % (
    count_two, count_two_no63, first_two_no63))
print("largest ratio among elements not divisible by 63: %s, first attained at %s" % (best_no63, arg_no63[:5]))
for m in (315, 1225, 3969, 117649):
    if m in lengths:
        print("   lengths(%d) = %s" % (m, sorted(lengths[m])))

# formula of the note: m = 3^a 7^b k, c = Omega(k); lengths = { c + a + v + x - w : u + 2v + w + 3x = b, 2w <= a, u <= c' }
def formula_lengths(m):
    a = b = 0
    k = m
    while k % 3 == 0:
        k //= 3
        a += 1
    while k % 7 == 0:
        k //= 7
        b += 1
    c = 0
    kk = k
    p = 5
    while kk > 1:
        while kk % p == 0:
            kk //= p
            c += 1
        p += 2
    res = set()
    for w in range(0, a // 2 + 1):
        for x in range(0, b // 3 + 1):
            for v in range(0, b // 2 + 1):
                u = b - 2 * v - w - 3 * x
                if 0 <= u <= c:
                    res.add(c + a + v + x - w)
    return res

bad = 0
for m, L in lengths.items():
    if m > 1 and set(L) != formula_lengths(m):
        bad += 1
print("length formula of the note (L = c + a + v + x - w with u + 2v + w + 3x = b, 2w <= a, u <= c) disagrees on %d elements" % bad)
