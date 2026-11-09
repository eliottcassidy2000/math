"""Ternary lane reduction (the base-3 analogue of THM-4094's lanes).

M = {v >= 1 : v = 1 (mod 3)} under multiplication; M3 = M minus {10, 31}.  M3 is closed under multiplication (10 and 31
have no factorisation into elements of M other than themselves: 10 = 2 * 5 with 2, 5 = 2 (mod 3), and 31 is prime).
H_3 is multiplicative over order joins (THM-4940 (1)), so Spec_3 is the multiplicative closure of the strong ternary
values, and Spec_3 is contained in M3 (THM-4940 (3), (4)).  Hence

    ternary Chan-Pak (Spec_3 = M3)  <=>  every irreducible element of M3 is the H_3 of a STRONG tournament.

This script enumerates the irreducibles of M3 below X by direct factorisation, classifies them by shape, and checks
the shape theorem:
  irreducibles of M3 = { p prime, p = 1 (mod 3), p != 31 }                      (lane P)
                     u { p q : p <= q primes = 2 (mod 3), p q != 10 }           (lane Q)
                     u { 10 p : p prime = 1 (mod 3), p != 31 }                  (lane 10P)
                     u { 31 p : p prime = 1 (mod 3) } (includes 961 = 31^2)     (lane 31P)
                     u { 31 p q : p <= q primes = 2 (mod 3), p q != 10 }        (lane 31Q)
                     u { 40, 250, 310, 9610, 29791 }    (sporadic: 2^3 5 = 4*10, 2 5^3 = 10*25, 10*31, 10*31^2, 31^3)
usage: python ternary_lanes.py X [attained_values_file]
"""
import sys


def sieve(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return s


def main():
    X = int(sys.argv[1])
    isp = sieve(X)
    inM3 = lambda v: v % 3 == 1 and v not in (10, 31)
    # reducible[v] = True if v = x * y with x, y in M3, x, y > 1
    red = bytearray(X + 1)
    elems = [v for v in range(4, X + 1) if inM3(v)]
    for i, x in enumerate(elems):
        if x * x > X:
            break
        for y in elems[i:]:
            v = x * y
            if v > X:
                break
            red[v] = 1          # x, y in M3 and x*y = v; v = 1 mod 3 automatically, and v != 10, 31
    irr = [v for v in elems if not red[v]]
    # classify by shape
    def factor(v):
        f = []
        d = 2
        while d * d <= v:
            while v % d == 0:
                f.append(d)
                v //= d
            d += 1
        if v > 1:
            f.append(v)
        return f
    def lane(v):
        f = factor(v)
        r1 = [p for p in f if p % 3 == 1]
        r2 = [p for p in f if p % 3 == 2]
        if len(f) == 1 and r1 and v != 31:
            return "P"
        if len(f) == 2 and len(r2) == 2:
            return "Q"
        if v in (40, 250, 310, 9610, 29791):
            return "sporadic"
        if len(f) == 3 and sorted(r2) == [2, 5] and len(r1) == 1 and r1[0] != 31:
            return "10P"
        if len(f) == 2 and f.count(31) >= 1 and len(r1) == 2:
            return "31P"
        if len(f) == 3 and f.count(31) == 1 and len(r2) == 2 and sorted(r2) != [2, 5]:
            return "31Q"
        return "OTHER"
    counts = {}
    other = []
    for v in irr:
        L = lane(v)
        counts[L] = counts.get(L, 0) + 1
        if L == "OTHER":
            other.append(v)
    print(f"irreducibles of M3 below {X}: {len(irr)}; by lane: {dict(sorted(counts.items()))}")
    print(f"irreducibles of any other shape: {other[:20]} ({len(other)})")
    # converse check: every element of the predicted lanes below X is irreducible
    pred = set()
    for v in range(4, X + 1):
        if not inM3(v):
            continue
        if lane(v) != "OTHER":
            pred.add(v)
    print(f"shape theorem check: irreducibles == predicted lane elements below {X}: {set(irr) == pred}")
    print(f"first irreducibles: {irr[:40]}")
    if len(sys.argv) > 2:
        att = set(int(x) for x in open(sys.argv[2]).read().split())
        miss = [v for v in irr if v not in att]
        print(f"irreducibles below {X} not attained by a strong piece in the supplied list: {len(miss)}; first {miss[:30]}")
        by = {}
        for v in miss:
            by[lane(v)] = by.get(lane(v), 0) + 1
        print(f"  missing by lane: {dict(sorted(by.items()))}")


if __name__ == "__main__":
    main()
