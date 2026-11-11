#!/usr/bin/env python3
"""Closed forms for the triangle and pentagon parity races (Theorem P of the note).

For a tournament on n vertices with adjacency matrix A over F_2 let s = A1 (odd scores) and t = As.  Then (mod 2)
    c_3 = C(|s|, 2) + floor(n/4),    c_3 + c_5 = C(|t|, 2) + floor(n/4),    c_5 = C(|s|, 2) + C(|t|, 2).
Joint law of (s, t) for a uniformly random labelled tournament (n >= 1):
    s is uniform on {|s| = C(n,2) mod 2};  given s with support O, 1 <= |O| <= n-1, t is uniform on
    {t : |t restricted to O| = C(|O|,2), |t restricted to the complement| = C(n-|O|,2) (mod 2)};
    t = 0 if O is empty, t = 1 (all ones) if O is everything.
Part 1 checks the three congruences and the joint law by exhaustive enumeration of all labelled tournaments (n <= NMAX).
Part 2 evaluates E_n[{3}], E_n[{5}], E_n[{3,5}] exactly from the joint law, compares them with the census values
(n <= 11) and with the closed forms stated in the note.
Closed forms (kron2(n) = Kronecker symbol (2/n)):
    E_n[{3}]   = 2^(-floor((n-1)/2));
    E_n[{3,5}] = 2^(-(n-2)/2) for n = 2 (mod 4);  2^(-(n-1)/2) (2 + (2/n) 2^(-(n-1)/2) - 2^(-(n-2))) for odd n;
                 the value at n - 1 for n = 0 (mod 4);
    E_n[{5}]   = 2^(-(n-2)) for n = 2 (mod 4);  2^(-(n-1)) (3 - (2/n) 2^(-(n-3)/2)) for odd n >= 3;
                 the value at n - 1 for n = 0 (mod 4).
Usage: python3 pentagon_closed_forms.py [NMAX exhaustive, default 6] [largest n printed, default 30] [largest n checked]"""
import sys
from fractions import Fraction as Fr
from math import comb
from itertools import combinations
from collections import Counter

def exhaustive(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    law = Counter(); bad = 0
    for bits in range(1 << len(pairs)):
        out = [0] * n
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1: out[i] |= 1 << j
            else: out[j] |= 1 << i
        s = 0
        for a in range(n):
            if bin(out[a]).count('1') & 1: s |= 1 << a
        t = 0
        for a in range(n):
            if bin(out[a] & s).count('1') & 1: t |= 1 << a
        law[(s, t)] += 1
        # cycle counts
        c3 = comb(n, 3) - sum(comb(bin(o).count('1'), 2) for o in out)
        c5 = 0
        for S in combinations(range(n), 5):
            v = S[0]; rest = S[1:]
            # 5-cycles on S through v: orderings of the rest
            from itertools import permutations
            for p in permutations(rest):
                seq = (v,) + p
                if all((out[seq[i]] >> seq[(i + 1) % 5]) & 1 for i in range(5)): c5 += 1
        q1 = bin(s).count('1'); q2 = bin(t).count('1'); kap = (n // 4) & 1
        if (comb(q1, 2) + kap - c3) % 2: bad += 1
        if (comb(q2, 2) + kap - c3 - c5) % 2: bad += 1
        if (comb(q1, 2) + comb(q2, 2) - c5) % 2: bad += 1
    # joint law
    total = 1 << len(pairs); badlaw = 0; support = 0
    for s in range(1 << n):
        q = bin(s).count('1')
        if (q - comb(n, 2)) % 2: continue
        ps = Fr(1, 2 ** (n - 1))
        for t in range(1 << n):
            if q == 0: ok = (t == 0); pt = Fr(1)
            elif q == n: ok = (t == (1 << n) - 1); pt = Fr(1)
            else:
                k1 = bin(t & s).count('1'); k2 = bin(t & ~s & ((1 << n) - 1)).count('1')
                ok = ((k1 - comb(q, 2)) % 2 == 0) and ((k2 - comb(n - q, 2)) % 2 == 0); pt = Fr(1, 2 ** (n - 2))
            pred = ps * pt * total if ok else 0
            if pred: support += 1
            if law.get((s, t), 0) != pred: badlaw += 1
    return bad, badlaw, support

def expect(n, f):
    """E f(|s|, |t|) from the joint law."""
    tot = Fr(0)
    for q in range(0, n + 1):
        if (q - comb(n, 2)) % 2: continue
        if q == 0: tot += f(0, 0); continue
        if q == n: tot += f(n, n); continue
        p1 = comb(q, 2) % 2; p2 = comb(n - q, 2) % 2; inner = Fr(0)
        for k1 in range(p1, q + 1, 2):
            for k2 in range(p2, n - q + 1, 2):
                inner += comb(q, k1) * comb(n - q, k2) * f(q, k1 + k2)
        tot += comb(n, q) * inner / 2 ** (n - 2)
    return tot / 2 ** (n - 1)

def E3(n):  return expect(n, lambda q, k: (-1) ** (comb(q, 2) + n // 4))
def E35(n): return expect(n, lambda q, k: (-1) ** (comb(k, 2) + n // 4))
def E5(n):  return expect(n, lambda q, k: (-1) ** (comb(q, 2) + comb(k, 2)))

def closed35(n):
    half = Fr(1, 2)
    if n % 4 == 2: return half ** ((n - 2) // 2)
    if n % 4 == 0: n -= 1
    h = (n - 1) // 2
    sign = (-1) ** ((n - 1) // 4) if n % 4 == 1 else -(-1) ** ((n - 3) // 4)
    return half ** h * (2 + sign * half ** h - half ** (n - 2))

def kron2(n):
    """Kronecker symbol (2/n) for odd n: +1 if n = +-1 (mod 8), -1 if n = +-3 (mod 8)."""
    return 1 if n % 8 in (1, 7) else -1

def closed35_odd(n):
    half = Fr(1, 2); h = (n - 1) // 2
    return half ** h * (2 + kron2(n) * half ** h - half ** (n - 2))

def closed5(n):
    half = Fr(1, 2)
    if n % 4 == 2: return half ** (n - 2)
    if n % 4 == 0: n -= 1
    return half ** (n - 1) * (3 - kron2(n) * half ** ((n - 3) // 2)) if n >= 3 else Fr(1)

CENSUS = {  # exact values from the censuses (cycrace.c for n <= 10, k35cube.c for n = 11)
 5: (Fr(1,4), Fr(7,32), Fr(13,32)), 6: (Fr(1,4), Fr(1,16), Fr(1,4)), 7: (Fr(1,8), Fr(11,256), Fr(67,256)),
 8: (Fr(1,8), Fr(11,256), Fr(67,256)), 9: (Fr(1,16), Fr(23,2048), Fr(263,2048)), 10: (Fr(1,16), Fr(1,256), Fr(1,16)),
 11: (Fr(1,32), Fr(49,16384), Fr(1007,16384)) }

if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    nprint = int(sys.argv[2]) if len(sys.argv) > 2 else 30
    print("Part 1: exhaustive check over all labelled tournaments")
    for n in range(1, nmax + 1):
        bad, badlaw, support = exhaustive(n)
        print("  n = %d: %d tournaments; congruence failures %d; joint-law mismatches %d (support %d pairs (s, t))" % (n, 2 ** comb(n, 2), bad, badlaw, support))
    print("Part 2: exact expectations from the joint law")
    print("  n   E[{3}]                 E[{5}]                           E[{3,5}]                         E[{3,5}]/E[{3}]   E[{5}]/E[{3}]^2")
    fail = 0
    for n in range(1, nprint + 1):
        e3, e5, e35 = E3(n), E5(n), E35(n)
        if e3 != Fr(1, 2 ** ((n - 1) // 2)): fail += 1
        if e35 != closed35(n): fail += 1
        if n % 2 == 1 and e35 != closed35_odd(n): fail += 1
        if e5 != closed5(n): fail += 1
        if n in CENSUS and (e3, e5, e35) != CENSUS[n]: fail += 1
        print("  %-3d %-22s %-32s %-32s %-17.10f %.10f" % (n, e3, e5, e35, float(e35 / e3), float(e5 / e3 ** 2)))
    print("mismatches with Theorem T, the closed forms for {3,5} and {5}, or the census values (n <= %d): %d" % (nprint, fail))
    nbig = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    if nbig > nprint:
        fail2 = 0
        for n in range(nprint + 1, nbig + 1):
            if E3(n) != Fr(1, 2) ** ((n - 1) // 2) or E35(n) != closed35(n) or E5(n) != closed5(n): fail2 += 1
        print("closed forms against the joint-law sums for %d < n <= %d (exact rational arithmetic): %d mismatches" % (nprint, nbig, fail2))
        print("(for n in a fixed class mod 8 both sides are combinations of at most five geometric sequences in n,")
        print(" so agreement at more than five consecutive members of each class proves the closed forms for every n)")
