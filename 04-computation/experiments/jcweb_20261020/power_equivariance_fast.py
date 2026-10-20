#!/usr/bin/env python3
"""Fast version (bitmasks over the (p-1)/2 classes {x,-x}) of the power-equivariance census, with the Gamma-even test.
For each prime p: Mult = multipliers, Gamma = {x -> a x^j}. Counts: Mult-types T, E (even circulants); Gamma-orbits on
T-sets and on E-sets; Gamma-even E-sets (no psi in Stab_Gamma(S) with an odd number of ι-flipped class-orbits inside S)
and their Gamma-orbits; the signed-set Burnside prediction (1/|Gamma|) sum_psi 2^{c(psi)} [psi has no flipped class-orbit]."""
import sys
from math import gcd
from sympy import isprime
PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 37
for p in [q for q in range(5, PMAX + 1) if isprime(q)]:
    h = (p - 1) // 2
    cls = list(range(1, h + 1)); idx = {}
    for i, x in enumerate(cls): idx[x] = i; idx[p - x] = i
    js = [j for j in range(1, p - 1) if gcd(j, p - 1) == 1]
    Gamma = [(a, j) for a in range(1, p) for j in js]
    def order(a):
        k = 1; x = a % p
        while x != 1: x = x * a % p; k += 1
        return k
    # class permutation and flip data for each psi
    perm = {}; flipsign = {}
    for (a, j) in Gamma:
        pr = [0] * h; fl = [0] * h   # fl[i] = 1 if psi maps class rep i to the NEGATIVE of class rep pr[i]
        for i, x in enumerate(cls):
            y = a * pow(x, j, p) % p
            pr[i] = idx[y]; fl[i] = 1 if y > h else 0
        perm[(a, j)] = pr; flipsign[(a, j)] = fl
    # orbits of psi on classes with flipped flag
    def class_orbits(key):
        pr, fl = perm[key], flipsign[key]
        seen = [False] * h; orbs = []
        for i in range(h):
            if seen[i]: continue
            O = []; x = i; s = 0
            while not seen[x]:
                seen[x] = True; O.append(x); s ^= fl[x]; x = pr[x]
            orbs.append((O, s))   # s = 1 iff psi^m(rep) = -rep (flipped orbit)
        return orbs
    corb = {key: class_orbits(key) for key in Gamma}
    # tournament sets as (orientation bitmask), graphs as class bitmask
    def act_T(mask, key):  # S -> psi(S): bit i = orientation of class i (1 means cls[i] in S)
        pr, fl = perm[key], flipsign[key]; out = 0
        for i in range(h):
            b = (mask >> i) & 1
            if b ^ fl[i]: out |= 1 << pr[i]
        return out
    def act_G(mask, key):
        pr = perm[key]; out = 0
        for i in range(h):
            if (mask >> i) & 1: out |= 1 << pr[i]
        return out
    N = 1 << h
    # Mult-even test via the criterion (THM-4721): for all even-order a in Stab_Mult(S): 2 ord(a) | |S| (|S| = 2 popcount)
    mult = [(a, 1) for a in range(1, p)]
    def mult_even(mask):
        c = bin(mask).count('1')
        if c == 0: return True
        if c == h: return False
        for key in mult:
            a = key[0]
            if act_G(mask, key) == mask:
                o = order(a)
                if o % 2 == 0 and (2 * c) % (2 * o) != 0: return False
        return True
    def gamma_even(mask):
        for key in Gamma:
            if act_G(mask, key) == mask:
                flips = sum(s for O, s in corb[key] if (mask >> O[0]) & 1)
                if flips % 2: return False
        return True
    evens = [m for m in range(N) if mult_even(m)]
    gevens = [m for m in evens if gamma_even(m)]
    def norbits(masks, act, group):
        S = set(masks); seen = set(); c = 0
        for m in masks:
            if m in seen: continue
            c += 1; stack = [m]; seen.add(m)
            while stack:
                x = stack.pop()
                for key in group:
                    y = act(x, key)
                    if y not in seen: seen.add(y); stack.append(y)
        return c
    tG = norbits(range(N), act_T, Gamma); eG = norbits(evens, act_G, Gamma); geG = norbits(gevens, act_G, Gamma)
    tM = norbits(range(N), act_T, mult); eM = norbits(evens, act_G, mult)
    # signed-set Burnside prediction for Gamma
    from fractions import Fraction
    pred = Fraction(0)
    for key in Gamma:
        orbs = corb[key]
        if all(s == 0 for O, s in orbs): pred += 2 ** len(orbs)
    pred /= len(Gamma)
    print(f"p = {p} ({'1' if p % 4 == 1 else '3'} mod 4): Mult types T {tM} E {eM} | Gamma-orbits T {tG}, E {eG}, Gamma-even E {geG} "
          f"(Gamma-even sets {len(gevens)} of {len(evens)} even) | Burnside prediction {pred} | "
          f"{'T = Gamma-even OK' if pred == tG == geG else 'MISMATCH'}; {'E has Gamma-odd members' if eG != geG else 'all even are Gamma-even'}")
