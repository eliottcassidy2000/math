#!/usr/bin/env python3
"""
The moment-count multiplier is prime-structured -- the "factor 2" is only p_1  (mac-mini-S150)
================================================================================
Owner's directive: numbers like 2, 3, 7 are not constants but the x-th member of the prime
family; when a prime appears in an equation, ask what the equation says for EVERY other prime.
No coincidences.

Applied to THM-1725's moment bound M* <= 2 * max_pairs (p+|n|)/gcd(p,|n|).  The factor 2 is
opus THM-1685's CT(m0)+CT(2m0), i.e. the PRIMITIVE level m0 plus its FIRST composite multiple
2m0.  By THM-415 (prime modulus = no nontrivial vanishing; composite = collisions) the extra
levels are needed EXACTLY to break collisions -- which happen at COMPOSITE multiples of m0.
So the honest object is the multiplier
    mu(pattern) = M* / m0,   m0 = max coprime charge-pair sum,
and the prediction is: mu is governed by the PRIME FACTORISATION of m0 and of the charge
lattice, NOT by the bare constant 2.

Tests:
  (A) single coprime straddle, m0 = q+1 across q -- does mu track primality of m0?
  (B) two independent coprime structures (charge sums a two distinct primes) -- does M* see
      the PRODUCT / LCM, i.e. BOTH primes, not just the larger?
  (C) charge-count growth -- does the multiplier bring in higher primes as k grows, refuting
      the universal '2'?
"""
import sympy as sp
from math import factorial, gcd
from sympy import primerange, factorint

def moments(monos, coeffs, M):
    P = {}
    for i, (a, b) in enumerate(monos): P[(a, b)] = P.get((a, b), 0) + coeffs[i]
    def mul(X, Y):
        o = {}
        for (a1, b1), c1 in X.items():
            for (a2, b2), c2 in Y.items():
                k = (a1+a2, b1+b2); o[k] = o.get(k, 0) + c1*c2
        return o
    Pm = {(0, 0): sp.Integer(1)}; g = []
    for m in range(1, M+1):
        Pm = mul(Pm, P)
        g.append(sp.expand(sum(c*factorial(a) for (a, b), c in Pm.items() if a == b)))
    return g

def minM(monos, coeffs, Mmax=26):
    pos = [coeffs[i] for i, (a, b) in enumerate(monos) if a > b]
    neg = [coeffs[i] for i, (a, b) in enumerate(monos) if a < b]
    if not pos or not neg: return 0
    w = sp.Symbol('w'); allg = moments(monos, coeffs, Mmax)
    for M in range(1, Mmax+1):
        ok = True
        for cp in pos:
            for cn in neg:
                G = sp.groebner(list(allg[:M]) + [1 - w*cp*cn], *(list(coeffs)+[w]), order='grevlex')
                if not any(x.is_number and x != 0 for x in G.exprs): ok = False; break
            if not ok: break
        if ok: return M
    return None

def m0_of(charges):
    P = [c for c in charges if c > 0]; N = [-c for c in charges if c < 0]
    return max((p+n)//gcd(p, n) for p in P for n in N)

print("=" * 78)
print("PART A -- single coprime straddle: does the multiplier track PRIMALITY of m0?")
print("=" * 78)
print("  Pattern [Z^q, W, ZW] : charges +q, -1, 0(rad2).  m0 = q+1 (coprime pair (q,1)).")
print("  Prediction (THM-415): m0 PRIME => primitive level alone saturates (mu=1);")
print("                        m0 COMPOSITE => a collision at the first prime factor forces")
print("                        an extra level (mu>1).")
print()
print(f"{'q':>3} {'m0=q+1':>7} {'factor(m0)':>14} {'M*':>4} {'mu=M*/m0':>9} {'m0 prime?':>10}")
for q in range(1, 9):
    monos = [(q, 0), (0, 1), (1, 1)]           # Z^q, W, ZW
    coeffs = sp.symbols('c0:3')
    Ms = minM(list(monos), list(coeffs))
    m0 = q+1
    fac = factorint(m0)
    prime = sp.isprime(m0)
    mu = (Ms/m0) if Ms else None
    print(f"{q:>3} {m0:>7} {str(dict(fac)):>14} {str(Ms):>4} "
          f"{(f'{mu:.3f}' if mu else 'n/a'):>9} {str(prime):>10}")

print()
print("=" * 78)
print("PART B -- two coprime structures: does M* see BOTH primes, not just the bigger?")
print("=" * 78)
print("  Pattern [Z^a, W^b, Zbar-ish] engineered so the pos-neg pairs give coprime sums p, q")
print("  at two DISTINCT primes.  If M* ~ max(p,q): only the visible prime matters.")
print("  If M* ~ p*q or lcm: the OTHER prime is really there (owner's 'no coincidences').")
print()
print(f"{'charges':>16} {'pair sums':>16} {'M*':>4} {'max':>5} {'sum':>5} {'lcm':>5} {'which fits':>12}")
for monos in ([(2, 0), (0, 3), (1, 1)],          # +2,-3,0 : pair (2,3) sum 5
              [(3, 0), (0, 2), (1, 1)],          # +3,-2,0 : (3,2) sum 5
              [(2, 0), (0, 1), (0, 3)],          # +2,-1,-3 : pairs (2,1)=3,(2,3)=5
              [(3, 0), (0, 1), (0, 2)],          # +3,-1,-2 : (3,1)=4,(3,2)=5
              [(4, 0), (0, 1), (0, 3)],          # +4,-1,-3 : (4,1)=5,(4,3)=7
              [(5, 0), (0, 1), (0, 2)]):         # +5,-1,-2 : (5,1)=6,(5,2)=7
    coeffs = sp.symbols('c0:%d' % len(monos))
    Ms = minM(list(monos), list(coeffs))
    charges = sorted({a-b for a, b in monos})
    P = [c for c in charges if c > 0]; N = [-c for c in charges if c < 0]
    sums = sorted({(p+n)//gcd(p, n) for p in P for n in N})
    mx = max(sums); sm = sum(sums)
    import math
    lc = sums[0]
    for x in sums[1:]: lc = lc*x//gcd(lc, x)
    fit = ("max" if Ms == mx else "sum" if Ms == sm else "lcm" if Ms == lc
           else "2*max" if Ms == 2*mx else "other")
    print(f"{str(charges):>16} {str(sums):>16} {str(Ms):>4} {mx:>5} {sm:>5} {lc:>5} {fit:>12}")

print()
print("=" * 78)
print("PART C -- charge-count growth: does the multiplier exceed 2 (refuting universal '2')?")
print("=" * 78)
print("  Nested coprime straddles.  If mu stays <= 2 the factor-2 bound survives; if it grows")
print("  the bound is only the FIRST term of a prime-indexed family.")
print()
print(f"{'monomials (a,b)':>34} {'charges':>16} {'m0':>4} {'M*':>4} {'mu':>6}")
for monos in ([(2, 0), (0, 1), (1, 1)],                          # k=3
              [(2, 0), (0, 1), (1, 1), (0, 0)],                  # k=4 + charge-0
              [(3, 0), (0, 2), (1, 1), (0, 0), (2, 2)],          # k=5 deep charge-0
              [(2, 0), (3, 0), (0, 1), (0, 2), (1, 1)]):         # k=5 two-sided rich
    coeffs = sp.symbols('c0:%d' % len(monos))
    Ms = minM(list(monos), list(coeffs), Mmax=20)
    charges = sorted({a-b for a, b in monos})
    m0 = m0_of(charges)
    mu = (Ms/m0) if Ms else None
    print(f"{str(monos):>34} {str(charges):>16} {m0:>4} {str(Ms):>4} "
          f"{(f'{mu:.3f}' if mu else 'n/a'):>6}")

print()
print("SUMMARY -- reading the multiplier through the prime family")
print("  The bare '2' in THM-1725 is p_1 = the first place a composite multiple of m0 can")
print("  collide (THM-415).  The honest statement is about the prime factorisation of m0 and")
print("  of the charge lattice, and the same equation, read at every prime, is one family.")
