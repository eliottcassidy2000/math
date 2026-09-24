#!/usr/bin/env python3
"""
procgen_atlas_20260924_primegame.py

Implication-atlas lane (session collatz-procgen-20260922, 2026-09-24): the logic layer.

  L1  Conway's PRIMEGAME as a generalized Collatz map: g(n) = f_j n for the first fraction f_j in
      (17/91, 78/85, 19/51, 23/38, 29/33, 77/29, 95/23, 77/19, 1/17, 11/13, 13/11, 15/14, 15/2, 55/1)
      with f_j n integral.  Conway's theorem: the powers of 2 in the orbit of 2 are exactly 2^p, p prime,
      in increasing order.  Hence the twin prime conjecture is equivalent to "the orbit of 2 contains
      infinitely many consecutive power-of-2 visits 2^p, 2^(p+2)" -- a Pi^0_2 statement about one
      explicit residue-class-wise linear map.  Here: the first primes are reproduced exactly.
  L2  arithmetical complexity of the atlas nodes (elementary classification; the note gives the reasons).
"""
import math

# Guy 1983 (Math. Mag. 56) variant, ending 15/14, 15/2, 55/1 [R]; the Wikipedia/OEIS list ends 15/2, 1/7, 55/1 [R].
FRACS_GUY = [(17, 91), (78, 85), (19, 51), (23, 38), (29, 33), (77, 29), (95, 23), (77, 19), (1, 17),
             (11, 13), (13, 11), (15, 14), (15, 2), (55, 1)]
FRACS_STD = [(17, 91), (78, 85), (19, 51), (23, 38), (29, 33), (77, 29), (95, 23), (77, 19), (1, 17),
             (11, 13), (13, 11), (15, 2), (1, 7), (55, 1)]
FRACS = FRACS_GUY
PR = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]

def fac(n):
    v = [0] * len(PR)
    for i, p in enumerate(PR):
        while n % p == 0:
            n //= p; v[i] += 1
    assert n == 1
    return v

def rules(fr):
    return [(fac(a), fac(b)) for a, b in fr]

def hdr(s):
    print()
    print("=" * 78)
    print(s)
    print("=" * 78)

def run(RULES, nprimes, budget):
    v = fac(2)
    steps = 0
    found = []
    stepsat = []
    while len(found) < nprimes and steps < budget:
        for num, den in RULES:
            if all(v[i] >= den[i] for i in range(len(PR))):
                for i in range(len(PR)):
                    v[i] += num[i] - den[i]
                break
        steps += 1
        if v[0] > 0 and all(x == 0 for x in v[1:]):
            found.append(v[0]); stepsat.append(steps)
    return found, stepsat

def L1(nprimes=15, budget=5 * 10 ** 6):
    hdr("L1  Conway's PRIMEGAME is a generalized Collatz map; the powers of 2 in the orbit of 2")
    lcm = 1
    for a, b in FRACS_GUY + FRACS_STD:
        lcm = lcm * b // math.gcd(lcm, b)
    print(f"  the fraction used depends only on n mod lcm(denominators) = {lcm}: a residue-class-wise linear map")
    primes = [p for p in range(2, 200) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
    for name, fr in (("Guy 1983 list (..., 15/14, 15/2, 55/1)", FRACS_GUY), ("standard list (..., 15/2, 1/7, 55/1)", FRACS_STD)):
        found, stepsat = run(rules(fr), nprimes, budget)
        ok = found == primes[:len(found)]
        print(f"  {name}: exponents of the powers of 2 in the orbit of 2: {found}")
        print(f"     steps: {stepsat}; equal to the first {len(found)} primes: {ok}")
        assert ok
    tw = [(a, b) for a, b in zip(found, found[1:]) if b - a == 2]
    print(f"  consecutive visits 2^p, 2^(p+2) so far (twin primes): {tw}")
    print("  Twin prime conjecture <=> infinitely many such consecutive visits (Conway's theorem, CITED).  The map is")
    print("  not 3x+1; this is the universality phenomenon behind the UNIFORM barrier, not a bridge to Collatz.")

def L2():
    hdr("L2  arithmetical complexity (elementary; used for the logic edges of the atlas)")
    rows = [
        ("no nontrivial Collatz cycle (NC)", "Pi^0_1", "each word/clock is a finite check; a counterexample is finite"),
        ("Terras sigma = tau (CST)", "Pi^0_1", "for all n, k: tau(n) = k => sigma(n) = k (sigma >= tau always)"),
        ("orbit of 8 under Collatz's permutation infinite", "Pi^0_1", "for all k >= 1: g^k(8) != 8"),
        ("Goldbach", "Pi^0_1", "bounded search for each even n"),
        ("Legendre", "Pi^0_1", "bounded search in (n^2,(n+1)^2)"),
        ("Riemann hypothesis", "Pi^0_1", "Lagarias 2002: sigma(n) <= H_n + exp(H_n) log H_n for all n"),
        ("graceful tree conjecture", "Pi^0_1", "finite search per tree"),
        ("lonely runner (all n)", "Pi^0_1", "for fixed n a finite check (Tao 2018)"),
        ("Erdos ternary problem", "Pi^0_1", "for each n > 8 inspect the ternary digits of 2^n"),
        ("Collatz conjecture", "Pi^0_2", "for all n there is k with T^k(n) = 1"),
        ("no divergence (T1)", "Pi^0_2", "for all n there are i < j with T^i(n) = T^j(n)"),
        ("Periodicity Conjecture", "Pi^0_2", "for all rationals x in Z_2 there are i < j with T^i(x) = T^j(x)"),
        ("E-SCC Q1, Q2; X_min", "Pi^0_2", "for all n there is a finite path / certificate"),
        ("HYP-9127 (2-adic irrationality)", "Pi^0_2", "for all a/b there is N with X != a/b mod 2^N"),
        ("twin prime conjecture", "Pi^0_2", "for all N there is p > N with p, p+2 prime"),
        ("abc conjecture", "Pi^0_3", "for all eps there is K with c <= K rad^(1+eps) for all triples"),
    ]
    for name, cls, why in rows:
        print(f"  {name:48s} {cls:7s} {why}")
    print("  Consequences: a Pi^0_1 statement that is independent of PA is true; each Pi^0_1 statement is")
    print("  equivalent to the non-halting of an explicit program, hence (Conway) to a statement about an explicit")
    print("  generalized Collatz map; Pi^0_2 statements reduce likewise (Kurtz-Simon).  None of these reductions")
    print("  lands on the 3x+1 map itself.")

if __name__ == "__main__":
    L1(); L2()
    print("\nDONE primegame")
