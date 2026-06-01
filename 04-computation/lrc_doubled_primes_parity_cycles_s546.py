#!/usr/bin/env python3
"""
lrc_doubled_primes_parity_cycles_s546.py    oracle-2026-06-01-S546o

Odd/even CYCLES, odd/even NUMBERS, primes, and the DOUBLED PRIME 2q.

THE PARITY DICTIONARY (additive prime basis):
  EVEN E = p + q        (Goldbach: odd + odd = even).
  ODD  O = p + 2q       (Lemoine/Levy: odd + EVEN = odd) -- the DOUBLED PRIME 2q is
                         the unique prime-structured EVEN summand that fixes parity.
So 2q is a PARITY HINGE: even numbers are covered by primes alone; odd numbers need
the doubled primes. The doubled primes are the PARITY-COMPLETION of the prime basis.

THE FACTOR OF 2 = the even prime = the DOUBLING / ANTIPODAL involution. It is the same
operator in three places:
  (numbers)  doubled prime 2q = parity hinge;
  (cycles)   even-length cycle = an antipodal-paired ('doubled') odd structure;
             Burnside for A000568: even cycles -> Fix=0; all-odd -> 2^{pairs}.
  (LRC)      n* = n/2 (even n) or n (odd n): the halving. n=2p (doubled prime) -> n*=p
             PRIME -> the CLEANEST channel structure (S533/S541). The even shadow of an
             odd prime.

This script: (1) Goldbach + Lemoine verification & rep-counts (the doubled-prime
summand); (2) A000568 via Burnside, showing EVEN cycles kill (Fix=0), ALL-ODD
contribute (the cycle-parity dictionary); (3) the LRC dimension table -- n -> n* ->
prime/prime-power/composite, marking DOUBLED-PRIME n=2p (clean), and why n=14 is the
'easiest' open even case while n=16,18 are 2x(prime power).
"""
from sympy import isprime, factorint
from itertools import permutations
from math import gcd, factorial
from collections import Counter

# ---------- (1) Goldbach + Lemoine ----------
def goldbach_reps(E):
    return [(p, E-p) for p in range(2, E//2+1) if isprime(p) and isprime(E-p)]
def lemoine_reps(O):
    # O = p + 2q, p,q prime
    return [(O-2*q, q) for q in range(2, O//2+1) if isprime(q) and isprime(O-2*q) and O-2*q>0]

def part1():
    print("="*70); print("(1) EVEN = p+q (Goldbach); ODD = p+2q (Lemoine, 2q=DOUBLED PRIME)"); print("="*70)
    print("  even E : #Goldbach (p+q)        odd O : #Lemoine (p+2q), the doubled prime 2q")
    for k in range(4, 30):
        if k % 2 == 0:
            r = goldbach_reps(k)
            print(f"    E={k:2d}: {len(r):2d}  e.g. {r[0] if r else None}")
        else:
            if k >= 7:
                r = lemoine_reps(k)
                ex = (r[0][0], 2*r[0][1]) if r else None
                print(f"    O={k:2d}: {len(r):2d}  e.g. p+2q = {ex[0]}+{ex[1]}" if ex else f"    O={k:2d}: 0")
    print("  => even covered by 2 odd primes; odd needs ONE doubled prime 2q (the even,")
    print("     parity-flipping summand). Doubled primes = parity-completion of the basis.")
    print()

# ---------- (2) A000568 Burnside: even cycles kill, all-odd contribute ----------
def cycle_types(n):
    def parts(n, mx):
        if n == 0: yield (); return
        for k in range(min(n, mx), 0, -1):
            for r in parts(n-k, k): yield (k,)+r
    yield from parts(n, n)
def perms_with_type(part):
    n = sum(part); c = Counter(part); r = factorial(n)
    for l, cnt in c.items(): r //= (l**cnt * factorial(cnt))
    return r
def fix_tournaments(part):
    """# tournaments fixed by a permutation of this cycle type. 0 if any even cycle;
    else 2^{(number of pair-orbits)}. Pair-orbits for all-odd cycle type:
    within a cycle of length l: (l-1)/2 orbits; between cycles a,b: gcd(a,b)."""
    if any(l % 2 == 0 for l in part):
        return 0
    pairs = 0
    L = list(part)
    for l in L: pairs += (l-1)//2
    for i in range(len(L)):
        for j in range(i+1, len(L)):
            pairs += gcd(L[i], L[j])
    return 2**pairs
def A000568(n):
    tot = 0
    for part in cycle_types(n):
        tot += perms_with_type(part) * fix_tournaments(part)
    return tot // factorial(n)

def part2():
    print("="*70); print("(2) A000568 via Burnside: EVEN cycles -> Fix=0; ALL-ODD -> 2^{pairs}"); print("="*70)
    print("  n : A000568(n)   (computed: only ALL-ODD cycle types contribute; even cycles kill)")
    for n in range(3, 9):
        print(f"    {n} : {A000568(n)}")
    print("  => the iso-class count is an 'all-odd-cycle' sum; the DOUBLING (an even cycle =")
    print("     antipodal-paired structure) annihilates Fix. Cycle-parity = the number-parity")
    print("     dictionary: odd cycles are the free/atomic contributors (like the odd primes).")
    print()

# ---------- (3) LRC dimension table ----------
def nstar(n): return n//2 if n % 2 == 0 else n
def classify(m):
    f = factorint(m)
    if m == 1: return "unit"
    if len(f) == 1:
        p, a = next(iter(f.items()))
        return "PRIME" if a == 1 else f"prime-power {p}^{a}"
    return "composite (CRT)"

def part3():
    print("="*70); print("(3) LRC DIMENSION TABLE: n -> n* -> channel type; DOUBLED PRIMES n=2p are clean"); print("="*70)
    print("  n  parity   n*   n* type              n structure        clean prime channels?")
    for n in range(3, 25):
        ns = nstar(n); typ = classify(ns)
        if n % 2 == 1:
            struct = "odd " + ("prime" if isprime(n) else f"comp {dict(factorint(n))}")
        else:
            half = n//2
            struct = f"2*{half} " + ("(DOUBLED PRIME)" if isprime(half) else f"(2*{dict(factorint(half))})")
        clean = "YES" if ("PRIME" in typ and "power" not in typ) else "no"
        mark = "  <<< open" if n in (14,16,18) else ""
        print(f"  {n:2d}  {'even' if n%2==0 else 'odd ':4s}  {ns:2d}   {typ:18s}  {struct:20s}  {clean}{mark}")
    print()
    print("  => clean PRIME channel modulus n* occurs iff n is an ODD PRIME or a DOUBLED PRIME")
    print("     n=2p. Among the open even cases: n=14=2*7 is a DOUBLED PRIME (n*=7 PRIME, clean,")
    print("     QR/Paley channels); n=16=2*8=2^4 and n=18=2*9=2*3^2 are 2x(prime power) ->")
    print("     n* a prime power -> FILTERED channels (S534, the messy ones). The doubled-prime")
    print("     structure of n controls the LRC channel cleanliness: n=14 is the easiest open even n.")

def main():
    part1(); part2(); part3()
    print()
    print("="*70)
    print("SYNTHESIS: the factor of 2 (even prime) = DOUBLING = ANTIPODAL involution = the")
    print("universal PARITY operator. Doubled prime 2q = parity hinge (Lemoine, completes the")
    print("additive basis for ODD numbers). Even cycle = doubled odd cycle (antipodal) -> kills")
    print("Burnside Fix (A000568 = all-odd sum). LRC: n=2p doubled-prime dimensions inherit the")
    print("clean PRIME channel structure of p (n*=p); n=14 (=2*7) is the clean open even case.")
    print("ODD = atomic/free/prime; EVEN = doubled/antipodal; 2q & n=2p carry odd structure across.")
    print("="*70)

if __name__ == "__main__":
    main()
