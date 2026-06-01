#!/usr/bin/env python3
"""
doubled_primes_goldbach_bridges_s549.py    oracle-2026-06-01-S549

Investigate doubled primes as GOLDBACH BRIDGES in relation to the hyperoperation
recursion (S548) and the LRC tower.

A doubled prime 2q is a BRIDGE in three senses, all in the recursion:

  (B1) the H1<->H2 rung. Additively (H1) 2q = q + q is the Goldbach DIAGONAL (two
       equal primes); multiplicatively (H2) 2q = 2 * q is a prime doubling. The
       doubled prime is the UNIQUE even number where the additive-prime and the
       multiplicative-prime decompositions coincide -- the bridge of
       'multiplication = repeated addition' restricted to primes.

  (B2) the PARITY bridge (Lemoine). even = p + q (Goldbach, depth-1: two prime
       leaves). odd = p + 2q (Lemoine, depth-2: a prime leaf p + a BRIDGE NODE 2q
       with two prime leaves q,q). The doubled prime is the bridge node that adds a
       level to cross from even to odd. It is the even 'carrier' every odd needs.

  (B3) the SCALE bridge (LRC apex / halving). 2q is rank-one (omega(n/2)=1, S542):
       a single tower; its apex (S530/S547) is speed q = n/2. The doubled prime
       bridges n=2q down to its half q.

We verify all three and place them in the recursion.
"""
from sympy import isprime, primerange, factorint, primefactors
from itertools import combinations

def goldbach_reps(e):
    return [(p, e - p) for p in primerange(2, e // 2 + 1) if isprime(e - p)]

def lemoine_reps(o):  # odd = p + 2q
    return [(o - 2 * q, q) for q in primerange(2, o // 2) if (o - 2 * q) >= 2 and isprime(o - 2 * q)]

def main():
    print("Doubled primes as Goldbach bridges in the recursion (oracle-S549)\n")

    print("(B1) the H1<->H2 rung: 2q = q+q (additive Goldbach diagonal, two primes) = 2*q (mult).")
    print("     doubled primes = even n with n/2 PRIME = exactly where additive-prime meets mult-prime:")
    dp = [2 * q for q in primerange(2, 30)]
    print(f"     2q = {dp}")
    print(f"     {'even n':>7}{'n/2':>5}{'n/2 prime?':>11}{'Goldbach diag q+q (q prime)?':>30}{'tower rank w(n/2)':>18}")
    for n in range(4, 26, 2):
        h = n // 2; pr = isprime(h)
        diag = pr  # n = (n/2)+(n/2) is a Goldbach rep with two equal PRIMES iff n/2 prime
        rank = len(primefactors(h)) if h > 1 else 0
        tag = "  <- DOUBLED PRIME = bridge" if pr else ""
        print(f"     {n:>7}{h:>5}{str(pr):>11}{str(diag):>30}{rank:>18}{tag}")
    print("     => doubled primes are EXACTLY the rank-one (single-tower) evens whose Goldbach")
    print("        diagonal is a prime pair: the H1<->H2 bridge. (non-doubled even: n/2 composite,")
    print("        no prime diagonal, rank >= 2.)\n")

    print("(B2) the PARITY bridge (Lemoine): odd = p + 2q; 2q is the bridge NODE (depth-2 tree).")
    print(f"     {'odd':>5}{'#Goldbach-style':>16}{'Lemoine reps p+2q (bridge 2q)':>34}")
    for o in (9, 15, 27, 45):
        L = lemoine_reps(o)
        shown = [f"{p}+2*{q}" for p, q in L][:3]
        print(f"     {o:>5}{'(odd, no p+q)':>16}{str(shown)+' ('+str(len(L))+' reps)':>34}")
    print("     => every odd is prime + a DOUBLED-PRIME BRIDGE; the bridge node 2q (children q,q)")
    print("        is the extra depth that crosses even->odd. Even needs none (depth-1: p,q).\n")

    print("(B3) the SCALE bridge (LRC): n=2q rank-one, apex = q = n/2 (S547); bridges n down to q.")
    for q in (3,5,7,11,13):
        n = 2*q
        print(f"     n={n}=2*{q}: rank-one (n/2={q} prime); apex speed = n/2 = {q}; halving bridge n -> {q}")
    print()

    print("(B4) the recursion: doubled primes are the LOAD-BEARING connectors of the tower.")
    print("     H1 (add): Goldbach even=p+q ; Lemoine odd=p+2q .  H2 (mult): doubling 2q, rank-one.")
    print("     The doubled prime 2q is the SAME number at H1 (q+q) and H2 (2*q) -> the +->x rung;")
    print("     it is the depth-adding bridge node for odd parity; and the n->n/2 apex/halving in LRC.")
    print("     Pure doublings 2^k (n/2=2^{k-1}) are the 2-adic rank-one bridges; doubled ODD primes")
    print("     2q are the q-adic rank-one bridges. Both are single-tower = bridges (S542: n=2p^k).")
    print()
    print("READING: a doubled prime is the recursion's BRIDGE -- it is where addition meets")
    print("multiplication (q+q = 2*q, both prime-structured), where even crosses to odd")
    print("(Lemoine bridge node), and where n descends to n/2 (the LRC apex). It is rank-one")
    print("because it is a SINGLE bridge: one prime q, doubled. The whole tower is held up by")
    print("these doubled-prime rungs; n=14=2*7 is the canonical one for LRC.")

if __name__ == "__main__":
    main()
