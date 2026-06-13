---
id: HYP-2044
status: PARTIALLY-TRUE
source: oracle-2026-06-01-S546
related:
  - HYP-2035
  - HYP-1994
  - HYP-2008
---

# HYP-2044: doubling is pairing -- the doubled primes 2q are the resonant diagonal across odd-cycles, Goldbach, and LRC

**Web (verified, `doubled_primes_odd_cycles_goldbach_s546.py`):**
1. Parity = parity of #odd-parts: even n = p+q (Goldbach, 2 odd primes); odd n = p+2q
   (Lemoine/Levy, prime + DOUBLED prime = 3 odd primes p,q,q). doubled prime 2q = q+q
   = the DIAGONAL of Goldbach (= 4,6,10,14,22,26,...).
2. A000568 = odd-cycle Burnside (verified n=3..8): (1/n!) sum_{all-odd lambda}
   |class| 2^{e}, e=sum(l_i-1)/2 + sum_{i<j} gcd(l_i,l_j); Fix=0 if ANY even cycle
   (even cycle reverses an antipodal edge). Even cycle = 0 = "even is not a single prime".
3. DOUBLED PRIME = equal odd-cycle pair (q,q): contributes gcd(q,q)=q = the MAXIMAL
   between-cycle term (distinct coprime -> 1). Verified (3,3)->3, (5,5)->5, (7,7)->7.
   So the doubled prime is the resonant diagonal where the Burnside symmetry concentrates.
4. LRC: n=2q (doubled prime) = rank-one p-adic tree (omega(n/2)=1, S542), the cleanest
   hard case, apex speed = n/2 = q. n=14=2*7 is a doubled prime -> the canonical target.

**PRINCIPLE:** doubling (x2) is PAIRING -- q -> 2q turns a single odd prime into a
resonant pair (q,q); equal things maximize gcd, make Goldbach diagonal, make the LRC
tower rank-one. The lone 2 (only even prime) is the parity-fixer: it converts the
even 2-object story into the odd 3-object story (p+q+q) and halves n to the apex.

**OPEN:** (A) does the doubled-prime cycle type (q,q) control the loneliest LRC
configuration at n=2q (apex speed = repeated cycle length q)? (B) Lemoine reps
(odd = p+2q) as a 'one odd cycle + one doubled cycle' decomposition -- do they track
an odd-n tournament invariant (vs Goldbach for even n)?

**Files:** `04-computation/doubled_primes_odd_cycles_goldbach_s546.py` (+.out).
Reflection: `07-reflections/doubled-primes-doubling-is-pairing-s546.md`.
