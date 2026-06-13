#!/usr/bin/env python3
"""pisano_s116l.py — Pisano periods and what they reveal.

The Fibonacci sequence mod p is periodic with period pi(p) (the Pisano period).
What are the Pisano periods of our fundamental primes?
And what do they tell us about the formal group?
"""
from math import sqrt, log, gcd

phi = (1+sqrt(5))/2

# Generate Fibonacci
fib = [0, 1]
for _ in range(500):
    fib.append(fib[-1] + fib[-2])

# Generate Lucas
luc = [2, 1]
for _ in range(500):
    luc.append(luc[-1] + luc[-2])

def pisano_period(m):
    """Find the Pisano period pi(m) = period of F_n mod m."""
    if m <= 1: return 1
    a, b = 0, 1
    for i in range(1, 6*m + 10):
        a, b = b, (a+b) % m
        if a == 0 and b == 1:
            return i
    return -1

def entry_point(m):
    """Find alpha(m) = smallest k > 0 with m | F_k."""
    for k in range(1, len(fib)):
        if fib[k] % m == 0:
            return k
    return -1

print()
print("  PISANO PERIODS AND THE FORMAL GROUP")
print()
print("="*70)
print()

# ============================================================
print("  I. PISANO PERIODS OF THE FUNDAMENTAL PRIMES")
print("  " + "-"*40)
print()
print(f"  {'p':>4s}  {'pi(p)':>6s}  {'alpha(p)':>8s}  {'pi/alpha':>9s}  {'pi(p) factors'}")
print("  " + "-"*55)

for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 42, 43, 47]:
    pp = pisano_period(p)
    ap = entry_point(p)
    ratio = pp // ap if ap > 0 else 0
    # Factor pp
    temp = pp
    factors = []
    d = 2
    while d*d <= temp:
        while temp % d == 0:
            factors.append(d)
            temp //= d
        d += 1
    if temp > 1: factors.append(temp)
    f_str = "*".join(str(f) for f in factors) if factors else "1"
    mark = ""
    if p in [2,3,5,7,11]: mark = " <-- fundamental"
    elif p == 42: mark = " <-- Hurwitz"
    print(f"  {p:4d}  {pp:6d}  {ap:8d}  {ratio:9d}  {f_str}{mark}")

print()
print("  The Pisano periods of {2,3,5,7,11}:")
print("  pi(2) = 3")
print("  pi(3) = 8")
print("  pi(5) = 20")
print("  pi(7) = 16")
print("  pi(11) = 10")
print()
print("  LOOK AT THAT:")
print("  pi(2) = 3 = the curvature quantum")
print("  pi(3) = 8 = the octonion dimension")
print("  pi(5) = 20 = 4*5 = period * pentagon")
print("  pi(7) = 16 = 2^4 = the sedenion dimension")
print("  pi(11) = 10 = T_4 = the reset cost = Petersen vertices")
print()

# ============================================================
print()
print("  II. THE PISANO PERIODS ENCODE THE CAYLEY-DICKSON CHAIN")
print("  " + "-"*40)
print()
print("  pi(2) = 3.  F_n mod 2 repeats every 3 steps.")
print("  pi(3) = 8 = 2^3 = dim(octonions).")
print("  pi(7) = 16 = 2^4 = dim(sedenions).")
print()
print("  For p = 2: pi = 3.")
print("  For p = 3: pi = 8 = 2^3.")
print("  For p = 7: pi = 16 = 2^4.")
print()
print("  The primes {2, 3, 7} give Pisano periods {3, 8, 16}.")
print("  8 and 16 are consecutive POWERS OF 2.")
print("  8 = dim(O), the last good algebra.")
print("  16 = dim(S), the first broken algebra.")
print()
print("  THE PISANO PERIOD OF 3 = THE OCTONION DIMENSION.")
print("  THE PISANO PERIOD OF 7 = THE SEDENION DIMENSION.")
print("  The curvature quantum and the forbidden threshold have")
print("  Pisano periods that are exactly the Cayley-Dickson dimensions!")
print()

# ============================================================
print()
print("  III. pi(11) = 10 AND THE PETERSEN GRAPH")
print("  " + "-"*40)
print()
print("  pi(11) = 10 = |V(Petersen)| = C(5,2) = T_4.")
print()
print("  The Fibonacci sequence mod 11 has period 10.")
print("  F_0, F_1, ..., F_9 mod 11:")
for k in range(10):
    print(f"  F_{k} = {fib[k]} = {fib[k] % 11} mod 11")
print()
print("  The residues: 0, 1, 1, 2, 3, 5, 8, 2, 10, 1.")
print("  Then F_10 = 55 = 0 mod 11 (since 11|55), F_11 = 89 = 1 mod 11.")
print("  Back to (0, 1). Period = 10. CHECK.")
print()
print("  The 10 residues mod 11: {0,1,1,2,3,5,8,2,10,1}.")
print("  Distinct values: {0, 1, 2, 3, 5, 8, 10}.")
print("  Missing: {4, 6, 7, 9}.")
print()
print("  The MISSING residues mod 11: {4, 6, 7, 9}.")
print("  = {4, 6, 7, 9}.")
print("  7 is missing! The forbidden number does NOT appear as F_k mod 11.")
print()
print("  The PRESENT residues: {0, 1, 2, 3, 5, 8, 10}.")
print("  These are the Fibonacci numbers mod 11.")
print("  The NON-Fibonacci residues mod 11: {4, 6, 7, 9}.")
print()

# What are these non-Fibonacci residues?
print("  Are the non-Fibonacci residues mod 11 meaningful?")
print("  {4, 6, 7, 9}:")
print("  4 = the period of Q.")
print("  6 = the flat number.")
print("  7 = the FORBIDDEN number.")
print("  9 = 3^2 = first H past the wall.")
print()
print("  The non-Fibonacci residues mod 11 ARE {4, 6, 7, 9}.")
print("  = {period, flat, forbidden, first-past-wall}.")
print("  These are EXACTLY the numbers from the 5-6-7 bridge!")
print()

# ============================================================
print()
print("  IV. THE ENTRY POINTS")
print("  " + "-"*40)
print()
print("  alpha(p) = smallest k with p | F_k.")
print("  alpha(2) = 3 (F_3 = 2). alpha(3) = 4 (F_4 = 3).")
print("  alpha(5) = 5 (F_5 = 5). alpha(7) = 8 (F_8 = 21 = 3*7, so 7|F_8).")
print("  alpha(11) = 10 (F_10 = 55 = 5*11).")
print()
print("  alpha(7) = 8. The forbidden prime enters Fibonacci at index 8 = dim(O)!")
print("  alpha(11) = 10. The Paley prime enters at 10 = T_4 = Petersen.")
print()
print("  And alpha(21) = ?")
a21 = entry_point(21)
print(f"  alpha(21) = {a21}. F_{a21} = {fib[a21]}. {fib[a21]}/21 = {fib[a21]//21}.")
print(f"  F_8 = 21 itself! So alpha(21) = 8.")
print(f"  The forbidden number 21 IS F_8. Its entry point is its own index.")
print()
print(f"  alpha(42) = ?")
a42 = entry_point(42)
print(f"  alpha(42) = {a42}. F_{a42} = {fib[a42]}. {fib[a42]}/42 = {fib[a42]//42}.")
print()

# ============================================================
print()
print("  V. THE LUCAS NUMBERS AT THE FUNDAMENTAL INDICES")
print("  " + "-"*40)
print()
print("  The entry points alpha give us KEY indices.")
print("  What are the LUCAS numbers at these indices?")
print()
print(f"  {'p':>4s}  {'alpha(p)':>8s}  {'F_alpha':>10s}  {'L_alpha':>10s}  {'L_alpha prime?':>14s}")
print("  " + "-"*55)
for p in [2, 3, 5, 7, 11, 13, 21, 42]:
    ap = entry_point(p)
    fa = fib[ap]
    la = luc[ap]
    # Check if la is prime
    la_prime = True
    if la < 2: la_prime = False
    else:
        for d in range(2, int(sqrt(la))+1):
            if la % d == 0:
                la_prime = False
                break
    print(f"  {p:4d}  {ap:8d}  {fa:10d}  {la:10d}  {'YES' if la_prime else 'no':>14s}")

print()
print("  L_3 = 4. L_4 = 7. L_5 = 11. L_8 = 47. L_10 = 123.")
print()
print("  At alpha(2)=3: L_3 = 4 = the Cayley period!")
print("  At alpha(3)=4: L_4 = 7 = the FORBIDDEN number!")
print("  At alpha(5)=5: L_5 = 11 = the Paley prime! (Wait, we said L_5=7 earlier?)")
print()
# Check
print(f"  L_5 = {luc[5]}. Actually L_0=2, L_1=1, L_2=3, L_3=4, L_4=7, L_5=11.")
print(f"  With 0-indexing: L_5 = {luc[5]}. YES, L_5 = 11.")
print()
print("  CORRECTION from earlier: I was using 1-indexed Lucas before.")
print("  With 0-indexing (standard): L_0=2, L_1=1, L_2=3, L_3=4, L_4=7, L_5=11.")
print()
print("  NOW THE PATTERN IS CRYSTAL CLEAR:")
print("  The prime p ENTERS Fibonacci at index alpha(p).")
print("  The Lucas number at that index GIVES THE NEXT FUNDAMENTAL NUMBER:")
print()
print("  p=2: enters at alpha=3. L_3 = 4 = the PERIOD.")
print("  p=3: enters at alpha=4. L_4 = 7 = the FORBIDDEN.")
print("  p=5: enters at alpha=5. L_5 = 11 = the PALEY PRIME.")
print("  p=7: enters at alpha=8. L_8 = 47 (prime).")
print("  p=11: enters at alpha=10. L_10 = 123 = 3*41.")
print()
print("  THE LUCAS NUMBER AT THE FIBONACCI ENTRY POINT OF p")
print("  IS THE 'NEXT LEVEL' OF THE THEORY!")
print("  2 -> 4 (period). 3 -> 7 (forbidden). 5 -> 11 (Paley).")
print()
print("  The SUCCESSOR in our framework:")
print("  2 --Lucas-at-entry--> 4")
print("  3 --Lucas-at-entry--> 7")
print("  5 --Lucas-at-entry--> 11")
print()
print("  This is a MAP from fundamental primes to fundamental numbers!")
print("  L_{alpha(p)} maps the theory forward one level.")
print()

# ============================================================
print()
print("  VI. THE ENTRY-LUCAS MAP")
print("  " + "-"*40)
print()
print("  Define EL(p) = L_{alpha(p)} = the Lucas number at the Fibonacci entry of p.")
print()
print(f"  {'p':>4s}  {'alpha(p)':>8s}  {'EL(p)':>8s}  {'meaning'}")
print("  " + "-"*50)
for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
    ap = entry_point(p)
    el = luc[ap]
    meaning = ""
    if p == 2: meaning = "4 = Cayley period"
    elif p == 3: meaning = "7 = FORBIDDEN"
    elif p == 5: meaning = "11 = Paley prime"
    elif p == 7: meaning = "47 = prime"
    elif p == 11: meaning = "123 = 3*41"
    elif p == 13: meaning = f"{el}"
    print(f"  {p:4d}  {ap:8d}  {el:8d}  {meaning}")

print()
print("  The chain: 2 -> 4 -> ... but 4 is not prime.")
print("  If we follow only PRIMES: 2 -> 3 -> 7 -> 47 -> ...")
print()
print("  Is there a chain?")
print("  Start with 2. alpha(2)=3. L_3=4. Not prime, but contains F entry alpha(3)=4.")
print("  Start with 3. alpha(3)=4. L_4=7. PRIME. Forbidden.")
print("  Start with 7. alpha(7)=8. L_8=47. PRIME.")
print("  Start with 47. alpha(47)=?")
a47 = entry_point(47)
print(f"  alpha(47) = {a47}. L_{a47} = {luc[a47]}.")
el47 = luc[a47]
# Is it prime?
el47_prime = True
if el47 < 2: el47_prime = False
else:
    for d in range(2, int(sqrt(el47))+1):
        if el47 % d == 0:
            el47_prime = False
            break
print(f"  L_{a47} = {el47}. Prime? {el47_prime}")
print()
print("  The chain through primes: 3 -> 7 -> 47 -> ...")
if el47_prime:
    a_next = entry_point(el47)
    print(f"  Next: alpha({el47}) = {a_next}. L_{a_next} = {luc[a_next]}.")

print()
print("  WHETHER OR NOT this chain continues,")
print("  the first three steps are remarkable:")
print("  3 -> 7 -> 47")
print("  curvature -> forbidden -> a prime that is L_8 = Lucas at the octonion index.")
print()

# ============================================================
print()
print("  VII. F_n * L_n = F_{2n}")
print("  " + "-"*40)
print()
print("  This identity means:")
print("  F_8 * L_8 = 21 * 47 = 987 = F_16.")
print(f"  Verify: F_16 = {fib[16]}. 21*47 = {21*47}. Match: {fib[16] == 21*47}")
print()
print("  The forbidden number (F_8=21) times the next-chain prime (L_8=47)")
print("  equals F_16 = 987.")
print(f"  987 = 3 * 7 * 47. Factor: {3}*{7}*{47}.")
print(f"  = curvature * forbidden * next-chain-prime.")
print()
print("  F_{2n} = F_n * L_n means DOUBLING the Fibonacci index")
print("  = MULTIPLYING the Fibonacci value by its Lucas companion.")
print()
print("  In rapidity: rapidity(F_{2n}) = rapidity(F_n) + rapidity(L_n).")
print("  Doubling the index ADDS the Lucas rapidity to the Fibonacci rapidity.")
print("  This is the FORMAL GROUP at work:")
print("  the index-doubling IS velocity composition in the golden helix.")
print()

# ============================================================
print()
print("  VIII. THE PISANO PERIODS AS CAYLEY-DICKSON DIMENSIONS")
print("  " + "-"*40)
print()
print("  Summary of what we found:")
print()
print("  pi(2) = 3 = curvature quantum")
print("  pi(3) = 8 = octonion dimension = 2^3")
print("  pi(5) = 20 = 4*5 = period * pentagon")
print("  pi(7) = 16 = sedenion dimension = 2^4")
print("  pi(11) = 10 = Petersen vertices = T_4 = reset cost")
print()
print("  The Pisano periods of the fundamental primes ARE the")
print("  structural constants of the Cayley-Dickson chain,")
print("  the Petersen graph, and the self-reference hierarchy.")
print()
print("  Fibonacci 'knows' about our entire theory through its mod-p behavior.")
print("  The periodicity of F_n mod p encodes the DIMENSION and STRUCTURE")
print("  that p governs in the Cayley framework.")
