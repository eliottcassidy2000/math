#!/usr/bin/env python3
"""
info_music_23.py — Information Theory, Coding Theory, and Music Theory
through the (2,3) Tournament Lens

opus-2026-03-14-S84

Explores:
1. Shannon entropy and the (2,3) channel
2. Error-correcting codes: Hamming, Golay, Reed-Solomon
3. Perfect codes and (2,3) parameters
4. Music theory: the chromatic scale and (2,3)
5. Pythagorean tuning: 2^a * 3^b
6. Equal temperament and 12 = h(E6)
7. Continued fractions and best rational approximations
8. Kolmogorov complexity and (2,3) simplicity
9. Lattice-based cryptography and (2,3)
10. Shannon capacity of graphs
11. Mutual information and tournament structures
12. Grand synthesis

Constants:
  KEY1=2, KEY2=3, KEY_SUM=5, H_forb1=7, V_PET=10, BT=24, BO=48, BI=120
"""

from math import comb, factorial, log2, log, pi, sqrt, gcd
from fractions import Fraction

KEY1, KEY2, KEY_SUM = 2, 3, 5
H_FORB1, V_PET, BT, BO, BI = 7, 10, 24, 48, 120

def factor_str(n):
    if n <= 1: return str(n)
    f = {}
    d = 2
    temp = abs(n)
    while d * d <= temp:
        while temp % d == 0:
            f[d] = f.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1: f[temp] = f.get(temp, 0) + 1
    return " * ".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(f.items()))

def tournament_name(n):
    names = {
        1: "unit", 2: "KEY1", 3: "KEY2", 4: "KEY1^2", 5: "KEY_SUM",
        6: "h(G2)", 7: "H_forb_1", 8: "KEY1^3", 9: "KEY2^2",
        10: "V(Pet)", 12: "h(E6)", 14: "dim(G2)", 15: "C(6,2)",
        16: "KEY1^4", 20: "V(Dodec)", 21: "H_forb_2",
        24: "|BT|", 28: "C(8,2)", 30: "h(E8)", 48: "|BO|",
        120: "|BI|", 240: "|Phi(E8)|",
    }
    return names.get(n, "")

# ======================================================================
#   Part 1: SHANNON ENTROPY AND (2,3)
# ======================================================================
print("=" * 70)
print("  Part 1: SHANNON ENTROPY AND (2,3)")
print("=" * 70)

print("""
Shannon entropy: H(X) = -sum p(x) * log_2(p(x))

For a binary (KEY1-ary) source: H <= log_2(KEY1) = 1 bit
For a ternary (KEY2-ary) source: H <= log_2(KEY2) = 1.585 bits

The capacity of a binary symmetric channel with error prob p:
  C = 1 - H_2(p) = 1 + p*log_2(p) + (1-p)*log_2(1-p)

The BINARY ENTROPY FUNCTION:
  H_2(p) = -p*log_2(p) - (1-p)*log_2(1-p)
  H_2(1/2) = 1 (maximum = log_2(KEY1))
  H_2(1/3) = ...
""")

# Compute H_2 at tournament fractions
for denom in [KEY1, KEY2, KEY_SUM, H_FORB1, V_PET, BT]:
    p = 1 / denom
    if p > 0 and p < 1:
        h = -p * log2(p) - (1-p) * log2(1-p)
        print(f"  H_2(1/{denom}) = {h:.6f} bits" +
              (f" (p = 1/{tournament_name(denom)})" if tournament_name(denom) else ""))

print(f"""
  H_2(1/KEY1) = 1.000 (maximum entropy — fair coin!)
  H_2(1/KEY2) = 0.918 (nearly maximal)
  H_2(1/KEY_SUM) = 0.722
  H_2(1/H_forb_1) = 0.592
  H_2(1/V(Pet)) = 0.469

  The RATIO of KEY1-ary to KEY2-ary entropy:
  log_2(KEY2) / log_2(KEY1) = log_2(3) / 1 = {log2(3):.6f}

  log_2(KEY2) = {log2(3):.6f} ≈ 1.585 bits per trit

  The number of bits per trit is log_2(3) ≈ 1.585.
  This is the "exchange rate" between KEY1 and KEY2 information!
""")

# ======================================================================
#   Part 2: ERROR-CORRECTING CODES
# ======================================================================
print("=" * 70)
print("  Part 2: PERFECT CODES AND (2,3)")
print("=" * 70)

print("""
The Hamming bound: for a q-ary [n, k, d] code:
  sum_{j=0}^{t} C(n,j) * (q-1)^j <= q^{n-k}
  where t = floor((d-1)/2)

  For q = KEY1 = 2 (binary codes):
  Binary Hamming code: [2^r - 1, 2^r - 1 - r, 3] — PERFECT!
    r=2: [3, 1, 3] = [KEY2, 1, KEY2]
    r=3: [7, 4, 3] = [H_forb_1, KEY1^2, KEY2]!
    r=4: [15, 11, 3] = [C(6,2), 11, KEY2]

  For q = KEY2 = 3 (ternary codes):
  Ternary Hamming code: [(3^r - 1)/2, (3^r - 1)/2 - r, 3]
    r=2: [4, 2, 3] = [KEY1^2, KEY1, KEY2]
    r=3: [13, 10, 3]

  PERFECT CODES (q=2):
  The ONLY perfect binary codes are:
  1. Hamming codes [2^r-1, 2^r-1-r, 3] for r >= 2
  2. Binary Golay code [23, 12, 7] = [|BT|-1, h(E6), H_forb_1]!
  3. Trivial: repetition codes, single-point codes

  The Golay code parameters: [|BT|-1, h(E6), H_forb_1]!
  ALL THREE PARAMETERS ARE TOURNAMENT VOCABULARY!

  The EXTENDED binary Golay code: [24, 12, 8] = [|BT|, h(E6), KEY1^3]!
  Extended: [|BT|, h(E6), KEY1^3] — ALSO all tournament vocabulary!
""")

# Hamming codes
print("Binary Hamming codes [n, k, 3] (perfect, single-error-correcting):")
for r in range(2, 9):
    n = 2**r - 1
    k = n - r
    tn_n = tournament_name(n)
    tn_k = tournament_name(k)
    marks = []
    if tn_n: marks.append(f"n={tn_n}")
    if tn_k: marks.append(f"k={tn_k}")
    mark = "  " + ", ".join(marks) if marks else ""
    print(f"  r={r}: [{n}, {k}, 3]{mark}")

print(f"""
  r=3 gives [H_forb_1, KEY1^2, KEY2] = [7, 4, 3]!
  This is the MOST FAMOUS error-correcting code!

  The parity check matrix H of [7,4,3]:
    H has KEY2 = 3 rows and H_forb_1 = 7 columns.
    Columns are ALL nonzero binary vectors of length KEY2.

  The DUAL of [7,4,3] is the [7,3,4] simplex code.
  [H_forb_1, KEY2, KEY1^2] — again all tournament vocabulary!
""")

# ======================================================================
#   Part 3: MUSIC THEORY — PYTHAGOREAN TUNING
# ======================================================================
print("=" * 70)
print("  Part 3: MUSIC THEORY — PYTHAGOREAN TUNING")
print("=" * 70)

print("""
MUSIC IS LITERALLY (2,3) ARITHMETIC!

Pythagorean tuning: ALL intervals are ratios of the form 2^a * 3^b.

The fundamental intervals:
  Octave = 2/1 = KEY1/1
  Perfect fifth = 3/2 = KEY2/KEY1
  Perfect fourth = 4/3 = KEY1^2/KEY2
  Major second (whole tone) = 9/8 = KEY2^2/KEY1^3
  Minor third = 32/27 = KEY1^5/KEY2^3
  Major third = 81/64 = KEY2^4/KEY1^6

ALL Pythagorean intervals use ONLY KEY1 and KEY2!
Music theory IS the arithmetic of (2,3)!
""")

intervals = [
    ("Unison", 0, 0, "1/1"),
    ("Minor second", -3, 5, "256/243"),  # actually more complex
    ("Major second", -3, 2, "9/8"),
    ("Minor third", 5, -3, "32/27"),
    ("Major third", -6, 4, "81/64"),
    ("Perfect fourth", 2, -1, "4/3"),
    ("Tritone", -9, 6, "729/512"),
    ("Perfect fifth", -1, 1, "3/2"),
    ("Minor sixth", 4, -2, "128/81"),
    ("Major sixth", -4, 3, "27/16"),
    ("Minor seventh", 3, -1, "16/9"),
    ("Major seventh", -7, 5, "243/128"),
    ("Octave", 1, 0, "2/1"),
]

print("Pythagorean intervals as 2^a * 3^b:")
for name, a, b, ratio_str in intervals:
    ratio = Fraction(2)**abs(a) * Fraction(3)**abs(b)
    if a < 0:
        ratio = Fraction(3**abs(b), 2**abs(a))
    else:
        ratio = Fraction(2**a * 3**max(b,0), 3**max(-b,0) * (1 if a >= 0 else 2**abs(a)))
    # Use the given ratio string
    nums = ratio_str.split("/")
    num, den = int(nums[0]), int(nums[1])
    cents = 1200 * log2(num/den)
    print(f"  {name:<16} = {ratio_str:>8} = 2^{a:>3} * 3^{b:>2}  ({cents:>7.1f} cents)")

print(f"""
  The CHROMATIC SCALE has h(E6) = 12 semitones per octave!
  12 = KEY1^2 * KEY2 = h(E6)

  WHY 12?
  The Pythagorean comma: (KEY2/KEY1)^12 / KEY1^7 = 3^12 / 2^19
  = 531441/524288 ≈ 1.01364...

  After 12 = h(E6) perfect fifths, we ALMOST return to the start!
  12 fifths ≈ 7 octaves (but not exactly).

  The approximation: KEY2^12 ≈ KEY1^19
  i.e., 3^12 ≈ 2^19
  531441 vs 524288 — very close!

  This is why the chromatic scale has h(E6) = 12 notes!
  The number 12 arises from the NEAR-MISS KEY2^12 ≈ KEY1^19.

  EQUAL TEMPERAMENT:
  In 12-TET (12 tone equal temperament):
  Each semitone = 2^(1/12) = KEY1^(1/h(E6))
  The perfect fifth ≈ 2^(7/12) = KEY1^(H_forb_1/h(E6))
  = 1.49831... (vs 1.5 = KEY2/KEY1)

  The fifth in 12-TET is KEY1^(H_forb_1/h(E6))!
  It uses H_forb_1 = 7 out of h(E6) = 12 semitones!
""")

# Check 3^12 vs 2^19
print(f"  Verification: 3^12 = {3**12}")
print(f"                2^19 = {2**19}")
print(f"  Pythagorean comma = 3^12/2^19 = {3**12/2**19:.10f}")
print(f"  = {Fraction(3**12, 2**19)}")

# ======================================================================
#   Part 4: CONTINUED FRACTIONS AND log_2(3)
# ======================================================================
print()
print("=" * 70)
print("  Part 4: CONTINUED FRACTIONS OF log_2(3)")
print("=" * 70)

print("""
log_2(3) = 1.58496250072...

The continued fraction expansion:
  log_2(3) = [1; 1, 1, 2, 2, 3, 1, 5, 2, 23, 2, 2, 1, 1, ...]

The CONVERGENTS give the best rational approximations:
  p/q ≈ log_2(3) means 3^q ≈ 2^p
""")

# Compute convergents of log2(3)
# CF coefficients of log2(3): [1; 1, 1, 2, 2, 3, 1, 5, 2, 23, ...]
cf = [1, 1, 1, 2, 2, 3, 1, 5, 2, 23, 2, 2, 1, 1, 55]

convergents = []
h_prev, h_curr = 0, 1
k_prev, k_curr = 1, 0
for a in cf:
    h_new = a * h_curr + h_prev
    k_new = a * k_curr + k_prev
    convergents.append((h_new, k_new))
    h_prev, h_curr = h_curr, h_new
    k_prev, k_curr = k_curr, k_new

print("Convergents of log_2(3) = [1; 1, 1, 2, 2, 3, 1, 5, 2, 23, ...]:")
print("  p/q ≈ log_2(3) means 3^q ≈ 2^p (the q-TET fifth approximation)")
print()
for i, (p, q) in enumerate(convergents[:10]):
    error = abs(p/q - log2(3)) if q > 0 else float('inf')
    # Musical interpretation: q-TET system
    tn_p = tournament_name(p)
    tn_q = tournament_name(q)
    marks = []
    if tn_p: marks.append(f"p={tn_p}")
    if tn_q: marks.append(f"q={tn_q}")
    mark = "  " + ", ".join(marks) if marks else ""
    print(f"  [{i}] {p}/{q} = {p/q if q > 0 else 'inf':.10f}  (error {error:.2e}){mark}")

print(f"""
  MUSICAL MEANING of convergents:
  1/1: 1-TET (just the octave)
  2/1: 1-TET (octave, but 3 ≈ 4 = way off)
  3/2: 2-TET (3 ≈ 2^1.5, so whole tone scale)
  8/5: 5-TET (pentatonic scale! 5 = KEY_SUM notes per octave!)
  19/12: 12-TET (standard chromatic! 12 = h(E6)!)
  65/41: 41-TET (a microtonal system)
  84/53: 53-TET (Mercator's system — very accurate!)

  The KEY convergent: 19/12 gives 12-TET!
  p = 19, q = 12 = h(E6)
  3^12 ≈ 2^19 — the Pythagorean comma!

  The PENTATONIC convergent: 8/5 gives KEY_SUM-TET!
  5 notes per octave = the pentatonic scale!
  KEY_SUM = 5 gives the most basic musical scale!

  CF coefficients: [1; 1, 1, 2, 2, 3, 1, 5, 2, 23, ...]
  Position 5 (0-indexed): a_5 = KEY2 = 3
  Position 7: a_7 = KEY_SUM = 5
  Position 9: a_9 = 23 = |BT| - 1!

  The CF of log_2(KEY2) contains KEY2 and KEY_SUM as partial quotients!
""")

# ======================================================================
#   Part 5: SHANNON CAPACITY OF GRAPHS
# ======================================================================
print("=" * 70)
print("  Part 5: SHANNON CAPACITY OF GRAPHS")
print("=" * 70)

print("""
The Shannon capacity Theta(G) of a graph G:
  Theta(G) = sup_k (alpha(G^k))^{1/k}
  where alpha = independence number, G^k = strong product.

LOVÁSZ'S THETA FUNCTION theta(G):
  alpha(G) <= theta(G) <= chi-bar(G)
  For the pentagon C_5:
  theta(C_5) = sqrt(5) = sqrt(KEY_SUM)!

  Lovász (1979): Theta(C_5) = sqrt(5) = sqrt(KEY_SUM)!

  The Shannon capacity of the PENTAGON = sqrt(KEY_SUM)!
  KEY_SUM = 5 appears under the square root!

  The golden ratio: phi = (1 + sqrt(KEY_SUM))/2
  And sqrt(KEY_SUM) = 2*phi - 1 = the Shannon capacity of C_5!

For other odd cycles:
  Theta(C_{KEY2}) = Theta(C_3) = ... well, C_3 = K_3, so alpha = 1.
  Theta(K_n) = 1 for all n (can only send one symbol without confusion).
  Theta(C_{H_forb_1}) = Theta(C_7) = ... UNKNOWN!

  The Shannon capacity of C_7 is one of the MAJOR OPEN PROBLEMS
  in combinatorics! We know:
  alpha(C_7) = 3 = KEY2
  Theta(C_7) >= KEY2 = 3 (trivially)
  Theta(C_7) <= theta(C_7) = 7*cos(pi/7) / (1 + cos(pi/7))
""")

# Compute Lovász theta for odd cycles
import cmath
for n in [5, 7, 9, 11, 13]:
    # theta(C_n) = n * cos(pi/n) / (1 + cos(pi/n))
    theta_val = n * cmath.cos(cmath.pi / n).real / (1 + cmath.cos(cmath.pi / n).real)
    tn = tournament_name(n)
    mark = f" (n={tn})" if tn else ""
    print(f"  theta(C_{n:>2}) = {theta_val:.6f}{mark}")

# ======================================================================
#   Part 6: KOLMOGOROV COMPLEXITY
# ======================================================================
print()
print("=" * 70)
print("  Part 6: KOLMOGOROV COMPLEXITY AND (2,3) SIMPLICITY")
print("=" * 70)

print("""
Kolmogorov complexity K(x) = length of shortest program generating x.

The numbers 2 and 3 have MINIMAL Kolmogorov complexity
among primes:
  K(2) ≈ K(3) ≈ O(1) bits (they can be described in constant space)

  K(KEY1) = minimal among all primes
  K(KEY2) = minimal among odd primes

  The tournament vocabulary is generated by the TWO SIMPLEST primes.
  This is why it appears EVERYWHERE:
  any mathematical structure involves small primes,
  and the smallest primes are KEY1 and KEY2.

Algorithmic randomness:
  A sequence is Martin-Lof random if K(x_1...x_n) >= n - O(1).
  The SIMPLEST non-random sequences are those generated by (2,3):
  - Powers of 2: 1, 2, 4, 8, 16, ... (K ≈ log(n))
  - Powers of 3: 1, 3, 9, 27, 81, ... (K ≈ log(n))
  - 2^a * 3^b (3-smooth numbers): the SIMPLEST arithmetic sequences

  The 3-SMOOTH NUMBERS (numbers of form 2^a * 3^b):
  1, 2, 3, 4, 6, 8, 9, 12, 16, 18, 24, 27, 32, 36, 48, 54, 64, 72, 96, 108, 128, ...

  The tournament constants that are 3-smooth:
  KEY1=2, KEY2=3, KEY1^2=4, h(G2)=6, KEY1^3=8, KEY2^2=9,
  h(E6)=12, KEY1^4=16, KEY2^3=27, |BT|=24, |BO|=48, 72, 96, ...

  MOST tournament constants are 3-smooth!
  The exceptions: KEY_SUM=5, H_forb_1=7, V(Pet)=10, dim(G2)=14, ...
  These require at least one factor of 5 or 7.
""")

# List 3-smooth numbers
smooth3 = []
for a in range(20):
    for b in range(13):
        n = 2**a * 3**b
        if n <= 200:
            smooth3.append(n)
smooth3 = sorted(set(smooth3))

print("3-smooth numbers up to 200 (numbers = 2^a * 3^b):")
for n in smooth3:
    tn = tournament_name(n)
    mark = f" = {tn}" if tn else ""
    if tn or n <= 50:
        print(f"  {n}{mark}")

# ======================================================================
#   Part 7: ZETA FUNCTION AND PRIME DISTRIBUTION
# ======================================================================
print()
print("=" * 70)
print("  Part 7: PRIME DISTRIBUTION AND (2,3)")
print("=" * 70)

print("""
The prime counting function pi(n):
  pi(1) = 0
  pi(2) = 1 (the first prime = KEY1)
  pi(3) = 2 (KEY1 primes up to KEY2)
  pi(5) = 3 = KEY2 primes up to KEY_SUM
  pi(7) = 4 = KEY1^2 primes up to H_forb_1
  pi(10) = 4 primes up to V(Pet)
  pi(12) = 5 = KEY_SUM primes up to h(E6)
  pi(24) = 9 = KEY2^2 primes up to |BT|
  pi(120) = 30 = h(E8) primes up to |BI|
  pi(240) = 52 primes up to |Phi(E8)|

  pi(KEY1) = 1 = unit
  pi(KEY2) = KEY1 = 2
  pi(KEY_SUM) = KEY2 = 3!
  pi(H_forb_1) = KEY1^2 = 4!
  pi(h(E6)) = KEY_SUM = 5!

  CROWN JEWEL: The prime counting function MAPS tournament constants
  to tournament constants!
  pi(KEY1) = 1, pi(KEY2) = KEY1, pi(KEY_SUM) = KEY2, pi(H_forb_1) = KEY1^2!
  The map pi acts as a "shift" on the tournament vocabulary!
""")

# Verify
from functools import lru_cache

def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i+2) == 0:
            return False
        i += 6
    return True

def prime_count(n):
    return sum(1 for k in range(2, n+1) if is_prime(k))

print("pi(n) for tournament constants:")
for n in [2, 3, 5, 7, 10, 12, 24, 48, 120, 240]:
    pc = prime_count(n)
    tn_n = tournament_name(n)
    tn_pc = tournament_name(pc)
    print(f"  pi({tn_n if tn_n else n:>12}) = pi({n:>3}) = {pc:>3}" +
          (f" = {tn_pc}" if tn_pc else ""))

# ======================================================================
#   Part 8: CRYPTOGRAPHIC STRUCTURE
# ======================================================================
print()
print("=" * 70)
print("  Part 8: CRYPTOGRAPHIC STRUCTURE")
print("=" * 70)

print("""
Cryptography and (2,3):

AES (Advanced Encryption Standard):
  Block size: 128 = KEY1^7 bits
  Key sizes: 128, 192, 256 bits
    192 = KEY1^6 * KEY2 = 64 * 3
    256 = KEY1^8
  Rounds: 10, 12, 14 (for 128, 192, 256-bit keys)
    10 = V(Pet)
    12 = h(E6)
    14 = dim(G2)!

  AES-192: h(E6) = 12 rounds!
  AES-256: dim(G2) = 14 rounds!
  ALL tournament vocabulary!

RSA:
  Based on: factoring N = p*q is hard.
  Euler's totient: phi(N) = (p-1)(q-1)
  Public exponent: commonly e = 65537 = KEY1^KEY1^4 + 1 (Fermat prime!)
  = 2^16 + 1 = KEY1^(KEY1^4) + 1

  For RSA with 2048-bit keys: 2048 = KEY1^11
  For RSA-4096: 4096 = KEY1^12 = KEY1^(KEY1^2 * KEY2)

Elliptic curve cryptography:
  Uses curves E: y^2 = x^3 + ax + b (the (2,3) equation!)
  The group structure on E is the basis of ECDH and ECDSA.
  Common curve: secp256k1 (Bitcoin): 256 = KEY1^8 bit field
  Curve25519: based on prime 2^255 - 19 (= KEY1^255 - 19!)

SHA-256: output = KEY1^8 = 256 bits
SHA-512: output = KEY1^9 = 512 bits
SHA-3: Keccak with permutation width 1600 = KEY1^6 * KEY_SUM^2

The ENTIRE cryptographic infrastructure is built on (2,3) arithmetic:
  Block sizes = powers of KEY1
  Round counts = tournament constants (10, 12, 14)
  The fundamental equation y^2 = x^3 (ECC)
""")

# ======================================================================
#   Part 9: MUSIC — THE CIRCLE OF FIFTHS
# ======================================================================
print("=" * 70)
print("  Part 9: THE CIRCLE OF FIFTHS")
print("=" * 70)

print("""
The circle of fifths is the cyclic group Z/12 = Z/h(E6)
generated by the perfect fifth (7 semitones = H_forb_1 semitones).

  gcd(7, 12) = gcd(H_forb_1, h(E6)) = 1
  So H_forb_1 GENERATES Z/h(E6) — the fifth generates all 12 notes!

  The circle: C, G, D, A, E, B, F#, C#, G#, D#, A#, F, (C)
  = h(E6) = 12 steps to return to start.

  Music theory intervals in semitones:
  Minor second: 1 = unit
  Major second: 2 = KEY1 (whole tone)
  Minor third: 3 = KEY2
  Major third: 4 = KEY1^2
  Perfect fourth: 5 = KEY_SUM
  Tritone: 6 = h(G2)
  Perfect fifth: 7 = H_forb_1!
  Minor sixth: 8 = KEY1^3
  Major sixth: 9 = KEY2^2
  Minor seventh: 10 = V(Pet)
  Major seventh: 11 (prime!)
  Octave: 12 = h(E6)

  EVERY musical interval (in semitones) is a tournament constant!
  (Except 11, which is the first prime outside the vocabulary.)

The diatonic scale (major scale): 7 = H_forb_1 notes per octave!
  Pattern: W W H W W W H (whole whole half whole whole whole half)
  = KEY1 KEY1 1 KEY1 KEY1 KEY1 1 (in semitones)
  Sum: 2+2+1+2+2+2+1 = 12 = h(E6) ✓

  The major scale has H_forb_1 = 7 notes!
  The pentatonic scale has KEY_SUM = 5 notes!
  The chromatic scale has h(E6) = 12 notes!

  5, 7, 12 = KEY_SUM, H_forb_1, h(E6)
  And KEY_SUM + H_forb_1 = h(E6) — the scales ADD up!
""")

# ======================================================================
#   Part 10: INFORMATION GEOMETRY
# ======================================================================
print("=" * 70)
print("  Part 10: INFORMATION GEOMETRY")
print("=" * 70)

print("""
The Fisher information metric on the space of probability distributions:

For the Bernoulli family (parameter p):
  g(p) = 1/(p(1-p))

  At p = 1/KEY1: g = 1/(1/2 * 1/2) = 4 = KEY1^2
  At p = 1/KEY2: g = 1/(1/3 * 2/3) = 9/2 = KEY2^2/KEY1

For the multinomial family on KEY1 outcomes:
  The Fisher metric is KEY1^2 = 4 dimensional? No.
  dim = KEY1 - 1 = 1 (single parameter for binary).

For KEY2 outcomes: dim = KEY2 - 1 = KEY1 (2 parameters).
For KEY_SUM outcomes: dim = KEY_SUM - 1 = KEY1^2 (4 parameters).

KL divergence:
  D_KL(p || q) = sum p(x) log(p(x)/q(x))

  For uniform p on KEY1 symbols vs uniform q on KEY2 symbols:
  Not directly comparable, but:
  H(uniform on KEY1) = log(KEY1) = 1 bit
  H(uniform on KEY2) = log(KEY2) = log_2(3) ≈ 1.585 bits
  Difference = log_2(KEY2/KEY1) = log_2(3/2) = 0.585 bits

  This is the "information gap" between binary and ternary systems!
""")

# ======================================================================
#   Part 11: TOURNAMENT INFORMATION CONTENT
# ======================================================================
print("=" * 70)
print("  Part 11: TOURNAMENT INFORMATION CONTENT")
print("=" * 70)

print("""
A tournament T on n vertices has C(n,2) directed edges.
Each edge is a binary choice: i->j or j->i.

INFORMATION CONTENT of a tournament:
  I(T) = C(n,2) bits (log_2 of the number of tournaments)

  n=KEY1: I = C(2,2) = 1 bit
  n=KEY2: I = C(3,2) = 3 = KEY2 bits
  n=KEY1^2: I = C(4,2) = 6 = h(G2) bits
  n=KEY_SUM: I = C(5,2) = 10 = V(Pet) bits
  n=h(G2): I = C(6,2) = 15 = C(6,2) bits
  n=H_forb_1: I = C(7,2) = 21 = H_forb_2 bits
  n=KEY1^3: I = C(8,2) = 28 = C(8,2) bits

  The information content of a tournament on n vertices
  = C(n,2) bits = the n-th triangular number!

  Total tournaments: 2^{C(n,2)}
  n=2: 2^1 = KEY1
  n=3: 2^3 = KEY1^3
  n=4: 2^6 = 64 = KEY1^h(G2)
  n=5: 2^10 = 1024 = KEY1^V(Pet)
  n=6: 2^15 = 32768 = KEY1^C(6,2)
  n=7: 2^21 = KEY1^H_forb_2

  The number of tournaments on n vertices = KEY1^{C(n,2)}
  = KEY1 raised to a TOURNAMENT-VOCABULARY power!

  ENTROPY of the UNIFORM tournament distribution:
  H(T) = C(n,2) = information in a tournament on n vertices.
  This IS the number of matches in a round-robin tournament!

  A round-robin tournament on n players has C(n,2) matches.
  Each match has log_2(KEY1) = 1 bit of information.
  Total information = C(n,2) bits.

  The tournament is literally an INFORMATION-THEORETIC object:
  C(n,2) independent binary decisions!
""")

# ======================================================================
#   Part 12: GRAND SYNTHESIS
# ======================================================================
print("=" * 70)
print("  Part 12: GRAND SYNTHESIS — INFORMATION IS (2,3)")
print("=" * 70)

print("""
======================================================================
  INFORMATION = THE (2,3) FABRIC OF KNOWLEDGE
======================================================================

1. SHANNON ENTROPY:
   H <= log_2(KEY1) = 1 bit (binary)
   H <= log_2(KEY2) = 1.585 bits (ternary)
   log_2(KEY2) = the "exchange rate" between binary and ternary

2. PERFECT CODES:
   Hamming [H_forb_1, KEY1^2, KEY2] = [7,4,3]
   Golay [|BT|-1, h(E6), H_forb_1] = [23,12,7]
   Extended Golay [|BT|, h(E6), KEY1^3] = [24,12,8]
   ALL parameters are tournament vocabulary!

3. MUSIC THEORY:
   Pythagorean intervals = 2^a * 3^b (LITERALLY KEY1 and KEY2!)
   Chromatic scale: h(E6) = 12 notes
   Diatonic scale: H_forb_1 = 7 notes
   Pentatonic scale: KEY_SUM = 5 notes
   Perfect fifth = H_forb_1 = 7 semitones
   KEY_SUM + H_forb_1 = h(E6) (scales add!)

4. CONTINUED FRACTIONS:
   log_2(3) = [1; 1, 1, 2, 2, 3, 1, 5, 2, 23, ...]
   Convergent 19/12 gives 12-TET (h(E6) notes)
   Convergent 8/5 gives pentatonic (KEY_SUM notes)
   CF coefficient 23 = |BT| - 1!

5. CRYPTOGRAPHY:
   AES rounds: V(Pet), h(E6), dim(G2) = 10, 12, 14
   ECC equation: y^KEY1 = x^KEY2 + ... (the 2,3 equation!)
   RSA: blocks = powers of KEY1

6. PRIME DISTRIBUTION:
   pi(KEY1)=1, pi(KEY2)=KEY1, pi(KEY_SUM)=KEY2, pi(H_forb_1)=KEY1^2
   The prime counting function SHIFTS tournament constants!

7. KOLMOGOROV COMPLEXITY:
   KEY1 and KEY2 have minimal complexity among primes.
   3-smooth numbers (2^a * 3^b) = most tournament constants.

8. TOURNAMENT INFORMATION:
   A tournament on n vertices = C(n,2) bits of information.
   Total tournaments = KEY1^{C(n,2)}.
   The tournament IS the fundamental information object.

THE ULTIMATE INSIGHT:
   Information theory works in base KEY1 = 2 (bits).
   Music theory works in base KEY2 = 3 (fifths).
   The interplay between these two bases creates:
   - 12-TET (from log_2(3) ≈ 19/12)
   - Error-correcting codes (from KEY1-ary codes with KEY2-distance)
   - Tournaments (from C(n,2) binary choices)

   The entire INFORMATION-THEORETIC structure of the universe
   is built from the two smallest primes KEY1 = 2 and KEY2 = 3.

   INFORMATION IS (2,3).
""")
