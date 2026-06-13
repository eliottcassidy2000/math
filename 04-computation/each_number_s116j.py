#!/usr/bin/env python3
"""each_number_s116j.py — Getting to know each natural number personally.

Not their properties. Their CHARACTER.
Who they are in the Cayley world.
"""
from math import log, sqrt, pi, exp, atanh, gcd, factorial
from fractions import Fraction
import cmath

phi = (1+sqrt(5))/2

def factorize(n):
    if n <= 1: return []
    f = []
    d = 2
    while d*d <= n:
        while n % d == 0:
            f.append(d)
            n //= d
        d += 1
    if n > 1: f.append(n)
    return f

def is_prime(n):
    return n > 1 and factorize(n) == [n]

print()
print("  GETTING TO KNOW EACH NUMBER")
print()
print("="*70)
print()

numbers = {
    0: {
        "address": "Q^{-1}(0) = -1 (the left boundary)",
        "rapidity": "-infinity",
        "character": [
            "The boundary. The left pole. The place where Q sends -1.",
            "Not nothing — the EDGE of nothing. The cliff before the void.",
            "Its Cayley address -1 is as far from center as +1 (infinity),",
            "but in the opposite direction. 0 and infinity are TWINS.",
            "Q(-x) = 1/Q(x): zero IS the reciprocal of infinity.",
            "In music: silence. Not the absence of sound — the BOUNDARY of sound.",
            "In probability: impossibility. The event that never happens.",
            "In physics: absolute zero. Unreachable but approachable.",
        ],
    },
    1: {
        "address": "0/2 = 0",
        "rapidity": "0",
        "pair": "(0, 1)",
        "character": [
            "The identity. The center. Rapidity zero. Address zero.",
            "Q(0) = 1. The Cayley transform of nothing is unity.",
            "The pair (0, 1) — nothing and something — sums to 1.",
            "The ONLY number that is neither comparison nor duplication.",
            "It is the REFERENCE. Everything else is measured from 1.",
            "1! = 1. 1^n = 1. log(1) = 0. The fixed point of everything.",
            "In music: unison. Two notes at the same pitch.",
            "In tournaments: the transitive tournament. H = 1. One path. ORDER.",
            "In probability: certainty of indifference. The fair coin, from the other side.",
        ],
    },
    2: {
        "address": "1/3",
        "rapidity": f"{log(2)/2:.6f} = ln(2)/2 = one octave",
        "pair": "(1, 2) — but 2 is even, no consecutive pair in the odd sense",
        "character": [
            "The FIRST. The only even prime. The atom of doubling.",
            "Address 1/3: the simplest non-trivial Cayley point.",
            "Q(1/3) = 2. One-third of the way from center to pole gives DOUBLING.",
            "Rapidity ln(2)/2 = the OCTAVE. The fundamental unit of rapidity.",
            "Every 'bit' is one copy of 2. Every doubling is one octave.",
            "2 is the base of binary, the base of computation, the base of life (cell division).",
            "The Cayley orbit of 1/3: {1/3, 2, -3, -1/2}.",
            "2's shadow is -3. 2's prerequisite is -1/2.",
            "The ONLY number whose shadow (-3) is the negative of the next prime.",
            "In music: the octave. The most consonant interval. The IDENTITY of pitch.",
            "In the formal group: the generator. Everything is built from 2.",
        ],
    },
    3: {
        "address": "1/2",
        "rapidity": f"{log(3)/2:.6f} = ln(3)/2",
        "pair": "(1, 2)",
        "character": [
            "The CURVATURE QUANTUM. The smallest odd cycle.",
            "Address 1/2: exactly HALF of the Cayley line.",
            "Q(1/2) = 3. Halfway to the pole gives the first odd prime.",
            "3 = 1 + 2. The comparison between nothing-something (1) and the first doubling (2).",
            "Rapidity ln(3)/2 = one TWELFTH in music (octave + fifth).",
            "The C-H stretch / C-C stretch ratio = 3. Chemistry knows 3.",
            "In tournaments: the 3-cycle. Rock-paper-scissors. The simplest contradiction.",
            "3 is the number of vertices of the triangle, the simplest polygon.",
            "The triangle IS the curvature. Flat space has no triangles wrapping around.",
            "In the formal group: F_h(1/3, 1/3) = (2/3)/(1+1/9) = (2/3)/(10/9) = 3/5 = addr_5.",
            "Composing 1/3 with itself gives the address of 5. TWO is THREE in the formal group.",
            "3 imaginary dimensions = quaternions. Where commutativity breaks.",
        ],
    },
    4: {
        "address": "3/5",
        "rapidity": f"{log(4)/2:.6f} = ln(4)/2 = ln(2) = TWO octaves",
        "pair": "(2, 2) — duplication",
        "character": [
            "The first SQUARE. 2^2. The duplication of the doubler.",
            "Address 3/5: three-fifths of the way to the pole.",
            "Rapidity = ln(2) = the Doppler octave velocity (v = 3/5 gives D = 2).",
            "4 = 2+2. Pure duplication. No comparison. Symmetric.",
            "4 = |Q^n = identity|. The PERIOD of the Cayley transform.",
            "Q^4 = Id. Four steps around the orbit. Four normed division algebras.",
            "4 dimensions = quaternions. 4 faces of every orbit.",
            "In the Cayley-Dickson chain: 4 = dim(H), the quaternion algebra.",
            "In music: Q^{-1}(4) = 3/5. The address of 4 IS the Doppler octave velocity.",
        ],
    },
    5: {
        "address": "2/3",
        "rapidity": f"{log(5)/2:.6f}",
        "pair": "(2, 3)",
        "character": [
            "The PENTAGON. Five vertices. Five Platonic solids. Five Fermat primes (known).",
            "5 = 2 + 3. The comparison between the doubler and the curvature quantum.",
            "The pair (2, 3): the two smallest primes, side by side.",
            "The ONLY number whose Cayley pair consists of two primes (consecutive primes 2,3).",
            "Address 2/3. Q(2/3) = 5. Two-thirds of the way to the pole.",
            "phi^2 = phi + 1 ~ 2.618. And sqrt(5) = 2*phi - 1.",
            "5 is the PRIME OF THE GOLDEN RATIO. phi = (1+sqrt(5))/2.",
            "In the Cayley-Dickson chain: 5 is the FIRST number past dim(H) = 4.",
            "In tournaments: the cyclic tournament on 5 has H = 15 = 3*5.",
            "H/|Aut| for Paley T_5... well, T_5 doesn't exist (5 = 1 mod 4, not 3 mod 4).",
            "5 exceptional Lie groups. 5+1 structure. The Kervaire pattern.",
        ],
    },
    6: {
        "address": "5/7",
        "rapidity": f"{log(6)/2:.6f}",
        "pair": "(3, 3) — duplication",
        "character": [
            "The first PERFECT NUMBER. 6 = 1 + 2 + 3 = 1*2*3.",
            "The duplication of 3. Even. 6 = 2*3.",
            "Address 5/7. Note: 5 and 7 are twin primes!",
            "The Cayley address of the first perfect number has TWIN PRIME components.",
            "6 = the number of edges of the tetrahedron K_4.",
            "6 = the number of faces of the cube.",
            "6 = FLAT. In the 5-6-7 bridge: 6 triangles per vertex = the flat plane.",
            "The transition between spherical (5) and hyperbolic (7).",
            "In our theory: n=6 is where the H-spectrum first has many achievable values.",
            "6 is the THRESHOLD of complexity in tournaments.",
        ],
    },
    7: {
        "address": "3/4",
        "rapidity": f"{log(7)/2:.6f}",
        "pair": "(3, 4)",
        "character": [
            "FORBIDDEN. H = 7 is impossible for any tournament.",
            "7 = 3 + 4. Geometry + symmetry. Curvature + period.",
            "Address 3/4 = Kleiber's allometric exponent.",
            "Q(3/4) = 7. The biological scaling law IS the forbidden threshold.",
            "7 = 2^3 - 1 = Mersenne prime M_3. The third Mersenne prime.",
            "7 = imaginary dimension of the octonions.",
            "7 = where COMMUTATIVITY BREAKS in the Cayley-Dickson chain.",
            "7 triangles per vertex = the hyperbolic plane. Negative curvature.",
            "The 7th harmonic is the first 'out of tune' in Western music.",
            "In the cascade: 3 pairwise-intersecting cycles OVERSHOOT to 4+.",
            "7 is the WALL. Below 7: the commutative world holds. At 7: it breaks.",
            "The pair (3, 4) encodes the crisis: 3 is the smallest geometry",
            "and 4 is the smallest symmetry. Their sum is the first impossibility.",
        ],
    },
    8: {
        "address": "7/9",
        "rapidity": f"{log(8)/2:.6f} = 3*ln(2)/2 = THREE octaves",
        "pair": "(4, 4) — duplication",
        "character": [
            "2^3. The cube of the doubler. The octonion dimension.",
            "8 = 2*4 = 2*2*2. Pure powers of 2. The binary cube.",
            "3 octaves of rapidity. The Bott periodicity period.",
            "Cayley-Dickson: dim(O) = 8. The last normed division algebra.",
            "In our theory: n=8 is where beta_4 first appears (path homology deepens).",
            "The period-8 of Tr(M^n) mod 8 in the transfer matrix.",
            "8 = number of vertices of the cube, edges of the octahedron.",
        ],
    },
    9: {
        "address": "4/5",
        "rapidity": f"{log(9)/2:.6f} = ln(3)",
        "pair": "(4, 5)",
        "character": [
            "3^2. The square of the curvature quantum.",
            "9 = 4 + 5. The comparison between symmetry (4) and pentagon (5).",
            "The FIRST achievable H-value past the wall (H=9 exists, H=7 doesn't).",
            "9/8 = the Pythagorean major second (the interval between 4th and 5th harmonics).",
            "In the 3-3-4 pattern: 9 ~ 10 (one decade). 3^2 ~ 10.",
            "Rapidity = ln(3) = the full twelfth as a single rapidity step.",
        ],
    },
    10: {
        "address": "9/11",
        "rapidity": f"{log(10)/2:.6f} = ln(10)/2 = one DECADE",
        "pair": "(5, 5) — duplication",
        "character": [
            "2*5. The product of the doubler and the pentagon.",
            "10 = number of fingers. The decimal base. The human base.",
            "Address 9/11: both components are significant numbers in the theory.",
            "Rapidity = ln(10)/2. One DECADE. The unit of precision in logarithmic scales.",
            "pH = -rapidity/1.151 uses ln(10)/2 as its conversion factor.",
            "10 = triangular number T_4. 10 = C(5,2). 10 = dim of the Paley QR set at p=23.",
            "10 is the HUMAN INTERFACE to the number system.",
            "It is not mathematically special. It is ANTHROPOLOGICALLY special.",
            "The 3-3-4 beat pattern exists because log10(2) ~ 3/10.",
        ],
    },
    11: {
        "address": "5/6",
        "rapidity": f"{log(11)/2:.6f}",
        "pair": "(5, 6)",
        "character": [
            "The PALEY PRIME. T_11 has H = 95095 = 5*7*11*13*19.",
            "11 = 5 + 6. Pentagon + flat. The comparison between golden and flat geometries.",
            "Address 5/6. The minor third in music. Q(1/11) = 6/5.",
            "H(T_11)/|Aut(T_11)| = 95095/55 = 1729 = the TAXICAB NUMBER.",
            "11 is the first number where the tournament H-value is truly large.",
            "The 11th harmonic in music is the 'undecimal tritone' — ambiguous, between.",
            "In the Cayley algebra: 1/11 (+) 1/19 = 1/7. Two Paley primes compose to forbidden.",
        ],
    },
    12: {
        "address": "11/13",
        "rapidity": f"{log(12)/2:.6f}",
        "pair": "(6, 6) — duplication",
        "character": [
            "2^2 * 3. The smallest number divisible by both 2 and 3.",
            "12 = the number of semitones in equal temperament.",
            "12 = number of faces of the dodecahedron, vertices of the icosahedron.",
            "12 fifths ~ 7 octaves. The Pythagorean comma lives at 12.",
            "12 is the MUSICAL number. The number of the chromatic scale.",
        ],
    },
    15: {
        "address": "7/8",
        "rapidity": f"{log(15)/2:.6f}",
        "pair": "(7, 8)",
        "character": [
            "3*5. Curvature times pentagon. The first product of two odd primes.",
            "15 = 7 + 8. Forbidden number + octonion dimension.",
            "Address 7/8 = 1 - 2^{-3}. On the Mersenne ladder. Q(7/8) = 15.",
            "H(cyclic tournament on 5) = 15. The maximum at n=5.",
            "15 = 2^4 - 1 = Mersenne M_4. NOT prime (3*5).",
            "The FIRST composite Mersenne number.",
            "The pair (7, 8): forbidden + octonion. Their sum is the first non-prime Mersenne.",
        ],
    },
    21: {
        "address": "10/11",
        "rapidity": f"{log(21)/2:.6f}",
        "pair": "(10, 11)",
        "character": [
            "FORBIDDEN. The second impossible H-value. 3*7 = curvature * threshold.",
            "21 = 10 + 11. Decade + Paley prime.",
            "21 = 10101 in binary. Alternating bits. The skip-depth pattern.",
            "Address 10/11. The product of Q(1/p) for primes p <= 19.",
            "21 = T_6 = the 6th triangular number = C(7,2) = edges of the complete graph K_7.",
            "Rapidity(21) - rapidity(7) = ln(3)/2 = ONE OCTAVE.",
            "The two forbidden numbers are separated by exactly one octave.",
            "21 = number of spots on a standard die (1+2+3+4+5+6).",
            "21 = |Aut(T_7)| = the automorphism group of the Paley tournament on 7 vertices.",
        ],
    },
    42: {
        "address": "41/43",
        "rapidity": f"{log(42)/2:.6f}",
        "pair": "(21, 21) — duplication",
        "character": [
            "2*3*7. The product of the three fundamental rapidity primes.",
            "42 = the DUPLICATION of the forbidden number 21.",
            "42 = the answer to the ultimate question (Adams, Hitchhiker's Guide).",
            "42 = C_5 (the 5th Catalan number). Counts binary trees on 5 nodes.",
            "42 = 2 * 21 = the doubling of the forbidden.",
            "Is 42 achievable as H? Yes — it would need to be odd, and 42 is even. Moot.",
            "But 42's rapidity = rapidity(2) + rapidity(3) + rapidity(7).",
            "= octave + twelfth + forbidden-rapidity.",
        ],
    },
}

for n in sorted(numbers.keys()):
    info = numbers[n]
    print(f"  {n}")
    print(f"  " + "-"*30)
    if "pair" in info:
        print(f"  Cayley pair: {info['pair']}")
    print(f"  Address: {info['address']}")
    print(f"  Rapidity: {info['rapidity']}")
    print()
    for line in info["character"]:
        print(f"  {line}")
    print()
    print()
