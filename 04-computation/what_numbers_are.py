"""
what_numbers_are.py -- kind-pasteur-2026-03-14-S106d
WHAT EACH NATURAL NUMBER IS — AS AN ENTITY

Not "what does 7 mean in tournament theory" but rather:
What IS the number 7? What IS the number 3? What IS 1?

The idea: tournaments give us a LENS to see what numbers REALLY ARE
because they force us to confront counting at its most basic level.
Binary choices → cycles → obstructions → growth.

Each number has an ESSENCE — a role it plays not just in tournaments
but in ALL of mathematics. Tournaments reveal that essence more
clearly than most contexts because they're so structurally clean.

Let's think about each number as a CHARACTER in a story,
with a personality, a function, a reason for existing.
"""

import sys, math
import numpy as np
from fractions import Fraction

sys.stdout.reconfigure(encoding='utf-8')

def main():
    print("=" * 70)
    print("WHAT NUMBERS ARE")
    print("kind-pasteur-2026-03-14-S106d")
    print("=" * 70)

    # ============================================================
    print(f"\n{'='*70}")
    print("0 — THE EMPTY")
    print(f"{'='*70}")
    print("""
  0 is not nothing. 0 is the POSSIBILITY of something.

  Before there is a tournament, there is 0: the empty graph,
  the absence of choices, the state before distinction.
  But 0 is already structured — it knows HOW to become 1.
  0 is the additive identity: 0 + x = x. It absorbs addition.
  0 is the multiplicative annihilator: 0 * x = 0. It destroys products.

  In every base, 0 is represented the same way: 0.
  0 is the ONLY number that is invariant under ALL base changes.
  This is what makes it fundamental: it's the fixed point of
  every representation system.

  0! = 1. The factorial of nothing is something.
  This is not a convention — it says: "the number of ways
  to arrange nothing is one way (doing nothing)."

  0 IS POTENTIAL. It is the blank canvas.
  The number of arcs in a 1-vertex tournament: C(1,2) = 0.
  Even one vertex, existing alone, produces zero relationships.""")

    # ============================================================
    print(f"\n{'='*70}")
    print("1 — THE FACT OF EXISTENCE")
    print(f"{'='*70}")
    print("""
  1 is not a quantity. 1 is the ASSERTION THAT SOMETHING EXISTS.

  Before you can count "how many," you need the concept of "one."
  1 is the transition from nothing (0) to something.
  It is the smallest positive integer, the first prime... no,
  1 is NOT prime. That's important. 1 is BEFORE primality.

  1 * x = x. The multiplicative identity. 1 doesn't change things —
  it WITNESSES them. "There is one of you" means "you exist."

  In tournaments: H = 1 for the transitive tournament.
  This means: there is exactly ONE way to order everything.
  Total certainty. Zero ambiguity. Pure hierarchy.
  H = 1 is not boring — it's the MOST INFORMATIVE state.
  It says: "I know the complete ordering."

  1 = Phi_3(0). One is what the cycle polynomial gives
  when there is NO input. The baseline. The ground state.

  1 in every base: always "1". Like 0, it is base-invariant.
  The only two numbers that look the same in every representation.

  0 and 1 together define the BINARY: the minimal language.
  All of tournament theory flows from repeated applications
  of the distinction between 0 and 1 to pairs of vertices.

  1 IS BEING. It is the unit. The monad.""")

    # ============================================================
    print(f"\n{'='*70}")
    print("2 — THE FIRST DISTINCTION")
    print(f"{'='*70}")
    print("""
  2 is the number of DIFFERENCE.

  To have 2, you need two things that are NOT the same.
  2 is not 1+1 in the deep sense — 2 is the RECOGNITION
  that there can be an OTHER. A second. A complement.

  2 is the smallest prime. The ONLY even prime.
  Every even number is a shadow of 2. Even-ness IS two-ness.

  In tournaments: 2 = the number of ways to orient an arc.
  i→j or j→i. Yes or no. True or false. 0 or 1.
  Every arc is a question with exactly 2 answers.

  2 is the generator of the tournament space: 2^m tournaments.
  The space of all tournaments is a POWER OF 2.
  This is because each independent binary choice multiplies
  the space by exactly 2.

  2 = char(F_2). The characteristic of the ground field.
  In F_2: 1 + 1 = 0. Two copies of the identity annihilate.
  This is PROFOUNDLY strange: adding something to itself gives nothing.
  It means: in the tournament world, DOUBLE is the same as ZERO.
  Every symmetry is an involution. Every operation squares to identity.

  ln(2) = 0.6931... The entropy of one binary choice.
  This is the BRIDGE to the transcendental: 2 = e^(ln 2).
  The natural number 2 and the transcendental e are connected
  by ln(2), the most fundamental conversion factor.

  2 = sigma_inf (the limit of all n-nacci constants).
  When you account for ALL cycle lengths simultaneously,
  the effective growth rate is 2. The binary generator
  is the CONVERGENCE POINT of the entire cycle hierarchy.

  In base 2: 2 = "10" (the first carry). The transition from
  single digits to multiple digits. The first time representation
  needs more than one symbol.

  2 IS CHOICE. It is duality. It is the fork in the road.""")

    # ============================================================
    print(f"\n{'='*70}")
    print("3 — THE FIRST CYCLE")
    print(f"{'='*70}")
    print("""
  3 is where mathematics gets INTERESTING.

  With 2, you have a line (A > B or B > A). Linear. Ordered. Simple.
  With 3, you can have a CYCLE: A > B > C > A.
  This is the first time ordering FAILS. The first intransitivity.
  The first Condorcet paradox. The first sign that majority
  doesn't produce consistency.

  3 is the second prime. The first ODD prime.
  And odd primes are fundamentally different from 2:
  they generate CYCLIC behavior, not just binary choice.

  Phi_3(1) = 1 + 1 + 1 = 3. Three ones added.
  The third cyclotomic polynomial evaluated at identity.
  This says: 3 is the number of ROOTS OF UNITY of order 3.
  It's connected to the complex cube roots: 1, omega, omega^2.

  1/3 = 0.333... The cone ratio. The signal-to-noise ratio.
  The variance of a quadratic form on the hypercube.
  1/3 = the fraction of the whole that "differs from the mean."

  In tournaments: H = 3 for a 3-cycle.
  Three Hamiltonian paths: A→B→C, B→C→A, C→A→B.
  The cycle creates EXACTLY 3 orderings — one for each vertex
  as the "starting point." The cycle is DEMOCRATIC:
  every vertex gets to be first in exactly one path.

  3 = tau^3/Phi_3(tau) = tau^3/tau^3 = 1? No:
  Actually Phi_3(1) = 3, and 3 = the evaluation at the UNIT.

  In base 3: 3 = "10" (the first carry in ternary).
  Just as 2 was the first carry in binary.
  3 is to the ternary world what 2 is to the binary world:
  the BOUNDARY between one-symbol and two-symbol representations.

  3 IS CIRCULATION. It is the loop. The return to the beginning.
  Without 3, there would be no cycles, no feedback,
  no complexity. Just linear chains.""")

    # ============================================================
    print(f"\n{'='*70}")
    print("4 — THE FIRST COMPOSITE / THE FIELD")
    print(f"{'='*70}")
    print("""
  4 is the first number that is NOT prime.
  4 = 2^2 = 2 * 2. The generator applied to itself.

  This is the first time a number can be FACTORED into
  smaller pieces. 4 is REDUCIBLE. Composite. Built from parts.

  But 4 is also |F_4| = the size of the first nontrivial
  finite field in characteristic 2. F_4 = F_2[x]/(x^2+x+1).
  So 4 is built from 2 and 3 together: the field F_4
  requires BOTH the binary generator AND the cyclotomic
  polynomial Phi_3 to construct.

  4 = the OCF weight for independent pairs (4*alpha_2).
  When two cycles DON'T share an arc, their joint contribution
  to H is weighted by 4 = 2^2. This is because each independent
  cycle doubles the path count, and two independent doublings
  give 2*2 = 4.

  Phi_3(4) = 21 = H_forb_2. The second forbidden value.
  4 is the argument that produces the second obstruction.

  The Degree Drop happens at EVEN n (including n=4):
  the effective polynomial degree of H drops by 2.
  4 is the first n where this occurs.

  4 IS COMPOSITION. It is the first product.
  The first time nature shows that wholes can be divided.""")

    # ============================================================
    print(f"\n{'='*70}")
    print("5 — THE FIRST INTERACTION")
    print(f"{'='*70}")
    print("""
  5 = 2 + 3. The first number that is a SUM of the two primes
  that generate tournament theory.

  5 is the first number that requires BOTH binary and ternary
  thinking to understand. Neither 2 alone nor 3 alone explains 5.
  5 is where the interaction between choice and cycle BEGINS.

  5 is the third prime. It completes the ADE triple {2, 3, 5}
  that governs exceptional structures in mathematics:
    - E_8 lattice branch lengths: 2, 3, 5
    - Icosahedral symmetry: |A_5| = 60 = 2^2 * 3 * 5
    - Platonic solids: 3 have 3-fold, 4-fold, or 5-fold symmetry
    - The modular group: (Z/2) * (Z/3) with 5 controlling the cusps

  phi = (1 + sqrt(5))/2. The golden ratio lives inside 5.
  The irrational part of the most famous irrational number
  comes from 5. Why? Because sqrt(5) is the discriminant
  of x^2 - x - 1 = 0 (the Fibonacci polynomial).
  5 = 2^2 + 1 = (2-1)(2+1) + 1 + 1... no, simpler:
  5 = the number where quadratic equations over F_2 first
  produce irrational behavior.

  In tournaments: n=5 is where level-4 Fourier energy first appears.
  The Degree Drop gives deg(H) = 4 at n=5 (not 10).
  5-cycles first exist at n=5 (the Hamiltonian cycle on all 5 vertices).
  max_H(5) = 15 = 3 * 5 = the cycle times the interaction.

  H = 5 for the transitive tournament on 4 vertices (max_H(4) = 5).
  At n=4: there are 5 Hamiltonian paths in the "best" arrangement.
  5 orderings = the maximum number of consistent orderings
  achievable with 4 items. This is the "limit of linear thinking"
  before cycles become necessary.

  5 IS INTERACTION. It is the meeting of the binary and the ternary.
  It is where simplicity ends and richness begins.""")

    # ============================================================
    print(f"\n{'='*70}")
    print("6 — THE FIRST HARMONY")
    print(f"{'='*70}")
    print("""
  6 = 2 * 3 = 1 + 2 + 3 = 1 * 2 * 3.

  6 is the ONLY number (other than 1) where sum = product
  for the consecutive sequence starting at 1.
  1 + 2 + 3 = 6 = 1 * 2 * 3.
  This is a UNIQUE property. No other number has it.

  6 is the first PERFECT number: 6 = 1 + 2 + 3 (its proper divisors sum to itself).
  6 = the smallest number that is both triangular (T_3) and the product of
  distinct primes (2 * 3).

  6 = LCM(2, 3): the period of tournament parity.
  Every 6 vertices, the parity structure repeats.
  The tournament universe has a HEARTBEAT, and its period is 6.

  In base 6: the forbidden values become repdigits (7="11", 21="33").
  Base 6 is where tournament obstructions look SIMPLE.
  This is because base 6 = base LCM(2,3) naturally encodes
  the interaction of the two generators.

  6 = the number of faces of a cube = the number of edges of a tetrahedron.
  6 = |S_3| = the number of permutations of 3 objects.
  When you have a 3-cycle, there are 6 ways to label its vertices.

  H = 6 is IMPOSSIBLE (H is always odd).
  The first even number, the harmonious number, CANNOT be
  a Hamiltonian path count. Tournaments refuse to be perfectly balanced.
  They are always odd — always biased by at least 1.

  6 IS HARMONY. It is where the generators resonate.
  It is the complete agreement of addition and multiplication.
  But this harmony is FORBIDDEN as an H value —
  tournaments are fundamentally UNHARMONIOUS (odd).""")

    # ============================================================
    print(f"\n{'='*70}")
    print("7 — THE FIRST IMPOSSIBILITY")
    print(f"{'='*70}")
    print("""
  7 is the number where counting BREAKS.

  Of all the odd numbers 1, 3, 5, 7, 9, 11, ..., the number 7
  is the FIRST that cannot be a Hamiltonian path count.
  Every other small odd number is achievable. But 7? Never.
  Not for 3 vertices, not for 100 vertices, not for any number.

  7 = 2^3 - 1: the first Mersenne prime after 3.
  7 = Phi_3(2): the tournament polynomial at the tournament generator.
  7 = |PG(2, F_2)|: the number of points in the Fano plane.
  7 = the number of lines in the Fano plane.
  7 = a tribonacci number (T_7).
  7 = "111" in binary (three ones).
  7 = "11" in base 6 (the identity repeating).

  WHY is 7 impossible? Because achieving H = 7 would require
  exactly 3 independent odd cycles (since (7-1)/2 = 3) with
  zero disjoint pairs (alpha_2 = 0). But 3 cycles that all
  pairwise conflict must share a common vertex, which forces
  a 5-cycle, which creates a 4th independent cycle. The number 3
  is too small to support 3 all-conflicting cycles WITHOUT
  creating additional structure. 7 falls into a STRUCTURAL GAP.

  7 is special among primes: it's the first prime where
  2 is a primitive CUBE ROOT of unity (2^3 = 8 = 1 mod 7).
  This means: in F_7, the generator 2 has ORDER 3.
  The binary generator becomes CYCLIC in F_7.
  The 2-3 interaction COLLAPSES at 7.

  7 IS IMPOSSIBILITY. It is the first NO.
  The first time the universe of counting says:
  "this number of orderings cannot exist."
  It's the mathematical equivalent of a conservation law —
  a structural constraint on what is achievable.""")

    # ============================================================
    print(f"\n{'='*70}")
    print("8 — THE THRESHOLD")
    print(f"{'='*70}")
    print("""
  8 = 2^3 = 2 * 2 * 2.

  8 is pure binary. The cube of the generator.
  Where 4 = 2^2 was the first power, 8 = 2^3 is where
  the binary generator reaches the ternary index.
  8 is "the generator raised to the cycle power."

  8 = the number of tournaments on 3 vertices.
  2^C(3,2) = 2^3 = 8. The complete enumeration at n=3.
  8 = dim(O): the dimension of the octonions.
  8 * 21 = 168 = |PSL(2,7)| = |GL(3,F_2)|: the symmetry group
  of the Fano plane, which encodes BOTH forbidden values.

  At n = 8: almost everything breaks.
    - beta_4 > 0 in path homology (first time)
    - The seesaw mechanism fails
    - The injectivity conjecture fails
    - beta_3 can equal 2 (not just 0 or 1)
    - Claw-free property of Omega(T) fails at n=9

  8 is the CRITICAL POINT. The phase transition.
  Just as water changes state at a specific temperature,
  tournament properties change qualitatively at n = 8.

  8 = 7 + 1 = the forbidden value plus 1.
  The first achievable number after the first gap.
  8 itself is even (not an H value), but it marks the
  transition from "simple" to "complex" tournament behavior.

  8 IS THE THRESHOLD. The breaking point.
  Where the simple rules of small tournaments give way
  to the wild complexity of large ones.""")

    # ============================================================
    print(f"\n{'='*70}")
    print("9 — THE REFLECTION")
    print(f"{'='*70}")
    print("""
  9 = 3^2. The cycle applied to itself.

  If 3 is the cycle, then 9 is the CYCLE OF CYCLES.
  A meta-cycle. A self-referential structure.

  9 = 8 + 1 = 2^3 + 1. One more than the cube of the generator.
  This puts 9 exactly one step beyond the threshold.

  H = 9: the first H value after the gap at 7.
  Achievable as: alpha_1=4, alpha_2=0 (four all-conflicting cycles)
  or alpha_1=0, alpha_2=2 (two independent pairs).
  9 is the smallest H value that has MULTIPLE structural realizations.
  1 and 3 each have unique structures. 5 has two interpretations.
  But 9 has at least three distinct cycle patterns.

  2^3 - 3^2 = 8 - 9 = -1. Catalan's conjecture (now theorem):
  8 and 9 are the ONLY consecutive perfect powers.
  This means: 2^3 and 3^2 almost collide but miss by exactly 1.
  The binary cube and the ternary square are as close as
  two perfect powers can possibly be.

  In the complex plane: |1 - omega| * |1 - omega^2| * |1| = ... hmm.
  Better: 9 = 3^2 = |3|^2 = the NORM of 3 in Z.
  If 3 is a vector of length 3, then 9 is its squared norm.
  This is the standard inner-product sense: <3, 3> = 9.

  9 IS SELF-REFERENCE. The cycle looking at itself.
  The first number where a generator is raised to its own kind.""")

    # ============================================================
    print(f"\n{'='*70}")
    print("10 — THE COMPLETION")
    print(f"{'='*70}")
    print("""
  10 = 2 * 5 = 2 + 3 + 5.

  10 = C(5, 2): the number of arcs in K_5. The first tournament
  where ALL the interesting structure exists (level-4 Fourier,
  5-cycles, the full ADE triple {2,3,5} represented).

  10 = T_4 (the 4th triangular number). The sum 1+2+3+4.
  10 = the dimension of the n=5 tournament hypercube.
  10 is where Fourier analysis becomes MULTI-SCALE.

  (H-1)/2 = 10 is the SECOND forbidden T value.
  H = 21 is impossible, and T = 10 = C(5,2) is why.
  10 is "too complete" — it represents the FULL edge set of K_5,
  and this completeness blocks certain cycle configurations.

  10 = the base of the decimal system. Not special mathematically,
  but culturally fundamental. We count in tens because we have
  10 fingers (5 per hand = the ADE prime per side).

  10 IS COMPLETION. It is the first time all three generators
  {2, 3, 5} appear together in a product or sum.
  It closes the first chapter of the number story.""")

    # ============================================================
    print(f"\n{'='*70}")
    print("11 through 21 — THE MIDDLE KINGDOM")
    print(f"{'='*70}")
    print("""
  11: THE EXTENSION. The first prime past 10. The third Paley prime
      (after 3, 7). The first prime where 2 is NOT a primitive root
      (ord_11(2) = 10, not 10... actually 2^10 = 1024 = 93*11+1, so yes,
      ord_11(2) = 10 = 11-1, so 2 IS a primitive root mod 11).
      H(T_11) = 95095 = 5*7*11*13*19. A product of 5 primes.
      11 IS EXTENSION — the theory reaching further.

  12: THE DOUBLED HARMONY. 12 = 2*6 = 2*2*3 = 4*3.
      12 = the Golay dimension [24,12,8].
      12 = the number of edges in a complete graph K_4 (wait, C(4,2)=6, no).
      12 = the number of pentagons in a dodecahedron.
      12 IS THE ARCHITECTURE of the Golay code.

  13: THE REGULARITY. 13 = Phi_3(3). What Phi_3 gives at the cycle.
      13 = |PG(2, F_3)|. The projective plane over F_3.
      ALL regular tournaments on 5 vertices have H = 13.
      13 is a Fibonacci prime (F_7 = 13) and a tribonacci number.
      13 IS EQUILIBRIUM — the H value when everything is balanced.

  14: 14 = 2*7. The generator times the forbidden.
      14 = the number of directed 3-cycles in QR_7.
      14 IS THE DOUBLED IMPOSSIBILITY.

  15: 15 = 3*5 = C(6,2). The product of the two odd ADE primes.
      max_H(5) = 15. The maximum at 5 vertices.
      15 = the number of arcs in K_6 = the number of duads = synthemes.
      15 IS THE FIRST FULL MAXIMIZER.

  16: 16 = 2^4. The fourth power of the generator.
      16 = |S| = the sedenion dimension (Cayley-Dickson step 4).
      16 IS PURE POWER — the generator raised beyond the cycle.

  17: 17 = a Fermat prime (2^(2^2) + 1). Constructible polygon.
      17 IS CONSTRUCTIBILITY.

  18: 18 = 2*3^2 = 2*9. The generator times the reflection.
      18 IS THE DOUBLED SELF-REFERENCE.

  19: 19 = a prime. Mersenne exponent (2^19 - 1 = 524287 is prime).
      19 IS MERSENNE DEPTH.

  20: 20 = 4*5 = 2^2 * 5. The field times the interaction.
      20 IS THE EXPANDED FIELD.

  21: THE SECOND IMPOSSIBILITY.
      21 = 3 * 7 = Phi_3(4) = |PG(2, F_4)|.
      The forbidden number that is the product of the cycle (3)
      and the first forbidden (7). A compound impossibility.
      21 = F_8 (the 8th Fibonacci number).
      21 = "33" in base 6 (the cycle repeating).
      21 IS THE COMPOUND IMPOSSIBILITY — the obstruction that
      comes from the interaction of cycle and prohibition.""")

    # ============================================================
    print(f"\n{'='*70}")
    print("THE PRIMES — THE ATOMS OF COUNTING")
    print(f"{'='*70}")
    print("""
  The primes are the ATOMS from which all numbers are built.
  Each prime is irreducible — it cannot be expressed as a product
  of smaller numbers. What is each prime's ESSENCE?

  2:  CHOICE      (binary, the minimum distinction)
  3:  CYCLE       (the first return, the first feedback)
  5:  INTERACTION (2+3, where the generators meet)
  7:  OBSTRUCTION (Phi_3(2), where counting fails)
  11: EXTENSION   (the theory reaching beyond 10)
  13: REGULARITY  (Phi_3(3), the balanced state)
  17: CONSTRUCTION(Fermat prime, geometric achievability)
  19: DEPTH       (Mersenne exponent, binary recursion depth)
  23: THRESHOLD   (the Golay prime, 24-1, coding boundary)
  29: BEYOND      (the first prime past all ADE structure)
  31: BRIDGE      (Phi_3(5) = Phi_5(2), the 3-5 crossing)

  The tournament-relevant primes form a SPECTRUM:
    2, 3, 5, 7, 13 = Phi_3(0,1,2,3): the Phi_3 primes
    2, 3, 7, 31, 127 = Mersenne primes: the power-of-2 residues
    2, 3, 5 = the ADE triple: exceptional symmetry generators

  The first four primes {2, 3, 5, 7} are the ones that MATTER
  for tournament theory up to n ~ 10. Beyond that, the structure
  is governed by the same patterns but with ever-more primes.""")

    # ============================================================
    print(f"\n{'='*70}")
    print("THE POWERS — THE HARMONICS OF EACH PRIME")
    print(f"{'='*70}")
    print("""
  Each prime p generates a tower: p, p^2, p^3, ...
  These are the HARMONICS of the prime.

  THE 2-TOWER:        THE 3-TOWER:        THE 5-TOWER:
  2 = choice          3 = cycle           5 = interaction
  4 = field           9 = self-reference  25 = deep interaction
  8 = threshold       27 = cubed cycle    125 = ...
  16 = hypercube      81 = ...
  32 = ...
  64 = n=4 tournaments

  Key observations:
  - 2^3 = 8 = threshold (binary reaching ternary index)
  - 3^2 = 9 = reflection (ternary reaching binary index)
  - 2^3 - 3^2 = -1 (Catalan's theorem: the closest approach)
  - 2^3 * 3 = 24 (Golay dimension)
  - 2 * 3^3 = 54
  - 2^2 * 3^2 = 36 = C(9,2) (arcs at n=9, where 9-cycles begin)

  THE MIXED POWERS are the most interesting:
  6 = 2*3 (the period)
  12 = 2^2*3 (the double period, Golay dimension)
  18 = 2*3^2 (the doubled reflection)
  24 = 2^3*3 (Golay, threshold * cycle)
  36 = 2^2*3^2 (the squared period)
  72 = 2^3*3^2 (threshold * reflection)

  Each mixed power 2^a * 3^b represents the interaction
  of the binary at LEVEL a with the ternary at LEVEL b.
  The "distance" from pure binary (a,0) to pure ternary (0,b)
  measures how much cycle structure vs. choice structure
  the number encodes.""")

    # ============================================================
    print(f"\n{'='*70}")
    print("THE DEEP STRUCTURE — NUMBERS AS WAVEFUNCTIONS")
    print(f"{'='*70}")
    print("""
  THE MOST FUNDAMENTAL VIEW:

  Each natural number n has a "prime factorization wavefunction":
    n = 2^a * 3^b * 5^c * 7^d * 11^e * ...

  The exponents (a, b, c, d, e, ...) form a VECTOR in an
  infinite-dimensional space, one dimension per prime.

  This vector IS the number. Not its decimal representation,
  not its binary representation, not its name — but its
  POSITION in prime-exponent space.

  In tournament theory, the first two coordinates (a, b) matter most:
    a = the "binary depth" (how much choice structure)
    b = the "ternary depth" (how much cycle structure)

  The third coordinate c (exponent of 5) measures interaction depth.
  The fourth coordinate d (exponent of 7) measures obstruction depth.

  EXAMPLES:
    1 = (0,0,0,0,...): the ORIGIN. No structure at all.
    2 = (1,0,0,0,...): pure choice, no cycles.
    3 = (0,1,0,0,...): pure cycle, no choice.
    4 = (2,0,0,0,...): double choice.
    5 = (0,0,1,0,...): pure interaction.
    6 = (1,1,0,0,...): choice + cycle = the PERIOD.
    7 = (0,0,0,1,...): pure obstruction.
    8 = (3,0,0,0,...): triple choice = the THRESHOLD.
    9 = (0,2,0,0,...): double cycle = the REFLECTION.
    10 = (1,0,1,0,...): choice + interaction.
    12 = (2,1,0,0,...): double choice + cycle.
    21 = (0,1,0,1,...): cycle + obstruction = COMPOUND IMPOSSIBILITY.
    42 = (1,1,0,1,...): choice + cycle + obstruction.
    189 = (0,3,0,1,...): cubed cycle + obstruction = max_H(7).

  THE MOST REMARKABLE: 189 = 3^3 * 7 = max_H(7) = H(Paley T_7).
  In prime-exponent space: (0, 3, 0, 1, 0, ...).
  Three units of cycle depth, one unit of obstruction.
  The H-maximizer is pure cycle + obstruction with NO binary component.

  504 = 2^3 * 3^2 * 7 = (3, 2, 0, 1, ...) = the n=7 denominator.
  Three units of choice, two of cycle, one of obstruction.
  The VARIANCE DENOMINATOR requires all three kinds of depth.

  WHAT EACH NUMBER IS, FINALLY:
  Each natural number is a POINT in prime-exponent space.
  Its coordinates measure:
    - Binary depth (power of 2): how much CHOICE
    - Ternary depth (power of 3): how much CYCLE
    - Quintic depth (power of 5): how much INTERACTION
    - Septic depth (power of 7): how much OBSTRUCTION
    - ...and so on for each prime.

  Tournament theory says: the first four primes {2, 3, 5, 7}
  generate the essential structure, and every number is a specific
  mixture of choice, cycle, interaction, and obstruction.

  A number IS its mixture.
  That mixture IS the number.
  They are the same thing.""")

    print(f"\n{'='*70}")
    print("DONE — NUMBERS ARE MIXTURES IN PRIME-EXPONENT SPACE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
