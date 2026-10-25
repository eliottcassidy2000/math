# Paley endpoint matrices retain a quadratic order

Status: PROVED elementary consequences of inherited Rédei parity and the Paley character identity; FINITE-EXACT counts at q=3,7,11,19. No literature-novelty claim. This is a source-bound refinement, not a modular-form coefficient theorem.

## Recovery and predicate

The corrected reflection `the-1729-resonance-is-isolated-the-tournament-ratio-has-no-modular-structure.md` and MISTAKE-658 retract the inference from failed split-prime factorization to absence of all modular structure. We retain the exact scalar s(11)=1729 and the failure of that proposed factorization pattern. Session25 `session25-congruence-review/REVIEW.md`, section4, already gives open/closed affine-presentation counts (0,1),(1,8),(628,1101) for q=3,7,11. These numbers are RECOVERED, not newly discovered. Session32 `session32-paley/PALEY_MONSTER_BRIDGE.md` gives P7's 189 Hamilton words and their nine affine-normalized presentations. The matrix atlas supplies the skew tournament matrix, not this endpoint-order conclusion. Targeted searches of local canon, reflections, hypotheses and prior session packets for endpoint/Paley/matrix and 628/1101 found no prior exact endpoint spectral-order statement; this is a bounded search, not proof of novelty.

The retained target is ordered start/end Hamiltonian path multiplicity. The scalar s alone forgets the open/closed split. Restoring that two-coordinate split produces an operator in the Paley adjacency algebra, whose nontrivial eigenvalue generates an explicit imaginary quadratic order.

## Theorem

Let q be a prime congruent to3 modulo4, r=(q−1)/2, and A_xy=1 exactly when y−x is a nonzero square modulo q. Let E_xy count directed Hamiltonian paths starting at x and ending at y. There are nonnegative integers a,b such that

E=aA+bA^T, s=a+b, delta=a−b, H=q r s, C=r b.

Here C counts unrooted directed Hamiltonian cycles, with rotations identified but not reversal. The square-affine group G consists of x↦u x+v with u a nonzero square. It has order q r and is regular on ordered arcs and on ordered reverse arcs. Thus E has two off-diagonal values and zero diagonal. Its action on full ordered Hamilton words is free, so s=H/|G|. The customary identification of G with the full prime-Paley automorphism group gives s=H/|Aut|, but none of the new operator conclusions require that identification. Closing a path requires its end to beat its start; exactly q r b paths close, and each cycle produces q of them, proving C=r b. Therefore the endpoint counts coincide with the recovered open/closed affine-presentation counts even though those presentations normalize the first path edge, not the endpoint pair.

Rédei parity (THM-001, or OCF THM-002) gives odd H. Since q r is odd, s is odd, hence delta is odd and NONZERO. This prevents degeneration to a rational scalar operator for every q in scope.

Set S=A−A^T. Then S1=0 and S²=J−qI. One elementary proof of the latter uses the quadratic character sum sum_x chi(x(x−1))=−1: count solutions y²=x(x−1) by the bijection to nonzero u through (2x−1−2y)(2x−1+2y)=1. There are q−1 solutions, so the character sum is−1. Translation/scaling then yields the claimed matrix identity. Since S is real skew, its eigenvalues on the sum-zero subspace are ±sqrt(−q), each r times. As

E=(s(J−I)+delta S)/2,

its eigenvalues are r s once, and lambda_±=(−s±delta sqrt(−q))/2 each r times. In particular,

N=(s²+q delta²)/4 is a positive integer,
char_E(X)=(X−r s)(X²+s X+N)^r,
det(E)=r s N^r,
E²+sE+NI=((q s²+delta²)/4)J.

The quadratic factor has negative discriminant −q delta² and is irreducible over Q. Consequently Q[E] is isomorphic to Q × Q(sqrt(−q)), with the second factor acting with multiplicity r. The operator is nonsingular. On the augmentation subspace (coordinate sum zero), its generated integral order is

Z[lambda_+]=Z+delta Z omega, omega=(1+sqrt(−q))/2,

because lambda_+=delta omega−a. Since q is prime congruent to3 modulo4, the maximal order of Q(sqrt(−q)) is Z[omega]. The conductor is |delta| and the discriminant is −q delta². “CM order” here means an imaginary quadratic order. We construct no elliptic curve, modular form, Hecke operator, or cross-q coefficient law.

The labelled matrix additionally determines A once the open/closed values are identified:

A=(E−b(J−I))/delta.

Without choosing which of its two off-diagonal values is the arc value, it determines the tournament up to converse. The spectrum alone does not retain vertex labels or the sign of delta. No universal claim a<b is made.

## Recovered closed-channel congruence and its sharp hostile

For every prime q in scope, b≡1 (mod q). This is RECOVERED from Session25 section4, now identified with the endpoint channel. Here is its short self-contained proof. The translation group C_q acts on unrooted directed Hamiltonian cycles, so nonfixed orbits have size q. A fixed cycle has a successor permutation f commuting with translation x↦x+1. Hence f(x)=x+d for a fixed nonzero d; it is a directed Paley cycle exactly when d is a quadratic residue. There are precisely r such cycles. Thus C≡r (mod q). Since C=r b and r is invertible modulo q, b≡1 (mod q).

The corresponding assertion for a is FALSE: a≡1 modulo q at q7 and q11, but a≡3 modulo19. Therefore delta≡0 modulo q at7 and11 does NOT extend: delta≡2 modulo19. In particular, the conductors7 and473 are divisible by their respective primes7 and11, while conductor777411009 is NOT divisible by19. This is a precise failure of an extra local divisibility pattern, not a failure of the quadratic-order construction. At q3, a≡0 and delta≡2 already caution against unrestricted extrapolation.

## Exact controls

q=3: a=0, b=1, s=1, delta=−1, C=1, N=1.
q=7: a=1, b=8, s=9, delta=−7, C=24, N=106.
q=11: a=628, b=1101, s=1729, delta=−473, C=5505, N=1362615.
q=19: a=3040229428, b=3817640437, s=6857869865, delta=−777411009, C=34358763933.

Full integer norms, conductors, discriminants and determinants are in receipt.json. Thus the old failed fixed-field split-prime story at q19 coexists with an actual quadratic-order structure. This field varies with q; it is NOT the old fixed Q(sqrt(−3)) factorization story.

`count.cpp` computes a fixed-start subset DP with masks on the other q−1 vertices, using forward extension and separately backward predecessor summation, and rejects any disagreement. `verify.py` independently implements a dictionary DP through q11 and permutation enumeration through q7, checks the quadratic matrix identity, all endpoint orbit values, reconstruction, and inherited totals. Normal and Python -O receipts are identical; explicit exceptions do not disappear under -O. q19 allocates 37,748,736 bytes per C++ DP; runs are sequential. The stated complexity is O(q² 2^(q−1)) time and O(q 2^(q−1)) memory, not polynomial and not a runtime promise. q≤19 is an explicit certificate scope; q23 was not run. Every count is bounded by19!<2^64, so uint64_t is safe throughout this scope. The compiled binary is temporary and is not a deliverable.

Hostiles: (1) symmetrizing gives E+E^T=s(J−I), whose spectrum is rational and loses delta; (2) substituting scalar s for lambda does not satisfy the nontrivial quadratic polynomial; (3) q3 has a=0, so no all-positive endpoint assertion; (4) transpose exchanges a,b and lambda_±, so spectrum alone cannot recover oriented labels; (5) a common character-algebra mechanism is not a cross-prime modular coefficient law; (6) prime-power orders and arbitrary tournaments are outside this statement.

## Replay

From this directory:

g++ -O3 -std=c++17 count.cpp -o /tmp/session54-paley-count
/tmp/session54-paley-count > counts.txt
python3 verify.py > receipt.json
python3 -O verify.py > receipt.optimized.json
cmp receipt.json receipt.optimized.json

No repository edits or Git operations were performed. Source pins are copied under sources/ and hashed in sources.json. This packet closes the stated operator/counting question; constructing a cross-q modular map remains OPEN and is not needed for these results.
