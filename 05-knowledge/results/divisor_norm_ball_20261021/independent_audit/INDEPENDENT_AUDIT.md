# Independent audit: exact divisor norm ball

Date: 2026-10-21. Verdict: **PASS, within the stated scope. No mathematical correction is required.**

The reviewed producer manifest is SHA-256
`7f00140b9b85d0b670231f8481f9cc78ec5f615975a60dde36a7e423f9072265`.
Its payload consists of `README.md`, `DIVISOR_NORM_BALL.md`, `verify_norm_ball.py`,
`verification.json`, `verification.optimized.json`, and `SOURCE_PROVENANCE.json`.
All payload hashes and the producer's `SHA256SUMS` were verified before and after
independent computation. This verdict applies to these bytes, not later edits.

## 1. Proof audit

### Exact carrier and labels

For distinct primes and positive integer k, the excluded divisors are explicitly
1 and p^kqr. This avoids any ambiguity in the phrase “proper divisors.” Unique
factorization makes the exponent triples faithful labels.

The completed-square inequality forces the coefficient of kω to be −1, 0, or 1.
On these three rows, the norm inequality gives respectively the integer intervals
[−k,0], [−k,k] without zero, and [0,k]. The four inverse cases in Theorem 1
recover exactly one allowed triple from every point. The two deleted corners are
the only collision in the full exponent box. Thus this is genuinely the entire
nonzero norm ball, not an injection supported only by matching cardinalities.

The substitutions for C and Q are correct. Complement becomes negation; exchanging
the two singly occurring prime labels becomes complex conjugation. The inherited
chart is exactly `(2 Re(z), 2 Im(z)/(k sqrt(3)))`. The fixed counts are 0, 2k, and
respectively 2 or 0 for C, Q, and CQ, with the CQ count depending on evenness.
The k=2 literal table, the two perpendicular-axis partitions, and all three norm
levels are correct.

### Full isometry group, including the small exception

Negation gives centroid zero, so a Euclidean isometry of the finite ball fixes
zero and is orthogonal. When k≥2, every nonreal order element has norm at least
3k²/4≥3. The unique shortest pair is therefore ±1. Preserving that pair forces
the real axis, leaving precisely identity, negation, conjugation, and negative
conjugation. All four preserve both order and ball. The lattice statement follows
from the same shortest-vector argument. At k=1 the unit hexagon has twelve
Euclidean symmetries, so the stated exception is necessary and is retained.

### Conductor and rotation gate

Multiplication by ω takes x+kyω to −ky+(x−ky)ω. It stays in O_k exactly when
k divides x. The resulting set is kO_1. Since O_1 has basis 1,ω, this is also
exactly the condition zO_1⊆O_k, and therefore the conductor ideal. Every
ω-stable additive sublattice of O_k lies in it.

Intersecting with the ball and dividing by k leaves positive integral norm at
most one, hence precisely the six Eisenstein units. The six decoded divisors
and both displayed three-cycles are correct. At k=1 this core is the whole ball;
for every k≥2 all other points leave the order after one rotation. The statement
about losing order preservation is read in this explicitly established k≥2
context. Localization removes the order constraint, not the original finite
carrier constraint.

### Even-conductor idoneality

The class group here is that of proper invertible ideal classes. The conductor
formula must divide by the unit index three; the proof does. The maximal order
has class number one, and each nonmaximal order has only the units ±1.

For f=2^a 3^b times factors at primes ℓ≥5, with a≥1, the inert-prime-2 factor
cancels the denominator three. The remaining class number is

    2^(a−1) 3^b ∏ ℓ^(e_ℓ−1)(ℓ−χ_(−3)(ℓ)).

Each prime ℓ≥5 supplies a factor divisible by three, in both residue classes
modulo three. Consequently every even conductor with an odd prime factor has
class number divisible by three and cannot have class-group exponent dividing
two. There is no overlooked cancellation of this obstruction.

At f=2,4,8 the listed primitive reduced forms are complete and all ambiguous.
Their counts equal the independently derived class numbers 1,2,4. At f=16,
(13,8,16) is primitive and strictly reduced, with discriminant −768. For
f=2^a, a≥5, the displayed form (16,8,1+3·4^(a−3)) has odd final coefficient,
discriminant −3·4^a, and strict inequalities 0<8<16<C. Its inverse is a distinct
reduced representative. This proves the entire higher dyadic tower, not merely
a checked initial range. The claimed classification 2,4,8 is unconditional
within the cited classical theory and does not require global 65-number
completeness.

## 2. Connection and non-identification audit

The inherited labelled reflection obstruction remains intact. Every involution
of at least two vertices fixes an unordered pair: a transposition itself gives
one, while the identity fixes every pair. Thus graph complement followed by
that vertex involution cannot fix a graph. An even-k divisor carrier has
CQ-fixed points and cannot map equivariantly to that target with the specified
operation labels, even by a many-to-one map. Swapping the two source axis
reflections does not remove the obstruction, because both then have fixed
points. The new norm-ball chart neither weakens this theorem nor replaces its
source-preserving resolution.

The native twelve-vector alphabet really has norms 1 and 3 only; the k=2 ball
has norms 1,3,4. Therefore neither norm preservation nor any nonzero common
scalar rescaling identifies the full carriers. The three norm levels cannot
become two under a fixed nonzero scaling.

The map is not multiplicative: p maps to 1 while p² maps to 2. A ball point is
not an ideal class. Prime labels, invertibility conditions, and the difference
between point cardinality and class number must remain explicit. No source ID,
terminal, braid clock, or cancellation certificate is manufactured by this
static chart. No dynamical, tournament-class, or Collatz conclusion is audited
or certified here. Historical novelty is not claimed.

## 3. Independent executable controls

`independent_audit.py` imports no producer functions. It uses:

- doubled real coordinates and an enlarged exact lattice rectangle for the ball;
- actual integer trial division for three distinct prime assignments;
- Gram-matrix constraints for all lattice isometries in the stated range;
- direct finite-ring unit quotients, rather than the producer's class-number product;
- factor-pair enumeration of AC=b²+d, rather than the producer's A,B loop;
- exact big-integer dyadic witnesses and involutive vertex-permutation checks.

For the finite-ring computation, a residue x+tω is a unit modulo f exactly when
its norm is coprime to f. The subgroup generated by scalar units and the six
Eisenstein units is enumerated explicitly. The standard conductor exact sequence,
with maximal-order class number one, identifies the class group with the quotient.
The quotient has exponent dividing two exactly when every unit square lies in
that subgroup. This independently tests the group exponent, not merely whether
its order is a power of two.

Both normal and optimized Python runs produced byte-identical receipts with
**480,652 live checks**. They additionally reproduced the producer's two
byte-identical **72,384-check** receipts. Explicit exception checks survive `-O`.
Declared independent ranges are:

- balls and reflections: k=1 through 256;
- Gram isometries: k=1 through 32;
- actual integer divisors: k=1 through 8, for three prime assignments;
- whole-order rotation membership: k=1 through 64, on stated finite boxes;
- finite-ring quotient groups: f=2 through 192;
- even-discriminant reduced forms: every even f through 192;
- dyadic witness exponents: a=4 through 1024;
- all involutive vertex permutations: 2 through 8 vertices.

These are finite controls, not substitutes for the uniform arguments in section 1.

Reproduce with the producer packet directory and, optionally, a repository clone
containing its pinned source commit:

    python3 -B independent_audit.py --packet PATH_TO_PACKET --repository PATH_TO_REPOSITORY
    python3 -B -O independent_audit.py --packet PATH_TO_PACKET --repository PATH_TO_REPOSITORY

Omitting `--repository` skips the five inherited-source hash checks, so the check
total and corresponding receipt field then differ by five. The mathematical
checks and frozen producer-payload validation still run. No repository writes
or network access are performed by this script.

## 4. Source and integrity verification

All five repository-source byte hashes match commit
`8730121c06121d741fa68fdb366bc1bbc75d00a4`, using read-only `git show`.
The inherited arithmetic chart, labelled reflection theorem, and native alphabet
were inspected directly. Their roles are represented correctly.

The following primary sources were independently accessed on 2026-10-21:

- [Keith Conrad, The conductor ideal of an order](https://kconrad.math.uconn.edu/blurbs/gradnumthy/conductor.pdf),
  Definition 1.1, Theorem 1.3, Example 2.1, and Theorem 5.3 with exact sequence (5.7).
  These support the conductor and class-number quotient, including the unit denominator.
- [Ernst Kani, Idoneal numbers and some generalizations](https://www.labmath.uqam.ca/~annales/volumes/35-2/PDF/197-227.pdf),
  Theorem 6, Remark 7(b), and Remark 18(a). These give the idoneality/group-exponent
  equivalence and the interior reduced-form obstruction. Theorem 12 supplies the
  broader classical context and confirms that the family is not a novelty claim.
- [MIT 18.785 Fall 2018, Problem Set 9](https://math.mit.edu/classes/18.785/2018fa/ProblemSet9.pdf),
  Problem 5, especially parts (b)–(c), records the reduction convention, inverse,
  and ambiguous boundaries, with an explicit extension to nonmaximal orders.

The repository's separate adopted qualification concerning global idoneal
completeness is preserved, not endorsed as a new proof or used as an input here.
The final audit verdict certifies only the frozen packet's stated elementary
bridge, conductor and isometry results, classical restricted classification,
finite receipts, provenance, and explicitly limited comparison claims.
