# Independent audit: niche recurrences and their boundaries

Source pin: `d3741fe6ed02f6072dbe6bfda2f3d04193d75ce0`.
Status: independent proof review and finite exact controls; not a literature-priority certificate.
The audit scripts import no producer code and use explicit exception checks, so
normal and optimized Python exercise the same tests.

## Reviewed mechanisms

### Translated poles and two different recurrence axes

The source convolution splits, past its finite span, into prefix sums, a
constant plateau, and suffix sums. Counting nonzero coefficients after signed
merging gives the stated affine growth law. The all-r moment recurrence has
one distinct factor per surviving site. A Vandermonde proof gives minimality,
and the leading Hankel determinant remains nonzero with arbitrary nonzero
signed weights. The site's being zero matters: it contributes an initial
spike, so an *eventual* recurrence can omit it. This convention was made
explicit during audit.

The all-time low-jet hostile is valid. The finite-difference weights annihilate
the first R+1 source moments. Translation preserves those vanishings by the
binomial formula at every time. The R+2 rational sites lie in distinct integer
translation orbits, so no cross-time collision removes their linearly growing
support. In the inherited fourth-order two-flow family, the unchanged low jet
does not decide the rational first-integral condition. The first differing
term is at degree R+3, consistent with the inverse coordinate change.

The h=4 rational-integrability equivalence is inherited from the cited
metabelian theorem, not newly proved or extended to unrelated flows here.

The finite cyclic-site corollary also holds: a time period must preserve the
spatial source, hence is divisible by its least spatial period d. Over Z/M,
one d-block adds its fundamental-period sum S_d to every coordinate. The
least return clock is therefore dM/gcd(M,S_d), not a clock derived from the
grand total alone. All3412 source vectors for M=2,3,4,5 were checked.

### Periodic valuation clocks and critical trees

For an h-residue event, probability1/d is possible exactly when d divides
2^h−1. The binary representation of (2^h−1)/d proves existence and uniqueness.
Its least possible period is ord_d(2). This classifies periodic valuation
predicates, not arbitrary randomness-extraction algorithms.

The ordered full d-ary tree formula has n=dm+1 and L=(d−1)m+1. Multiplying
the m internal and L leaf valuation series gives the joint size/division
clock coefficient; there is no extra tree-shape variance because those counts
are fixed. The simplified conditional mean and variance follow by total
expectation/variance, retaining the boundary correction. Tree extinction's
natural-density-one arithmetic consequence follows from finite cylinder lower
bounds. It is not a native Collatz descent theorem.

Two branching observers driven by the same valuation have correlated source
increments. Their Mersenne internal-event intersection uses the lcm clock,
not a product probability. The audited covariance is−2/3 for clocks2,3 and2
for clocks2,4. The synchronous frontier update must subtract one for each
active observer and freeze a stopped observer. This is not a statement that
two recursively observed trees stop at aligned depths.

### Arithmetic boundary recurrence and named multipliers

The zero-product formula correctly uses all-residue CRT counts first and
then subtracts tuples with a globally zero coordinate. It does not assume
that globally nonzero residues have independent CRT components. The spectral
bases are distinct by the largest-prime argument after cancellation. Their
nonzero leading polynomial coefficients prove the claimed minimal scalar
recurrence for composite n. Prime n is correctly excluded.

For arbitrary fixed iid weights, the truncated-valuation algebra splits into
local nilpotent components. Its displayed product is an annihilator, not
necessarily a minimal polynomial. The independent tests retain signed source
weights at zero and correlated CRT coordinates.

The36-modulus invisible-mode certificate hM=12h and h_36=0 proves equality
for every uniformly averaged future length; it is stronger than a finite
prefix match. Named multiplication by2 reveals a difference2. Divisor probes
P_d=sum_(d|g)m_g invert by the ordinary divisor Möbius formula. With total mass
known, tau(n)−1 probes are sufficient and dimensionally necessary among fixed
linear banks on arbitrary normalized real mass vectors. Adaptive discrete
state identification is outside that lower bound, as stated.

### Schur boundary carrier for cyclic substitution

The Cauchy expansion must retain the paired Schur boundary indices. The q
variable is the total nonnegative edge flow, not a signed flow or a class
weight. At fixed q cutoff every cross-edge partition size is bounded, so the
formal sums are coefficientwise finite. In C3 substitution, total balance
forces a common cross-block degree. Root-block rotation has its own fixed-flow
trace, evaluated with q cubed; it is not a full-automorphism quotient.

An independent literal enumeration of every nonnegative edge-flow assignment
of total at most6 on C3[T2,T2,T2] gives coefficients

    raw:                 1,0,0,8,12,6,41
    fixed root rotation: 1,0,0,2,0,0,5.

These match separately written two-variable Schur trace formulas. The full
LR carrier closure follows by splitting the external Schur polynomials and
multiplying the local factors; the mixed trace ordering M1 M3 M2 is correct.
The later stronger hostile uses a five-vertex tournament and its converse,
both with scalar series1/(1−q³), but boundary q³ coefficients878 and874.
The independent audit directly enumerates all degree-three assignments on its
ten arcs and uses explicit s_(3),s_(2,1) monomial multiplicities.

The actual attachment contexts C3[T,K1,T2] and C3[T^op,K1,T2] have identical
coefficients through q11 but q12 coefficients9470 and9462. A separate audit
uses no Schur or LR computation: it enumerates all cross-block matrices of
common degree at most4 and matches their divergences against literal local
edge-flow kernels. It evaluates240850 cross triples and obtains both full
prefixes independently. In
particular the8 degree-three flows expose the failure of the naive scalar
child-flow substitution.

### Attachment orbit-mass moments

The audit constructs all81 permutations of the depth-two cyclic wreath
 action and its three-element reflection centralizer directly. For binary
 iid sources with probabilities1/2,1/3,1/5, literal pushforward orbit masses
 agree with the recursive moments of orders0 through4, including the paired
 palette centralizer formula. Global converse acts freely on these attachment
 orbits, but biased sources have unequal paired masses. Orbit-count halving
 does not authorize halving each original source mass. For three nonidentical
child laws, the inclusion-exclusion correction is exactly a six-way overlap.
The32-symbol hostile has identical nine pair overlaps but collision
probabilities1534/16^6 and1536/16^6. Independent weighted orbit dictionaries
verify this difference. The entropy formulas were also reviewed: they subtract
expected log orbit size, not a fixed-width coding cost.

## Incoming operation-type controls

The fresh source repaired a previous operation name: the fixed-parent
one-bit selector is for **global converse**, not Hamilton-path-native
complement. The explicit strong source S=(2,12,9,1), path P=(3,0,1,2), has a
cyclic three-vertex parent, but its actual path-native complement makes that
parent transitive. The independent audit reproduces that distinction.

Incoming PORT_LAWS also gives equal ordinary cycle inventories with different
labelled transporters. Thus complete cycle counts cannot be silently promoted
to a complete realized joint-action carrier.

## Reproduction

Run these scripts normally and with `python -O`:

- `check_independent.py`:16703 checks,162 signed source families,3412 cyclic
  source vectors, literal tree
  clock compositions, periodic event selectors, conditional clock moments,
  and the operation-type hostile.
- `check_flow_boundary.py`:24 checks, literal degree-six flow universe and
  separately implemented Schur traces.
- `check_attachment_mass.py`:45 checks, explicit permutation actions and
  source-mass pushforwards and higher-overlap hostile.
- `check_arithmetic_boundary.py`:2219 checks, residue-based dynamic programming
  through modulus60, signed weighted recurrence tests, the all-order36
  certificate and divisor-probe inversion.
- `check_flow_context.py`:2 aggregate prefix checks,240850 evaluated
  cross-block triples, independently confirming the actual q12 separator.

The stated counts describe finite test universes. Universal claims depend on
the displayed proofs. All five final proof documents, including their bounded cyclic-clock and
probe-restoration extensions, have been reviewed. The audit is frozen with
18993 independently implemented checks in each mode. No mathematical blocker
was found; no unrestricted open-problem resolution is inferred.

## Synthesis review

The root finite-dimensional observation theorem is correctly restricted to
linear operations and linear encodings. The annihilated-observation kernel is
invariant under the generated unital algebra. Its quotient is a lower bound
for an exact linear encoding. A scalar sequence only needs to obey the
observable operator's minimal polynomial; the report does not assert that
every initial vector attains that order. Nonlinear mass hierarchies and
infinite stopping grammars are explicitly kept outside this finite-dimensional
statement.
