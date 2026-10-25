# Session60 final independent audit: exact Picard rank and Boolean reflection

**Gate: PASS.** For the specified hourglass arrangement, the outer six-root K3
has geometric Picard number **19**, and its seven-root minimal resolution has
Picard number **58**. The modified Boolean complement is realized by the
central base half-turn; it preserves no Fano marking, pairs the thirty quotient
kernels, and its quotient has geometric genus **4**.

No genericity assertion, conjectural Tate input, primitive-lattice claim,
Fano-incidence promotion, or arithmetic-square construction is used. The
pre-existing Session60 deck/Hodge/central-Kummer audit remains an explicit
inherited dependency; this audit independently closes the new lower and upper
rank arguments and the incoming reflection connection.

## 1. Audited inputs and exact controls

Inputs are `session60-outer-picard/{LOWER_BOUND.md,check_lower_bound.py}`,
`session60-picard-upper/{INDEPENDENT_UPPER_BOUND.md,check_counts.py}`, and
`session60-reflection/{BOOLEAN_REFLECTION.md,check.py}`. The inherited inputs
are `session60-deck-k3/HODGE_DECOMPOSITION.md` and
`session60-hourglass-symmetries/RESULT.md`, already reviewed in
`session60-audit-deck-kummer/AUDIT.md`. Exact source bytes are pinned in
`FROZEN_SOURCE_SHA256SUMS`; the published paths there are relative to this audit directory.
Original capture paths are retained in the publication source receipt.

`independent_check.py` does not import any packet checker. It executes **5,316
explicit exception gates** in both ordinary and optimized Python, with
byte-identical JSON. It independently checks:

* rational and mod-5/7/11/13 line incidences, and the rank-18 Gram determinant;
* all four required finite-field counts, using quadratic-norm characters in
  F_(p²), also verified against direct field-square membership;
* quartic coefficients, factorizations, nonintegral normalized root sums,
  and discriminant square classes (plus extension-degree 1–8 controls);
* the elliptic coordinate change, section, independent invariant formulas,
  all singular-place valuations including infinity, and denominator growth;
* all 30 labelled Fano incidence structures, all 15 reflection pairs, and
  all 168 elements / 21 nonidentity involutions of GL(3,F2);
* the explicit Boolean marking and canonical residue determinant ratio.

The original lower, upper, and reflection checkers were also rerun in normal
and `-O` modes. Each output pair agrees. The lower source has 87 active gates.
These are exact arithmetic certificates, not machine-checked proofs of the
geometric theorems. A deliberately corrupted count expectation is rejected
in both modes; see `negative-control.normal.txt` and
`negative-control.optimized.txt`.

## 2. Good reduction and the eighteen classes: PASS

Let Z be the minimal resolution of

    w² = E(E²−x²)(E²−z²)(E+x+z).

The exact incidence calculation gives nine nodes and two ordinary triples.
Every nonzero three-line determinant is a unit over R=Z[1/6], so intersections
and tangent directions do not collide there. All relevant centers are rational
sections. Blow up the two triples and include their exceptional curves in the
branch; then blow up the nine original branch nodes and the six new ones,
excluding these fifteen exceptional curves from the branch. The final branch
is eight disjoint smooth relative rational curves. Its double cover is smooth
and projective over R because 2 is invertible. This constructs a simultaneous
smooth model explicitly, rather than presuming good reduction from a singular
plane equation or from stable point counts.

The branch line bundle transforms as L' = pi*L−floor(m/2)E. For m=2,3,
K'+L'=pi*(K+L)=0. The covering formula and
f_*O=O plus O(K_base) show trivial canonical bundle and h¹(O)=0 on each fiber:
the base is a blowup of P² and Serre duality applies. The nonempty reduced
branch makes the degree-two function-field extension geometrically connected.
Thus the fibers are K3 surfaces, including at 7 and 13. This argument excludes
2 and 3; it does not assert anything about those primes.

A base exceptional curve at a branch node lifts to a genus-zero curve branched
at two rational points, so it is P¹ over the field itself, without requiring a
square root of an unspecified unit. Each triple contributes a rational central
ramification P¹ and three such rational arms, all intersection points rational:
a split D4 configuration. Adjunction on the K3 gives self-intersection −2.
The resolution fibers therefore contribute q and 4q extra points respectively,
for a total correction 17q. In particular, the D4 components are not being
silently assumed split from their graph alone.

The pullback h of a plane line has h²=2 and is orthogonal to all exceptional
curves. Their Gram matrix is <2> direct-sum A1(−1)^9 direct-sum D4(−1)^2,
with determinant **−2^14**, independently recomputed. Hence all eighteen
classes are independent and rational. No saturation, primitiveness, or claim
that they exhaust the Néron–Severi lattice is needed. Their nondegenerate
18-dimensional cohomological span suffices to split off (T−p)^18.

## 3. Frobenius, functional sign, and roots of unity: PASS

The point formula is

    #Z(F_q)=q²+18q+1+sum_(x,z) chi((1−x²)(1−z²)(1+x+z)).

The affine chart accounts for q² plus the character sum; the branch line at
infinity contributes q+1; the preceding resolution contributes 17q. The
independent counts are:

| p | #Z(F_p) | #Z(F_(p²)) | residual trace | residual square trace |
|---|---:|---:|---:|---:|
| 7 | 180 | 3384 | 4 | 100 |
| 13 | 408 | 32088 | 4 | 484 |

Here Frobenius is the convention with untwisted H² eigenvalues of complex
absolute value p, giving the standard point-count trace formula. On the
orthogonal complement of the eighteen classes, F/p is an isometry. If its
determinant were −1, degree-four reciprocity would force the middle
coefficient to vanish. Newton's identity gives b=(a²−A)/2 equal to −42 and
−234, so the sign is **forced positive**, not guessed. Thus

    P_7(T)  = (T−7)^20 (T²+10T+49),
    P_13(T) = (T−13)^20(T²+22T+169).

For the remaining normalized roots beta,beta^−1, their sums are −10/7 and
−22/13. A root of unity and its inverse have algebraic-integer sum, and a
rational algebraic integer is integral. These nonintegral sums exclude roots
of unity exactly, not numerically. No power of either remaining root can
become a root of unity. Therefore the geometric algebraic-eigenvalue upper
bound is 20 over every finite extension. Counts at 5 and 11 are not used to
infer a functional-equation sign from a vanishing middle coefficient.

## 4. Tate / Artin–Tate / Brauer / specialization: PASS

**Tate applicability.** Z_7 and Z_13 are smooth projective geometrically
connected K3 surfaces in characteristic at least 5, so Charles's Corollary 2
applies without a degree-of-polarization restriction. It identifies divisor
classes with Frobenius invariants. There is no unproved semisimplicity step:
NS over the finite field is nondegenerate by Hodge index, hence splits off
orthogonally. Any eigenvalue p on the complement would give an actual
p-eigenvector there, contradicting Tate. Consequently multiplicity 20 really
is the arithmetic Picard rank. The geometric rank is also 20 by the root-of-
unity upper bound. The finite-order action of Frobenius on geometric NS has
all eigenvalues 1, so it is the identity over Q. The upper packet now states
this orthogonal-splitting argument explicitly.

**Artin–Tate normalization.** Milne supplies the implication from Tate to the
full Artin–Tate formula and finiteness of Br. For a K3, chi(O)=2, Pic^0=0,
and NS is torsion free. If Phi is the characteristic polynomial of F/p,

    |disc NS| #Br = p lim_(T→1) Phi(T)/(T−1)^rho.

The exponent is one, not zero or two. At p=7 and 13 the limits are 24/7 and
48/13, giving absolute discriminants times Brauer order 24 and 48. Hodge
index gives sign negative in rank 20.

**Square-order issue.** Milne's original skew-symmetric-pairing argument alone
would give square or twice square and is insufficient here: the obstruction
itself is a factor of two. This audit checked Liu–Lorenzini–Raynaud Theorem 1
and their 2018 corrigendum. The correction inserts an index **square** in an
intermediate formula and retains the square-order conclusion. Thus the
geometric reduction discriminant classes really are **−6 and −3**.

**Finite extensions.** Arithmetic NS is full rank in geometric NS; even
without using integral equality, its finite index changes discriminant only
by a rational square. No new divisor eigenvalues arise after extension,
because the two omitted normalized roots are not roots of unity. Therefore
these are geometric discriminant classes, not values that disappear on a
splitting extension. The independent recurrence controls degrees 1–8 as a
check only; the preceding eigenvalue argument proves the all-extension claim.

**Specialization and primitiveness.** For the explicit smooth proper model,
specialization of geometric NS is an intersection-preserving injection. One
may first extend the characteristic-zero field to define a chosen NS basis;
the geometric special fiber is unchanged. If the characteristic-zero rank
were 20, both images would be finite-index sublattices of rank-20 lattices.
For an index m sublattice, det(L')=m² det(L). Thus both reduction square classes
would equal the characteristic-zero one, contrary to their ratio 2 being
nonsquare. **A primitive embedding is not required or asserted.** Unknown
indices, including possible p-primary indices, cannot repair this discrepancy.
This proves rho(Z_Qbar)<=19. Picard rank is unchanged on extending the
algebraically closed characteristic-zero base from Qbar to C.

### Primary-source pins inspected on 2026-10-25

* [Charles, *The Tate conjecture for K3 surfaces over finite fields*, Corollary 2,
  printed p.1](https://arxiv.org/pdf/1206.4002): unrestricted K3 Tate theorem
  for finite fields of characteristic at least 5. Do not substitute the more
  restrictive polarization-based Theorem 4 for that corollary.
* [Milne's author clarification of *On a conjecture of Artin and Tate*](https://www.jmilne.org/math/articles/1975a.html):
  Theorems 4.1/6.1, full Brauer finiteness and Artin–Tate implication, plus the
  warning about skew-symmetric versus alternating pairings. The original
  scanned paper link was opened but its screenshot rendering was unavailable;
  this audit does not claim a successful page-image inspection.
* [Elsenhans–Jahnel, *On the computation of the Picard group for K3 surfaces*,
  §4.2–4.3, printed pp.4–5](https://wwwuser.gwdguser.de/~jjahnel/Arbeiten/pic_dst_cam2a.pdf):
  the exact normalized K3 formula above; §2.2 also states the
  intersection-preserving geometric specialization injection.
* [Liu–Lorenzini–Raynaud, *On the Brauer group of a surface*, Theorem 1,
  printed p.673](https://dinolorenzini.franklinresearch.uga.edu/sites/default/files/inline-files/papers/LLR2InventionesPublished.pdf),
  together with [the 2018 corrigendum, corrected Corollary 3, printed p.2](https://arxiv.org/pdf/1804.11158):
  the required square-order theorem with its correction accounted for.
* [van Luijk, *K3 surfaces with Picard number one and infinitely many rational
  points*, §2 and Proposition 2.2 / Corollary 2.3](https://arxiv.org/pdf/math/0506416):
  the specialization/eigenvalue bound and finite-index discriminant identity.
  The present application is to rank 20, not a misuse of a rank-one-specific
  conclusion.

## 5. Nineteenth class and transfer to 58: PASS

The birational substitution X=−(1−t²)z, Y=(1−t²)w gives

    Y²=(X−d)(X+d)(X−d(t+1)), d=1−t²,
    P=(dt,d²).

Both identities were independently verified. The model has I2 fibers at
0,−2, I0* at 1,−1, and I2* at infinity. The infinity substitution has weights
4 and 6; its valuations (2,3,8) are essential. At every singular place
v(Delta)<12 certifies minimality; there are no additional discriminant zeros.
The trivial-lattice rank is 2+1+1+4+4+6=18. The Kodaira and Shioda–Tate inputs
are precisely [Schuett–Shioda, *Elliptic Surfaces*, §4.5 Table 1 and Corollary
6.13](https://arxiv.org/pdf/0907.0298); the characteristic-zero hypotheses hold.

At smooth t=2 an integral translation gives y²=u³−36u and Q=(−3,9), whose
double is (25/4,−35/8). The exact duplication identity yields
v_2(u(2R))=v_2(u(R))−2 whenever v_2(u(R)) is a negative even integer.
Starting at 2Q, this proves infinitely many distinct doubles. It is not merely
an experimental denominator list or an assumption that a specialization map
is injective. Torsion would specialize to torsion on the smooth fiber, so P
is nontorsion and rho(Z)>=19. This meets the upper bound exactly.

The inherited actual degree-16 symplectic quotient correspondence identifies
rational transcendental Hodge structures of Z and the outer six-root K3 S,
with pairing scaled by 16, hence rho(S)=19. The already audited seven-summand
Hodge decomposition, central rank20, and D12 transitivity give

    rho(Y_hourglass)=20+6*19−76=58.

Neither an integral lattice isometry nor Kummer status is transferred by that
argument. The area-six specialization yields the signed 3–4–5 triangle and,
after doubling, the fixed-gap triple (1,25,49); it does not solve a magic-square
or three-centers problem and does not infer global CM from j=1728 at one fiber.

## 6. Incoming Boolean reflection and quotient genus: PASS

The incoming source was read at exact commit
`82691d391b4beef9682ef956526391bef5f91948`, path
`05-knowledge/results/reflection_bridge_20261025/seven.md`, not at an unpinned
later branch. The captured source is `incoming-seven.at-82691d39.md`, SHA256
`93234b5f14b6f364a90a2e618dbf0a1cf77e87293e76886326be48ab86186e72`.
Its §§2–4 assert exactly the boundary contraction and modified complement
used by the parent. Its §1 explicitly rejects the second-level Keller
reflection, and that boundary is retained.

Mark boundary -> E; singleton masks 1,2,4 -> E+x,E+z,E+x+z; complements
6,5,3 -> E−x,E−z,E−x−z. Modified complement is then exactly the branch
permutation (0)(1 2)(3 4)(5 6). The projective transformation
(E,x,z)->(E,−x,−z) permutes normalized forms literally, so the pure root
permutation is an honest involution and lifts to the minimal resolution.
It is the central half-turn (determinant +1 on the base), not one of the
orientation-reversing D12 mirrors with determinant −1.

For any nonidentity order-two T in GL(3,F2), N=T−I satisfies N²=0 and
rank N=1 in dimension3. Its fixed space has dimension2, giving exactly three
nonzero fixed vectors. A permutation with one fixed point therefore cannot
be a collineation under **any** Fano marking. Kernel stabilization is equivalent
to preserving the induced quotient incidence: two surjections to F2³ with the
same kernel differ by GL(3,F2). The independent enumeration of all thirty
labelled Fano planes finds no fixed one and exactly fifteen pairs. The actual
geometric automorphism conjugates the corresponding kernels, so paired
quotients are isomorphic; it does not descend to an automorphism of a fixed
marked quotient whose kernel it moves.

The determinant twist is independently verified from the seven form matrix
and its four-dimensional relation kernel: the root-permutation determinant
and relation-space determinant are both −1, so their residue quotient is +1.
The canonical representation has four +1 and three −1 eigenvalues.
Consequently the quotient resolution has p_g=4. This last step can be checked
locally without a blanket rational-singularity citation: an involution is
analytically linearizable; a fixed curve has local action (u,v)->(−u,v), where
invariance of f du∧dv makes f odd in u and the form descends regularly through
U=u². An isolated fixed point has action (−u,−v), giving an A1 quotient
singularity with crepant resolution. Invariant canonical forms therefore
correspond exactly to canonical forms on the resolved quotient.

This is an isomorphism of explicitly marked **sets with involution**, not of
Boolean orders, divisor arithmetic, Fano incidence, or the historical
second-level Keller tree. The packet makes no claim about the quotient's K²,
minimal birational type, or the full automorphism group. The contrast with
invariant genera 1,3,5 for deck involutions is valid because this involution
permutes branches instead of changing only root signs.

## 7. Reproduction and release scope

Run `python3 independent_check.py` and `python3 -O independent_check.py`;
Python 3.10+ and SymPy are the only runtime requirements. No source files,
Git checkout, network calls, or absolute input paths are needed by the checker.
Compare the resulting JSON with the two archived outputs. Source replays are
separate evidence; in particular the lower checker writes its own result JSON
when rerun. No bare assert statements are used for any audit gate.

The exact-rank and reflection conclusions are cleared for integration with the
already audited Session60 packets. The audit does not create a new claim of
historical novelty or erase inherited-dependency boundaries.
