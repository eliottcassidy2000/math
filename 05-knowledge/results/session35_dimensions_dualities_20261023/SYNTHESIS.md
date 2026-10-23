# Session35: dimensions, dualities, and what survives a change of representation

2026-10-23. This session repeatedly alternated repository recovery, proposed bridges, exact tests, hostile examples, and independent audits. Results below distinguish elementary proved applications, inherited/cited theorems, and finite certificates. No discovery-priority claim is made merely because a connection was newly recovered in this session.

The unifying question is not whether the same integer occurs twice. It is: **what map preserves the operations and the target observable, and which discarded coordinate prevents the next inference?**

## 1. Direct progress on the blue/black tournament lens

[Paley two-design balance](recursions/PALEY-TWO-DESIGN-BALANCE.md) gives an exact new-to-this-session characterization of the nine Hamiltonian-path orbits of P7.

The two cyclic Fano planes partition its fourteen directed triangles. For a marked Hamiltonian path, let u+ and u− count lines avoiding every path edge. The complete nine-orbit certificate gives:

- (1,1): one orbit,21 paths, BLUE.
- (2,2): two orbits,42 paths, BLUE.
- (2,3): three orbits,63 paths, BLACK.
- (3,2): three orbits,63 paths, BLACK.

Thus BLUE is exactly balance between the two available path-preserving line-flip cube dimensions. Equivalently, BLUE paths admit2 or4 directed-triangle reversals preserving the path, while BLACK paths admit5. All such moves preserve regularity. All1932 path/plane/cube-state instances were checked; every nonorigin in each unused-line subcube leaves P7 for the H171 regular class.

This is not a claim that a full Fano cube has only one Paley state: it has16. Nor are these triangle moves the original off-path-complement operation. A same-order H171 example has both BLUE and BLACK paths of degree5, proving the criterion is P7-specific. A uniform survivor is proved for Paley paths reversed by an affine reflection: their admissible-triangle degree is even. The converse already fails at Paley11.

## 2. The +2 and doubling operations really do share a carrier

[Recursion formulas](recursions/RESULT.md) put source–sink extension E and skew-Sylvester doubling D on the same marked skew spectral data

    p_M(z)=det(zI−M),  F_M(z)=1ᵀ(zI−M)^−1 1.

E acts by a Möbius transformation on F; a Cayley change of coordinate linearizes its iteration. D acts by the quadratic substitution w²=(z²+1)/2 with an explicit prefactor. The operations do not commute.

[Context and ports](recursions/CONTEXT-AND-PORTS.md) identifies this pair exactly with ordinary adjacency characteristic data. Hence adjacency cospectrality is preserved by every word in E,D. This is a closure theorem, not a stronger-than-adjacency spectral invariant.

A minimal six-vertex hostile has identical data but Hamiltonian counts29 and25. Both have six directed triangles and six five-cycles; only one has a disjoint triangle pair. That single compatibility difference contributes4 through the inherited odd-cycle formula. The exact square-zero, vertex-labelled port resolvent restores the missing repeated-vertex exclusion and recovers all induced Hamiltonian counts, at the cost of exponential label state.

## 3. A genuine theta lift, and its limit

[Zeta/theta construction](zeta/RESULTS.md): Q=I−M² is positive. Finite spectral zeta has positive values at negative integers, so multiplying by Gamma creates poles, not Riemann-style trivial-zero cancellation.

Adding the actual infinite lattice gives

    Theta_D(M)(t)=Theta_M(2t)^2,
    Theta_Q(t)=det(Q)^−1/2 t^−n/2 Theta_(Q^−1)(1/t).

The completed Epstein function reflects s↔n−s; the uncompleted Epstein zeta has simple zeros at negative even integers. The inverse quadratic form, not tournament reversal, supplies the dual. Normalizing the center to1/2 changes the negative-even spacing to2/n. Only dimension1 recovers exactly the Riemann normalization.

The infinite doubling tower is a different object: its zeta is Z_Q(s)/(1−2^(1−s)), with an unavoidable pole at1 and log-periodic heat asymptotics. Neither construction establishes RH.

[Conference collapse](zeta/CONFERENCE_COLLAPSE.md) is a strong reconstruction hostile: every switch of a skew-conference matrix has exactly Q=nI. At order8 this gives identical full theta functions and spectra but H189 or333. The difference144 splits into three contributions of48 from seven-cycles, disjoint triangle pairs, and triangle–pentagon pairs. The older189/333 and conference classification are explicitly credited; the added lens is the full entrywise Q/theta collapse.

## 4. Dimension57 leads to a precise E8 completion theorem

[Intermediate algebra](lie/RESEARCH_NOTE.md): E7½ is E7⋉H57, with dimensions133+56+1=190. Every invariant symmetric bilinear form on it factors through E7, so the Heisenberg radical is invisible to the form. A naive ordinary inverse-metric Sugawara construction is unavailable.

The initial [graded completion](lie/COMPLETION_ADDENDUM.md) and [equality rigidity](lie/EQUALITY_RIGIDITY.md) led to the stronger [ungraded theorem](lie/UNGRADED_METRIC_COMPLETION.md):

Over C, among ALL finite-dimensional Lie algebras with a nondegenerate invariant symmetric bilinear form and an injective Lie embedding of E7⋉H57, the minimum dimension is248; every equality ambient is simple E8. No grading is assumed. The centralizer/multiplicity argument excludes247 and forces a centralizer sl2 at equality; the original contact grading is then recovered as a conclusion. Adding just the grading coordinate gives a nonunimodular191-dimensional algebra; cotangent completion gives a larger380-dimensional positive control.

The general ungraded lower bound is dim(s)+2dim(V)+3 for simple complex s and a nontrivial irreducible symplectic module V. Equality forces simplicity and the contact grading when it exists. An explicit symplectic family attains the bound for every odd Heisenberg dimension at least3, but arbitrary Levi/module pairs need not attain it. The [split-real corollary](lie/SPLIT_REAL_COROLLARY.md) identifies the equality ambient as split E8(8) for the specified split E7 input. The invariant bilinear form may be indefinite; it is not a Hilbert-space positivity certificate or a mass-gap estimate. No marked-embedding uniqueness or Yang–Mills continuum conclusion is claimed.

The invariant-form obstruction also excludes E7½ as the full weight-one algebra of a simple strong-CFT-type VOA. It does not exclude an embedded subalgebra or a construction with different hypotheses.

## 5.57 really bridges arithmetic, projective geometry and Heisenberg actions

[Arithmetic](arithmetic/ARITHMETIC.md) starts from57=3·19=2^5+5². Quadratic norm filters constrain Leyland–Blum candidates; the bounded census below10^12 contains111 parameter pairs and12 Blum values. CRT principal-root choices give explicit4/2/1 sign-torsor normal forms compatible with the literal marked deletion-bit moves. This is a fibrewise map, not a global metagraph identification.

[Singer geometry](arithmetic/SINGER_LINE_FLIPS.md) constructs the cyclic (57,8,1) difference set and its CRT transport. For odd projective order q, line incidence over F2 has kernel the all-ones line mask. Consequently a regular starting tournament's full-line-flip family has exactly two regular members, itself and converse. The Fano q2 regular cube does not extrapolate to q7.

For the directed-triangle observable, every quadratic line-character survives when q is odd. At q7 the exact deterministic continuation count is2^56 and the linear translate-span dimension is1597. These are not Hamiltonian-count statements.

[Actual tournament Heisenberg lift](ACTUAL_TOURNAMENT_HEISENBERG_LIFT.md): every odd tournament has a nondegenerate alternating form on its sum-zero lattice over Q. There are1024 labelled C57⋊C3-invariant Cayley orientation rules compatible with the chosen multiplier7. Each gives a Heisenberg57 action.

[The E7 obstruction](SINGER_HEISENBERG_E7_OBSTRUCTION.md): the distinguished cyclic augmentation action has order3 trace−1, whereas E7's minuscule56 permits only−25,−7,2,20 (or56 for identity). Thus this action cannot extend to E7, despite the matching dimension and symplectic structure. Forgetting the finite action removes that obstruction.

[Fourth-order boundary](FOURTH_ORDER_E7_BOUNDARY.md) makes the distinction sharper: the standard E7 quartic-stabilizer theorem excludes any GL-equivalent Freudenthal quartic invariant under this particular C57 action, even without fixing a compatible symplectic form. Invariant quartics of other types are not excluded.

## 6. Integer versus formal dimension has an exact positivity horizon

[Signed Burnside and positivity](DIMENSION_POSITIVITY_HORIZON.md): the cycle trace t^(number of permutation cycles) has signed average binomial(t,k). Its S_k Gram matrix is positive semidefinite exactly when t is a nonnegative integer or t≥k−1 (for real t≥0).

A noninteger t≥0 first fails at order floor(t)+2. The formal dimension7/2 passes every fourth-order permutation-trace check but fails at fifth order. Rational, algebraic irrational and transcendental parameters in the same interval share this horizon. This does not prohibit noninteger quantum dimensions in braided categories, where the symmetric-group trace assumptions change.

[Augmentation lattice](recursions/AUGMENTATION-LATTICE.md) retains another arithmetic layer: det(BᵀMB)=det(I+2A), the inherited Pfaffian square realized as an integral alternating volume. It is not the signed Rédei count. Even the full integral Smith profile fails to determine H, already at five vertices with H11 versus5. Through odd order7 the census found no strict Smith refinement beyond the Pfaffian product; it explicitly does not claim one.

## 7. Reciprocal Virasoro parameters and mirror bases meet exactly

[Coulomb/mirror duality](lie/COULOMB_MIRROR_DUALITY.md) uses c(t)=13−6(t+t^−1). Distinct positive mirror bases x^(1/x)=y^(1/y) have

    x=t^(1/(t−1)),  y=t^(t/(t−1)),  t=y/x.

Exchanging the bases is exactly t↔1/t, so c is an actual symmetric coordinate, not a verbal analogy. For rational c<1, square discriminant (c−1)(c−25) gives two algebraic bases; nonsquare discriminant gives two transcendental bases by Gelfond–Schneider. The coalesced t1 endpoint is separately excluded.

[Powering](lie/COULOMB_POWERING_ADDENDUM.md) descends to Chebyshev polynomials and preserves this discriminant square class. At t=phi²,

    c_m=13−6 Lucas_(2m),
    discriminant(c_m)=180 Fibonacci_(2m)^2.

The recovered cat-map trace and path rooted-forest determinant give these two quantities exactly. The charges−5,−29,−95,−269,... all produce transcendental mirror pairs.

[Elliptic charge dynamics](charge-dynamics/ELLIPTIC_CHARGE_DYNAMICS.md) supplies a contrasting exact orbit at the formal parameter38/5: the reciprocal parameter is algebraic on the unit circle, its angle divided by pi is transcendental, and all positive integer powers produce distinct rational charges equidistributed with the arcsine law on(1,25). Their reduced denominators are10^m/2. This does not prove density for the dyadic subsequence, and does not construct VOAs at those charges.

Rational central charge is not rational VOA: simple Virasoro L(−5,0) is nonrational, but the tensor of the (2,5) and (3,5) minimal models is rational at the same charge−5. Ising provides a rational VOA with irrational quantum dimension sqrt2. No equivalence of these arithmetic predicates is asserted.

## 8. Repository recovery repaired a generating function

[Good-cut recovery](arithmetic/GOOD_CUT_RECOVERY.md) found the proposed SCC formula already proved in THM340/349. Recovery prevented a false novelty claim and corrected an unpublished guessed strong count: s6=903.

It also exposed a genuine THM341 algebraic error: the positive-good-cut bivariate GF omits1/(1−x), already wrong at n4,k2. The missing factor restores the arbitrary leading singleton-SCC prefix. The repaired formal three-marker series tracks vertices, good cuts and nontrivial SCC count; fixed-k columns have degree floor(k/2) and a proved conditional asymptotic morphology. The series is formal, not convergent analytic near generic positive arguments. Black path-orbit deficit and good-cut count are different observables.

A separate narrow THM1926 correction changes Paley adjacency eigenvalue modulus to sqrt(p+1)/2; sqrt(p) belongs to the skew spectrum. Neither repair invalidates the surrounding structural toolkit.

## 9. What the unification does and does not say

[Observable-completion framework](OBSERVABLE_COMPLETION_SYNTHESIS.md) states the exact common principle: a sufficient quotient must preserve the chosen observable under every allowed continuation. Its Walsh specialization computes deterministic state counts and linear span dimensions separately.

The session produced explicit maps, minimal completion results, sharp counterexamples, and qualified mathematical connections. It did not prove RH, a4D Yang–Mills mass gap, Collatz, general lonely-runner, or complete metagraph realizability. The next substantial frontier is to retain enough path/port/incidence information to make new wiring or enumeration statements without paying for an unnecessarily full labelled state space.

All substantive lanes have independent mathematical audits and replayable scripts. Finite universes and untested generalizations are stated in their own notes. Publication receipts separately record actual Git state; this synthesis does not itself certify a push.
