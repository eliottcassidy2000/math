# Geometry, arithmetic content, and actual-source constraints

20 October 2026. Session map for the simplex/spheroid/cuboid investigation and the next Collatz and planar Jacobian targets.

**Outcome:** exact scoped mechanisms and explicit remaining gates. Neither the Collatz conjecture nor the planar Jacobian conjecture is proved. Novelty here means new synthesis within this investigation; historical priority is not asserted.

## Start with these proof packets

1. [Simplex, sphere and cube completion](../simplex_sphere_cube_completion_20261020/SIMPLEX_SPHERE_CUBE_COMPLETION.md): invertible geometric charts, exact observation-repair rank, full Boolean moment completion, and the failure of squared orthogonal kernels to compose coherently. [Independent audit](../simplex_sphere_cube_completion_20261020/review/AUDIT.md).
2. [Collatz rotation-defect ideal](../collatz_rotation_defect_ideal_20261020/LOCAL_ROTATION_DEFECT_IDEAL.md): local differences between cyclic rotations recover the exact cycle-denominator content. [One-swap extension](../collatz_rotation_defect_ideal_20261020/ONE_SWAP_TRINOMIAL_GATE.md), [audit](../collatz_rotation_defect_ideal_20261020/review/SCOPED_AUDIT.md).
3. [Inverse trace lattice](../inverse_trace_lattice_20261020/TRACE_HYPERBOLIC_PACKET.md): all Smith slopes, sharp integral repair index, and finite-sheet approximation precision. [Independent audit](../inverse_trace_lattice_20261020/review/AUDIT.md).
4. [Global selector and Euler gate](../inverse_trace_lattice_20261020/FINITE_SELECTOR_EULER_GATE.md): actual finite graph geometry converts a local selector into an exact omitted-point locus and intersection budget. [Independent audit](../inverse_trace_lattice_20261020/selector-review/AUDIT.md).
5. [Every source-line specialization](../inverse_trace_lattice_20261020/SOURCE_LINE_AND_IRREDUCIBLE_BOUNDARY.md): reduced irreducibility, immersed affine-line normalization, and compulsory collisions outside the nonproperness curve. [Independent audit](../inverse_trace_lattice_20261020/source-line-review/AUDIT.md).

6. [Collision-gap and singular-locus gate](../inverse_trace_lattice_20261020/COLLISION_GAP_AND_SINGULAR_LOCUS.md): an explicit unit certificate for the separation of colliding source parameters, and a finite map onto the graph singularities. [Independent audit](../inverse_trace_lattice_20261020/collision-review/AUDIT.md).

7. [Observable-algebra completion](../observable_algebra_completion_20261020/OBSERVABLE_ALGEBRA_COMPLETION.md): the smallest product-and-transition closure is exactly the function algebra of the observational quotient; a register can be state-observable while its separate output means fail to identify a probability law. [Independent audit](../observable_algebra_completion_20261020/review/AUDIT.md).

The original author notes sometimes retain their dated “review pending” handoff text. The linked separate audits and publication-review records give the completed review status without silently rewriting the reviewed source.

## 1. What the three shapes actually contribute

A centered regular simplex is polar-self-dual up to similarity. A Euclidean unit ball is self-polar for its fixed pairing. A cube is polar to a crosspolytope; its useful different self-duality is order reversal on the coordinate box. These are distinct operations.

There are nevertheless exact, productive transitions:

- A marked positive ray has simplex, spherical and cubical-section representatives. Keeping the marker preserves scale and makes the chart changes invertible. Discarding it can lose arithmetic valuation and integrality.
- A probability simplex maps bijectively to a spherical orthant by p_i ↦ 2√p_i. This carries the Fisher metric to the sphere, but does not turn every moment projection into a bijection.
- The cube's sign characters provide a complete coordinate algebra for laws modulo global sign. Their multiplication rule turns a moment matrix into exact atom-positivity conditions.
- An m-gon of complex root phases matches the circle's holomorphic moments through order m−1. The mth moment detects the packet. A triangle can therefore look circular to two observations while carrying a pole invisible to them.

There is no universal complexity ordering of these bodies. The relevant complexity is the number and kind of observables required by the intended consumer.

## 2. Collatz: replace a picture by its exact integer ideal

For a chronological shortcut parity word w of length q with p ones, let B(w) be its affine carry and D=2^q−3^p. Its rational cycle starts at B/D and has reduced denominator |D|/gcd(D,B).

Choose any connected graph on the distinct cyclic rotations. For every edge, strip any common outside prefix and suffix, leaving an equal-shape core difference δ. The proved identity is

    gcd of all edge |δ| = gcd(D,B).

A spanning tree suffices. Thus sparse local comparisons recover the full denominator obstruction; Bézout coefficients give independently checkable certificates. This is stronger than knowing that the cycle lies near a particular geometric shape.

For one adjacent swap of a zero-intercept balanced mechanical word, all nonexceptional differences become explicit three-term expressions. This gives a concrete near-balanced target. The strongest tempting guess is already false: some one-swap words have nonunit common divisor. What remains useful is to prove that the common divisor is a proper divisor of D whenever D>1 in an appropriate nontrivial positive-cycle family.

**Precise open conjecture (one-swap positive-cycle exclusion).** Let q≥4, 0<p<q, gcd(p,q)=1 and D=2^q−3^p>1. Set c_k=floor((k+1)p/q)−floor(kp/q). For every word w obtained by swapping cyclically adjacent unequal bits of c, conjecture gcd(D,B(w))<D. The probe supports every qualifying case through q=80; its broader universe contained 52,792 swap occurrences for 3≤q≤80. Rotations, duplicates and swaps returning to the balanced necklace require no extra exclusions. This remains an all-length conjecture, not a consequence of the finite probe.

This is a cycle exclusion program. Excluding cycles would still leave divergent trajectories. The separate all-orbit program needs genuine stopping/entry control; a finite bank of uniformly comparable homogeneous norms cannot force descent along arbitrarily long lawful expanding return blocks.

## 3. Jacobian: preserve the source, not just a root polynomial

Assume the stated local complete-DVR hypotheses: m escaping inverse sheets, one retained sheet, and integral first m−1 trace moments. If the first pole has order a, every Smith exponent of the rank-(m−1) trace lattice is −a. The sharp integral repair index is

    ceil(a(m−1)/2).

The unique finite root is approximated by a coefficient ratio, with Newton error orders at least (2^(n+1)−1)a. A trace-unimodular lattice containing the source coordinate can still fail multiplication closure, so integrality cannot be inferred from the repaired bilinear form alone.

For an actual polynomial Keller source in the primitive monic coordinate gauge, impose the explicit global low-trace hypothesis. Clearing denominators gives a finite graph hypersurface

    R = H Z^d + Hc₁ Z^(d−1) + … + Hc_(d−2) Z² + A Z + B.

With h=rad(H) and C=V(h), the exact omitted set and topological budget are

    O = V(h,A),
    #O = (d−1)(1−χ_c(C)) ≤ deg(h) deg(A).

Here B is a unit at the finite intersection scheme, and −B/A is the actual finite-sheet selector along C away from O. These conclusions use the actual source and étaleness, not merely a formally admissible polynomial R.

For every a, R(u,v,a) must itself be reduced irreducible, with finite immersed A¹ normalization. Noninvertibility forces a collision on every source line, and those collisions lie outside C. If C is irreducible, the cited self-intersection theorem gives #O≥d−1. In the selected degree-six low-trace sector this requires at least five distinct intersections, while the local rank-four trace lattice requires repair index 2a at each boundary component.

The last algebraic step makes those collisions concrete. Write U=(P(a,s)−P(a,t))/(s−t), V=(Q(a,s)−Q(a,t))/(s−t). Constant nonzero Jacobian J gives a polynomial Ω with

    J=P_x V−Q_x U−(s−t)Ω.

Thus s−t has explicit inverse −Ω/J in the collision quotient. This quotient maps finitely onto Sing(W), preserving the parameter a=Z. Moreover Sing(W)=V(R,R_Z) set-theoretically, and it lies entirely off H=0. Noninvertibility therefore requires this singular curve to meet every plane Z=a. Dominance alone does not suffice. Even surjectivity with a unit gap and finite fibers does not force properness, as the explicit punctured-line hostile demonstrates.

**Next exclusion target:** obtain an intersection upper bound, a forbidden specialization, or a source-compatible singular-locus obstruction contradicting these necessary conditions. No such universal contradiction is established here. The low-trace premise itself is not automatic for arbitrary Keller maps.

## 4. Broader mechanisms worth carrying forward

- **Consumer-relative completion:** add the missing row-space directions needed for a specified observable, rather than assuming every low-order relaxation reconstructs the state.
- **Arithmetic content:** retain the integer ideal or Smith data under a change of representation. An invertible real map need not be unimodular.
- **Multiplication closure:** additive/trace compatibility does not imply that an integral lattice is an algebra.
- **Composition data:** entrywise squaring of an orthogonal map preserves a one-step stochastic coupling but loses cross terms required for coherent composition.
- **Actual-source realizability:** local root packets and formal hypersurfaces must still arise from the specified polynomial source. Requiring every source-coordinate slice is a much stronger test than one Newton polygon.
- **Boundary accounting:** local poles, finite surviving sheets, omitted points and normalization branch counts have to agree globally. Delta invariants cannot replace branch counts in the Euler calculation.

These are precise transfer tests, not assertions that the underlying open problems are isomorphic.

## 5. An additional bridge: registers, character products and probability laws

For a finite deterministic state space, start with declared observables and close under transition pullbacks and pointwise products. The resulting algebra consists exactly of functions constant on states with identical complete observation records. This identifies the coarsest quotient supporting the specified observations and dynamics.

For binary affine transitions and individually observed Walsh characters, compute the smallest transpose-stable binary label space L. If its rank is r, the quotient has 2^r states and its full probability algebra has real dimension 2^r. But if K is the orbit of observed labels before additive closure, the linear output-statistics space has dimension only |K|. Exactly 2^r−|K| further independent expectations are required to reconstruct every law on the quotient.

The three-bit cyclic register shows the distinction sharply. Three successive bit observations identify a single initial state. Yet the separate timewise means, together with the constant, span only four functions on eight states. Uniform even-parity and odd-parity laws have identical output means at every time. Products of outputs recover the missing joint information.

This gives a precise extension of the shape program: a cube character quotient determines a simplex of observable laws, whose square-root representation is a spherical orthant. The binary rank, linear statistic count and full-law dimension are different quantities. This finite-state theorem provides no termination bound for the infinite Collatz state space.

## 6. Scoped repairs made during the session

- A finite positive Collatz cycle has rational odd-step ratio p/q strictly below log₃2. The latter is only a conditional limit as the minimum tends to infinity. [Product proof and correction](../collatz_rotation_defect_ideal_20261020/SCOPED_CRITICAL_SLOPE_CORRECTION.md).
- Generic/function-field degree four for complex polynomial Keller maps is already excluded by [Domrina's 2000 theorem](https://www.mathnet.ru/eng/im273). Quartic A4/S4 packets remain local method controls, not globally open cases. No exact first open generic degree is claimed.
- The source-line and irreducible-boundary corollaries explicitly use the theorems recorded in [Jelonek's versioned primary preprint](https://arxiv.org/html/2011.03472v3).

All research checkpoints are ordinary commits directly to main. Each packet includes readable proofs, finite controls, independent review and exact scope. Computation supports the displayed arguments; it does not substitute for an all-order proof.
