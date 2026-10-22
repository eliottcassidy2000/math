# Session30: complete formal near-regular loop corrections

ACTIVE checkpoint,2026-10-22. The mathematical session has not yet ended.

## Main answer

The linear-forest expansion extends beyond regular tournaments by retaining
terminal leaves. For almost-regular even tournaments the correct grading is
weight=excess+number of terminal leaves. We classified and implemented all
five weight-one and26 weight-two cores, with no fitted parameters.

The missing first-order term in the prior path-only extension is the lollipop:

    L(z)=−2 diag(H−G)^T G1.

It vanishes for regular tournaments but can genuinely contribute at first order
in almost-regular families. An explicit deleted non-vertex-transitive blow-up
has a nonzero rational limit for nL; this is proved for the formal term without
assuming a full Hamiltonian-count asymptotic.

The complete first/second predictors combine all5/26 cores with the first two
inverse-factorial operators. Their formulas and runnable implementation are in
[the engine report](session30-decorated-engine/README.md); full derivation is in
[the independent theory report](session30-nonregular/nonregular_forest_extension.md).

What is proved:
- Exact forest identities for arbitrary tournaments and exact coefficients throughz4.
- Complete reduced-core inventories, automorphism factors and rational kernels.
- Uniform fixed-coefficient remainder O_L(n^-3) for each fixed forest degreeL.
- Uniform resolvent bounds and score-flow/regular-completion identities.

What is NOT proved: a uniform bound after evaluating the infinite expansion at
z=1, or a universal Hamiltonian-count error/ranking guarantee. That remaining
long-forest bound is a real gap, not a fitted constant.

## Complete finite validation and failure boundaries

An independently mass-certified census contains85 almost-regular order8 classes
and3230080 labelled tournaments. All85 H values are independently checked via
signed forests, in addition to subset DP. The second-order predictor:

- Improves absolute y error for73 of85 classes, worsens12.
- Reduces unweighted RMSE from0.0964332 to0.0478943.
- Reduces reversed distinct-H pairs from636 to465 out of3263 compared with
  first order; it is still NOT an exact ranking algorithm.

The second order also improves all17 selected larger/structured controls,
including nine fresh held-out cases, but worsens the tiny orders2 and4.
These universes are stated separately; no aggregate hides those failures.
See [complete census](session30-almost-regular/CENSUS8_REPORT.md) and
[held-out report](session30-almost-regular/SECOND_ORDER_REPORT.md).

## Structural information that spectra lose

[Two order6 witnesses](session30-root/ISOSPECTRAL_METAGRAPH_OBSTRUCTION.md) share
skew spectrum and blue count9, yet have H41 and45, and fundamental black degrees
16 and3. The paired switching construction also preserves their ordered score
vector. At order8, equal spectrum AND equal path-tree scalar still allow
H653 versus651. The next mixed invariants distinguish that pair but do not yet
rank it correctly. This is a concrete invariant hierarchy, not a numerical
analogy between unrelated exceptional numbers.

The [score-flow decomposition](session30-root/CENTERING_AND_REGULAR_COMPLETION.md)
is an exact complete-graph divergence/gradient decomposition. In almost-regular
tournaments, centering changesrho onlyO(n^-2), although endpoint corrections
enter the forest expansion atO(n^-1).

For unrestricted independent orientations, [exact random controls](session30-root/SCORE_SCALE_AND_RANDOM_CONTROLS.md)
prove Eκ=1 while Eρ tends toe. Thus omitting score corrections cannot have
vanishing mean absolute error on that ensemble. This does not assert thatκ
itself concentrates or that the regular-tournament limit is universallye.

## Collatz: a genuinely uniform relation bank, but not uniform coverage

[The arrival-boundary theorem](session30-collatz/REPORT.md) gives a completed
relation energy−2x²+3xy−3y² and an exact rank-one contraction with sharp factor4/9.
Two additional congruence families double the aligned-scale relation supply.
All original source IDs remain disjoint, with literal mod15/mod35 registers and
full clocks. Current-boundary and exterior interactions remain explicit in a
Schur sufficient gate. Negative pair credits are not spent twice or erased
from their surroundings.

Independent audit verifies811022 explicit checks and the all-scale algebra.
Uniform153N still requires a bound on exterior channels plus the current arrival
boundary. Numerous relations alone do not settle that obligation.

## Incoming deficit-six and composed-bound claims

The current synced repository did not contain the subsequent certificates for
the user-reported R6=z5+z6+z7 closure or composed bounds520738,1039055,3848613.
They are retained as user-reported, not refuted and not silently promoted.

[The integration note](session30-deficit-six/INTEGRATION.md) derives exact refined
known-family series and shows that, conditional on the reported closure, the
canonical degree-three-node recurrence has minimal recurrence order18 and sharp
onset19. “Degree-three” describes the counted metagraph vertices, not recurrence
order3. Almost-regular tournaments of even order≥4 are strong, so the SCC
padding mechanisms in those series are absent in this sampling family.

## Independent audit scope

[Loop audit](session30-loop-audit/AUDIT.md): independent general low-degree proof,
all labelled tournaments throughn6, all31 core types and automorphism factors,
connected-partition controls through8arcs, exact symbolic tensor agreement,
score-flow identities, and complete order8 census certification.

[Collatz audit](session30-collatz-audit/AUDIT.md): independent symbolic algebra,
source disjointness, actual stopped arrival seams and guards. No uniform-coverage
or arbitrary-source-entry claim is made.

All runnable packets use sibling directories as named here. Compiled binaries
are omitted; source and reproduction instructions are included. Mathematical
proof, finite evidence, inherited source and open obligations remain separate.
