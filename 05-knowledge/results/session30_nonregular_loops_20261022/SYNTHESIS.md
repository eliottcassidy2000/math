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
16 and3. Within an explicit paired-switching family, isospectral Hamiltonian-count
ambiguity persists even with the full ordered score vector held fixed. At order8, equal spectrum AND equal path-tree scalar still allow
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


## Further structural results from the next proof/test cycle

### Exact flat-family constants, not fitted coefficients

For a tournament obtained by deleting one vertex from a doubly regular
n+1 tournament, [exact core contraction](session30-deleted-drt/DELETED_DRT_FORMAL_CONSTANTS.md)
gives

    κ1/e=1−5/(2n)+20/(3n²)+O(n^-3),
    κ2/e=1−5/(2n)+65/(8n²)+O(n^-5/2),
    y2=3/2−139/(24n)+O(n^-3/2).

These are proved asymptotics for the finite predictors, NOT the exact H count.
The previously numerically inferred vacuum constant29 is proved here by exact
integer kernels; the additional leaf core contributes4, giving33. An independent
[audit](session30-pointing-audit/deleted-drt/AUDIT.md) verifies the spectral
remainder bound and Euler-operator arithmetic.

### Converse symmetry and an exact blue-boundary map

[Unique completion and anti-involution conjugacy](session30-root/SELFCONVERSE_POINTING_IDENTITY.md)
prove that the number of self-converse almost-regular classes at even n equals
the number of self-converse regular classes at n+1. At n8 there are85 classes,
15 self-converse, hence50 converse-merged nodes within this restricted family.
This is not closure under metagraph wiring.

The [completion boundary flux](session30-root/COMPLETION_BOUNDARY_FLUX.md)
compares blue counts using the SAME signed half-path transfer states and a
changed endpoint acceptance condition. Completion does not preserveblue count:
of the15 self-converse n8 classes, four change by−2, seven by0, four by+2.
This is an explicit bridge to the previous half-path work, with centralizers
and endpoint information retained.

### Parity reduction and honest performance accounting

[Even/odd propagators](session30-parity-basis/README.md) reduce2364 symbolic terms
to553, or263 after signed graph-automorphism contraction equivalences. Palette
storage drops31396→12970 tensor slots. Those are exact finite reductions.
Naive sparse evaluation is slower than the existing optimized engine; compressed
dense timings are planner/cache-dependent, not a guaranteed4.27-fold speedup.

### New incoming work preserved, not silently audited

A concurrent source update at commit57dcd8e5 adds a proposed four-loop engine
and DRT23 data. It is preserved in Git. Its111-core and related claims are not
included in this session's independent audit scope or claimed as our findings.

### A larger sufficient flatness class

[The universality extension](session30-flat-universality/FLAT_ALMOSTREGULAR_UNIVERSALITY.md)
retains the same formal65/8 coefficient whenever S is almost-regular,
off-diagonal entries ofS² stay bounded, and Sd stays bounded, where d=S1.
The [independent audit](session30-flat-audit/AUDIT.md) further weakens the last
condition to the energy bound||Sd||²=O(n). The operator-norm hypothesis follows
from flatness rather than being an extra assumption.

The class is preserved by bounded-per-vertex score-preserving reversals.
Maximal disjoint triangle packings leave onlyO(√n) vertices, so this supplies
exponentially many distinct labelled perturbations with the same predictor
constants, including cases with linearly many reversed arcs. No exponential
isomorphism-class count or exact-H asymptotic follows from that statement.

### Runtime and scope correction from independent benchmarking

An [equal-planner-cap audit](session30-parity-audit/HYBRID_ADDENDUM.md) finds the
compressed engine slightly slower atn10 and20, about10% faster atn40 in its
controlled run. Tensor storage and symbolic reductions are exact; there is no
uniform runtime-speedup claim. Earlier larger apparent gains reflected cache
or contraction-planner behavior and are not credited to the algebra.

## Next precise obligations

1. Bound the omitted long-forest/large-core tail atz=1. Fixed-degree bounds and
   successful finite predictors do not discharge this.
2. Determine which additional higher-weight invariants repair the recorded
   order8 ranking failures, without fitting away counterexamples.
3. Prove a uniform exterior-plus-current-arrival-boundary bound for the
   source-owned Collatz bank. Its supply/contraction alone is insufficient.
4. Recover the missing six-order deficit-six and composed-bound certificates
   before promoting their reported closure in dependent canon.

The common structural lesson is boundary retention: eliminating internal
variables leaves an endpoint vector, a rooted class, or an exterior/boundary
quadratic term. Dropping that surviving data causes the concrete failures above.
This is a comparison of explicit elimination maps, not a claimed equivalence
between tournament counting and Collatz dynamics.
