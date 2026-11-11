# Careful reading of the Wilson Yang–Mills source

Repository: yang-mills-group/OpenAI-YM-20261028. Exact commit: 3456e5e75cf8d4cd558955a5c7638fba2d68f464.

This is mathematical reading of trusted results, not a build, proof audit, or reverification. Twenty-three source files were read; exact blob identities and immutable links are in SOURCE_PIN_MANIFEST.json. External source code is not vendored in this packet. Proposed bridges below are deductions or research directions, not assertions already proved by the source. The accompanying lanes connect these interfaces to earlier repository mechanisms.

## 1. Exact object and quantifier map

`YangMills.JointExistence` in ChallengeSpec.lean has the order

∀ gauge group g with simple nonabelian real Lie algebra, ∃ one Wilson family C, [regularity C and ∀ strictly increasing volume sequence p, ∃ strictly increasing refinement q, ∃ hierarchy T, T.Properties(p∘q)].

C is fixed before the volume sequence. T, its field quotients, selected character, OPE reader, and phase measure are inside the refinement quantifier. This is substantial simultaneous existence, but not identification of different refinements. Do not exchange ∀p∃q∃T for ∃T∀p∃q.

`YangMills.constructed` in Construction.lean additionally fixes an arbitrary positive blocking-kernel penalty b before constructing C. It then quantifies over an extraction order with mesh and period coordinates and supplies one table S, one Final F, one every-phase Almanac L, one Reader Q, contact devices M, and scaling X. All retained conclusions refer to these same witnesses. In particular M.basis = Q.port. The ordered limit retained has full mesh sequence, not merely the mesh component of the requested refinement. The arbitrary positive penalty is an internal construction parameter, not a theorem that theories at different penalties are equivalent or have a common physical mass.

The gauge group is a closed subgroup of a finite unitary matrix group. Its Lie algebra is specified exactly by the one-parameter exponential condition. Simplicity concerns this real Lie algebra and includes nonabelianness. No connectedness hypothesis is separately written into `GaugeGroup`.

## 2. What converges and in which order

`WilsonFamily.Converges` (PublicWilson.lean): for each arity, finite token tuple, and compact smooth test supported in the pairwise-separated configuration locus, suitable large volumes eventually contain its support. There is a row limit at every suitable fixed period p:

lim(N→∞) smear(N,p,v,t) = row(p), then lim(k→∞) row(volumeOrder(k)) = ∫ t(x) unselected(v,x) dx.

The mesh is denominator(N)⁻¹, the physical period is period(p), and the discrete side length is their integer product. This is not a diagonal limit N=N(p), not arbitrary joint convergence, and not yet convergence on tests meeting collision diagonals. Row is existential after the test tuple, so no row-level probability law is included in this statement.

`WilsonRegularity.running` gives a positive slope before ε, then an N threshold uniform over every p: β(N,p)/log(mesh(N)⁻¹)→slope. `HierarchyBridge.regular` identifies its constructed value as 2·Pert.Asymptotics.b0(g.core). Each token has one finite stencil in lattice units before N,p,x. The observable is continuous and gauge invariant. Its shrinking physical support is valuable for local composite-field interpretations, but its N,p-dependent normalization is not a boundedness claim.

`FullVolume.PositionsUV` (YM/FullVolume/Profiles.lean), retained by constructed with Sequence.full, is stronger at fixed period: a smooth function on the periodic no-collision locus, limits of suitable smears, and convergence of every fixed finite difference/derivative list through the interpolation statements Cluster and Meets. Thresholds may depend on period, arity, fields, derivative list, compact set, and tests. It is not a uniform estimate over an unbounded field or test family.

`Extraction.Order.unique` proves uniqueness of two limits for the SAME order in a Hausdorff target. `retained` proves stability under further refinement. Neither compares disjoint subsequences. Two increasing sequences need not possess an infinite common refinement, so a common-refinement argument cannot silently prove full-volume uniqueness.

## 3. Phase mixture versus reconstructed theory

Public `Hierarchy.correlation` is indexed by Option Character. `none` is the unselected lattice-limit hierarchy; `some χ` is a character phase. `Reconstructed.mixture` averages the latter with a full-support probability measure to obtain the former at separated configurations. Character space is compact metric and separated correlations depend continuously on the character.

The selected theory's holomorphic continuation roots at `some selected`, not at `none`. Its complex-time boundary, after product Schwartz smearing, equals vacuum expectations of its actual fields. Public W0 provides a forward-cone spectral condition and uniqueness of translation-invariant vectors up to vacuum scalar. W1 supplies common dense invariant domain, cyclic vacuum, tempered field action, and adjoint identity. W2 is full proper orthochronous affine Poincaré covariance, not merely coordinate subgroup covariance. W3 is commutation for spacelike-separated smeared supports. Faithfulness and coverage identify the field labels with the graded local field union; the lattice token projection is onto the graded fields.

Deeper `YM.Nonperturbative.Headline.Reconstruction` reconstructs every character, not just selected, and its average identity holds for joint boundary distributions. `YM.OS.Spread.Resolution.Out` (YM/Phase/Mixture.lean) makes this concrete: the character space bijects with the characters of the commutative closed star algebra generated by whole-word insertions compressed to the Euclidean heat-ground subspace. It is not an arbitrary latent variable or character assigned separately to each insertion. `Selection.Certificate.retains` says every phase has exactly the unselected null subspace at every grade. Consequently all phases share physical field quotients and differential structure.

`YM.EveryPhase.Completed` (YM/EveryPhase/Extension.lean) adds:

- exact zero-energy sector = scalar vacuum in each character;
- for each fixed finite SumOfWords A, an M≥0 uniform over characters with νletter(i,A)(lowE(u)) ≤ M u² for 0<u≤1;
- for the corresponding word vectors, 0≤⟨x,(e^(−sH)−P0)x⟩≤M s⁻² uniformly over characters for s≥1;
- scalar weak commutant for actual unbounded fields, without requiring the bounded commuting operator to preserve the common domain;
- continued stress tensor in every character, and continuous character dependence of boundary moments.

These are real quantitative infrared inputs, but M depends on A and there is no positive spectral exclusion interval. Low-energy tails of order u² and inverse-square heat decay are compatible with gapless spectra. Irreducibility in each phase is not uniqueness of the phase or of the extraction limit.

## 4. OPE and contact scope

Public `OperatorProductExpansion` uses one coefficient system for all truncation degrees and all viewers (mixture plus phases). For fixed insertions, spectators, geometry constants and ε>0, and each D, one small-radius threshold works for all allowed positions and viewers. The remainder is at most r^(D+1−sum of input grades−ε). Spectators remain a fixed positive distance from the cluster center and from each other. Coefficient bounds are componentwise with j and channel preceding the radius threshold. Finite channels permit compact finite-dimensional norm estimates; this does not give a threshold uniform over every grade.

The all-orders statement means every finite truncation has its asymptotic estimate. It does NOT assert convergence as D→∞ at fixed r. This distinction matters for bootstrap or Banach-algebra completion analogies.

`Headline.SeparatedOPE.profiles`/`Independent.Core.joint` put the radius threshold BEFORE arbitrary integrable complex densities, admissible measurable positions, and measurable varying character choices, provided total variation is ≤B. Thus geometries and phase labels may vary with r and be integrated without new extraction. `Headline.NestedOPE` extends this to admissible finite trees, nodewise truncations/scales, buffer conditions, and bounded-total-variation joint laws. This is stronger than pointwise OPE.

`Hierarchy.Contacts` states existence of a continuous distribution on the full configuration space, agreeing with the ordinary correlation on tests whose support avoids diagonals. `MultipleLight` adds collision-forest estimates and subtraction-dependent thresholds. It does not establish uniqueness of contact extensions, raw lattice convergence on diagonals, reflection positivity on overlapping supports, arbitrary-subtraction contact Ward identities, or real-time time-ordered products. The source explicitly says so. Do not infer positivity of an extended full moment functional merely from positivity on the original separated/reflection geometry.

## 5. Six exploitable structural bridges and their exact extra inputs

### A. Extraction uniqueness as a stability/observability problem

The source already removes the ultraviolet subsequence at fixed volume. Thus the remaining uniqueness problem can be isolated to the sequence of continuum finite-volume rows. Analogy: compactness plus a separating family of conserved observables, or attractor uniqueness for a dynamical system. A concrete target is a volume comparison modulus for a countable convergence-determining family of separated test moments. If all cluster limits agree on that family, continuity/density identifies their whole separated hierarchy; relative compactness then promotes subsequential convergence to full convergence. Missing inputs: a common topology/state space for comparisons, relative compactness there, a separating family, and a boundary-condition/volume stability bound or a rigidity theorem making every cluster point agree. Order.unique alone does not provide them. Phase uniqueness is a second, separate target.

### B. Ground compression as conditional expectation / Choquet coordinates

The character algebra is generated by compressed whole words, and retains the same graded null ideals. This resembles a central decomposition in which local algebraic relations survive in every fiber while long-range order lives in the base. An exploitable identity is variance decomposition: whenever phasewise connected correlations decay, the mixed connected correlation can retain Cov_μ(phase one-point values). Thus failure of clustering in the mixture need not mean failure of clustering in a physical phase. Missing for a rigorous use: precise connected-observable identifications and large-separation limits; enough domination to pass the limit under μ. Collapse of the character algebra to scalars, or constancy of a separating set of its character evaluations, would force a single phase. Scalar weak commutants in each fiber do not force this collapse.

Full support and continuity make a useful upgrade available: a continuous nonnegative defect with zero phase average vanishes at EVERY character, not only almost everywhere. This could transfer an averaged positivity/rigidity identity to all phases without selecting a new subsequence. The defect must be genuinely nonnegative and continuous; equality of complex averages alone is insufficient.

### C. Polynomial infrared decay to a gap via a norm-controlled contraction

The u² spectral-tail bound and s⁻² heat estimate form a Tauberian pair. They provide a concrete baseline for a bootstrap, rather than a blank mass-gap problem. One sharp bridge target is, at some physical time t0>0, ||e^(−t0 H)(1−P0)||≤q<1 uniformly in the target family; then spectral calculus gives gap ≥−log(q)/t0. Equivalently, require an exponential quadratic-form bound with the SAME positive exponent on a dense set of centered word vectors (constants may depend on vector); positivity of spectral measures then removes all spectral weight below that exponent and density finishes the argument. Missing input: exponential decay or an operator-norm/coercivity estimate. The existing word-dependent M cannot be optimized over a unit sphere to manufacture q<1. For a uniform continuum gap, scale conversion and uniformity over cutoff/volume/phase must be specified before taking limits.

### D. Bounded-total-variation OPE as a module over phase-valued measures

Because the same OPE coefficients act in every phase and the threshold precedes measurable phase/geometry choices, the OPE behaves like a uniform fiberwise finite-rank approximation stable under signed integration. This is analogous to robust quadrature, operator interpolation, or a module over a commutative parameter algebra. It supports inserting randomized or adaptively chosen geometries and phase weights without losing the estimate. Nested OPE supplies a tree-composition version. Missing for global associativity or convergent factorization: control uniform in truncation depth/degree, summability of channel growth, and compatibility across all competing collision trees. Existing finite-tree asymptotics alone are not a convergent vertex algebra or a global factorization theorem.

### E. Contact extension as renormalization ambiguity, not probability existence

The same field basis for OPE and contacts makes local counterterm comparison tractable. Two extensions agreeing off diagonals differ by a distribution supported on the diagonals; under finite-order and scaling-degree hypotheses, such differences can be organized as transverse delta derivatives with coefficients on the diagonals and compatibility along their intersections. Analogy: extension of singular kernels, relative cohomology, or renormalization torsors. A viable measure bridge must select contact data satisfying a positive full moment functional and suitable growth/determinacy estimates. For a characteristic-functional route, require a normalized continuous positive-definite functional on a specified nuclear real test space and identify its derivatives with these moments. For a lattice-measure route, require a common measurable field/distribution space, tightness, compatible pushforwards, and uniform integrability. Neither the finite lattice probability measures nor the contact distributions alone supply these missing hypotheses. Reflection positivity is distinct from ordinary probabilistic moment positivity.

### F. UV non-freeness is compatible with an IR continuum of masses

`TraceTest.shortDistance` fixes a positive coefficient for r⁸ log²(1/r) times the centered trace two-point function. `TraceTest.excludes` excludes finite sums of specified mass-shell Wick profiles; the detailed TraceExclusions adds specific finite-Wick-polynomial and Maxwell-Fock comparison countertests. The logarithmic UV signature is a potentially useful rigidity observable for comparing extraction limits or normalizations. It is not a universal infrared classifier: a continuous-mass generalized free field can have logarithmic UV behavior, so the finite-profile exclusion does not prove a particle or a gap. A particle bridge requires nonzero atomic/isolated mass-shell spectral weight for some physical observable; stable scattering additionally needs suitable isolated-shell and localization/energy-bound hypotheses. A uniform gap still permits purely continuous spectrum above threshold.

## 6. Wilson loops and area law: what must be added

The word Wilson in the public files identifies the plaquette action and finite-volume measure. The theorem's tokens are local finite-stencil insertions. A macroscopic loop of fixed physical size typically uses a number of links growing like mesh⁻¹ and is not automatically one such fixed-stencil token. No rectangular-loop area law or string-tension bound is present in these statements.

An actual loop bridge requires definitions of loop holonomy/trace, control under continuum loop approximation, perimeter/cusp renormalization where needed, and large-loop estimates with uniform physical constants. A possible transfer-matrix route would identify charged/flux-sector energies and prove E(R)≥σR−c, with overlap and boundary control sufficient to estimate rectangular loops. The vacuum gap alone does not establish this flux-sector coercivity; an area law alone does not identify ordinary particle poles. Strong-coupling bounds at fixed lattice spacing also cannot be carried to this running β→∞ regime without a scale-uniform argument.

## 7. Source links and declarations

All links use the pinned commit. Base URL:
https://github.com/yang-mills-group/OpenAI-YM-20261028/blob/3456e5e75cf8d4cd558955a5c7638fba2d68f464/

- ChallengeSpec.lean: YangMills.WilsonRegularity; Hierarchy.Reconstructed; Hierarchy.OperatorProductExpansion; Hierarchy.TraceTest; Hierarchy.Properties; JointExistence.
- PublicWilson.lean: wilsonLaw; WilsonFamily.finiteMoment; WilsonFamily.smear; WilsonFamily.Converges.
- PublicTheory.lean: GaugeGroup.IsSimple; FieldSystem.DomainAxioms; FieldSystem.LocalCommutativity; WightmanTheory.W0–W3; WightmanTheory.Axioms.
- Construction.lean: YangMills.constructed.
- YM/Extraction/Orders.lean: Extraction.Order.limits, retained, unique, from_late.
- YM/FullVolume/Profiles.lean: FullVolume.PositionsUV.
- YM/Phase/Mixture.lean: OS.Spread.atom, algebra, Resolution.Out, Resolution.receives, Resolution.barycenter; OS.Selection.Certificate.
- YM/EveryPhase/Extension.lean: EveryPhase.Completed; all_characters; recognizes; exactly_fixed.
- YM/EveryPhase/Irreducible.lean: EveryPhase.weaklyCentral; commutant; precise; irreducible.
- YM/Nonperturbative/HeadlineSpec.lean: Headline.Reconstruction, SeparatedOPE, NestedOPE, TraceExclusions, WilsonExistence.
- YM/Independent/Core.lean: Independent.Core.joint; Independent.Core.point.
- ContactSpec.lean: YangMills.Hierarchy.Contacts.
- YM/Contacts/Solstice/Result.lean: Contacts.Solstice.Resolution, Sunshine, MultipleLight.
- HierarchyBridge.lean and ContactBridge.lean: explicit identification of the public objects with these deeper constructions.

The read source README itself identifies positive mass gap, general continuous-mass generalized-free-field exclusion, all-orders MS-bar matching, uniqueness across contact prescriptions, and raw collision-diagonal cutoff convergence as outside the stated scope. These boundaries are used here as directions for new research, not objections to the trusted construction.

## Spectral-definition addendum

YM/EveryPhase/Light.lean defines lowE(u) as time⁻¹((0,u]), not an invariant-mass interval. νletter(L,i,A) is the diagonal joint energy–momentum spectral measure of stands(actors(i).operators,A). Projection.lean defines zeroE={p:time(p)=0}, aboveE={p:time(p)>0}, boltz(s,p)=exp(−s time(p)), and bHeat(J,s)=J* E_OS(s) J. Its bHeat_pair identifies the diagonal pairing with the Laplace integral of the boundary spectral measure. residual_pair removes the full zeroE projection and gives exactly the integral over aboveE. Decay.decay uses positive tail integration and Gamma(3) to turn the u² bound into s⁻² decay. These defining files are now staged and pinned.
