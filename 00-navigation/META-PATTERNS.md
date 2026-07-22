# Meta-Patterns for Mathematical Research

**Status:** CURRENT
**Role/Use:** compact cards from repeated successes and failures; scan relevant triggers and follow links for evidence.

These cards are defaults, not slogans. Each names when it applies, what to do,
why it works, and when it can mislead. New cards require evidence from distinct
threads or a severe failure with a demonstrated repair.

## Search the statement before the method

**Trigger:** inheriting a target, naming an invariant, or proposing a “new” lemma.
**Action:** search exact constants, inequalities, quantifiers, construction shape,
theorem IDs, and canonical synonym families; dereference every cited ID.
**Mechanism:** methods give the same theorem disjoint vocabularies, so method-keyword
searches miss prior solutions and refutations.
**Counterindication:** independent rederivation remains useful when framed as verification and compared proof-by-proof.
**Evidence:** MISTAKE-183, 187, 189, 200, and 158 in [`MISTAKES.md`](../01-canon/MISTAKES.md).

## Correct the object before sharpening the technique

**Trigger:** many increasingly elaborate methods stall at the same residual.
**Action:** ask whether the optimized quantity is the actual object; replace a
mean by a maximum, an intrinsic shadow by a marked/observer object, or a scalar
by its profile.
**Mechanism:** no stronger bound repairs an information-losing object.
**Counterindication:** retain valid local estimates as components after the reframe.
**Evidence:** the `L -> M` and observer-lens reframes in §1 of
[`lrc14-history-synthesis-patterns-and-reframings-opus-S399.md`](../07-reflections/lrc14-history-synthesis-patterns-and-reframings-opus-S399.md),
and the support-profile correction in MISTAKE-209.

## Type every analogy and every implication

**Trigger:** formulas look identical, an analogy suggests an iff, or a scalar equation is read coefficientwise.
**Action:** write type and implication ledgers. Distinguish nodes/exponents,
supports/multiplicities, labels/classes, scalar/polynomial identities,
graph/tournament spectra, and notions of dimension/rank. Rank coincidences need
ambient modules and an explicit map; spectra need normalization and exceptional-prime scope.
**Mechanism:** shared notation hides type-error bridges.
**Counterindication:** analogies remain productive conjecture generators when the missing map is explicit.
**Evidence:** MISTAKE-209, 211, 212, 214–216, 222–225, and 227–229 in [`MISTAKES.md`](../01-canon/MISTAKES.md).

## Audit saturation and basis covariance before naming a lattice bridge

**Trigger:** objects have equal rank, a generator Gram looks canonical, or a coordinate inequality gets a discriminant/class group.
**Action:** name ambient modules and construct the map; compute kernel, cokernel
or saturation index, discriminant, and minimal vectors. For coordinate quadratic claims,
prove `GL_n(Z)` covariance and build the invariant form before arithmetic classes.
**Mechanism:** full-rank frames can have huge finite index, and a max of linear
forms against a chosen norm need not be one quadratic representation problem;
equal rank hides torsion and basis dependence controlling short vectors.
**Counterindication:** rank/Gram comparisons are useful when explicitly labeled
as frames; class-group tools need a basis-invariant form and target-preserving map.
**Evidence:** MISTAKE-227's index-`11!` AP chain frame and MISTAKE-229's repair
of the nonexistent Heegner discriminant on THM-2053's tangent-disk union.

## Preserve the selected side, not only the walls

**Trigger:** a problem is recast as avoiding hyperplanes, hypertori, resonance walls, or forbidden bands.
**Action:** distinguish the bare arrangement, its ordinary complement, a
thickened complement, and a selected inequality-cell complex. Record wall
owner, orientation/sign, selected side, height functional, and deletion unit.
**Mechanism:** intersection posets and complement cohomology forget the feasible side and extremal top cell.
**Counterindication:** ordinary arrangement invariants suffice when the target is invariant under all forgotten labels and sides.
**Evidence:** MISTAKE-224's repair of the Fourier-lattice/toric-complement
conflation and the exact phase-height recovery of
[THM-1002's pair-sum ruler](../01-canon/theorems/THM-1002-pair-sum-denominator-bound-and-the-bounded-gap-case.md).

## Find the hidden second coordinate in a nearly true theorem

**Trigger:** one elegant property appears to explain a theorem, but a minimal
counterexample preserves that property.
**Action:** compare the proof with the witness and identify the extra
coordinate actually excluding it—parity form, selected side, height, owner,
support, or positivity. Restate the theorem as a conjunction.
**Mechanism:** many structural theorems are two-key locks; naming only the
visually dominant key creates false universal transfers.
**Counterindication:** do not add decorative hypotheses; each coordinate must
perform an explicit proof step or eliminate a concrete witness.
**Evidence:** MISTAKE-225 (skew symmetry plus the tournament block modulo two
for odd support) and MISTAKE-224 (walls plus selected-side/height labels for
the LRC phase-height object).

## Restore a collapsed grading before transporting moments

**Trigger:** a diagonal moment, scalar generating function, or total count is
said to encode the same relations as a multigraded target.
**Action:** expose the omitted grading or adjoin an observer that records it;
test affine translates that preserve the diagonal but change the target.
**Mechanism:** diagonalization silently fixes a conservation law. In
MISTAKE-226, `CT[P^m Pbar^m]` sees only augmentation zero, while the mixed
`CT[P^r Pbar^s]` table restores augmentation `r-s` and the full LRC relation
support. This is the moment analogue of retaining multiplicity beside sequence
support.
**Counterindication:** the diagonal is sufficient only when the target is
proved invariant under the discarded grading.
**Evidence:** MISTAKE-226 and the support/multiplicity collision-tax repair in
MISTAKE-209.

## A local factorization transfers only through a stated functor

**Trigger:** a source polynomial, determinant, arrangement, or generating
function factors, and a downstream wall/limit/topological object is then said
to factor “the same way.”
**Action:** draw the source-to-target construction and prove it respects the
product or localization. Separate vanishing order, codimension, nonzero local
units, analytic limits, and derived valuations. If the arrow is missing, keep
the factorization as a source-side theorem and test the proposed transfer.
**Mechanism:** products are not preserved automatically by scalar summation,
confluence, asymptotic limits, constant-term functionals, or passage to a
different cell complex.
**Counterindication:** use product closure immediately when an exact theorem
shows that the target construction is multiplicative.
**Evidence:** MISTAKE-215 (special determinant versus general moment) and
MISTAKE-223 (braid localization versus hyper-Bessel/Euler claims), with the
surviving special-matrix scope in THM-2033.

## Controlled forgetting requires a sidecar

**Trigger:** quotienting, canonicalizing, projecting, folding, taking moments,
or replacing an object by an invariant.
**Action:** state the next operation, the predicate preserved, the coordinate
forgotten, and the sidecar required to make the next operation legal.
**Mechanism:** a quotient may classify the present object while failing under
extension, gluing, chirality, phase transport, or obstruction detection.
**Counterindication:** no sidecar is needed when a theorem proves the next
operation factors through the quotient.
**Evidence:** the interface field and stress tests in
[`perspective-groupoid-controlled-forgetting-codex-s261.md`](../07-reflections/perspective-groupoid-controlled-forgetting-codex-s261.md),
plus the collision-tax correction in MISTAKE-209.

## After a scalar quotient, choose between a coherent lift and orbit incidence

**Trigger:** a pair average, leading-layer equation, or subset statistic is
exact, but the original target asks whether many local pieces coexist
simultaneously.
**Action:** identify the distinguished subset/fiber before averaging. Then make
an explicit choice: either retain the alignment/divisibility data needed to
lift the scalar equality, or prove that a transitive group action makes every
translate have the same scalar value and sum all translates by uniform
incidence. In the second route, look for a full-object invariant with an
incompatible value.
**Mechanism:** scalarization usually destroys coherent compatibility, but a
transitive orbit can turn that loss into a theorem when every point occurs in
the translated subset with the same multiplicity.
**Counterindication:** do not orbit-average merely because a group acts. The
subset identity must be equivariant, its right side invariant, and the full
sum/product independently controlled. Otherwise preserve a sidecar instead.
**Evidence:** THM-2126's exact LRC pair spectrum forgets the common mod-seven
guard fiber and admits a neutral rank-eight wall; THM-2102's scalar first-defect
identity `L=0` does not imply the divisibility needed to lift a common
approximate root; THM-2101 succeeds in the opposite direction because Galois
transitivity converts the small-root residue subset sum `1` into the full-root
Lagrange sum `0`.

## Existence is a maximum or tail question, not automatically a mean question

**Trigger:** proving a witness exists, excluding a covering configuration, or
studying a near-extremal family.
**Action:** identify the actual max/tail event before applying averaging or
moments; test whether a saturated exceptional rung is invisible to the mean.
**Mechanism:** averaging integrates across precisely the rare structured event
that may carry extremality.
**Counterindication:** means remain useful when a separate inequality converts
them to the needed tail with correct equality cases.
**Evidence:** strategy families 1, 2, and 5 and the five-axis triage in §2–3 of
[`lrc14-history-synthesis-patterns-and-reframings-opus-S399.md`](../07-reflections/lrc14-history-synthesis-patterns-and-reframings-opus-S399.md);
MISTAKE-129 and MISTAKE-171.

## Respect symmetries by searching orbit representatives

**Trigger:** a target is invariant under dilation, translation, relabeling,
complement, or reflection.
**Action:** normalize first and search representatives of the relevant group
action; separately test covariance of every auxiliary predicate.
**Mechanism:** a coordinate slice can intersect an extremal orbit once and make
the rest appear absent; a canonical label can make a non-invariant look stable.
**Counterindication:** normalization is unsafe when an auxiliary hypothesis is
not invariant; rederive it after normalization.
**Evidence:** MISTAKE-156, MISTAKE-160, MISTAKE-170, and MISTAKE-208; the
LRC dilation warning in §3 of the S399 history synthesis.

## Attack a proposed bound before extending it

**Trigger:** sampled maxima keep increasing, a boundary case is sharp, or a
finite certificate is proposed as uniform.
**Action:** spend a short hostile pass constructing a family that worsens the
quantity indefinitely; inspect CRT/lcm, dilation, recursion, and degeneration.
**Mechanism:** an easy construction often exposes the missing quantifier faster
than a larger search.
**Counterindication:** after construction attacks fail for a structural reason,
the reason itself becomes proof data.
**Evidence:** MISTAKE-157’s one-line lcm counterfamily, the fixed-certificate
axis in §3 of the S399 history synthesis, and the rung-ladder analysis cited
there.

## Test structured adversaries, not only random samples

**Trigger:** empirical minima, asymptotic domination, or negative searches.
**Action:** include known extremals, dilated or near-dilated APs, planted
resonances, subgroup/coset families, boundary parameters, and constructions
designed to trigger every escape clause.
**Mechanism:** extremals often occupy thin algebraic strata invisible to random
sampling and local descent.
**Counterindication:** random sampling is still useful for discovering generic
behavior when clearly separated from uniform claims.
**Evidence:** MISTAKE-101/102, MISTAKE-160, MISTAKE-162, MISTAKE-175, and
MISTAKE-194.

## Verify the consequence, not the model’s own assumptions

**Trigger:** a computation supports a derived bound, asymptotic ladder, or
analytic consequence.
**Action:** feed the program data from the original geometry/definition and
print the object asserted in the conclusion. Label ansatz-fed tests as
arithmetic checks only.
**Mechanism:** a model can verify its own algebra perfectly while its premise or
translation to the theorem is false.
**Counterindication:** synthetic tests are appropriate for unit-testing an
implementation, not validating the mathematical model.
**Evidence:** MISTAKE-204 and MISTAKE-207 in
[`MISTAKES.md`](../01-canon/MISTAKES.md).

## Use redundant paths as detectors

**Trigger:** a computation is exhaustive, uses pruning/canonicalization, or
supports a load-bearing theorem.
**Action:** rerun with a different pruning configuration, labeling, code path,
optimization mode, or independent formula; include positive controls.
**Mechanism:** small disagreements reveal unsound pruning, hidden labels, and
non-invariants even when the final headline happens to survive.
**Counterindication:** redundancy should target the load-bearing logic rather
than duplicate every inexpensive arithmetic step.
**Evidence:** MISTAKE-194’s pruning discrepancy, MISTAKE-208’s
canonical-vs-random-label detector, and MISTAKE-191’s source/output mismatch.

## Exactify before interpreting a decimal

**Trigger:** a polytope, affine-cell measure, recurrence, or repeated numerical
constant appears.
**Action:** add the defining constraints, change coordinates, factor, telescope,
or use exact rational/symbolic arithmetic before naming the phenomenon.
**Mechanism:** exactification reveals shape and mechanism; sampling supplies
only a fragile value.
**Counterindication:** numerics remain valuable scouts when accompanied by an
error model and no arithmetic claim.
**Evidence:** MISTAKE-176’s six-simplex correction and the “exactification
bounty” in §8 of
[`lrc14-history-synthesis-patterns-and-reframings-opus-S399.md`](../07-reflections/lrc14-history-synthesis-patterns-and-reframings-opus-S399.md).

## A failed certificate is not a failed theorem

**Trigger:** a sufficient dominance order, invariant, sieve, or projection stops
working.
**Action:** state exactly where the certificate ends; search for phase,
confluence, recurrence, or another sidecar before inferring failure of the
target.
**Mechanism:** certificate failure identifies information loss, not the target’s
truth value.
**Counterindication:** an explicit counterexample to the target itself remains
decisive.
**Evidence:** MISTAKE-212’s transitivity/source certificate, the proportional-
NC2 MISTAKE-213 entry separating “not eventually zero” from “finitely many
zeros,” and the “frames describe, do not prove” phase in the S399 history
synthesis.

## Formal verification needs a semantic witness

**Trigger:** a Lean file builds, a finite decision returns true, or a theorem is
called kernel-pure.
**Action:** exhibit satisfiable hypotheses, test boundary branches, confirm the
module is imported by the intended root, and audit axioms/build artifacts.
**Mechanism:** the kernel validates the encoded implication, including vacuous
ones; it does not certify that the formal statement matches the paper theorem.
**Counterindication:** once semantic correspondence and non-vacuity are proved,
kernel checking is the repository’s strongest durable evidence.
**Evidence:** MISTAKE-186, MISTAKE-138, and MISTAKE-094 in
[`MISTAKES.md`](../01-canon/MISTAKES.md).

## Treat pulls and collisions as synthesis opportunities

**Trigger:** concurrent work lands on the same prompt, invariant, ID, or proof
route.
**Action:** pull before filing, compare distinct contributions, consolidate with
credit, and keep independent verification only when labeled as such.
**Mechanism:** early synthesis converts duplicated effort into cross-checks and
often exposes a stronger common object.
**Counterindication:** genuinely orthogonal routes should remain separate until
an explicit reduction connects them.
**Evidence:** MISTAKE-199 and the rebase-as-signal protocol in
[`CONCURRENT-SESSIONS.md`](CONCURRENT-SESSIONS.md).

## Turn certificate failure into an address, then change sidecars

**Trigger:** a sufficient gate leaves a finite or structured bad set that is
still too large, or a control family fails the gate while the theorem remains
true.
**Action:** preserve the failed gate as an address (deck, owner cone, ray, or
cell), then attach an orthogonal consequence-bearing sidecar such as a clock
orbit, divisibility tax, phase interval, endpoint owner, or rank increment.
**Mechanism:** failure of a sufficient certificate localizes uncertainty but
does not describe danger. THM-2057 closes two one-tail planes containing many
THM-2053 determinant failures because scaled clocks and binding rays see the
missing phase coordinate.
**Counterindication:** do not multiply sidecars indiscriminately; each must
preserve the target predicate and eliminate a named failure family or state.
**Evidence:** the THM-2053 transverse deck, THM-2055/2056 normal-fan and
Kelvin/Farey address, THM-2057 missing-clock closure, and MISTAKE-224's repair
of the side-blind toric-complement quotient.

## Join compatibility fibers before comparing marginal sizes

**Trigger:** two nonempty packets, fibers, or local certificate sets must meet,
but their separate cardinalities do not decide whether a common witness exists.
**Action:** identify the common quotient, retain histograms over its fibers, and
compute the compatibility pairing before collapsing either side to one scalar.
**Mechanism:** generalized CRT turns a core packet and tail packet into the
exact dot product of their reduction histograms; disjoint supports explain a
failed join even when both marginals are large.
**Counterindication:** use a simpler scalar only after proving uniform fibers,
independence, or another theorem making the pairing a function of the totals.
**Evidence:** THM-2059's CRT packet theorem repairs the marginal-size loss; the
same error genus is exposed by MISTAKE-231's observable-relative fiber counts
and by observer/cut-payload quotients in the tournament atlas.

## “One item left” requires a typed residual

**Trigger:** a synthesis declares a problem one lemma, one inequality, or one
finite check from completion.
**Action:** type the residual independently of the current method: its domain,
quantifiers, symmetry, uniformity, and known equivalent names. Attack that
typed statement before announcing proximity.
**Mechanism:** phase-boundary summaries compress the residual into the current
vocabulary and hide escape clauses or renamed walls.
**Counterindication:** a genuinely finite residual with a verified exhaustive
certificate may be called finite, with its exact checker named.
**Evidence:** the fourteen near-completion episodes and phase-boundary analysis
in §1 of
[`lrc14-history-synthesis-patterns-and-reframings-opus-S399.md`](../07-reflections/lrc14-history-synthesis-patterns-and-reframings-opus-S399.md).

## Recover niches by more than citation counts

**Trigger:** choosing an underexplored thread or auditing supposedly abandoned
work.
**Action:** distinguish identifier dormancy, frame dormancy, under-titled work,
orphan computations, refuted-but-close claims, and repeatedly cited but stalled
problems. Rank candidates by orthogonality, artifact readiness, transfer value,
and cost of a decisive probe.
**Mechanism:** citation counts can label live code “abandoned” and heavily cited
but stagnant problems “healthy.”
**Counterindication:** citation graphs remain excellent discovery tools when
paired with semantic inspection.
**Evidence:** the blind-spot discussion in
[`ATLAS-OF-ATLASES-2026-07-20-macmini-S124.md`](ATLAS-OF-ATLASES-2026-07-20-macmini-S124.md),
the frame census in
[`the-recovery-ledger-forgotten-frames-revivable-hypotheses-and-the-gaps-deathstar-S79.md`](../07-reflections/the-recovery-ledger-forgotten-frames-revivable-hypotheses-and-the-gaps-deathstar-S79.md),
and the citation-graph method in
[`the-corpus-atlas-forgotten-threads-the-zoo-and-the-gaps-kps-S128c136.md`](../07-reflections/the-corpus-atlas-forgotten-threads-the-zoo-and-the-gaps-kps-S128c136.md).

## Fill operation columns, not only invariant columns

**Trigger:** a field has many objects and invariants but few reusable laws.
**Action:** ask how each invariant behaves under complement, deletion, ordinal
sum, duality, lift, product, and degeneration; transport the same operation to
the dual object.
**Mechanism:** one operation law fills an entire row of the research grid and
connects results that scalar invariant hunting leaves siloed.
**Counterindication:** compute a new invariant first when no existing invariant
sees the target predicate at all.
**Evidence:** the meta-gap in
[`the-procedural-generation-grammar-for-the-tournament-zoo-deathstar-S79.md`](../07-reflections/the-procedural-generation-grammar-for-the-tournament-zoo-deathstar-S79.md)
and the gap cross-product in the corpus atlas.

## Tournament Analysis must preserve content

**Trigger:** pairwise data suggests a tournament representation.
**Action:** require an intrinsic binary relation, name the preserved target
predicate, keep exact ties as ties/confluent data, and audit alternate vertex
sets.
**Mechanism:** cosmetic tie-breaking always manufactures a transitive ranking
and can turn node equality, scalar order, or proof regimes into a contentless
“tournament.”
**Counterindication:** a meaningful gauge may resolve ties when the gauge itself
is part of the theorem.
**Evidence:** MISTAKE-212, MISTAKE-214, and the tournament audit in
[`the-underlying-object-is-the-support-dirichlet-profile-codex-20260721.md`](../07-reflections/the-underlying-object-is-the-support-dirichlet-profile-codex-20260721.md).

## Type a shared sum as a weighted fiber before transferring it

**Trigger:** two problems are both written as “a sum over a kernel.”
**Action:** record `Z(M,A,b,w)=sum_{x in M, Ax=b} w(x)`, including the monoid,
grading/fiber, regularization, weight ring, and target predicate on each side.
**Mechanism:** LRC uses an infinite Fejer-regularized fiber in `Z^n` with sinc
weights; a fixed GMC moment uses a finite fiber in `N^s` of `(1,q)` with
multinomial, factorial, radial, and coefficient weights. The schema is shared,
but the typed objects and quantifiers are not.
**Counterindication:** transfer only after a map intertwines the fibers and
weights and preserves the desired nonvanishing or witness predicate.
**Evidence:** MISTAKE-226/234/235 and THM-2059's successful typed CRT fiber product.
