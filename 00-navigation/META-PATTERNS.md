# Meta-Patterns for Mathematical Research

**Status:** CURRENT
**Role:** compact method cards distilled from repeated successes and failures
**Use:** scan the triggers relevant to the current session; follow links for the full evidence.

These cards are defaults, not slogans. Each names when it applies, what to do,
why it works, and when it can mislead. New cards require evidence from distinct
threads or a severe failure with a demonstrated repair.

## Search the statement before the method

**Trigger:** inheriting a target, naming an invariant, or proposing a “new”
lemma.
**Action:** search exact constants, inequalities, quantifiers, construction
shape, theorem IDs, and canonical synonym families. Dereference every cited
ID.
**Mechanism:** different methods give the same theorem disjoint vocabularies;
method-keyword searches therefore miss prior solutions and refutations.
**Counterindication:** independent rederivation remains valuable when explicitly
framed as verification and compared proof-by-proof.
**Evidence:** MISTAKE-183, MISTAKE-187, MISTAKE-189, MISTAKE-200, and
MISTAKE-158 in [`MISTAKES.md`](../01-canon/MISTAKES.md).

## Correct the object before sharpening the technique

**Trigger:** many increasingly elaborate methods stall at the same residual.
**Action:** ask whether the optimized quantity is the theorem's actual object;
replace a mean by a maximum, an intrinsic shadow by a marked/observer object,
or a scalar by its profile when necessary.
**Mechanism:** an information-losing object cannot be repaired by a stronger
bound on that object.
**Counterindication:** retain successful local estimates as components after
the reframe; an object correction does not erase valid lemmas.
**Evidence:** the `L -> M` and observer-lens reframes in §1 of
[`lrc14-history-synthesis-patterns-and-reframings-opus-S399.md`](../07-reflections/lrc14-history-synthesis-patterns-and-reframings-opus-S399.md),
and the support-profile correction in MISTAKE-209.

## Type every analogy and every implication

**Trigger:** two formulas look identical, an analogy suggests an iff, or a
scalar equation is read coefficientwise.
**Action:** write a type ledger and implication ledger before transporting the
claim. Distinguish nodes/exponents, supports/multiplicities, labels/classes,
scalar/polynomial identities, and different notions of dimension or rank.
**Mechanism:** many compelling false bridges are type errors hidden by shared
notation.
**Counterindication:** analogies remain productive as conjecture generators if
their missing map is stated openly.
**Evidence:** MISTAKE-209, MISTAKE-211, MISTAKE-212, MISTAKE-214,
MISTAKE-215, and MISTAKE-216 in [`MISTAKES.md`](../01-canon/MISTAKES.md).

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
