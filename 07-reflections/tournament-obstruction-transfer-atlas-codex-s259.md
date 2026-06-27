# Tournament Obstruction-Transfer Atlas

*codex-2026-06-27-S259b.  Companion to HYP-3105, LTI-243, and LTT-141.*

The user's prompt points at an important evolution in the tournament program.
The H=7 and H=21 arguments are not just small-value curiosities.  They are
examples of a proof-by-transfer pattern:

```text
creative subproblem -> conflict carrier -> forbidden spectrum or forced
expansion -> contradiction
```

The hard part is the first arrow.  A tournament analogy has no proof value
until the transfer functor says which predicate survives and which coordinate
has been destroyed.  This is exactly the lesson of the incoming Lean
observer-gluing frontier: `coarseWinding_degenerate_not_terminal` formally
demotes coarse mod-14 winding H to a shadow, while
`ObserverGluingCoverage` names the real target.

## What H=7 and H=21 Teach

THM-029 says the H=7-looking component cannot remain isolated: the true
tournament/OCF structure forces a common vertex, a 5-cycle, or a larger
independent-set contribution.  THM-079 and THM-115 do the same for H=21:
component factorization, P4 components, K6-minus-two-edges, K8-minus-one-edge,
and K10 skeletons all look numerically possible only until the full carrier
forces extra payload.

That pattern is more valuable than the numbers.  The reusable move is:

```text
assume the clean forbidden skeleton;
prove legal closure adds a sidecar, cycle, overlap term, or component factor;
watch the invariant jump out of the forbidden spectrum.
```

This is why "the subproblem has H=7" is a weak statement, while "the transfer
to a connected OCF carrier preserves the LRC predicate and forces the H=7
skeleton" can become proof currency.

## New Tournament Uses

The atlas generated several ways to use tournaments beyond direct
contradiction:

- Forbidden-spectrum transfer: import THM-029/079/115 only after a faithful
  carrier map.
- Forced-expansion closure: show an attractive skeleton cannot be closed under
  legal packet operations.
- Edge-flip stress disproof: perturb sidecars while preserving the predicate;
  if a proposed scalar ranking flips, it is telemetry.
- Median-center legality: test whether route triples acquire unique legal
  centers after sidecars are attached.
- Cycle-class observability: use cycles as basis atoms for missing proof
  columns, not as decoration.
- Typed OCF evaluation: replace a scalar H by a vector of component,
  inclusion-exclusion, overlap, and exit types.
- Valuation-residue audit: use p-adic or residue analogies only when the lift
  fiber and metric/discrepancy payload are explicit.
- Redei/Landau certificates: encode a count or score profile only if the
  transfer is faithful, then use parity or score inequalities as hard
  obstructions.
- Cycle-census holes: treat missing `(c3,c5,...)` fibers as their own
  spectral obstruction layer, separate from H.
- Improvement nontransitivity: use local minima in an exchange tournament as
  a finite-check scheduler when greedy exchange fails.

The Tournament Analysis vertices are technique applications and proof
obligations, not runners.  The pairwise observable is retained payload versus
transfer burden and destroyed-coordinate safety.

## Incoming S31ah/S65

The concurrent tournament certificate engine makes the atlas more concrete.
It self-tests Paley `T_7` with `H=189`, confirms H parity, and treats H=7,
H=21, and even H as certificate failures.  More importantly, it turns
"find another H=7" into a generator: enumerate achieved spectra for multiple
tournament invariants, then promote persistent holes to candidate
impossibility certificates.

The follow-up KPS spectrum pass validates that generator on the original
frontier.  It rediscovers H=7 and H=21 as persistent H gaps and identifies the
single-component form: clique `Omega=K_m` gives `H=1+2m`, with `K_3` and
`K_10` as the missing clique sizes behind H=7 and H=21.  This is the cleanest
current obstruction-transfer target: prove the construction lands in a
connected clique-like `Omega` component, not just a scalar H shadow.

The S65 cap exchange scout is the cleanest example of a tournament technique
failing productively.  It verifies that cap minimizers are bounded over
`{1..13}`, `{1..16}`, and `{1..20}`.  Greedy descent reaches the known
minimizer for `j=2,3,4`, but gets stuck for `j=5`; the improvement tournament
has spurious local minima.  So the exchange proof is not transitive enough,
but the failure is structured: optimality reduces to a finite local-minimum
audit.

The apex-7 bridge audit blocks an attractive false transfer.  The seven at
the apex is the number of antipodal diameter pairs in a 14-clock perfect
matching, not the value `I(K3,2)=7`.  Winding tournaments avoiding H=7 is
vacuous because H=7 is forbidden for every tournament.  Any future apex-to-H
argument must build a different functor.

The baby-Hodge crosscheck adds a second obstruction layer.  The `(c3,c5) =
(8,10)` hole is a `c5`/power-sum spectral hole inside the regular score class,
not an H=7/H=21 alpha-vector event.  This justifies using cycle-census holes
as distinct certificates rather than forcing every gap into H language.

## Applications

For LRC14 observer gluing, the best tournament has vertices like arithmetic
lift, normalized arc, Pascal scissors, moment/Perron, branch/K33, and formal
witness chart.  A contradiction would say: no chart supplies a terminal exit,
yet the overlap ledger realizes a forbidden obstruction-transfer skeleton.
That is a sharper target than finding one universal scalar.

For THM-577, the tournament object is the inclusion-exclusion order vector.
The cap value is the shadow; the proof question is whether the `j=4,5` dip
must be produced by finite forced overlap terms or discharged by a
moment/Perron certificate.  The S65 exchange tournament adds a practical
route: prove boundedness first, then eliminate non-global local minima.

For HYP-3094, K33 state lift should be tested by median-center legality rather
than positive mass.  Covering rows and K33 rows are both positive-open; what
matters is which sidecars create legal terminal exits.

For HYP-2963, edge-flip stress is a practical disproof machine.  Any new
scalar quotient should be shaken against endpoint owners, exact M/Farey
height, topology, primitive deck, active binders, and route/state-lift fields
before it is promoted.

For q-cusp, sixth-power, and speculative p-adic/Roth/Steiner routes, the
tournament is mostly a rejection filter.  Keep finite principal parts,
support-six lane rank, valuation fibers, discrepancy bins, and metric payloads
or name the analogy as scalar shadow only.

## Next Work

Build an obstruction-transfer ledger over HYP-2963 and the S258/S259 observer
gluing rows.  Each row should declare:

```text
surrogate vertex set
transfer functor
preserved predicate
target H or typed OCF vector
minimal forbidden skeleton
forced expansion payload
destroyed coordinate
required sidecar
edge-flip stress result
certificate invariant family
single-component H gap
clique Omega realizability
Omega sparsity
cycle-count fiber
improvement-tournament local minima
apex-tie matching status
terminal exit or named debt
```

The creative reframe is that tournaments are not one object in the proof.
They are a proof-language for controlled forgetting.  They can be contradiction
engines, stress tests, route schedulers, or analogy rejectors depending on
what the transfer functor preserves.
