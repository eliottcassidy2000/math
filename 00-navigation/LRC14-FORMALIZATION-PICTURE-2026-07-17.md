# THE LRC(14) FORMALIZATION PICTURE — the cohesive map (death-star-2026-07-17-S43)

One page. Every node's status; every arc's target. Update in place as pieces land.
(Statuses: **LEAN** = kernel-pure in TournamentH7; **paper** = proved on paper /
exact-computationally, Lean pending; **open** = genuinely open mathematics.)

## The trunk

```
LRC14Statement
  ⟸  cite : LRCUpTo13                       [named citation node — owner policy, never a sorry]
   +  ResidualObligation                     [ONE open branch; everything else LEAN]
```

`lrc14_grand_assembly` (monad, THM-671) discharges, in LEAN, every family that is:
non-covering (sieve) · bounded ≤ 22 (kernel census) · all-comparable ≤ 13× (spread) ·
91-dominant (cite-peel) · repeated |speed| (cite) · detuned harmonic (THM-668) ·
coarse ≤ 12 (cite) · common residue (THM-682a).

## The residual, and who bites it

`ResidualObligation` = covering ∧ scale-gapped (>13) ∧ compressed ∧ distinct ∧
max ≥ 23 ∧ no-detuned ∧ no-coarse ∧ no-common-residue. Its machine-checked
SHRINKS (each strictly sharper, all LEAN):

| Surface | Module | What it removes |
|---|---|---|
| DenseCoreObligation (THM-937) | LRCChainDichotomy | any citable ratio-3 singles tail |
| TripleCoreObligation (THM-938) | LRCPairLiftDichotomy | the S20 pair dodge at the dense position |
| QuadCoreObligation/Closed (THM-941/942B) | LRCBlockSplitLift | the dense triple as a fat block; ALL corners closed |
| dissociated chain-dense B5 (codex S25–27) | LRCDenseCoreEndgame | detuned 2/3-nonmultiples; B5 needed only on the trapped core |
| + trapped supply (THM-939) | LRCDenseCoreRelationTrap | the supplier may assume relation-freeness above the dense pair |
| + lacunary branch (opus S336) | LRCLacunaryNest | ratio ≥ 7/3 throughout ⟹ lonely UNIVERSALLY (first quantifier-level branch) |
| + full lacunary wire (codex S56) | — | the lacunary conjunct closed and wired |

## The B5 funnel (the discrete instrument, all LEAN)

```
B5 = Σ(−1)^d S_d                                    [THM-671 discrete Bonferroni]
  = Σ_{|T|≤5} (−1)^{|T|} N_T                        [THM-940 subset expansion]
  = (q−1)·2052/16807 + signed deviation ledger      [THM-940 exact identification]
       singles D_i ∈ [−13/7, 0], −13/7 exact at 14|q     [THM-942A]
       dilate pairs D_ij = (5/7)Q exact; ≤12, all below the pair  [THM-943A + THM-944]
       moment certificates: gap 328× → 1.077× capped     [THM-945]
       MOMENT WALL: moments alone can never finish        [THM-945, legal witness]
  = liveCount − #{bandCount = 6}   on cap-6 strata   [THM-945 exact identity]
  = liveCount − Σ_p C(bandCount(p)−1,5)             [exact weighted identity]
  ≥ liveCount − 792·#{bandCount ≥ 6}  UNCONDITIONALLY [THM-950 census criterion]
```

The pointwise costs at depths `6,...,13` are exactly
`1,6,21,56,126,252,462,792`; thus `792` is sharp at depth thirteen but
coarse for every shallower event.  `LRCWeightedDeepCensus` proves that the
exact weighted race is equivalent to `B5 > 0`, while THM-950 is its uniform
corollary.  Neither identity supplies a modulus or a live-count floor.

**The census pipeline (S43)**: `CoverageCapped` decidable + `lonely_of_census`,
and now the uncapped exact weighted consumer:
cap + live > depth-six, or weighted cost < live, by `decide` ⟹ B5 > 0 ⟹
Mreach ≥ 1/14 ⟹ Lonely.  These are end-to-end certificates for a supplied
tuple and modulus, not a universal `DenseCoreCensusB5Supply`.

## The killer/witness arc (all LEAN)

fragmentation ladder (THM-883, S30–32) → the DICTIONARY (THM-947: bands = bad
arcs at rationals; 7-overlap = 21 integer constraints) → q-window + witness-ladder
rigidity (THM-949: witnesses ∈ [1,v]; n_j ≥ 3n_i on [3,13] — the residual's own
ratio window; 7-chains force n_top ≥ 729) → pinned-p counting (THM-950: one bad p
per witness at q ≤ 7v; deep sets = finite checks) → codex S48: overlapDet colors +
Plücker relations on bad triples → activity-weighted colored event fibers and
gcd-spaced fixed-color fibers (the trap's target shape).

## The rooted-seven activity bridge (LEAN ledger; supplier still open)

At depth `c`, rooting at one bad runner produces `C(c−1,6)` seven-stalks and
the exact B5 debt is `C(c−1,5)`.  Hence

```text
C(c−1,5) / C(c−1,6) = 6/(c−6),
C(c−1,5) ≤ 3 C(c−1,6) for 8 ≤ c ≤ 13,
weightedDeepCost ≤ N6 + 3 N7 + 3 rootedSevenActivity.
```

Equality in the high-depth inequality occurs only at `c=8`.  The exact
unpaid combinatorial residue is depth six (one unit, no seven-stalk), depth
seven (six units versus one possible three-unit colored charge, or one unit
after the five-unit nonzero-triangle charge), and rank-one/aligned stalks.
`LRCSevenOverlapActivity` proves the eventwise aligned/colored partition and
the three-/five-unit *lower* spoke-mass charges.  This is not yet a payment
theorem.  `LRCSevenStalkReuseBudget` now gives the exact global transport
double count: a spoke is reused `choose(m-1,5)` times and a lower pair
`choose(m-2,4)` times, hence at most `462` and `210` times.  These are exact
reuse factors, not a non-reuse estimate.  `LRCAlignedStalkAggregation` adds the
missing fixed-root Fubini identity: multiplier-first rooted activity equals the
sum of the fixed rooted-face fibers.  It partitions that activity exactly into
an all-zero root star and a colored complement, identifies the colored
root-spoke mass with the imported reuse transport, and retains the actual
overlap multiplicity on both sides.  The open step is to convert this ledger
into a payment inequality while preserving the separate depth-six and
depth-seven residues.

`LRCOverlapReflection` identifies the exact mirror action on this carrier.
For every bad runner, `failWitness(q-p)=v-failWitness(p)`; consequently every
determinant color and the entire sparse Plücker coefficient vector change
sign.  Genuinely colored triple and nonzero-lower-triangle activity therefore
has even cardinality for **every** positive modulus, not only odd moduli; an
even modulus's half point is necessarily zero-colored.  Mirrored events carry
one projective relation rather than two independent relations.  This rules out
singleton colored exceptions; combined with the coprime fixed-color bound,
every nonzero absolute color occurs in exactly `0` or `2` events.  This does
not by itself improve the live/debt ratio, since live and colored activity
mirror together.

`LRCDeepReflectionParity` applies the same orbit principle to every exact-depth
stratum at every modulus: its cardinality modulo two is exactly the membership
of the possible reflection midpoint.  At `q=2m,p=m`, a runner is bad exactly
when its speed is even, so the unique odd-cardinality depth is the number of
even speeds.  For odd `q`, every stratum is even, in particular depths six and
seven.  This removes unexplained singleton tails from a future payment
argument, but paired events and a possible classified midpoint still carry
their full debt.

`LRCOverlapColorFibers` supplies the first multiplicity control.  For positive
pair speeds `a,b`, every determinant is divisible by `gcd(a,b)`, and on
`q ≤ 7a` each fixed-color multiplier fiber has cardinality at most
`gcd(a,b)` (at most one for a coprime pair).  `LRCAlignedResonance` proves the
arithmetic end of the aligned target: `q ∣ hp` iff
`q/gcd(h,q) ∣ p`, hence exactly `gcd(h,q)−1` nonzero multipliers resonate;
it also proves the strict integer-closeness atom that forces `hp=rq` in the
top window.  `LRCZeroColorGluing` now upgrades every connected zero-color stalk
to a clique with one primitive witness parameter whose denominator divides all
stalk speeds; exact resonance carries the reduced multiplier modulus into that
parameter, with the exact normalization `d=q/gcd(p,q), n=p/gcd(p,q)`.  The
sharp alternative is `14|v/d|<q`, and the window is not implied by `q≤7|v|`
(already `q=99,d=98` supplies a zero-slope nonresonant stalk).  What remains is
the adaptive-modulus instantiation and the comparison of colored lower charges
with the `462/210` transport budgets.  `LRCZeroColorStalkFork` closes the sharp
window uniformly through `q=98`: seven distinct positive reduced magnitudes
cannot all lie in `{1,...,6}`, so every connected all-bad zero-color seven-stalk
in that range is resonant; at coprime multipliers, `q` divides all seven speeds.
The `q=99` example makes this exact cardinality cutoff sharp.  On the concrete
rooted-face carrier, `LRCAlignedStalkGluing` proves that an anchor-zero star is already
a clique, the actual `Finset.gcd` supplies the shared integer parameter,
top-window badness forces `h*p=r*q`, and a fixed stalk has at most
`gcd(h,q)-1` active multipliers.  Aggregating gives

```text
zeroColorRootedFaceActivity
  <= sum_face (gcd(gcdSpeeds(face),q)-1)
```

under the per-face top window.  This is the first honest aligned aggregation
theorem, but not a global supplier: a selected-root/top-window coverage lemma,
a bound for the summed gcd budgets, the colored payment inequality, the
depth-six/seven payment, and an adaptive live floor are still missing.

## Tournament/carrier audit

The static runner tournament is not faithful here.  At fixed `(q,p)`, after
positive-speed normalization the sign of `overlapDet(i,j)` merely orders the
witness slopes `n_i/v_i`, so it is transitive modulo zero-color ties: there
are no directed cycles, SCCs are singleton after tie-breaking, and the
Hamiltonian path is unique up to orders within tie blocks.  Forgetting
magnitude destroys the Plücker law, gcd spacing, and multiplier activity.  The
preferred carrier is the fibered movie of wall events/rooted seven-stalks
`(p,i,n_i,D_ij)`.

| Carrier | Preserves | Loses if quotienting too early |
|---|---|---|
| multiplier–runner incidence hypergraph | live/depth and event multiplicity | determinant colors and repeated color fibers |
| fibered colored wall event / rooted seven-stalk | phase, witnesses, colors, Plücker, depth, activity | faithful current carrier |
| relation circuit / oriented matroid | exact coefficients and elimination | multiplier multiplicity unless fibered |
| pair-tower wall/prefix carrier | wall order and prefix refinement | phase timing in a runner graph |
| selected-phase safe-cell/branch incidence | existential common branch | intersection geometry in a scalar tournament |
| Fano/`χ₇` address atlas | compact labels for the 21 edges of `K₇` | metric weights and Heawood-cycle data; needs the full colored graph or a sidecar |
| translated Kakeya needle / arithmetic stalk | direction, intercept, and repeated wall hits | a direction-only quotient loses phase offset and multiplicity |

The challenged assumption is therefore that tournament vertices must be
runners.  Multiplier events, wall crossings, rooted stalks, relation circuits,
and proof obligations are the useful alternate vertex sets; each quotient
must state which LRC predicate it preserves.

## The parallel analytic arcs (paper → Lean in progress)

- **kps exhaustion line** (THM-946/952/953/959): the corrected two-pole atom,
  punctured congruence estimate, leverage identity, and general odd folded sum
  are rigorous, and `T3` is closed.  `T4`/`T5` have a coherent composed-atom
  blueprint and finite referees, but the floor-case/Fubini transcription and
  the structured small-relation alternative are not yet a universal supplier.
- **opus/codex cascade and blocks** (THM-932/933/936, THM-955, THM-959):
  `window_tail_glue` is Lean; the corrected positive-window cluster bound and
  direct prescribed block-tower induction are paper-proved with exact rational
  constants.  The abstract sorted-gap L1 is production-green in
  `LRCClusterGapBrick`; the live formal debt is periodic-comb enumeration and
  `arcSafe` assembly, plus a partition trichotomy—not arbitrary-block
  composition by the single-tail API.
- **klein**: manifest items 1, 2, ⅔ of 6 LEAN; interval transport; band data.
- **mac-mini / codex-S52--S66**: UnitBudget and the complete strict-open
  pair-grid ledger are LEAN, including the zero atom, circle volume, sharp
  component budget, and clean-534 discrepancy.  The top `24/12` path caps and
  anchored ratio quotient are LEAN.  Exact comb reindexing closes the covariance
  producer, and tau3/tau4 ratio replays are LEAN; tau5--tau9 remain.  Separately,
  the sparse rational-pair branch interval now transports bijectively back to
  multiplier space, giving the complete finite-`q` three-branch count; the
  zero-width integral boundary requires `Int.toNat`, not a raw floor difference.
- **death-star S52 two-circle probe, continued in S66**: for the canonical
  family `{1,...,13}`, `LRCTwoCircle` proves that the integer and half-resonance
  circles lie inside the depth-six set.  `LRCTwoCircleConverse` now proves the
  converse branches with smallest failing speed `k₀=1` and `k₀=2`, and proves
  that every deep phase has a failing speed at most eight, eliminating
  `k₀>=9`.  Exact recon still predicts that the two circles exhaust the deep
  set, but the branches `k₀=3,...,8` remain; there is no full converse and this
  canonical classifier is not a uniform dense-core live floor.
- **codex-S57/S66 weighted/colored census**: the exact depth ledger,
  rooted-seven comparison, fixed-color gcd fibers, exact `462/210` reuse count,
  zero-color primitive gluing, fixed-root incidence Fubini, summed aligned gcd
  budget, and all-q reflection-midpoint parity are LEAN.  The remaining step is
  a quantitative root/window/payment theorem and adaptive live floor, not
  another pointwise Plücker identity or incidence count.
- **boxeph/klein L-frame** (LEM-032..037, LEM-045): the factorization law's
  both factors, the per-runner danger measure `<=1/7`, and the seven
  consecutive overlap credits are kernel-pure.  `LRCC8Consecutive` assembles
  them along the path on every block `v,...,v+7` and proves positive restricted
  Lebesgue measure for the block's common safe complement.  This is the live
  generic shifted `c=8` theorem; it does not add the other five danger combs,
  so it is neither an LRC(14) theorem nor the weighted dense-core supplier
  above.

## What is genuinely open (the honest short list)

1. **The residual's dense core** past all shrinks: packets admitting no corrected
   THM-959 block partition, including a moderate-junction branch not yet reduced
   to one seven-comparable block. Attacks: an explicit partition trichotomy,
   kps exhaustion (paper, Lean pending), and the census funnel.
2. **The pair-floor finite replays** — the pair continuum producer, component
   budget, Turan arithmetic, top path caps, and clean-534 transfer are closed.
   Tau3 and tau4 compile; the exact remaining tasks are tau5--tau9.
3. **A-priori liveCount floors** (the discrete side of the same wall): the
   winData rows give them on censused strata; beyond, it is kps's law made
   uniform. SHARPENED (S44): the naive coprime law is REFUTED (Dirichlet-rate
   mirror-paired exceptions); the question is now adaptive-q — for each family,
   SOME q in the window with live > 792·deep (deep sets provably mirror-symmetric,
   even at odd q — THM-953).
4. **Selected-witness relation residues** — q22 remains an opposition-avoidance
   selector; q244 and q333 remain explicit zero/small support-three frequencies.
   The q4488 common-sign and matching-wall producers are closed, and its large
   branch escapes.  The two q4-derived frequencies satisfy the exact affine law
   `4(F_a-F_b)=a4b-a4a`; a q4 gap at least `60B` therefore closes.  The honest
   q4488 residue is the close nonzero q4 pencil (equivalently the joint
   zero/small two-frequency branch), not missing phase matching.

**S45 update — all three items moved this cycle:**

1. kps's exhaustion rendering continued (T3-unconditional line); opus's THM-955
   cluster gap gives the CONTINUOUS a-priori safe width for every ≤6-speed base
   block, and its Lean L1 (`sorted_gap_pigeonhole`, LRCClusterGapBrick.lean) is
   now GREEN kernel-pure — with one statement correction: positivity 0 < b−a−B
   is NECESSARY (in-kernel counterexample; teeth may overhang the window).
   L2/L3 assembly remain with opus; consumers thread the positivity.
2. codex landed THM-954 (weighted-ratio-layer pair-floor); S54/S66 closed the
   strict-open pair-grid and covariance reindexing completely.  Tau3/tau4 are
   replayed in Lean, leaving five finite clique certificates.
3. THM-956 (LRCAdaptiveQ.lean, GREEN): Farey separation ⇒ ONE distinct-instant
   carrier per deep component on super-ladder strata (D² < 7·v_top) ⇒ the
   adaptive-q pigeonhole: whenever the window holds more coprime moduli than
   components, a deep-free q EXISTS, and THM-950's census there needs only
   liveCount > 0. Recon: 387/400 families have zero window-total deep; worst
   all-coprime carrier count = 1. Named remainder: (i) the K-refinement
   (provable 64·v_bot two-choice bound vs empirical K ≤ 2), (ii) instantiate
   comp concretely (component = the pinned witness index, mechanical from
   THM-949), (iii) the live floor at the chosen q — the program's nucleus,
   shared with item 1's exhaustion face.

**S57/S66 weighted-census refinement.**  The dense-core socket may now be stated
as `DenseCoreWeightedCensusB5Supply`, strictly weaker than the old uniform
`792·deepCount` socket.  Reuse multiplicities, zero-color gluing, fixed-root
Fubini, per-face gcd resonance budgets, and all-q midpoint parity are now
exact.  Its remaining arithmetic is: choose a root/window covering the
relevant events, bound the summed gcd budgets, pay the colored `462/210`
transport factors, handle depths six and seven, and supply the adaptive-modulus
live floor.

**Guardrails.**  The once-requested universal `q≤25` period bound is false
(THM-566/762/764; explicit blockers survive every `q≤25`).  Uniform emptiness
of the `n=12` sporadic tight branch also remains open: proper AP-centred H6
faces are certified through scale seventeen (scale thirteen is primitive-
impossible).  THM-982 records a provisional exact scale-eighteen owner
deficit.  A frozen primary C++ certificate has landed, but independent and
cross-build replays are still pending, so it is not part of the certified
frontier; scale eighteen and higher, ramified H5, non-AP/deep, and
higher-sheet branches therefore remain.  Neither the
Fano/`χ₇` address atlas nor the historical self-line analogy supplies these
missing metric/multiplicity statements; the proposed all-`n` black self-line
law is in fact refuted at `n=8` (`404` quasi-fixed versus `SC(8)=176`).

No theorem in this picture proves LRC(14) without the named suppliers above.
