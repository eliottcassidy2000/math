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
depth-six/seven payment, and a positive safe-period interval-table supplier at
the adaptive modulus are still missing.  The table-to-liveCount composition is
now kernel-checked in `LRCLiveFloorSampling`.

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
  producer, and tau3--tau5 ratio replays are LEAN; tau6--tau9 remain.  Separately,
  the sparse rational-pair branch interval now transports bijectively back to
  multiplier space, giving the complete finite-`q` three-branch count; the
  zero-width integral boundary requires `Int.toNat`, not a raw floor difference.
- **death-star S52/S53 two-circle theorem (THM-985)**: for the canonical family
  `{1,...,13}`, `LRCTwoCircle` proves that the integer and half-resonance
  circles lie inside the depth-six set, and `LRCTwoCircleII.deep_iff_circles`
  proves the converse.  The smallest-failing-speed branches `k₀=1,2` are
  explicit; `compat_card_le` collapses every middle case `k₀=3,...,8` to one
  finite kernel decision, and `k₀>=9` is impossible by cardinality.  Thus the
  canonical deep set is exactly the two circles in-kernel.  THM-987 now counts
  their disjoint union for every `q>=2`, including the half-circle parity term,
  and proves the canonical tight family lonely through the `q=14` census
  (`deep=1`, `live=6`).  THM-991 is now exact in `LRCLiveLawExact`: for every
  positive `q`, the canonical live count is six when `14|q` and zero otherwise;
  on `q=14m` the live set consists exactly of `m` times the six units modulo
  fourteen.  `LRCTwoCircleCount` exposes the same structure as
  literal low/high/half Finsets and proves their disjointness, parity and
  compact card formulas, and `canonicalDeepMultipliers 14 = {7}`.  These are canonical-family closures, not the uniform
  weighted/trapped-core supplier.
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
  Lebesgue measure for the block's common safe complement.  The same module
  now proves the `c=7` wall theorem: at exactly unit union budget, one positive
  adjacent-pair credit already forces a positive safe complement.
  `LRCC8ConsecutiveWitness` extracts literal `Fin 7` and `Fin 8` `Lonely 14`
  witnesses from the two positive-measure statements.
  `LRCC8ConsecutiveWitness` extracts a concrete point and packages a literal
  `Lonely 14` witness for that `Fin 8` family without duplicating the measure
  proof.  This is the live generic shifted `c=8` theorem; it does not add the
  other five danger combs, so it is neither the full thirteen-speed LRC(14)
  theorem nor the weighted dense-core supplier above.


## THE 7-WALL EXISTENCE CHAIN (opus S346–S351) — all LEAN, kernel-pure

Seven modules, every one at `[propext, Classical.choice, Quot.sound]`, no
`sorry`, no `native_decide`. Together they take the 7-wall's EXISTENCE
conclusion from end to end without the sawtooth identity:

| module | content |
|---|---|
| `LRCFoldedIdentity` | THM-965 `muNum_folded`: `14·muNum a b = 4ab + fold((a+b)%14) − fold((b−a)%14)` |
| `LRCFloorTable` | per-class floor `14·muNum ≥ 4ab − 49` (fold ∈ [0,49] via (r−7)²≥0); `overlap_floor_rat` |
| `LRCFloorAllocation` | the allocation law `r_j^(j+1) ≤ ∏ r ≤ C` for antitone ratios |
| `LRCWindowAverage` | the Fubini position step: inner integral = L, Tonelli window-average, `live_window_exists` |
| `LRCHunterAssembly` | uncovered ≥ Σ consecutive path-tree overlaps (tree-Hunter + `measure_compl`); n=7 capstone |
| `LRCSevenWallExistence` | positive uncovered ⟹ a LONELY POINT; the `∑ μ ≤ 1` weakening |
| `LRCPairOverlapFloor` | pair-overlap LOWER bound by containment (+ the gcd strengthening via the common period 1/g) |
| `LRCCombUpperBound` | the SHARP single-comb bound `≤ 2λ` on a half-cell-shifted unit window |
| `LRCArcCounting` | THE COUNTING LEMMA: `m` half-cell-aligned cells hold `m` whole arcs (engine of THM-1012's sharp nesting floor) |
| `LRCPairIndependence` | **THM-1012 LANDED**: the pair overlap reaches the independence constant `4λ²` up to a linear defect — general consecutive-cells induction + alignment finder + assembly |
| `LRCSharpWallBound` | **THE WALL WIRED**: `sharp_wall_bound` / `seven_comb_wall` — comb upper bounds + pair lower bounds ⟹ a lonely point INSIDE the window (bridge = `volume.restrict W`) |

**Two architectural findings recorded along the way.** (1) LRC's conclusion
needs a lonely TIME, not a positive-measure window — so the circle/line
reconciliation is OFF the existence critical path (window machinery is for
NESTING only). (2) The assembly's `∑ μ(A i) = 1` weakens to `≤ 1`, so only an
UPPER single-comb bound is needed. Together these removed the exact sawtooth
identity (THM-965/856) from the existence path entirely; it is now reserved
for the SHARP nesting floor (THM-1012, paper-proved by period counting —
also without the sawtooth).

## What is genuinely open (the honest short list)

1. **The residual's dense core** past all shrinks: packets admitting no corrected
   THM-959 block partition, including a moderate-junction branch not yet reduced
   to one seven-comparable block. Attacks: an explicit partition trichotomy,
   kps exhaustion (paper, Lean pending), and the census funnel.
2. **The pair-floor finite replays** — the pair continuum producer, component
   budget, Turan arithmetic, top path caps, and clean-534 transfer are closed.
   Tau3 through tau5 compile; the exact remaining tasks are tau6--tau9.
3. **A-priori liveCount floors** (the discrete side of the same wall): the
   winData rows give them on censused strata; beyond, it is kps's law made
   uniform. SHARPENED (S44): the naive coprime law is REFUTED (Dirichlet-rate
   mirror-paired exceptions); the question is now adaptive-q — for each family,
   SOME q in the window with live > 792·deep (deep sets provably mirror-symmetric,
   even at odd q — THM-953).  `LRCLiveCountLonely` now records the simpler
   terminal fact explicitly: once `liveCount(v,q)>0` is known, no cap or census
   premise remains—one live multiplier is already an inclusive rational
   `Lonely 14` witness.  Producing that positivity is still the crux.
4. **Selected-witness relation residues** — q22 remains an opposition-avoidance
   selector.  `LRCSelectedWitnessRelationRouter` now closes the exact-zero
   q333, q244, and q4488 alternatives: each becomes a signed unit support-three
   relation on the original speeds, and a factor-three dense ladder forces its
   top position to `lastDensePair+1` or lower.  The q4488 common-sign and
   matching-wall producers are closed, and its large branch escapes.  The two
   q4-derived frequencies satisfy `4(F_a-F_b)=a4b-a4a`; a q4 gap at least `60B`
   closes.  The honest residue is genuine small nonzero frequency plus the
   missing common-phase/payment input—shifted wall phases cannot be silently
   fed to the real relation lock.

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

**S67--S71 modulus-window and live-floor audit.**  Connected all-bad zero-color seven-stalks
are forced into exact reduced-modulus resonance for every `q <= 98`; the
published `q=99,d=98` example makes this cardinality cutoff sharp.  THM-979 and
THM-984 propose large explicit sampling moduli from continuous `B5` or safe-set
floors, respectively.  These advances are complementary, not composable by
size alone: for thirteen distinct positive magnitudes `E=2*sum|v_i|>=182`, so
THM-984's `ceil(2E/mu0)` is at least `364`, already outside the stalk window.
`LRCGridSampling` proves the abstract separated-interval lower bound
`card >= qL-n`, and `LRCLiveFloorSampling` now supplies the LRC application
layer: safe interval-table samples inject into `liveCount`, yielding
`liveCount >= q*mu0-error`, the strict explicit
`q0=ceil((error+1)/mu0)`, and direct/capped-five loneliness consumers.
`LRCGridCount` independently supplies the rational one-interval kernel.  The
actual strict-safe-set interval-table producer, its component/measure budget,
`T_s` assembly, and exhaustive trapped-core reduction remain.  Large-q
small-reduced-speed stalks and colored `462/210` payment are still honest
sockets.

The complementary margin route is also exact on the consumer side:
`LRCPairOverlapArcs.good_interval_of_margin` and `good_measure_of_margin`
turn an instant with `M>1/14` into an explicit safe interval and positive
measure floor.  The trapped-core analytic question is therefore the rigidity
of the tight equality case `M=1/14` plus a quantitative strict margin on every
non-tight residual; the Lipschitz converter itself is no longer debt.

**Guardrails.**  The once-requested universal `q≤25` period bound is false
(THM-566/762/764; explicit blockers survive every `q≤25`).  Uniform emptiness
of the `n=12` sporadic tight branch also remains open: proper AP-centred H6
faces are now certified through scale thirty-one (scale thirteen and its
multiples are primitive-impossible, and prime scales at least nineteen are
uniformly excluded).  THM-982's scale-eighteen owner deficit has independent C++
and Python exact replay, and THM-983 uniformly excludes every prime common
scale `p>=19`.  THM-986's composite scale-twenty owner deficit now has
independent C++ and Python exact replay.  THM-988 independently closes `c=21`;
its two all-order rows are scalar-tight at every owner but have exact local
maximum union `20/21`; a new cubic-character argument explains that projection
debt symbolically.  THM-989 independently
closes `c=22`: literal-CRT/C++ and algebraic-CRT/Python agree on ten hardened
certificate banks, and every row has at least five empty owners.  `c=23` is
prime.  THM-990 independently closes `c=24` (66,984 scalar rows, at least two
empty owners) and hashes all 101,961,528 reachable masks.  THM-992's
owner-normalized `Z/25 -> Z/5` fibre proof reduces `c=25` to 36 rows and makes
every owner miss at least three sheets; a primary and two independent referees
agree.  Scale 26 is multiple-of-thirteen excluded.  THM-993/994 close `c=27`
in three implementations: nested `Z/3,Z/9` relaxations kill every row, and a
separate nine-fibre flag proof makes every order-27 owner impossible; heredity
supplies at least two.  THM-1072 independently closes `c=28`: a
`Z/28 -> Z/4` fibre-capacity relaxation makes every one of 3,170 scalar
survivors impossible at at least four of its six owners, and three
implementations agree on all 6,628,500 exact reachable-mask incidences.
`LRCNestedFibreRelaxation.lean` kernel-checks the generic anchor/deviation
upper relaxation used here and at `c=27`; the large concrete banks remain
external.  Scale `c=29` is prime-excluded.  THM-1090 then closes `c=30` with
complementary mod-six and mod-ten thick-fibre relaxations: the common
live-owner histogram is `0:45998,1:6852,2:1200`, so no row reaches all six
owners.  An independent C++ referee reproduces the complete scalar bank and
  checks 18,874,368 literal residual words against both bounds.  Scale `c=31`
is prime-excluded.  THM-1096 closes `c=32` by a `Z/8` thick-fibre
relaxation: independent Python and C++ certificates agree on the 3,450 scalar
survivors and live-owner histogram `0:2802,1:456,2:192`.  THM-1124 closes
`c=33` by a `Z/3` anchor relaxation over the complete squarefree grammar, with
independent Python and C++ agreement.  The next untreated composite common
scale was `c=34`; THM-1125 now closes it by an exact `Z/2` anchor relaxation
whose 3,312 surviving owner bounds are all at most `29<34`, independently
replayed in Python and C++.  The next untreated common scale is `c=35`.
The clustered large-speed branch has also sharpened independently of this
ramified-scale ledger.  THM-1094 proves the uniform two-comb component theorem
and THM-1097 proves the uniform three-comb theorem, so clustered `r=3,4` are
closed at all scales.  `LRCSharpCombArithmetic.lean` kernel-checks their ratio
tails, final-killer inequality, exceptional rectangle, toothpick phase
transition, THM-1126 gap-energy arithmetic, and THM-1128's thirteen-grid scale
inequalities without `sorry` or `native_decide`.  It now also checks
THM-1129's needle-tail arithmetic, THM-1133's three-range dispatch, and
THM-1134's five-grid cone, step-two, and separated-ratio scale inequalities;
the interval discrepancy,
66/220-core atlases, and large exact banks are still external certificate
producers.  The phase quotient fails at four removals, and THM-1101 now has an
exact covering tuple above 235 missed by both sides of its former split.
THM-1123 nevertheless proves the sharp component target on a complete
45,238,050-row bounded bank.  THM-1126 replaces scalar mass by the truncated
autocovariogram/gap energy; THM-1127 proves the `(0,4,5,9)` fixed ray at every
scale; and THM-1128 proves every offset shape in the uniform Kakeya cone
`k1>=max(1272,26(k4-k1+1))`.  THM-1133 now closes every offset shape of span
at most 30 for every legal scale: THM-1123 bottom, a 3,539,936-row exact
middle, and THM-1129 individual tails meet without gaps.  Thus uniform
clustered `r=5` is reduced to spans above 30 outside the centred cone and the
other analytic gates.

At `r=6`, THM-1134 removes the apparent five-centre obstruction by choosing
the thirteenth-grid multiplier: ten affine orbits prove every five-residue
pattern has some chart with gap at least `6/13`.  This yields the cone
`B>=17 max(A,80)` and closes the entire step-two family by a 792-core rectangle
atlas plus 12,771-row finite complement.  Its complementary exact `Q5` gate
handles sufficiently separated ratios.  Uniform r6 outside those pieces
remains open.
Ramified H5, non-AP/deep, and higher-sheet branches
remain.  Neither the
Fano/`χ₇` address atlas nor the historical self-line analogy supplies these
missing metric/multiplicity statements; the proposed all-`n` black self-line
law is in fact refuted at `n=8` (`404` quasi-fixed versus `SC(8)=176`).

The former short capstone through `INVcov` is kernel-valid only as a
conditional implication.  THM-1158 refutes its literal premise even after
`Covering(2..14)` was added: the doubled AP `2*{1,...,13}` is Covering, has no
`Lonely13` time, and has no 13-dominant speed.  `LRCMSplit.lean` now labels
that interface historical and root-imports a formal counterexample.  The
sound exact capstone is `ResidualINV`, with
`residualINV_iff_LRC14`; this is an equivalence/restatement under the AP
bridge, not progress on the remaining mathematics.  A live noncircular route
must normalize gcd, rederive Covering after reduction, or use the explicit
Easy/Compact producers in `LRC14DispatchAssembly.lean`.

No theorem in this picture proves LRC(14) without the named suppliers above.

## 2026-07-18 boundary update — sieve, modular charts, and the AP-core consumer

The top-level divisibility dispatch is now **LEAN**.  `LRCSieveDispatch` defines

```text
Covering(v) := every n with 2 <= n <= 14 divides at least one speed,
```

and proves `not Covering(v) => exists t, Lonely 14 v t`, using the literal
witness `t=1/n` for a missing modulus.  Its dichotomy and capstone give

```text
CoveringCase [OPEN]  =>  LRC(14),
```

for positive thirteen-speed rows.  All four printed dispatch theorems audit to
`[propext, Classical.choice, Quot.sound]`.  Thus non-covering is proved; the
positive covering case, not the dispatch, is the open mathematical crux.

The representation-level chart law is now **LEAN** and exact:

```text
∃ k : ℤ, d ∣ u + qk    ⟺    gcd(q,d) ∣ u.
```

`LRCRationalScaleGuardrails` proves this for `q,d : ℕ`, `u : ℤ`, including
zero edge cases, by an explicit Bézout lift.  It also separates the two critical-
window inequalities (`q < 14 val` and `13 val < q`) and proves that a coprime
multiple/affine range is all of `ZMod q`, not a proper coset.  Its affine-lift
and coprime-range theorems audit to `[propext, Quot.sound]`; the combined ordered-
field window theorem audits to the standard foundational trio.

The compact/AP-core implication is also **LEAN**, while its supplier remains
**open**:

```text
CoveringCase [OPEN]
       └─ hard compact / M<1/14 class
          └─ INV Compact [OPEN: this chosen class supplies 13-fold dominance]
       + LRCUpTo13 [explicit named citation]
       ⇒ ap_core_bridge / ap_core_bridge_of_shape [LEAN]
       ⇒ ∃ t, Lonely 14 v t.
```

Here `Compact` is deliberately abstract.  `INV` is not the false blanket claim
that every covering row has a 13-dominant speed (`{2,...,14}` is an immediate
counterexample); it is the still-missing inverse theorem for the isolated hard
class.

The density/far-extension route is a **LEAN conditional consumer**, not another
independent branch:

```text
frame 1/13-lonely at supplied t0
 + |v_i| <= V for the frame, V>0
 + 91V <= v_far
 ⇒ frame stays 1/14-lonely on radius 1/(182V)
 ⇒ that interval contains a half-integer phase for v_far
 ⇒ Lonely 14.
```

`density_far_bridge` obtains the frame time from the named `LRCUpTo13` citation,
but it still assumes the far index, the bound `V`, and 91-fold separation.  It
does not supply those hypotheses, a universal density estimate, `CoveringCase`,
or `INV`.  Moreover `91V <= v_far` is strictly stronger than the 13-fold
dominance used by `ap_core_bridge` when the positive frame is bounded by `V`.
`LRCDensityDischarge` is therefore a genuinely different interval-completion
proof of an already covered subcase, not extra family coverage.

### Corrected top-level picture

The former global premise “no Lonely-13 time forces a 13-dominant speed” is
**FALSE**, not the open crux: `{1,...,13}` has exact maximum `1/14`, no
Lonely-13 time, and no such dominant speed.  Its old Lean consumer was a valid
implication from a false hypothesis.  The corrected live graph is:

```text
                         no Lonely14
                              |
                    sieve_dispatch [LEAN]
                              |
                           Covering
                              |
 ResidualINV [OPEN: positive + Covering + no Lonely14 => dominance]
                              |
             Lonely14 OR 13-fold dominance [LEAN]
                       /              \
                 immediate     LRCUpTo13 + ap_core_bridge
                       \              /
                LRC14_of_residual_INV [LEAN]
                              |
             lonely_fract + Finset enumeration [LEAN]
                              |
          LRC14_finset_of_residual_INV : LRC14.LRC14 [LEAN]
```

`M_split` and `crux_of_dominance` remain correct pedagogical/scoped lemmas, but
they are not dependencies of this capstone.  `ResidualINV` is an exact
counterexample-structural target and, with the cited bridge, is logically
equivalent to LRC(14); it is not identified with Tao's `n=12` inverse without
another theorem.

The noncircular route remains `LRC14DispatchAssembly`: explicitly choose
predicates `Easy` and `Compact`, prove `CoveringSplit Easy Compact`, supply an
`EasyCase Easy`, and prove `INV Compact`.  This generic chain is honest because
the easy/already-lonely branch is a visible input rather than silently forced
into dominance.  Density discharge and rational guardrails remain useful side
consumers, not capstone dependencies.

The root-imported modules use no proof-hole `sorry` and no `native_decide`.
The genuine mathematical gap is a predicate-preserving residual
classification—not a representation bridge and not the refuted universal
no-Lonely13 dominance statement.

## 2026-07-18 S74 frontier picture -- three residual objects, not one equivalence

The compact, twelve-speed, and four-comb lines now meet as follows:

```text
FOUR-COMB r=5
  span<=30 [THM-1133]
  OR four-residue cone [THM-1148]
  OR Q4>0 [THM-1148]
  OR exact Phi transfer (9/5 corollary) [THM-1148]
       |
       +-- residual begins with m(3,4,5,6) [OPEN method gap]

TWELVE-SPEED SHALLOW
  full residues
       |
  A12 mechanical word [THM-1143; arithmetic LEAN]
       |
  owner-stalk lowering law [THM-1150]
       |
  PEHD13 = central colour elimination + oriented fringe [OPEN]
       |
  height<=12 + no final AP lift [THM-770/795]

TWELVE-SPEED DEEP s>=2 -------------------------------- [OPEN]

COMPACT M<1/13
       |
  tight deletion OR all-loose essential crown [THM-1149]
       |                              |
       |                              +-- Cover14 crown collapse [OPEN]
       v
  n=12 equality classification [HYP-4382, OPEN]
       |
  regenerated d[12] core
       |
  13d|v [THM-1149]
       |
  primitive + 14-carrier + rho<13 contradiction [LEAN]
```

This corrects the S113/S114 identification.  The compact `1/13` floor is a
stronger sufficient route than the LRC(14) `1/14` target.  The n=12 equality
probe classifies a tight deletion but does not extract one.  Crown collapse
and shallow/deep equality rigidity are separately typed suppliers.

The new root-imported `LRCCompactEssentialCrown` module checks the finite
private-mass algebra and the post-regeneration divisibility/ratio
contradiction.  Analytic crown extraction, Farey regeneration, PEHD13, and
the deep branch are not formalized as proved theorems.  The A12 owner-stalk
functional law is exact in the research artifact but remains outside Lean.

Carrier choice is now explicit.  Four-comb proofs need labelled cyclic gaps,
cluster count, endpoint owners, and metric interval widths.  Shallow descent
needs moving two-sheet needles, private stalks, strip index, and an oriented
fringe bit.  Compact crown collapse needs the private/`>=3` overlap chamber
complex plus shared Cover14 lift congruences.  The corresponding naked
tournaments are transitive shadows and do not preserve the target predicate.
