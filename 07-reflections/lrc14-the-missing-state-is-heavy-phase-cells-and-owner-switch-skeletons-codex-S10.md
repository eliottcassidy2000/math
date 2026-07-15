---
date: 2026-07-14
source: codex-2026-07-14-S10
status: rigorous frontier synthesis; no claim that LRC(14) or the n=12 sporadic branch is closed
tags:
  - lonely-runner
  - lrc14
  - sporadic-branch
  - phase-cells
  - sheet-incidence
  - endpoint-owners
  - christoffel-words
  - tournament-analysis
  - controlled-forgetting
related:
  - THM-762
  - THM-769
  - THM-770
  - THM-771
  - THM-772
  - THM-773
  - THM-774
  - THM-775
  - THM-777
  - THM-778
  - THM-779
  - THM-780
  - THM-782
  - THM-783
  - THM-784
  - THM-786
  - THM-788
  - THM-789
  - THM-790
  - HYP-6820
  - HYP-6835
  - HYP-6840
---

# The missing state is heavy phase cells and owner-switch skeletons

## Honest status

The standard fourteen-runner case remains open.  The primitive tight
twelve-speed sporadic branch used by the current reduction is uniformly finite
in principle, but it has not been classified and is not proved empty.  This
session nevertheless changes the frontier in two material ways:

1. THM-780 proves a uniform positive safe-measure floor from the already known
   lower-dimensional strict margin.  The floor is explicit,

   ```text
   |G'_P| >= 182^(-12)
   ```

   for every twelve-speed core.  Height can no longer make the normalized
   regime disappear merely by sending its safe measure to zero.
2. THM-779's tempting empirical bound on the number of covered walls is false.
   An exact persistent seven-owner stalk supports arbitrarily long covered wall
   runs.  The corrected dynamic coordinate is the owner-switch skeleton after
   persistent stalks have been factored, not raw wall count.

The two changes point in the same direction.  The base of the problem is a
measured set in joint phase space; the obstruction over that base is a labelled
incidence-and-transport system.  Scalar period bounds, component counts,
unlabelled tournaments, and raw event counts are shadows of that object.

This is not an attempt to catalogue every historical analogy in the repository.
It isolates the chain of ideas that now has proved maps between consecutive
stages, and records exactly what is lost at each map.

## 1. What the successive failures taught us

The recent history is unusually coherent when read as a sequence of failed
compressions.

| proposed compression | exact content it really saw | first datum it forgot |
|---|---|---|
| a uniformly good denominator `q<=25` | zero owners and signed unit-pair blockers for each fixed `15<=q<=28` | scale and the adaptive modulus |
| residues modulo `13` | the shallow full-residue packet and its protected points | binding-scale depth `s` and lift ownership |
| number of occupied sheets | capacity necessary conditions | which labelled owner occupies which sheet, especially at equality |
| measure of a folded obstruction | a sharp necessary cap | component-by-component containment and endpoint colours |
| safe-component count or raw winding | a coarse amount of boundary | component widths, reduced winding, and endpoint-owner incidence |
| a static prime-seven token state | exact cover at one phase | endpoint order, inverse steps, and global carry |
| a transitive next-event tournament | pairwise chronological order | the labelled Hamiltonian-path movie and simultaneous blocks |
| number of covered walls | finite-sample difficulty | a persistent exact seven-owner stalk that makes those walls redundant |

None of these negative results says the corresponding viewpoint was useless.
Each viewpoint became useful when its lost coordinate was restored as a
sidecar.

### 1.1 Why the `q<=25` intuition failed

THM-762/764 give the exact small-period statement.  For `15<=q<=28`, a
covering family has a witness at denominator `q` exactly when it has no
`q`-divisible speed and its signed unit-pair blocker deck is incomplete.  This
is a clean finite theorem.  It does not imply a universal terminal denominator.

The coherent family

```text
26*{1,...,12} union {339}
```

has no witness through denominator `25` and first witnesses at denominator
`27`.  The gcd-incoherent family recorded in THM-762/764 first witnesses at
`26`.  Thus neither a large common-factor explanation nor its negation repairs
the bound.  The lesson is structural: raw denominators are not scale-normal.
The sheet constructions of THM-760/761 can close coherent families by choosing
an adaptive lift even when every bounded denominator bank is blind to them.

The correct finite object at a prescribed modulus is the obligation
hypergraph

```text
(multiplier, modulus) <- zero owner or signed-pair owner,
```

not a tournament ranking the moduli.  To globalize this object one must either
produce an adaptive missing obligation or classify simultaneous completeness
of the decks with scale retained.

### 1.2 Binding scales exposed the sheet fibre

THM-769 replaces a prime-residue snapshot by the reduced binding scale
`p/(13s)`.  The on-sheet speeds are multiples of `s`; the off-sheet speeds must
cover every lift of the loose set of the quotient core.  At `s=2`, a primitive
two-tightener packet is exactly

```text
A=2U union {x,y},       |U|=10,       x,y odd,
```

with persistent opposite parity ownership of the two lifts.  The three-owner
equality edge similarly produces a three-colour cover.  This is the point at
which “the runners have certain residues” becomes too weak: the theorem is a
cover predicate on an entire family of fibres, not at one rational time.

THM-772 then shows that the quotient core in the two-sheet packet is itself
primitive and contains a multiple of every `2,...,12`; its odd exceptions are
bounded relative to the quotient height.  In the three-sheet equality packet
there is an analogous divisor transfer.  These are genuine recursive
constraints, not merely residue telemetry.

### 1.3 Endpoint incidence repaired sheet counts

At seven exceptions, THM-771 replaces a Kirchhoff analogy by the exact defect
identity

```text
free sheets = capacity slack + overlap debt - ramification surplus.
```

In the unramified equality layer, a strict endpoint event deletes an owner at
equality and tears an exact seven-owner cover.  An entering danger arc cannot be
treated as handing off ownership at that same point: equality is safe.  This
strict endpoint convention is what turns a static sheet count into a proof.

The shallow branch exhibits the same lesson in a different form.  THM-770's
height-twelve classification is not proved by the residue tournament, whose
burdens all tie.  It is proved by an atomic endpoint-cell/owner incidence
hypergraph.  That finite exact result leaves higher lifts open; it does not
license a global extrapolation.

### 1.4 Folding made the two-sheet metric exact

For the two odd exceptions put

```text
a=(x+y)/2,       b=(x-y)/2.
```

THM-774 proves that eligibility and opposite sheet colours are equivalent to

```text
||a tau||+||b tau|| >= 11/13.                            (FD)
```

This is a valuable conjugacy: two labelled odd runners become one folded
diamond in half-frequency coordinates.  It gives the sharp universal measure
cap `8/117`, and it closes every quotient core in `[1,19]` even with unbounded
odd exceptions.  But the cap is only necessary.  Some loose cores already
have smaller measure, so the live predicate remains the pointwise containment
of every loose component in `(FD)`, including component position and width.

THM-775 explains the only possible imprimitive deletion of a quotient core.
It is dyadic: one unique odd guard sits above a primitive child, the first seam
is an exact `Z/4` ownership tiling of sizes `2+1+1`, and iteration gives a
binary sheet tree ending at a hereditarily primitive divisor-complete base.
The dyadic tower is a well-founded reduction of the gcd defect.  It is not yet
a proof that the terminal base or the tower is impossible.

### 1.5 Measure does not vanish: the heavy phase cell

THM-777 connected the normalized peel parameter to safe measure, but the
possibility of a sequence of taller primitive cores with safe measure tending
to zero remained the qualitative obstruction.  THM-780 removes it.

For `d` speeds, a `beta`-deep witness, and `alpha<beta`, partition the joint
phase torus into

```text
M^d,       M=ceil(1/(beta-alpha)),
```

half-open phase cells.  Some cell has orbit-preimage mass at least `M^(-d)`.
Subtracting one anchored orbit point from all other points in that cell gives a
simultaneous return of coordinate radius `<1/M`.  Translating that return set
to the deeper witness gives an `alpha`-safe set of the same mass.  With
`d=12`, `beta=1/13`, and `alpha=1/14`, this is `182^(-12)`.

The independent limiting proof identifies the correct height-free state as
the entire stabilized character relation lattice `L<=Z^d`; its annihilator
`L^perp` is the limiting compact subgroup, of whatever rank occurs.  This
supersedes attempts to guess in advance that only one- or two-dimensional
limit tori matter.  The deeper witness survives in the subgroup and lies in
the relative interior of the safe region, so Haar mass is positive.

The heavy-cell theorem supplies substrate, not targeting.  It says that safe
mass exists uniformly; it does not say which safe components meet a blocker
wall, which tooth owns their endpoints, or which continued-fraction convergent
lands in them.  Those are incidence questions on top of the measured base.

THM-782 specializes this substrate to the ten-speed quotient core of a
two-sheet packet.  Using a `1/11`-deep time and `72` phase cells per coordinate,
it produces an anchored return packet of measure at least `72^(-10)` whose
translate is strictly `1/13`-safe by an additional `1/10296`.  It also forces a
safe component of length at least `72^(-10)/(20 max(U))`.  Tightness places the
entire packet and component inside the folded diamond.  This is the precise
structured-noncontainment target that scalar comparison with `8/117` misses.

### 1.6 Endpoint schedules are centered mechanical words

THM-778 makes the event clock exact.  The danger endpoints of speed `u` are

```text
(2i+1)/(2u).
```

For a pair `u,v`, centered ceiling ranks give its complete merge word.  The
reduced continued fraction controls the Euclidean substitutions, but the gcd
and a transported parity bit must also be retained.  Simultaneous midpoint
walls occur exactly when the reduced pair is odd/odd, equivalently when the
two speeds have the same `2`-adic valuation; their positions form an exact
`gcd(u,v)` midpoint grid.

For many owners, centered Beatty ranks reconstruct the global wall order and
all simultaneous blocks.  This is a major gain: chronology no longer needs a
floating-point sort or an unstructured wall list.  It is still only the base
word of a skew product.  The same word acts differently with different
inverse steps, token assignments, carries, and core-safe intervals.

### 1.7 The prime-seven fibre and the corrected `r=8` object

THM-773 makes every unramified owner at the prime-seven lens a labelled token
in `F_7`; at its own wall the token is absent.  Coverage is the exact
polynomial divisibility predicate

```text
X^7-X divides product_a (X-k_a).
```

At eight owners, a covered generic state has one duplicated root:

```text
product_a (X-k_a)=(X^7-X)(X-z).
```

THM-779 normalizes by `z` and obtains the multiset
`{0,0,1,2,3,4,5,6}`.  Its event-word supportability equation is exact; the
rooted states form an `A_8` torsor of size `20,160`.  The normalized transition
graph is strongly connected, so no state-only potential can force an exit.
The mechanical schedule and metric core base are necessary.

This is where a finite-sample conclusion had to be withdrawn.  The observed
maximum of five consecutive covered walls is not universal.  On

```text
P={1,2,11,12,13},
W_0={1,4,5,6,8,9,10},
I=[25/182,27/182],
```

the seven owners in `W_0` form a fixed exact stalk.  The family

```text
7P union W_0 union {182m+1}
```

has `2m` consecutive covered walls in `I` and is divisor-complete.  These
extra walls do not measure new covering complexity: at each one, the
persistent seven-stalk already covers every sheet.

THM-783 supplies exact laws on the part that can carry genuine complexity.
Its wall-phi recurrence and period-sum identity force every complete
fast-owner period either to have no visitors or a zero-sum visitor cluster;
one visitor is impossible, fixed companion pairs de-phase, and absence of a
balanced companion gives a metric extent bound.  The persistent-stalk family
clarifies the scope rather than contradicting those laws: all of its fast
periods are visitor-free.  Thus THM-783 advances the owner-switch/visitor
coordinate while its sampled raw maximum `6`, like the earlier `5`, has no
universal status.  THM-784 independently proves the same raw-bound failure
with a simpler `21N`-wall family; the divisor-complete family above adds the
arithmetically stronger core-safe witness.

The surviving problem is therefore a first-return or owner-switch problem:
after contracting intervals supported by a persistent exact stalk, can the
centered mechanical schedule keep agreeing with the redundancy cocycle across
every relevant part of the core-safe set?  THM-778 already settles one useful
stratum: sufficiently dense equal-valuation double walls pierce because two
owners disappear simultaneously.

### 1.8 The even-maximum collar is now a clocked incidence problem

THM-790 starts from the exact remaining two-sheet max-peel collar rather than
asserting that the collar is empty.  If `R=max(U)`, `max(A)=2R`, and deleting
`2R` raises the margin above `1/12`, then an eleven-speed maximizer has reduced
denominator

```text
q<=4R-2.
```

It lies strictly inside a `2R` danger tooth, and either `q|2R` or the effective
order of `2R` modulo `q` is at least `14`.  The center case pins a missing
divisor and meets THM-775's dyadic seam when the deletion core is imprimitive.

Applying THM-780 one level down gives strict safe mass at least `156^(-11)`.
Tightness forces this mass into the top danger comb, hence into many occupied
top teeth.  Every occupied tooth has disjoint nonempty left- and right-flank
blocker sets, so an ordered flank-owner type repeats.  This converts uniform
measure into repeated boundary incidence.

In the deep branch, a forced odd exception `w=13c` gives a second finite
carrier.  When `w` is ineligible, the ten quotient runners form labelled
Cayley two-edges covering `Z/13Z`.  Their total overlap excess is seven.  A
simple owner event slides one edge and must transport one overlap chip.  The
faithful local state is therefore

```text
(mechanical event word, ten labelled Cayley edges, seven overlap chips).
```

This is a reduction, not an exclusion theorem.  In particular, an integer-grid
edge event should not be silently identified with a midpoint wall clock; the
appropriate mechanical word and phase convention must be carried explicitly
when this automaton is joined to THM-778.

The first bounded automaton tear is also exact.  If the forced exception is
`w=13` and `U subset [1,24]`, the event word on the ineligible interval
`2/13<s<11/13` eliminates every ten-core.  Of the
`binom(23,10)=1,144,066` possible cores, `101,850` cover the initial chamber
and zero cover the whole 117-group event word.  The primitive
divisor-complete subatlas likewise has `20,604` initial static covers and zero
survivors.  This is finite-exact for the stated height and forced exception;
it is not a uniform collar theorem.  Its conceptual value is that static edge
covers are plentiful while dynamic covers vanish in the first complete bounded
atlas: the moving labels and event order, not initial degree counts, do the
separating.

## 2. A candidate underlying object

The most economical object consistent with all of the proved corrections is a
**measured, stalked event system** over the quotient circle.  For a core `P`
and an exception deck `W`, its fields are

```text
base:
  the strict core-safe set G_P, with component endpoints and measure;

phase substrate:
  a heavy joint-phase cell or, in a limit, the character relation lattice;

event address:
  exact wall positions, centered/global ranks, simultaneous blocks,
  and the Euclidean phase cocycle;

fibre:
  labelled owner-to-sheet incidences, strict endpoint deletion,
  inverse steps, and global carry;

stalk decomposition:
  inclusion-minimal owner sets that already cover throughout an interval;

essential dynamics:
  changes of those minimal stalks, together with redundancy roots,
  flank owners, or overlap chips.
```

This package is larger than a tournament and much smaller than the original
speed vector.  It forgets absolute height when the relation lattice or
normalized clock suffices, but it does not forget the coordinates that decide
coverage.

One precise proposed definition is useful for the `r=8` route.  For a
wall-free phase `x`, let

```text
B(x)={B subset W : B is inclusion-minimal and its strict bad-sheet
                    incidences cover every sheet at x}.
```

A **persistent stalk on `J`** is a labelled `B` that covers every phase of
`J`, with the strict wall convention included.  A wall is **stalk-redundant**
if some persistent stalk not using its event owner covers a neighbourhood of
that wall.  Contract maximal runs of stalk-redundant walls.  The remaining
marked sequence, with the minimal-cover sets on both sides, is the
**owner-switch skeleton**.

This definition deliberately makes the THM-779 family simple: its `2m` walls
are redundant over the fixed `W_0` stalk.  It also explains why the definition
must remain over the metric base.  A stalk may persist on one core-safe
component and fail on another even when the abstract token state and owner
word are isomorphic.

The proposed skeleton is not yet proved to have bounded complexity.  It is a
better theorem target because its complexity cannot be inflated by adding a
redundant high-frequency owner to an already exact stalk.

## 3. Tournament Analysis after challenging the vertices

The repository's tournament viewpoint remains useful, but chiefly as a test of
whether a quotient has forgotten the proof predicate.  The standard endpoint
tournament uses owner clocks as vertices, orients `u->v` when `u` has the
earlier next endpoint, switches between chronological order and the Euclidean
shear presentation, and resolves a simultaneous block by the fixed owner-order
Hamiltonian path.  It is transitive in every chamber.  Its score histogram,
SCCs, directed cycles, and Hamiltonian-path count therefore stay essentially
constant while the actual sheet cover changes.

The following vertex audit is more informative than another runner
tournament.

| candidate vertices | relation / switch | what the quotient preserves | what it destroys |
|---|---|---|---|
| runners or exception owners | earlier next wall; chronology versus reflected/Euclidean gauge | labelled pairwise order and edge flips | joint phase, metric scale, simultaneous ownership, sheet cover, and core-safe targeting |
| moduli `q` | smaller blocker cost; pair-first versus compression-first gauge | fingerprints of a fixed small-period bank | multiplier identity and joint zero/signed-pair ownership; the tournament verdict can stay fixed while most edges flip |
| witness obligations `(q,a)` | blocked by a named runner; switch multiplier or sign | the exact small-period predicate when kept as a hypergraph | scale-adaptive witnesses if collapsed to pairwise comparisons |
| joint phase cells | same-cell equivalence; anchor subtraction | the simultaneous return and its measure | endpoint chronology and local blocker ownership; the relation is an equivalence, not naturally a tournament |
| character obligations `m in Z^d` | eventual vanishing under a height sequence | the limiting subgroup and all stabilized linear relations | which safe component meets which wall; subgroup closure is again not antisymmetric |
| safe components or top teeth | cyclic order, width, ordered flank owners | repeated collar incidence and tooth occupancy | exact sheet transport if reduced to an order tournament |
| sheets | owner incidence or Cayley-edge incidence; circle reflection as gauge | the fixed-phase cover predicate and overlap degree | endpoint order, inverse step, event owner, and future continuation |
| endpoint events | centered rank; chronological versus Euclidean decoding | the complete wall/tie schedule | current token assignment, persistent stalk, and metric restriction to `G_P` |
| rooted redundancy states | collision-hop transitions; choice of duplicated root | the local `A_8` cocycle and supportability of an event word | absolute root, global carry, metric schedule, and the core-safe base; strong connectivity prevents a state-only exit proof |
| minimal cover stalks | persistence across a wall; before/after minimal bases | genuine owner switches after redundant walls are contracted | raw event multiplicity by design; component location must remain as a sidecar |
| proof obligations `(event, chip)` or `(tooth side, owner)` | transports or discharges the obligation | exactly the boundary incidence used by THM-790 | raw runner geometry unless the event clock and edge label remain attached |

The challenged conclusion is that there is no single canonical tournament
vertex set.  Different theorems live on cells, characters, events, sheets,
stalks, or proof obligations.  When the natural relation is equivalence,
subgroup closure, hypergraph incidence, or labelled transport, forcing it into
an antisymmetric tournament discards information without adding leverage.

Tournament fingerprints are still valuable telemetry:

- edge flips reveal sensitivity to the chosen gauge;
- directed cycles can witness genuine nontransitive transport when they occur;
- SCCs of a finite transducer distinguish trapped from escaping state regions;
- Hamiltonian-path counts measure how much labelled chronology an isomorphism
  node has collapsed.

But the tie Hamiltonian path must stay attached to the event schedule, and a
fingerprint is not a substitute for the owner-incidence sidecar.  THM-778's
transitive event tournaments and THM-779's strongly connected `A_8` state
graph are complementary warnings: the former is too small, while the latter
is dynamically rich but detached from the metric base.

## 4. The exact remaining mathematical pieces

### 4.1 Uniform emptiness of the primitive tight twelve-speed sporadic branch

The branch is finite in principle, but the following uniform arguments are
still missing.

1. **Shallow full-residue rigidity beyond the bounded lift box.**  THM-770
   closes lift height twelve and the corresponding primitive range, but a
   scale-free endpoint-owner coherence theorem or a complete finite execution
   of the global bound is absent.
2. **The two-sheet folded branch.**  One must exclude every primitive
   divisor-complete quotient core, not only `max(U)<=19` or the bounded-height
   bank.  A proof may show that the folded containment `(FD)` cannot persist,
   or that THM-775's binary ownership tree cannot reach a terminal
   hereditarily primitive base.
3. **The even-maximum collar.**  THM-790 gives a bounded rational clock,
   repeated ordered flank types, and—on the forced `13`-multiple subbranch—a
   moving edge-cover automaton.  It excludes `w=13` through quotient height
   `24`, but what is missing is a uniform tear or congruence obstruction
   linking those three structures beyond the bounded atlas.
4. **Higher-sheet packets.**  The three-sheet equality edge and general deep
   `s`-fibres still require a classification of persistent colour covers with
   effective orders, ramification, and omit-one gcd descent retained.
5. **Assembly across deletion choices.**  A tight twelve-set must satisfy the
   leave-one-out and binding-scale constraints simultaneously.  Existing
   theorems constrain each packet sharply, but no theorem yet proves their
   global incompatibility.

### 4.2 The `r=8` prime-lens exit

The corrected alternatives are precise.

1. Prove that no deck blocks the **entire** closed core-safe set, or classify
   the exceptional whole-core covers.
2. Define and bound the owner-switch skeleton after persistent exact
   seven-stalks are contracted.  A raw wall bound is impossible.
3. Use THM-783's period-sum, single-visitor, cluster-balance, and de-phasing
   laws to bound the genuine visitor-cluster word or its metric extent after
   stalk contraction.
4. Combine THM-779's supportability equation with THM-778's centered schedule
   to prove non-synchronization on at least one core-safe component.
5. Complete the simultaneous-wall analysis below the current large-gcd pierce,
   including unequal `2`-adic valuations and core intervals shorter than the
   double-wall mesh.
6. Generalize the redundancy polynomial/stalk language to `r>8` and to
   ramified or non-prime effective lenses without losing strict endpoint data.

### 4.3 Making the global normalized atlas usable

THM-780 proves that the normalized regime is bounded, but its explicit cutoff
is enormous.  Three different improvements would be meaningful and should not
be conflated:

1. prove the conjectural sharp safe-measure floor `7/858` and its uniqueness;
2. obtain any practically sized uniform floor or structural dichotomy;
3. exploit the existing crude floor without enumerating the entire resulting
   atlas, by using repeated phase cells, teeth, or owner types recursively.

The third option may be the most robust.  The phase-pigeonhole proof produces
a heavy return set whether or not the sharp extremizer is understood.  A
successful argument could pass that mass to an incidence repetition, as
THM-790 begins to do, and then force a finite automaton tear.

### 4.4 LRC(14) assembly

Even a complete sporadic-branch theorem must be inserted into the current
covering/non-covering and scale-normal decomposition with its hypotheses
checked.  The live global residual still includes multi-exception sheets,
ramification/s-threshold decks, gcd descent, and families without an obvious
large common-factor core.  The proved bounded and sheet regimes are not yet an
assembly theorem for LRC(14).

## 5. A recursive research program suggested by the synthesis

The history suggests a disciplined recursion rather than a new scalar.

1. **Find mass.**  Use a deeper lower-dimensional witness and the heavy-cell
   argument to obtain a positive measured base.
2. **Expose boundary incidence.**  Decompose that base into exact components or
   top teeth and retain ordered endpoint owners.
3. **Choose the natural fibre.**  According to the binding scale, use `s`
   sheets, a prime-field token deck, a dyadic sheet tree, or the `Z/13` Cayley
   edge cover.
4. **Decode chronology.**  Use centered/global ranks with gcd, parity phase,
   simultaneous blocks, and metric positions retained.
5. **Factor persistent stalks.**  Remove event multiplicity already absorbed by
   a fixed exact cover.
6. **Study the essential skeleton.**  Apply the redundancy cocycle, overlap-chip
   transport, divisor shells, or endpoint-owner constraints only to genuine
   changes of minimal cover basis.
7. **Descend or tear.**  Either exhibit a free sheet/time, force a dyadic or gcd
   descent, or return a smaller normalized packet with all sidecars intact.

The key invariant is not “how complicated is the speed vector?”  It is “how
many genuinely different minimal cover stalks can the measured core-safe base
support under its exact mechanical schedule?”  That question survives the
counterexamples that defeat bounded denominators, raw component counts, and
raw wall counts.

## Closing reframe

The accumulated work has not found one magical tournament on runners.  It has
found a hierarchy of controlled forgettings:

```text
speed vector
  -> joint phase orbit / relation lattice
  -> measured core-safe base
  -> endpoint-ranked mechanical word
  -> labelled sheet-incidence fibre
  -> persistent cover stalks
  -> owner-switch skeleton
  -> lonely witness or a smaller residual packet.
```

The first arrow explains why height need not destroy safe mass.  The middle
arrows explain why boundary equality, labels, gcd, parity, and carry repeatedly
reappear after an overaggressive quotient.  The final arrow is the open one.

The sharpest current view is therefore not that LRC(14) is waiting for a
larger computation.  It is waiting for a theorem that couples a heavy measured
phase base to a finite labelled switch skeleton.  THM-780 supplies the former;
THM-778/779/783/786/788 and THM-790 supply complementary increasingly exact
versions of the latter.
The remaining work is to prove that the coupling cannot stay covered without
descending to one of the already controlled equality packets.
