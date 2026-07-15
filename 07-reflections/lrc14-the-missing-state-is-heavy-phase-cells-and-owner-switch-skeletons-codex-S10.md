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
  - packet-holonomy
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
  - THM-791
  - THM-792
  - THM-794
  - THM-795
  - HYP-6820
  - HYP-6835
  - HYP-6840
---

# The missing state is heavy phase cells, owner-switch skeletons, and packet holonomy

## Honest status

The standard fourteen-runner case remains open.  The primitive tight
twelve-speed sporadic branch used by the current reduction is uniformly finite
in principle, but it has not been classified and is not proved empty.  This
session nevertheless changes the frontier in four material ways:

1. THM-780 proves a uniform positive safe-measure floor from the already known
   lower-dimensional strict margin.  The floor is explicit,

   ```text
   |G'_P| >= 182^(-12)
   ```

   for every twelve-speed core.  Height can no longer make the normalized
   regime disappear merely by sending its safe measure to zero.
2. THM-779's tempting empirical bound on the number of covered walls is false.
   An exact persistent seven-owner stalk supports arbitrarily long covered wall
   runs.  Persistent stalks must be factored; raw wall count is not complexity.
3. THM-794 proves that counting the remaining active fastest periods is still
   too fine.  A full seven-visitor packet can repeat arbitrarily often at
   `ceil(f/g)=2` while its return is only diagonal sheet translation.  The
   corrected coordinate is the ordered packet path after both persistent
   stalks and zero reduced-holonomy loops have been quotiented.
4. THM-795 closes the entire arbitrary-height Hamming-one star around every
   shallow AP dilation.  Any residual shallow sporadic packet must use at least
   two replacement colours, so the next object is a sheet--tooth incidence
   graph rather than another one-coordinate lift bound.

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
| one-coordinate shallow lift | one replacement tooth against a complete missing-owner splice grid | two or more replacement colours sharing the sheet obligations |
| number of occupied sheets | capacity necessary conditions | which labelled owner occupies which sheet, especially at equality |
| measure of a folded obstruction | a sharp necessary cap | component-by-component containment and endpoint colours |
| safe-component count or raw winding | a coarse amount of boundary | component widths, reduced winding, and endpoint-owner incidence |
| a static prime-seven token state | exact cover at one phase | endpoint order, inverse steps, and global carry |
| a transitive next-event tournament | pairwise chronological order | the labelled Hamiltonian-path movie and simultaneous blocks |
| number of covered walls | finite-sample difficulty | a persistent exact seven-owner stalk that makes those walls redundant |
| number of active fastest periods | nonempty zero-sum visitor incidence | a full-support prefix-legal loop with zero reduced deck holonomy |

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
hypergraph.  THM-795 now extends one complete direction uniformly: a congruent
one-coordinate lift of any unit AP dilation is tight only when unchanged, at
every height.  Its proof uses the `c` missing-owner splice germs as vertices;
one replacement danger arc cannot contain their nontrivial deck grid.  Thus
the live higher-lift problem begins with at least two replacement colours that
can divide the splice obligations.  The next faithful object is their
bipartite sheet--tooth incidence graph with deck orders and one-sided germ
orientation, not an extrapolation from the height-twelve residue census.

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

THM-789 sharpens both the positive statement and the missing quantifier.  The
symmetric difference packet `A-A` doubles the safe-measure floor to
`2*72^(-10)` and improves the component-width floor to
`72^(-10)/(5 max(U))`.  Bohr erosion also turns depth into a pointwise tax on
each odd exception.  At the global level, if `g_i` are the cyclic gaps of the
deep set and `delta=2/(143 max(U))`, hypothetical tightness forces

```text
mu(E)+sum_i min(g_i,2delta) <= mu(H),
mu(E)+mu(R) <= mu(H),
mu(R) >= max(2delta,2*72^(-10)).
```

This is a real component-layout tax, but it does not yet choose the escaping
component.  The exact core

```text
U_0={1,2,3,5,7,8,9,10,11,12},   t_0=4/17,   (x,y)=(13,9)
```

traps the full natural return set, its symmetric packet, the Lipschitz
interval, and every literal refinement by the odd/folded coordinates inside
the diamond at `t_0`; the same core escapes globally at `14/19`.  Refining one
good anchor can therefore be monotone trapping.  The live theorem is not
“make the local packet finer,” but **select a deep component whose eroded
packet is not contained**.

Even that selector has a controlled-forgetting boundary.  On the same core,
the odd pairs `(13,9)` and `(17,13)` have identical raw margin order, identical
raw escape signs, identical transitive tournament fingerprints, identical
diamond measure `8/169`, and identical total eroded measure `212/5577`.
Nevertheless the return-thickened middle components escape only for `(13,9)`.
At `4/17`, the pairs `(13,9)` and `(43,13)` even share the unsigned odd-error
multiset, folded margin, parity opposition, and sharp determinant while their
signed local slopes make the same return interval behave differently.  The
component tournament is useful because its vertices preserve global
alternatives, but the faithful vertex must carry the signed affine tooth
address and eroded margin.  Order, signs, and scalar erosion are still only
telemetry.

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

The subsequent THM-786 audit isolates the next exact scale.  A fixed companion
serving `L` consecutive second-fastest-owner walls satisfies
`L<1+2gc/(f(g-c))`: compare only the first and last companion walls, so skipped
indices cause no loss.  Consecutive visitor sets satisfy a signed entrant-
leaver balance.  The density argument is now residue-sensitive: make the
possible balanced companion clusters into a hypergraph and fractionally hit
them with weights `lambda_c`.  If the speed-weight
`sum lambda_c c` is below `g`, grid counting gives a finite run bound.  This
strictly exceeds the ultra-sparse class: for `g=9,C={1,...,6}`, total companion
speed is 21 while the optimal transversal `{1,2,5}` has weight 8.  What remains
is not generic “co-landing.”  A second marginal looks between consecutive
`g`-walls: allowed `(extra f-wall, companion set)` packets form a finite
polytope, and lattice discrepancy bounds a run whenever the natural frequency
vector lies outside it.  This closes a dense `(f,g)=(65,64)` profile that the
transversal misses.  But `(f,g)=(69,29)` with companions
`(4,5,12,13,16,27)` lies inside both marginal polytopes and respects every
fixed-span constraint abstractly.  The missing state is therefore the common
centered-Beatty event order and carry coupling the `f`-period cluster word to
the `g`-period packet word, not either frequency vector separately.  The
originally reported factor-
one span and `sum c<f` threshold do not survive exact audit, and the stored
S304 script does not regenerate its `0.589` census.

THM-794 then stress-tests the repaired picture in the opposite direction.
For

```text
F=49H+1,   f=F,   g=F-7,   C={F-14,...,F-49},
```

every fastest period contains all seven slower owners, in the same order, and
`H-1` such periods are covered.  This refutes the universal extent conjecture
and makes active-period count unbounded even at `ceil(f/g)=2`.  Exact
cross-audit explains why neither new marginal theorem contradicts it.  The
balanced hypergraph has only the full edge `C`, with

```text
W_*=F-49,       1-W_*/g=42/g,
```

while the complete-`g` packet target has exact polytope distance

```text
dist_infinity(r,P)=49/(2g).
```

Both positive margins tend to zero, so their per-tuple bounds grow like `H`
and admit the construction.  More importantly, the repeated owner word

```text
(f,w_7,w_6,...,w_1)
```

is prefix-legal in the collision automaton and translates every token by the
same `-1`.  Its return is zero in `F_7^8/Delta`.  The two polytopes see only
the abelian occupation vector; they forget the ordered path and its return
map.  Thus the missing coupling is sharper than “use both marginals at once”:
it is a centered-Beatty packet path groupoid carrying reduced holonomy and a
metric/core-incidence lift.

The surviving problem is therefore a first-return or owner-switch problem:
after contracting intervals supported by a persistent exact stalk and legal
zero-holonomy packet loops, can the centered mechanical schedule keep agreeing
with the redundancy cocycle across every relevant part of the core-safe set?
THM-778 already settles one useful stratum: sufficiently dense equal-valuation
double walls pierce because two owners disappear simultaneously.

### 1.8 The even-maximum collar is now a clocked incidence problem

THM-792 starts from the exact remaining two-sheet max-peel collar rather than
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

The dynamic carrier now has an exact minimal form.  With initial excess
`e^0=d^0-1` and cumulative entry-minus-departure current `C`, every grouped
event is a sum of `A_12` roots and

```text
e=e^0+C,       coverage iff e>=0 coordinatewise.
```

The safe state space is therefore the 50,388 integer points of the seven-chip
simplex plus a dead state.  The energy is only its quadratic shadow,

```text
K=K_0+<e^0,C>+||C||^2/2.
```

All 14,184 initial `K=0` covers still tear.  Exact enumeration shortens the
bounded proof too: every initial cover tears by `3/8` (38 event groups), and
every divisor-complete cover by `4/11` (36 groups).  Once the labelled clocks
generate the root word, the edge multiset can be compiled out of predicate
testing; the root word itself cannot.  The uniform problem becomes an
intersection of two languages: a regular safe-current language on the chip
simplex and the much thinner arithmetic language of root words realizable by
divisor-complete rational clocks.

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

packet transport:
  ordered prefix-legal collision paths, reduced deck holonomy, and the
  metric translation of the centered-Beatty base cell;

fibre:
  labelled owner-to-sheet incidences, strict endpoint deletion,
  inverse steps, and global carry;

stalk decomposition:
  inclusion-minimal owner sets that already cover throughout an interval;

essential dynamics:
  changes of those minimal stalks and reduced packet-return classes,
  together with redundancy roots, flank owners, or overlap chips.
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

THM-794 shows that this is still not the terminal quotient if every retained
period is marked merely “active.”  For an ordered packet path `P`, also retain
its return class

```text
Hol(P)=[(-N_a(P) w_a^(-1))_a] in F_7^8/Delta.
```

Contract a maximal repetition only when the same prefix-legal collision path
returns by the diagonal class and remains in the same centered-Beatty base
cell type.  The result, with its nonzero metric translation still attached,
is the **holonomy-reduced packet skeleton**.  The metric clause prevents a
closed fibre loop from being mistaken for a closed orbit on the core-safe
base.

This definition deliberately makes the THM-779 family simple: its `2m` walls
are redundant over the fixed `W_0` stalk.  It also explains why the definition
must remain over the metric base.  A stalk may persist on one core-safe
component and fail on another even when the abstract token state and owner
word are isomorphic.

The proposed reduced skeleton is not yet proved to have bounded complexity.
It is a better theorem target because its complexity cannot be inflated either
by adding a redundant high-frequency owner to an exact stalk or by repeating
THM-794's diagonal full-support packet loop.

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
| active packet occurrences | chronology; switch to reduced deck-return class | visitor support and abelian owner frequencies | prefix legality, repeated zero holonomy, centered-Beatty cell, and metric translation |
| rooted redundancy states | collision-hop transitions; choice of duplicated root | the local `A_8` cocycle and supportability of an event word | absolute root, global carry, metric schedule, and the core-safe base; strong connectivity prevents a state-only exit proof |
| minimal cover stalks | persistence across a wall; before/after minimal bases | genuine owner switches after redundant walls are contracted | raw event multiplicity by design; component location must remain as a sidecar |
| proof obligations `(event, chip)` or `(tooth side, owner)` | transports or discharges the obligation | exactly the boundary incidence used by THM-792 | raw runner geometry unless the event clock and edge label remain attached |

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

The collar carrier nevertheless has one exact tournament-flow shadow.  If
`d_j` are its thirteen sheet degrees, then

```text
K=sum_j binom(d_j-1,2)
```

satisfies `Delta(-K)=d_departure-d_entry-1` at a simple edge slide, exactly the
endpoint-defect flux of THM-785; `8K` has THM-787's step-eight quantization.
This stratifies the automaton by overlap-chip concentration and suggests
measuring current across `K`-cuts.  It is telemetry, not a proof quotient:
the normalized `r=8` stalk has the same excess-degree pattern in all 20,160
states, and even the collar atlas has two cores with the same labelled degree
vector and `8K=80` but different first tears.  THM-791's Hamiltonian-path axis
is orthogonal in the metagraph, but the collar chronology itself is transitive
with one path.  Event labels, transition instances, and the metric base remain
the irreducible sidecar.

## 4. The exact remaining mathematical pieces

### 4.1 Uniform emptiness of the primitive tight twelve-speed sporadic branch

The branch is finite in principle, but the following uniform arguments are
still missing.

1. **Shallow full-residue rigidity from Hamming radius two onward.**  THM-770
   closes lift height twelve, and THM-795 closes every arbitrary-height
   Hamming-one star around an AP dilation.  What remains is a scale-free
   endpoint-owner coherence theorem when two or more replacement teeth share
   the missing-owner splice germs: prove a Hall defect, forced common divisor,
   or complete the finite global bound.
2. **The two-sheet folded branch.**  One must exclude every primitive
   divisor-complete quotient core, not only `max(U)<=19` or the bounded-height
   bank.  A proof may show that the folded containment `(FD)` cannot persist,
   or that THM-775's binary ownership tree cannot reach a terminal
   hereditarily primitive base.
3. **The even-maximum collar.**  THM-792 gives a bounded rational clock,
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
2. Define and bound the holonomy-reduced packet skeleton after persistent
   exact seven-stalks and repeated diagonal-return packet loops are contracted.
   Raw wall and active-period bounds are both impossible.
3. Use THM-783's period-sum, single-visitor, cluster-balance, and de-phasing
   laws to bound genuine transitions between reduced packet-return classes or
   their metric extent after both contractions.
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
THM-792 begins to do, and then force a finite automaton tear.

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
6. **Factor diagonal packet holonomy.**  Contract repeated prefix-legal packet
   loops whose return is zero in the reduced deck, while retaining their
   metric base translation.
7. **Study the essential skeleton.**  Apply the redundancy cocycle, overlap-chip
   transport, divisor shells, or endpoint-owner constraints only to genuine
   changes of minimal cover basis or reduced return class.
8. **Descend or tear.**  Either exhibit a free sheet/time, force a dyadic or gcd
   descent, or return a smaller normalized packet with all sidecars intact.

The key invariant is not “how complicated is the speed vector?”  It is “how
many genuinely different minimal cover stalks and reduced packet-return classes
can the measured core-safe base support under its exact mechanical schedule?”
That question survives the counterexamples that defeat bounded denominators,
raw component counts, raw wall counts, and active-period counts.

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
  -> ordered packet paths / reduced deck holonomy
  -> holonomy-reduced owner-switch skeleton
  -> lonely witness or a smaller residual packet.
```

The first arrow explains why height need not destroy safe mass.  The middle
arrows explain why boundary equality, labels, gcd, parity, and carry repeatedly
reappear after an overaggressive quotient.  The final arrow is the open one.

The sharpest current view is therefore not that LRC(14) is waiting for a
larger computation.  It is waiting for a theorem that couples a heavy measured
phase base to a finite labelled switch skeleton.  THM-780 supplies the former;
THM-778/779/783/786/788/790/791 and THM-792 supply complementary increasingly exact
versions of the latter.
The remaining work is to prove that the coupling cannot stay covered without
descending to one of the already controlled equality packets.
