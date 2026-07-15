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
  - THM-796
  - THM-797
  - THM-798
  - THM-799
  - THM-800
  - THM-801
  - THM-802
  - THM-803
  - THM-804
  - THM-806
  - THM-807
  - THM-808
  - THM-809
  - THM-810
  - THM-811
  - THM-812
  - HYP-6820
  - HYP-6835
  - HYP-6840
  - HYP-6865
---

# The missing state is heavy phase cells, owner-switch skeletons, and packet holonomy

## Honest status

The standard fourteen-runner case remains open.  The primitive tight
twelve-speed sporadic branch used by the current reduction is uniformly finite
in principle, but it has not been classified and is not proved empty.  This
session nevertheless changes the frontier in eight material ways:

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
4. THM-802 proves that THM-794 is one generator of a larger affine-pumpable
   language.  A phase-realizable prefix-legal loop with constant
   `d_a w_a^(-1)` modulo seven inflates to linearly many diagonal-return
   packets. More strongly, in every prescribed nonempty open five-core-safe
   interval it realizes each fundamental unequal one-fast class
   `(1,...,1,k)`, `k=2,...,6`. Multiplicity and reduced holonomy do not decide
   legality, so the unreduced insertion word and exact phase cell must remain.
5. THM-795 closes the entire arbitrary-height Hamming-one star around every
   shallow AP dilation.  Any residual shallow sporadic packet must use at least
   two replacement colours, so the next object is a sheet--tooth incidence
   graph rather than another one-coordinate lift bound.
6. THM-800 closes that two-colour stratum.  Oriented half-open splice-deck
   capacity forces both replacements to the common AP scale at exact
   tightness, after which every proper normalized double lift has the sharp
   floor `M>=2/25`. THM-804 proves the corresponding common-scale descent for
   three replacements, and THM-806 closes the scale-one triple chart by an
   exact collar/component argument. The whole Hamming-three star is therefore
   loose. THM-810 splits the first live radius-four chart into common scale and
   one quartic order-three coset packet that is already an `s=3` deep branch.
7. THM-797 turns odd exception divisors and their one-sided walls into global
   selectors.  At `q=13`, double-13 and full support are impossible; the sole
   survivor is the exact signed complement of `+/-y`, with a `2/13` core
   witness and sharp exception-speed caps.  A post-wall row silences every
   exception-divisor grid but escapes at denominator 17, proving that all deep
   components and core-generated denominators—not merely argmax or divisors—
   must remain. THM-803 then proves the half-grid ladder, full parity-twisted
   support, and an exact quadratic-size all-component selector. Its sharp row
   passes every grid and both maximizers but escapes at the nonmaximal
   component `7/22`; the selector is exact, while its uniform failure remains
   open.
8. THM-798/799 separate raw fragmentation from conservative state transport.
   Component count can grow transversely with no coherent divisor scale, while
   the two-number state `(safe mass, topological components)` still yields
   exact terminal cones for sufficiently lacunary far flags.

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
| raw positive-length component count | transverse tooth fragmentation | safe mass relative to the eventual peel frequency |
| only global deep maximizers | the sharpest quotient-core clearance values | escaping threshold components and their endpoint owners |
| exception-divisor and universal anti-grids | THM-797's exact shell plus THM-803's half-grid/parity support at `26,52,78` | non-grid component endpoints/cusps such as the sharp survivor's `7/22` |
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

THM-797 supplies the first uniform global selector on top of this corrected
state.  For every odd divisor `q` of an exception, it compares two exact
folded unit-class sets:

```text
D_q(U) = deep classes for the quotient core,
A_q(x,y) = exception-eligible, opposite-colour classes.
```

Any class in `D_q(U)\A_q(x,y)` is already a point of `E_U\H`, hence proves
the full erosion failure.  At the mandatory divisor `q=13`, hypothetical
containment forces the core's folded residue support to be either all six
classes or the five-class support missing exactly the off-divisor exception's
class.  This eliminates every support of size at most four, every misaligned
five-class row, and every non-full double-13 row without a height bound.

The one-sided boundary germs then recover the signs discarded by folding.  At
the inverse grid point of any rejected class, all minimum owners must occur
with both residue signs, or `G_U` leaks into the open diamond complement.  Ten
core runners cannot supply both signs in all six classes, or in five rejected
classes plus an accepted class.  Hence double-13 and full support are
impossible, and the sole survivor of the prime-grid signed-wall gate satisfies

```text
U mod 13=(Z/13Z)^* minus {+y,-y}
```

with every remaining signed residue exactly once.  At `p=y^(-1)` this gives
`phi_U(p/13)=2/13`; the pointwise tooth tax sharpens the exception bounds to
`x<=2B-1`, `y<=B-1`, and `13B+2xy<=2B(x+y)`.

The original aligned row

```text
U=(1,2,3,4,7,9,10,11,12,16),       (x,y)=(13,5),
```

both global maximizers `5/13,8/13` lie inside the diamond, while the threshold
endpoint `7/33` has `phi_U=1/11` and lies far outside it.  Thus “choose a
deepest point” is another lossy compression.  The faithful selector is a
bipartite incidence object between folded divisor obligations and **all**
owner-labelled closed deep components, with signed escape margin attached.
The signed walls now reject this row, but the stronger survivor

```text
U=(1,2,4,6,7,9,10,11,12,16),       (x,y)=(13,5)
```

passes the exact signed complement and every scalar tax.  Its only exception-
divisor grids are `q=5` (empty) and `q=13` (trapped), yet `t=6/17` has
`phi_U=2/17` and folded value `10/17`.  Thus exception-divisor grids themselves
are another lossy compression.  The faithful carrier joins arithmetic walls
to geometric component/tooth incidence and retains core-generated
denominators, instead of asking either quotient or a privileged divisor list
to finish alone.

THM-803 identifies the next exact layer. Both gate-sharp rows above are caught
by the quarter anti-grid `11/52`. Every residual must have full
parity-twisted support and silence on the universal `26,52,78` ladder. More
importantly, with `K_U=E_U+closure(R_U)`, the entire erosion predicate is
equivalent to its values at owner-labelled component endpoints and folded
cusps, at most `200B^2+22B` points. The sharp replacement row

```text
U=(2,4,6,7,9,10,11,12,14,16),       (x,y)=(13,5)
```

passes every grid and both global maximizers but escapes on the nonmaximal
singleton `7/22` (or, after return thickening, at `365/1144`). Thus the exact
all-component selector is no longer missing; the open step is a uniform
negative-margin or incompatibility theorem on its scale-dependent obligations.

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

THM-802 shows that this is a language phenomenon rather than a single
once-per-owner family.  Its affine lifting lemma pumps every phase-realizable,
prefix-legal loop with `d_a w_a^(-1)` constant modulo seven.  In every chosen
nonempty open five-core-safe interval, an analytic density construction
realizes all fundamental unequal one-fast classes `(1,...,1,k)`,
`k=2,...,6`; the offsets and onset scale are existential, not effectively
bounded. The explicit `k=2` word still supplies the quantitative family. Yet
only three of the `181,440` words with that multiplicity are legal from the
displayed state, and only one of the eight repeated-owner insertion positions
sharing the same first-owner order works.  Thus owner multiplicity, reduced
return, and collision-SCC membership together remain telemetry.  The exact
carrier is the centered-mechanical owner word with prefix state, phase-cell
wall coordinates, metric translation, and core incidence.  Here “noncentral”
describes the word shape; the holonomy itself is central and zero after the
diagonal quotient.

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

This is a reduction before quotient height is bounded, not by itself a uniform
exclusion theorem.  In particular, an integer-grid
edge event should not be silently identified with a midpoint wall clock; the
appropriate mechanical word and phase convention must be carried explicitly
when this automaton is joined to THM-778.

The first bounded automaton tear is also exact and uniform in the multiplier.
For `w=13c` and `U subset [1,24]`, if `c>=5` is odd then immediately right of
`2/(13c)` every edge is `{0,-u^(-1)}` and ten edges meet at most eleven sheets.
For `c=1`, the event word on `2/13<s<11/13` eliminates every ten-core.  Of the
`binom(23,10)=1,144,066` possible cores, `101,850` cover the initial chamber
and zero cover the whole 117-group event word.  The primitive
divisor-complete subatlas likewise has `20,604` initial static covers and zero
survivors.  At `c=3`, all 2,528 initial covers tear by `1/7`.  Thus every odd
`c` is excluded at this quotient height.  This is uniform in `c` but not in
`max(U)`.  Its conceptual value is that static edge
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

There is a smaller exact compositional state between the full word and the
energy.  For a grouped root word `W`, let `c_W` be its endpoint current and
`b_W` its coordinatewise minimum over all prefix currents.  Then

```text
T(UV)=(c_U+c_V,min(b_U,c_U+b_V)),
e^0 survives W iff e^0+b_W>=0.
```

Equivalently, if `J_s` is the oriented prefix edge-chain on `K_13`, then
`C_s=partial J_s` and
`<p,C_s>=sum_(r<=s)(p(entry_r)-p(departure_r))`.  The singleton nodal
potentials decide current coverage, while the ordered cycle lift still decides
arithmetic realizability.  This is the exact cut/current bridge suggested by
the Smith-diagram thread; it is not an identification with HYP-6865's harmonic
tournament potential.

The existential deficit `D(W)=-sum b_W` tells whether **some** seven-chip
allocation survives.  In the full height-24 atlas every initially covered
word has `10<=D<=20`; none of the 280 endpoint-zero words is abstractly safe.
But at the actual first tear only 294 of 101,850 words have `D>=8`, so most
tears depend on the specific initial chip allocation.  Thus `(e^0,T(W))` is a
predicate-exact block transfer; `D` alone is another one-sided shadow.

### 1.9 Transverse far flags need state transport, not a fragmentation bound

THM-798 closes another historical compression.  The four-far family

```text
P_N={1,...,9,15,110,N},       V_N=P_N union {1092N}
```

has no nontrivial divisor packet of size seven, yet the number of positive-
length good components grows at least as `N/1540-8/7`.  Raw fragmentation is
therefore independent of the coherent scale coordinate even in the first open
four-far chart.  This is not a hard family: its proportional top peel closes
uniformly.

THM-799 identifies the conservative state that explains that closure.  If a
base has safe mass `mu` and full topological component count `r`, adjoining a
frequency `N` gives

```text
mu_N >= 6mu/7-2r/(7N),       r_N<=r+N.                   (GS)
```

The state composes.  An adjacent-exchange identity proves that inserting
distinct frequencies in increasing order uniquely maximizes this affine lower
certificate, so the `k!` order gauges collapse to one canonical path per peel.
Combined with the heavy-cell floor and THM-755, (GS) proves terminal cones for
every fully lacunary far-count flag; in the literal four-far chart the uniform
factor is `412`, improved to `19` by the exact 2,002-core floor.

This does not make `(mu,r)` a lossless LRC quotient.  It deliberately forgets
endpoint owners and correlations, and may give a nonpositive lower bound on
comparable-scale faces.  Its lesson is subtler: different proof stages need
different controlled resolutions.  The scalar pair `(mu,r)` suffices for a
one-sided terminal cap, while folded containment, collar current, and packet
holonomy require the richer incidence stalk.  The unresolved infinity is now
concentrated on clustered comparable-scale faces, not fully lacunary ones.

## 2. A candidate underlying object

The most economical object consistent with all of the proved corrections is a
**measured, stalked event system** over the quotient circle.  For a core `P`
and an exception deck `W`, its fields are

```text
base:
  the strict core-safe set G_P, with component endpoints and measure;

terminal-cap shadow:
  the conservative pair (safe mass, full component count), used only when
  transverse insertion plus a proportional peel makes its inequalities fire;

phase substrate:
  a heavy joint-phase cell or, in a limit, the character relation lattice;

arithmetic selectors:
  odd-divisor and half-grid deep classes, signed/parity-twisted walls,
  exception acceptance shells, and THM-803's owner-labelled endpoint/cusp
  selector on K_U=E_U+closure(R_U);

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

### 2.1 Quotients, certificate shadows, and telemetry

The accumulated counterexamples support a useful formal distinction.  Let
`X` be a space of configurations or event histories, `P:X->{0,1}` the proof
predicate, `T_l` the transition with label `l`, and `pi:X->Y` a compression.

1. `pi` is **predicate-exact** when `P` is constant on every fibre of `pi`.
2. It is **transition-exact** (strongly lumpable for the labelled dynamics)
   when `pi(T_l x)` is determined by `(pi(x),l)` whenever the transition is
   defined.
3. It is a **certificate shadow** when a condition on `Y` proves one side of
   `P` uniformly over its fibre, even though neither exactness property need
   hold.
4. It is **telemetry** when it supplies an order, energy, or diagnostic but no
   fibrewise implication for `P`.

A recursive proof state needs the first two properties, possibly after a
stalk or phase sidecar is restored.  Static separation is not enough.  THM-796
gives an exact tournament-metagraph warning: its normalized face row separates
all 272 merged nodes at `n=7`, yet 1,206 of the 1,312 nonzero parent/child
blocks have nonconstant lift counts.  The aggregate node address identifies
the parent but is not a Markov state.  The exact repair is the joint
node/line/face coupling tensor together with the hidden endpoint-phase bit.

THM-801 adds the missing gluing lesson.  Its three staircase faces form an
exact Čech descent: compatible face tilings glue uniquely, while bare
complement lines require a relative phase cocycle on the overlaps.  At `n=5`
pairwise-compatible triples split evenly into zero- and nonzero-holonomy
classes because the triple overlap is empty.  The expanded `Omega+S2` address
is injective through `n=7`, but the theorem explicitly does not promote that
finite separator to an all-size Markov state.  For LRC this is a design rule,
not a transfer theorem: component, divisor, deletion, and packet charts may be
locally exact yet fail to glue unless their shared owner/phase data agrees on
overlaps.

THM-809 now extends the finite separation result to every one of the
`2^20` complement lines at `n=8`: a lower-first recursive address becomes
injective by the `tau=7` layer. THM-808 simultaneously gives the needed
dynamic warning. A centered-Christoffel owner-count block acts affinely on the
prime-sheet redundancy root, but the same literal mask and same block can have
different target masks; distinct source roots repair the displayed collision.
Static reconstruction, transported root, and continuation completeness are
therefore three different levels. The LRC carrier needs the fibre product of
the face/mask address with the owner-labelled sheet root, event word, carry,
and metric component—not a concatenation of independent tournament scores.

THM-811 locates Möbius/Smith curvature one layer lower in this ledger.  Its
factorized polynomial is an exact all-size law for flux, linear current, and
positional curvature, but `(C3,H_x)` is only node-injective through `n=7` and
curvature endpoint data reaches only `7248/8064` black reflection orbits.
Orbit-symmetrized `(B2,B3)` reaches all `8064`; the sixteen residual literal
collisions are precisely mirror pairs.  Thus curvature is faithful telemetry
on a reflection orbit, `(B2,B3)` is the stronger positional orbit address, and
literal continuation still needs mirror/path orientation.  Neither carries
the LRC metric witness, owners, threshold side, wall schedule, root, or carry.

THM-812 now proves the first direct centered-Christoffel action on the
metagraph stalk.  Its `(1,2)` coordinate replication embeds `X_5` in `X_6`,
commutes with complement/reflection, and transports all twenty projected
coloured edge cells injectively even though bare nodes spread.  Boolean
coefficients transform by subset-image pushforward; the unique fixed-core
collision surviving degree at most three is first separated by three explicit quartic
cross-leg/apex coefficients.  This is a genuine path/core transport theorem,
not an LRC transport theorem.  The missing LRC state is therefore even more
specific: couple THM-812's path/core action to THM-808's owner/root action and
retain metric component, threshold side, wall schedule, gcd, and carry.

THM-813 then separates action from observation.  Reflection-orbit line cells
`Q_n=E_n/<sigma>` transport injectively under every equivariant embedding;
colour plus unordered endpoint nodes is only an observation
`pi_n:Q_n->P_n`.  The first centered step works because `pi_5` is bijective.
At `n=6`, `pi_6` merges 52 independent obligations, and the next step splits
51 of those fibres.  A recursive state must therefore retain `Q`, not merely
the projected coloured edge cell.  Likewise, a low source-degree Möbius
sector pulls back from the `rho`-saturated target packet, which may contain
high target degrees.  Static edge identity and bounded-degree truncation are
both observations of the transported stalk, not the stalk itself.

The LRC evidence falls cleanly into the same ledger:

| level | present examples | legitimate use |
|---|---|---|
| predicate- and transition-exact after labels | collar `(e^0,root word,cursor)` or block transfer `(e^0,T(W))`; owner-labelled deep components with signed tooth address; ordered packet path with reduced holonomy and metric base | recursive tear, escape, or descent |
| certificate shadow | `(mu,r_top)` transverse transport; `W_*<g` balanced-transversal bound; positive packet-polytope distance; global erosion budgets | close a terminal cone or one tuple without reconstructing every event |
| telemetry | collision energy `K`; raw component/wall/period counts; argmax-only selectors; black Möbius/Smith current-curvature; bare tournament fingerprints | choose a gauge, locate a bottleneck, or falsify an overcompressed conjecture |

THM-799 makes the middle distinction exact.  Two bases can have the same safe
mass, component count, and entire component-length multiset, yet adjoining the
same frequency `59` changes their new safe masses by `1/826`.  Thus `(mu,r)`
is a valid one-sided transport shadow and provably not an exact transition
state.

This also clarifies the honest import of the Smith-diagram tournament thread.
The S307 metagraph entry currently sharing the colliding ID `HYP-6865` has
exact voltage/score concordance through `n=6`. Its numerically refined `n=7`
solve finds 132 well-separated discordant adjacent-level pairs, and the `n=8`
audit finds 65,921 discordances reaching four score levels: coarse pairwise
concordance is decisively refuted numerically, though not promoted here to an
exact theorem. Level-mean and level-median monotonicity are verified through
`n=7` (the `n=8` check is pending), so the score axis is at most a mean-field
harmonic coordinate. The separate `E_n` network is
rationally exact through `n=7`, and already at `n=4` the proposed electrical
reciprocity fails: `(3/7)(2/5)=6/35`, so the cut/cycle relation is algebraic
rather than a planar Smith duality. THM-796 meanwhile proves that actual current and
recursive continuation live on weighted line fibres.  This is a rigorous
potential-versus-current separation inside the tournament model, not a direct
LRC theorem.  In the collar the analogous exact current is the full `A_12`
cut-current vector; its quadratic energy is only a potential-like shadow.

The next computational object for the `r=8` route can now be stated more
sharply than “couple the marginals.”  Form a finite joint incidence tensor

```text
Xi_LRC(core-safe component, initial stalk;
       f-packet type, g-packet type;
       ordered event block, reduced holonomy;
       recursive face/overlap phase address).
```

Its packet-type marginals recover the two THM-786 frequency polytopes.  Its
conditional fibres test whether the next reduced stalk/holonomy class, or the
tear predicate, is constant.  If not, the differing event block is a literal
non-lumpability witness; if yes, that tensor cell is a valid state quotient.
This turns the missing common Beatty order/carry field into a finite strong-
lumpability and Čech-gluing question, modelled on THM-796's exact coupling and
THM-801's overlap-phase descent rather than on another scalar separator.

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
period is marked merely “active,” and THM-802 extends the obstruction to an
affine isotropy family containing all five fundamental unequal one-fast
classes in every open core-safe interval. For an ordered packet path `P`, also retain
its return class

```text
Hol(P)=[(-N_a(P) w_a^(-1))_a] in F_7^8/Delta.
```

Contract a maximal repetition only when the same exact centered-mechanical
owner word is prefix-legal, returns by the diagonal class, and remains in the
same decorated metric/core phase-cell type.  This includes all affine
phase-cell diagonal-return loops, not only THM-794's once-per-owner generator.
The result, with its nonzero metric translation still attached, is the
**holonomy-reduced packet skeleton**.  The metric clause prevents a closed
fibre loop from being mistaken for a closed orbit on the core-safe base.

This definition deliberately makes the THM-779 family simple: its `2m` walls
are redundant over the fixed `W_0` stalk.  It also explains why the definition
must remain over the metric base.  A stalk may persist on one core-safe
component and fail on another even when the abstract token state and owner
word are isomorphic.

The proposed reduced skeleton is not yet proved to have bounded complexity.
It is a better theorem target because its complexity cannot be inflated either
by adding a redundant high-frequency owner to an exact stalk or by repeating
any recognized affine prefix-legal diagonal-return loop.  A quotient keyed
only by multiplicity or bare collision SCC is still too coarse by THM-802.

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
| folded/signed/anti-grid classes | raw and parity-twisted multiplicity; inversion as gauge | THM-797's signed gate and THM-803's universal `26,52,78` ladder | non-grid endpoints/cusps and return-thickened component margins |
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
- SCCs of a predicate- and transition-exact decorated transducer distinguish
  trapped from escaping regions; bare collision SCCs are only telemetry;
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

1. **Shallow full-residue rigidity at Hamming radius four and beyond.** THM-770
   closes lift height twelve, THM-795 closes every arbitrary-height
   Hamming-one star around an AP dilation, and THM-800 closes the full
   Hamming-two star.  Its proof first forces both replacement deck orders to
   one at exact tightness, then gives the sharp normalized floor `2/25`.
   THM-804 proves the three-replacement common-scale descent by an exact
   half-open capacity/residue-ratio argument. THM-806 then closes every proper
   scale-one triple lift: a uniform collar forces one replacement into
   `[14,24]`, two-comb geometry sharpens the remaining box to `v<=262`,
   `w<=12v`, and a larger 5,713,539-row superset replay has zero tight rows.
   Thus radius at least four is the first shallow residual. THM-810 proves its
   first deck classification: scalar capacity leaves only common scale or the
   all-order-three labels in a coset of `<5>={1,5,8,12}`. Exact sheet overlap
   leaves four parity patterns, and division by the common gcd turns the latter
   into an `s=3` deep packet. The scale-one quadruple chart and arbitrary lifts
   of the coset interface remain; four least-CRT base rows are loose but do not
   give uniform radius-four closure.  Every coset has eight lift-invariant
   `q=39` equality points pinned by an opposite pair of core speeds.  Since
   each is a strict local cusp, the missing field is not lift height on that
   clock but the identity and exception-comb incidence of a different
   core-safe component.
2. **The two-sheet folded branch.**  One must exclude every primitive
   divisor-complete quotient core, not only `max(U)<=19` or the bounded-height
   bank.  THM-797's q=13 signed walls leave only the exact signed complement of
   `+/-y`, with `M(U)>=2/13` and strong exception-speed caps.  On that survivor,
   THM-803 further forces full parity-twisted support and the complete
   `26,52,78` anti-grid ladder, then decides the full erosion predicate on an
   exact owner-labelled endpoint/cusp set of size at most `200B^2+22B`. Prove
   that some selector point escapes uniformly, or show that the joint selector
   obligations are incompatible with THM-775's terminal ownership tree.
   THM-807 proves a linear central-return selector of size at most `42B-2`,
   exact in the connected-return branch. It isolates disconnected return
   satellites as the source of the quadratic component interaction rather
   than removing the uniform obligation. Its exact method-boundary rows also
   silence every multiplier grid through `d=7`, or every even grid through
   `d=18`, before escaping through the component selector.
3. **The even-maximum collar.**  THM-792 gives a bounded rational clock,
   repeated ordered flank types, and—on the forced `13`-multiple subbranch—a
   moving edge-cover automaton.  It excludes every `w=13c` through quotient
   height `24`: `c=1,3` tear exactly, while `c>=5` fails at the first boundary.
   What is missing is a uniform-in-quotient-height tear or congruence
   obstruction linking those three structures beyond the bounded atlas.
   The sharp pumpability test is to exclude arithmetic-realizable prefix-safe
   tropical return blocks up to `Z/13` sheet rotation; a literal zero-current
   word is only the simplest case.  Retain `e^0`, exact root insertion order,
   phase coset, and ineligibility-window incidence in that search.
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
2. Identify the affine diagonal-isotropy subgroupoid after persistent exact
   seven-stalks, and quotient every decorated prefix-legal diagonal-return
   loop—including all fundamental unequal one-fast classes—while retaining
   its unreduced word, phase-cell germ, metric translation, and continuation
   labels. Raw wall, active-period, multiplicity, and bare-SCC bounds are all
   insufficient; arbitrary open core-safe localization cannot force exit.
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
large common-factor core.  THM-798/799 remove raw fragmentation and every
fully lacunary far-count flag from this list; unresolved infinity is on
clustered comparable-scale faces.  The proved bounded, state-transport, and
sheet regimes are not yet an assembly theorem for LRC(14).

## 5. A recursive research program suggested by the synthesis

The history suggests a disciplined recursion rather than a new scalar.

1. **Find mass.**  Use a deeper lower-dimensional witness and the heavy-cell
   argument to obtain a positive measured base.
2. **Try the conservative terminal cap.**  Transport `(safe mass,component
   count)` across genuinely transverse frequencies and peel proportionally;
   if the inequality does not fire, restore the richer endpoint state.
3. **Fire arithmetic selectors.**  Test odd-divisor deep grids and their
   signed one-sided walls against the exception shells; retain every support
   class and core-generated denominator that survives.
4. **Expose boundary incidence.**  Decompose the whole remaining superlevel set into exact components or
   top teeth and retain ordered endpoint owners.
5. **Choose the natural fibre.**  According to the binding scale, use `s`
   sheets, a prime-field token deck, a dyadic sheet tree, or the `Z/13` Cayley
   edge cover.
6. **Decode chronology.**  Use centered/global ranks with gcd, parity phase,
   simultaneous blocks, and metric positions retained.
7. **Factor persistent stalks.**  Remove event multiplicity already absorbed by
   a fixed exact cover.
8. **Factor diagonal packet holonomy.**  Contract repeated prefix-legal packet
   loops whose return is zero in the reduced deck, while retaining their
   metric base translation.
9. **Study the essential skeleton.**  Apply the redundancy cocycle, overlap-chip
   transport, divisor shells, or endpoint-owner constraints only to genuine
   changes of minimal cover basis or reduced return class.
10. **Descend or tear.**  Either exhibit a free sheet/time, force a dyadic or gcd
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
