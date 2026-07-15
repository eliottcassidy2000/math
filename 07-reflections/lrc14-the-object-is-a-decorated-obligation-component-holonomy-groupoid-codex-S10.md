# LRC14: the object is a decorated component-obligation holonomy groupoid

**Session:** codex-2026-07-14-S10
**Status:** holistic frontier synthesis and proposed proof interface; not a proof
of LRC14 or of the uniform twelve-runner sporadic branch.

## 0. Two corrections before any synthesis

The historical conjecture that every covering thirteen-speed row has a good
period with denominator `q<=25` is **false**.  THM-762 and THM-764 give exact
counterexamples and the corrected signed-pair deck criterion.  In particular,

```text
26*{1,...,12} union {339}
```

has its first rational witness at `2/27`, not at a denominator at most `25`.
The failure is structural: a bounded denominator deck is an observer of one
residue fibre, not a scale-uniform model of the row.

The statement that the `n=12` sporadic branch is empty **uniformly is still
open in this repository**.  THM-795 removes the entire Hamming-one star around
every shallow arithmetic-progression dilation, at every lift height.  THM-797
proves a general odd-divisor grid obstruction and a strong `q=13` support
gate in the two-sheet branch.  Neither theorem closes Hamming radius at least
two, every binding scale, or every owner-labelled deep component.  A finite
census, a shallow rigidity theorem, and a uniform branch theorem are three
different statements.

Those corrections do not merely narrow the frontier.  They identify the
kind of object the frontier requires.

## 1. The recurring lesson of the repository

The repository has repeatedly found a useful quotient, found two packets in
one quotient fibre with different lonely-runner behaviour, and then promoted
the missing coordinate to a sidecar.  The history is therefore cumulative:

```text
residue             needed scale;
scale               needed sheets and ramification;
sheet counts        needed owners and overlaps;
owners              needed endpoint side and exact margin;
one component       needed all components;
component lists     needed placement relative to the next tooth grid;
event counts        needed order and carry;
order               needed return holonomy;
reduced holonomy    needed metric translation and core incidence;
tournament class    needed its observer/Hamiltonian-path fibre.
```

There is no evidence for a single scalar master invariant.  There is strong
evidence for a small *typed state* whose fields are retained or discharged
according to the next proof operation.

### A preservation ledger

| Viewpoint | Natural carrier or vertices | What it preserves | What it destroys when used alone | Lesson at the current frontier |
|---|---|---|---|---|
| Torus orbit and danger-arc circle | phase line, danger arcs, safe cells | the exact threshold predicate and simultaneous phases | arithmetic provenance and recursive scale | This is the faithful total space, but not a useful finite quotient by itself. |
| Farey cells, pair-sum rulers, continued fractions, centered Beatty words | rational cusps and endpoint events | exact candidate maxima, event order, carry, and local clock changes | simultaneous owner obligations if only denominators are retained | The denominator is an event address; the common event word, not its marginal counts, is the object. |
| Divisor covering and cheap rational periods | moduli `2,...,14` and missing divisors | explicit witnesses and the covering reduction | geometry inside the covering class | Covering is an arithmetic gate, not a geometric classification. |
| Small-period residue decks | zero bit plus signed unit-pair classes modulo `q` | exact existence of an `a/q` witness in the proved range | coherent scale, cross-modulus transport, component geometry | THM-762/764 replace the false `q<=25` slogan by a correct deck predicate. |
| Scale, sheets, and the binding time `p/(13s)` | quotient core, sheet cycle `Z/sZ`, exceptions | scale covariance, lift multiplicity, gcd ramification, shallow versus deep binding | persistent incidence if sheets are reduced to a capacity count | A prime residue fibre is only the height-one chart of a ramified cover. |
| Component-tooth containment | safe/deep components and danger teeth | the exact inequality `||wc||+wh<=alpha`, endpoint side, and escape margin | global alternatives if only one selected component is inspected | The quantified object is all owner-labelled components, not a preferred maximizer. |
| Folded erosion and the two-sheet diamond | folded deep components, return set, erosion signs | two-sheet symmetry, signed return geometry, global gap budgets | residue provenance, endpoint owners, and other divisor grids if aggregated | A local trap is not a global trap; signed component data must be transported. |
| Prime-seven token events and collision cocycles | wall events, token sheets, packet words | local legality, inverse residues, event order, carry, and deck holonomy | exact metric displacement and core incidence if normalized too soon | THM-794 shows that a legal full-support loop can repeat indefinitely at fixed fastest ratio. |
| Gcd, parity, and dyadic descent | deletion cores, effective sheet orders, colour classes | forced common factors, ramification, parity transport, a smaller quotient problem | primitive-core geometry and odd-prime component obligations | Descent is earned only after a simultaneous sheet obligation forces the divisor. |
| Safe measure, Fourier/Bernoulli, and capped envelopes | mass, Fourier modes, endpoint pairs, peel rate | strict-margin mass, correlation, exact fixed-core tails, and sound cap inequalities | isolated equality witnesses; owner placement if only coefficients or mass are kept | THM-799 gives a powerful enclosure state, not a lossless insertion state. |
| Safe-mass/component-count transverse insertion | `(mu,r,N,a)` plus a named peel | a composable lower mass bound and upper component load | phase placement relative to the new `N`-tooth grid, even if all component lengths are known | Exact insertion still requires owner-coloured endpoint/grid incidence. |
| Projective packets, affine suspensions, and cluster trees | normalized clusters, slopes, offsets, bounded-ratio children | scale-normal recurrence and coherent families | residue, owner, and threshold data unless decorated | Fully lacunary flags are terminal by THM-799; residual infinity is a recursive bounded-ratio cluster tree. |
| Safe peels and deletion | a row, a chosen complement, and a witness | exact reduction of `M` when the peel is safe | canonicity and completeness of the reduction tree | A peel is a morphism with a certificate, not an automatic normalization. |
| Finite exact and formal certificates | rational interval complexes, CSP cells, theorem lemmas | exact counterexamples, complete bounded atlases, replayable witnesses | uniformity outside the quantified parameter set | Computation should enumerate obligations and falsify quotients; a uniform theorem needs a transport or descent argument. |
| Tournaments, tilings, and metagraphs | pairwise orientations, Hamiltonian-path fibres, transition nodes | cyclic order, switches, SCCs, converse symmetry, observer multiplicity | lengths, signs, owners, equality flags, and continuation data | A tournament is a diagnostic two-skeleton or chart label.  Its metric and incidence stalk carries the theorem. |

Several old apparent disagreements disappear in this ledger.  Exact circle
geometry and Fourier analysis are not competitors: one is the primal endpoint
complex and the other is a weighted projection of it.  Residue decks and sheet
covers are not competitors: a residue deck is a fibre of the cover at a chosen
observer.  Tournaments and component words are not competitors: a tournament
records selected pair comparisons inside a richer labelled transition object.

## 2. Proposed underlying object

The smallest useful organizing object suggested by THM-794, THM-795,
THM-797, THM-798, and THM-799 is a **decorated
component-obligation holonomy groupoid**.  This is a proposed proof interface,
not a completed theorem or a claim of categorical minimality.

### 2.1 Objects

At a fixed threshold `alpha` an object is a packet

```text
X=(P, C, O, I; A, M, W, K),
```

with the following typed fields.

- `P` is a normalized row or quotient core, together with the chosen scale,
  deletion, fold, or affine chart.  Normalization never erases the data needed
  to reconstruct covering and primitivity.
- `C` is a connected safe/deep component, a closed component endpoint, or a
  one-sided endpoint germ.  When the proof is global, `C` ranges over the
  complete owner-labelled component family rather than a selected maximizer.
- `O` is the finite set or hypergraph of active obligations over `C`: sheets
  that must be burned, residue classes that must be supported, teeth that must
  contain the component, divisor-grid points that must be trapped, packet
  prefixes that must remain legal, or proof obligations that must be
  discharged.
- `I` is the labelled incidence relation between components and obligations.
  It records which owner/tooth/grid/event meets which component or germ.
- `A` is the arithmetic decoration: exact residues, scale, gcds, effective
  orders, inverse steps, parity or colour, ramification, divisor support, and
  the binding address such as `p/(13s)`.
- `M` is the metric decoration: rational endpoints, widths, exact escape or
  erosion margins, strict/closed flags, slopes, safe mass, component count,
  and the named peel rate where relevant.
- `W` is the transport decoration: centered Beatty or Christoffel word,
  simultaneous tie blocks, owner order, sheet/token assignment, local carry,
  marked observer cut, and return word.
- `K` is the proof decoration: the predicate being preserved, the certificate
  already obtained, the next legal operation, and every residual debt created
  by a quotient.

The vocabulary “component-obligation” is deliberate.  Sometimes the natural
vertices are safe components and inserted teeth (THM-798/799).  Sometimes they
are splice sheets and replacement teeth (THM-795).  Sometimes they are odd
divisor classes and owner-labelled deep components (THM-797).  Sometimes they
are collision events and packet prefixes (THM-794).  Runners are labels on
these incidences, not universally the vertices.

### 2.2 Reversible arrows and holonomy

Within one exact stratum, arrows are reversible transports that preserve the
decorations required by `K`.  Examples are:

- circle translation or reflection with its orientation bit;
- unit dilation and residue relabelling;
- change of marked observer cut or Hamiltonian-path gauge;
- sheet lift and its inverse after the sheet index is retained;
- wall transport through a chamber with the whole tie block retained;
- folded reflection between the two sheets;
- component return after one exact event packet.

A loop has at least two holonomy coordinates:

```text
h_red    = combinatorial/sheet/token return modulo declared gauge,
Delta_M  = metric translation of the base component and endpoint phases.
```

They must not be conflated.  In the THM-794 family, one full packet has raw
token map

```text
k -> k-(1,1,...,1),
```

so its reduced holonomy in `F_7^8/Delta` is zero.  Nevertheless the metric
base advances by `1/F`.  The packet is a closed loop in the normalized
collision automaton and a nonclosed translation in the metric skew product.
Moreover, THM-786 §3g shows that both exact marginal margins are
`Theta(1/H)` along this family.  Thus “positive for every fixed tuple” does
not imply a uniform compactness gap, and “zero reduced holonomy” does not mean
“nothing happened.”

### 2.3 Descent is not an invertible arrow

Gcd quotient, peeling, contraction of empty refinement, forgetting a metric
sidecar after a dual bound, and branching from one component to all components
are generally not invertible.  It would be inaccurate to put them into a
groupoid as if they were symmetries.

The intended structure is therefore a **groupoid-valued recursion**:

```text
each exact stratum       -> a decorated holonomy groupoid;
between exact strata     -> certified descent/branching functors;
forgotten information    -> an annihilation proof or named residual debt.
```

The full proof machine is a category whose reversible core is the holonomy
groupoid and whose condensation graph is well founded only after a genuine
complexity decrease is proved.  This separation is important.  THM-795 earns
a divisor descent.  THM-797 earns a branch to another component.  THM-794
earns only a gauge quotient of the diagonal loop; it does not by itself earn a
smaller row.

### 2.4 The observation functor

Every proposed quotient must name its observation.  For the local geometric
branches the basic observation is the sign vector

```text
Obs(X)=(escape/containment margin for every component-obligation incidence).
```

For a global row it also records whether at least one closed safe point
survives.  A quotient is theorem-safe only if one of the following is proved:

1. `Obs` is constant on every quotient fibre;
2. the forgotten coordinate is reconstructed by retained transport;
3. a metric, Fourier, Hall, or dual certificate annihilates its possible
   effect;
4. every legal continuation has the same terminal witness outcome; or
5. the ambiguous fibre is emitted as a named residual obligation.

The continuation clause is the dynamic version of the metagraph “stalk”
lesson.  Two states that agree now but respond differently to a future tooth
grid, wall crossing, peel, or divisor lift are not the same proof state.

## 3. THM-794, THM-795, and THM-797 are one recursion in three regimes

The three newest structural theorems look different only after their common
interface has been projected away.

```text
normalize a core and threshold
        |
enumerate every relevant component or one-sided germ
        |
build its finite obligation fibre and owner-labelled incidence
        |
transport obligations through one exact event/deck/grid word
        |
compute return holonomy and every metric escape margin
        |
        +-- an obligation misses a component --> witness / erosion escape
        |
        +-- simultaneous coverage forces a divisor --> descend and recurse
        |
        +-- a legal loop is pure gauge --> quotient the loop, retain Delta_M
        |
        +-- several components remain --> branch over all of them
```

### 3.1 THM-794: legal recurrence without descent

For `F=49H+1` and `w_j=F-7j`, every one of the seven visitors occurs in every
fastest period.  The packet word is exact, every inverse step is one, and the
same normalized collision state returns after eight events.  Active-period
count, raw wall count, and universal extent are therefore not decreasing
complexities.  The correct action is:

1. recognize the full-support diagonal packet as holonomy;
2. quotient its diagonal token translation;
3. retain the nonzero metric translation and the core-safe component it is
   moving through; and
4. study transitions among normalized collision SCCs only after that quotient.

This is also the exact stress test for the two marginal reductions in
THM-786.  Both remain valid per tuple, but their margins decay like `1/H`.
An argument that only sees the cluster occupation vector or packet polytope
has abelianized away prefix legality, return state, and metric residence.

### 3.2 THM-795: simultaneous sheet obligations force descent

For a one-coordinate lift of the shallow full-residue packet, component-tooth
geometry supplies a strict safe germ.  Lift it through the scale sheets.  If
the replacement is dangerous on every germ, the sheet-deck order must be one:

```text
c/gcd(c,w)=1.
```

Thus the simultaneous obligation does not merely count coverage.  It forces
`c|w`, permits quotient by `c`, and returns to the scale-one statement, where
the displaced residue is loose.  This is the ideal descent pattern:

```text
all sheets blocked
  -> trivial effective monodromy
  -> common divisor
  -> smaller normalized packet
  -> contradiction from the base theorem.
```

The remaining shallow frontier begins at Hamming radius at least two.  Its
natural carrier is not a tournament on speeds but the bipartite incidence of
missing-owner splice sheets with replacement teeth, decorated by deck order,
one-sided germ orientation, exact margin, and owner.  A Hall, Smith-normal-form,
or common-divisor lemma on this carrier would generalize the successful
Hamming-one descent.

### 3.3 THM-797: a silent grid means “look elsewhere,” not “trapped”

In the two-sheet branch, every odd divisor `q` gives a finite set of folded
deep obligations.  The exception shell says exactly which of them the two
exceptions can trap.  If one deep class falls outside that shell, global
erosion fails immediately.  At `q=13` this forces the folded support to be all
six classes, except for one sharply aligned five-class possibility.

The aligned survivor is the decisive correction.  Its deepest `q=13` grid is
fully trapped, but a different deep component has an endpoint at `7/33` that
escapes.  Therefore:

- one chosen global maximizer is not the component family;
- one prime grid is not the obligation family;
- aggregate folded support is not the signed erosion predicate; and
- a component theorem must retain endpoint owners and test every closed deep
  component against every relevant odd-divisor shell.

THM-797 is the branching regime of the same recursion.  When one obligation
chart is silent, transport to another divisor grid or another component before
declaring a trap.

## 4. The transverse-tooth state and the shape of infinity

THM-798 and THM-799 add a second axis to the proposed object.

THM-798 constructs a four-far family with no high-support divisor scale and
unbounded positive-length component count.  Raw fragmentation is therefore
not controlled by the largest divisor packet.  Yet a proportional top peel
closes the entire family.  The proof-facing quantity is not `r` but a
peel-relative load such as

```text
r_+(P)/(v*|G'_P|).
```

THM-799 proves the sound enclosure transition

```text
(mu,r) --insert N--> (6mu/7-2r/(7N), r+N)
```

and couples it to a named peel rate.  This gives a composable, scale-free
terminal calculus.  It also proves exactly where that compression ceases to
be faithful: two bases can have the same mass, component count, and entire
component-length multiset, yet respond differently to the same inserted
frequency.  The missing coordinate is phase placement relative to the pending
tooth grid, most naturally an owner-coloured endpoint/grid incidence word.

The uniform lacunary consequences change the global picture.  Every fully
lacunary far-count flag is terminal.  Consequently an unresolved unbounded
family cannot run to infinity by arbitrary successive scale separation.  It
must contain a bounded-ratio cluster at some level.  Recursively, the residual
infinity has the form

```text
root packet
  -> bounded-ratio cluster children
       -> bounded-ratio cluster grandchildren
            -> ...,
```

with a terminal cone whenever all remaining children become sufficiently
lacunary.  This is the correct outer base for the component-obligation
groupoid: a bounded-ratio cluster tree whose nodes carry residue, owner,
component, and holonomy stalks.  “Finite height” should now mean finite
normalized cluster types plus certified terminal cones, not a raw maximum
speed cutoff.

## 5. Corrected frontier targets

### Target A: uniform `n=12` sporadic emptiness

This remains the central imported obligation for several LRC14 reductions.
A credible proof must combine, rather than substitute among:

- the binding-scale representation `p/(13s)`;
- hereditary gcd and dyadic/parity descent;
- every owner-labelled component of the quotient core;
- exact component-tooth margins and one-sided endpoint germs;
- sheet ramification and effective deck orders;
- odd divisor-grid shells, beginning with but not ending at `q=13`; and
- a Hamming-radius-at-least-two incidence theorem extending THM-795.

The next local theorem should say that a fully saturated sheet-tooth incidence
either exposes a negative margin, forces a common divisor and descends, or
belongs to a finite explicitly named primitive residual.  “The sporadic branch
is empty” must not be written before that trichotomy is uniform in the scale.

### Target B: finish the two-sheet branch after the `q=13` gate

For full six-class support and the aligned five-class survivor, enumerate all
closed deep components, retain their endpoint owners and signed folded slopes,
and impose the complete odd-divisor obligation family.  The desired theorem is
a global matching obstruction between divisor-grid obligations and components,
not another statement about a selected maximum.  The sharp THM-797 survivor
should be a mandatory unit test because it defeats “deepest `q=13` point” while
still escaping globally.

### Target C: quotient the prime-seven diagonal packet correctly

For `r=8`, contract empty fastest refinements and full-support diagonal packet
loops.  The quotient state must still contain:

- the centered-Beatty base cell and exact ordered packet word;
- the reduced collision state and return map;
- the nonzero metric translation;
- factor-two fixed-owner spans and signed entrant/leaver balance;
- token-supportability; and
- incidence with the relevant core-safe component.

Then seek a theorem that every remaining long normalized collision-SCC path
either changes component incidence, forces a divisor/scale descent, or misses
the core.  Raw active-period boundedness and the old universal extent target
are refuted.  THM-786 §3g also rules out treating positive marginal distance
as a uniform constant.

### Target D: exact transverse insertion over a recursive cluster tree

Use `(safe mass, component count, new wall rate, peel rate)` as a sound
enclosure state, but attach the owner-coloured endpoint/grid incidence whenever
another exact insertion is needed.  The research question is whether that
incidence admits a smaller continuation-equivalent quotient on a bounded-ratio
cluster node.  Candidate coordinates include endpoint residues modulo `N`,
cyclic owner word, signed distances to nearest teeth, and their gcd pushforward.

### Target E: prove a well-founded global descent

The proposed global complexity is lexicographic and proof-relative, not a raw
speed statistic.  Possible decreasing coordinates are:

```text
number of unresolved bounded-ratio cluster nodes;
number of unsatisfied component-obligation incidences;
nontrivial effective sheet order;
normalized collision SCC after gauge loops are contracted;
Hamming radius from a classified rigid packet;
finite residual address.
```

Every move must state which coordinate decreases.  A holonomy loop is not a
decrease; a new observer cut is not a decrease; a finite calculation at one
height is not a decrease.  THM-795 supplies a model divisor descent, THM-797
supplies a model component branch, and THM-799 supplies terminal lacunary
cones.

## 6. Tournament analysis belongs inside the groupoid, not above it

The tournament thread remains useful after the assumption challenge, but its
role is now sharper.

### 6.1 Known fingerprints

- A scalar-order tournament is transitive: score histogram
  `{0:1,1:1,...,m-1:1}`, no directed triangles, singleton SCCs, and one
  Hamiltonian path.  The sheet-margin tournament in THM-795 has this form.
  It sorts margins but cannot encode the theorem predicate “all margins are
  nonnegative.”
- The residue-obligation and deep-component tournaments in THM-797 are also
  transitive under their scalar gauges.  The recorded observer change flips
  edges while leaving the coarse fingerprint transitive.  The missing class,
  exception shell, erosion sign, and endpoint owners decide the theorem.
- The unmarked seven-vertex return tournament in THM-794 is the regular
  tournament `R_7`: every score is `3`, it has `14` directed triangles, one
  SCC, and `175` Hamiltonian paths.  A marked-cut change produces six edge
  flips.  Global token translation does not change these fingerprints, even
  though the metric packet advances.

Thus cycles can signal holonomy and edge flips can signal observer dependence,
but neither is automatically an obstruction or a witness.

### 6.2 A better local tournament

On a fixed component, take **atomic obligations**, not runners, as vertices.
For two obligations orient

```text
o_i -> o_j
```

when the exact transport word discharges or crosses `o_i` first.  The switch
is change of component orientation or marked cut.  Simultaneous events are
kept as tie blocks; the declared tie Hamiltonian path is the exact owner order
in the stored event word.  Attach to every edge:

```text
(owner_i, owner_j, exact time gap, signs, sheet/grid classes, metric margins).
```

The unweighted tournament is then a checksum of the ordered word.  The
decorated tournament is one local chart of the obligation groupoid.  If a
pair is genuinely incomparable because the carrier is bipartite or
hypergraphic, forcing an orientation is information loss; retain the coloured
incidence graph instead.

### 6.3 The transition metagraph

A more faithful global combinatorial object has:

- nodes: isomorphism classes of fully decorated component-obligation states;
- reversible edges: wall, cut, lift, reflection, and return transports;
- directed noninvertible edges: certified descent, peel, or component branch;
- SCCs: recurrent holonomy classes before metric lift;
- condensation edges: candidates for well-founded proof descent; and
- node stalks: exact metric incidence and Hamiltonian-path/observer fibres.

This reframes earlier metagraph work.  A tournament class or converse-merged
node is a base-chart address.  THM-781's Hamiltonian-path quotient explains
the observer fibre.  Neither supplies the endpoint/grid stalk.  The lonely
runner predicate lives on a section of the decorated metagraph, not on its
unweighted node label.

Tournament fingerprints worth recording remain:

```text
score histogram;
directed triangle and longer-cycle counts;
SCC decomposition;
edge flips under every declared gauge;
Hamiltonian-path count and automorphism quotient;
which metric/owner decorations vary inside one fingerprint fibre.
```

The last line is the important one.  A fingerprint becomes mathematically
valuable when its liar fibre names the next sidecar.

## 7. A recursive working protocol

For a new hard family, the following protocol reflects the present frontier.

1. Normalize by a legal dilation, fold, deletion, or affine chart, retaining
   covering and primitivity data.
2. Build the bounded-ratio cluster tree.  Close every fully lacunary branch by
   THM-799 before doing finer enumeration.
3. Enumerate all safe/deep components and one-sided germs exactly, including
   isolated equality points when the closed threshold matters.
4. Choose obligations appropriate to the operation: divisor classes, sheets,
   teeth, token events, packet prefixes, Fourier modes, or proof debts.
5. Build owner-coloured incidence and exact margins.  Do not replace this by
   component count, safe mass, or a tournament unless fibre-purity is proved.
6. Compute the exact transport word and both reduced and metric holonomy.
7. Apply the trichotomy: escape, forced descent, or gauge loop.  If one chart
   is silent, branch over every component and relevant divisor grid.
8. Use tournament fingerprints to compare gauges and expose lost labels, with
   an explicit tie Hamiltonian path.
9. Emit every ambiguous fibre as a finite residual with a replayable exact
   certificate rather than turning empirical exhaustion into a universal
   statement.
10. Record the decreasing complexity or the terminal theorem for every
    recursive edge.

## 8. What the problem now appears to be

The underlying object is not a thirteen-element speed set viewed once.  It is
a section through a recursive arithmetic-geometric stack:

```text
bounded-ratio cluster tree
  + owner-labelled safe/deep components
  + residue/sheet/tooth/divisor obligation fibres
  + exact endpoint and margin incidence
  + centered event words and carry
  + reduced holonomy and metric translation
  + certified descent or residual proof debt.
```

The deepest common phenomenon is *controlled recurrence*.  Coherent packets
reappear at new scales; sheet decks return after a common divisor; collision
states return after diagonal translation; folded components return under a
two-sheet symmetry.  A proof must distinguish three outcomes that older
quotients merged:

```text
return as a harmless gauge symmetry,
return that forces arithmetic descent,
return to the same combinatorics at a different metric location.
```

THM-794 exhibits the third.  THM-795 turns the second into a contradiction.
THM-797 shows why the first observation chart may be silent while another
component escapes.  THM-798/799 show that transverse fragmentation can grow
without creating a new global obstruction and that arbitrary lacunary escape
to infinity is terminal.

This is meaningful progress in the description of the frontier, but it is not
the missing uniform theorem.  The two honest boundary statements remain:

```text
the universal good-period bound q<=25 is false;
uniform emptiness of the n=12 sporadic branch is open.
```

The next proof should be judged by whether it acts faithfully on the decorated
component-obligation state, whether its quotient keeps both reduced holonomy
and metric translation, and whether its recursion genuinely decreases after
all components—not merely a preferred one—have been tested.
