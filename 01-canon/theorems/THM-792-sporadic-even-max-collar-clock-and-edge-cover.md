---
id: THM-792
title: The sporadic even-maximum collar has a bounded rational blocker clock, repeated top-tooth flank types, a Z/13 moving edge-cover/root-current carrier, and no w=13 packet through quotient height 24
status: PROVED (quantifier-exact reductions and cut-current/energy identities) + VERIFIED (finite-exact w=13 closure for U subset [1,24], by 3/8 in general and 4/11 in the divisor-complete subatlas; does not exclude the unbounded collar)
source: codex-2026-07-14-S10
depends_on:
  - THM-668-pair-sum-ruler-witness-structure
  - THM-772
  - THM-775
  - THM-780
related: [THM-774, THM-776, THM-778, THM-782, THM-785, THM-787, THM-790, THM-791, HYP-6820]
verification:
  - 04-computation/lrc13_w13_sheet_edge_cover_h24_codex_S10.py
  - 05-knowledge/results/lrc13_w13_sheet_edge_cover_h24_codex_S10.out
---

# THM-792 — Sporadic even-maximum collar reductions

Let a primitive tight twelve-speed set lie in the two-sheet packet

```text
A=2U union {x,y},       |U|=10,       x,y odd,
R=max(U),               max(A)=2R.
```

Put

```text
V=U\{R},
P=A\{2R}=2V union {x,y},
H=max(P)<=2R-1,
gamma=M(P).
```

Assume the maximum deletion is sporadic, `gamma>1/12`.  THM-774(6) expresses
this as the nine-core collar.  The present theorem adds three exact structures:
a bounded rational center clock, a quantitative family of occupied top teeth
with disjoint flank owners, and, in the global deep branch, a moving edge cover
on thirteen sheets.  It then proves that this edge automaton always tears when
the forced odd exception is `13` and `U subset [1,24]`.

## 1. Bounded blocker clock and the center/high-order dichotomy

There is a maximizing time `t*=p/q` in lowest terms such that

```text
q<=2H<=4R-2,       gamma=k/q,       12k>q.                (1)
```

At that time the deleted maximum is **strictly** inside a danger tooth:

```text
||2R t*||<1/13.                                           (2)
```

If `r` is the least absolute residue of `2Rp modulo q`, then `13r<q`.
Consequently exactly one of the following holds:

```text
q divides 2R,                                             (3a)

D=q/gcd(q,2R)>=14.                                       (3b)
```

The sporadic margin is also arithmetically separated:

```text
gamma-1/12>=1/(12q)>=1/(48R-24).                         (4)
```

### Proof

The pair-sum/self-cusp theorem in
`THM-668-pair-sum-ruler-witness-structure.md` gives a global maximin point of
the eleven-speed set `P` with reduced denominator at most `2H`.  Every
clearance there has denominator `q`, giving `gamma=k/q` and (1).

All `P` runners are strictly clearer than `1/12` at `t*`, hence remain clearer
than `1/13` in a neighborhood.  Tightness of `A` puts that neighborhood inside
the closed `2R`-danger comb.  Equality in (2) would put `t*` on a regular tooth
boundary; moving to its outward side while staying `P`-safe would make all of
`A` strictly `1/13`-safe.  Thus (2) is strict.

If `r>0`, then `r` is a nonzero multiple of `g=gcd(q,2R)` because `p` is a
unit modulo `q`.  Hence `r>=g` and `q>13r>=13g`, proving (3b).  If `r=0`,
coprimality of `p,q` gives (3a).  Finally `12k-q` is a positive integer, so
`gamma-1/12=(12k-q)/(12q)>=1/(12q)`; use (1) for (4).  ∎

The center branch (3a) is itself pinned.  Put

```text
m=q/gcd(q,2).
```

Then `m|R`, `13` does not divide `q`, and no member of `V` is divisible by
`m`: otherwise its doubled speed in `P` would be divisible by `q` and have
zero clearance at `p/q`.  If `2<=m<=12`, THM-772's divisor completeness makes
`R` the unique `m`-multiple in `U`.  If `V` is imprimitive, THM-775 forces the
dyadic seam: `gcd(V)=2`, `R` is odd, and `q=m` or `2m` for an odd divisor
`m|R`.  Thus the “center” alternative exposes a deleted divisor pin rather
than an unstructured small clock.

## 2. An isolation tax and a uniform number of occupied top teeth

The `H`-Lipschitz interval around `t*` and tightness give

```text
H>=26R(gamma-1/13)
 >=R/6+13R/(24R-12).                                    (5)
```

Moreover, define the strict safe set of the eleven retained speeds

```text
G_P={t:min_(v in P)||vt||>1/13}.
```

Then

```text
Leb(G_P)>=156^(-11),                                     (6)
```

and `G_P` meets at least

```text
N=ceil(13R/156^11)                                       (7)
```

distinct danger teeth of the deleted speed `2R`.

### Proof

The interval

```text
|t-t*|<(gamma-1/13)/H
```

is contained in `G_P`, hence by tightness in one `2R`-danger tooth.  Such a
tooth has length `1/(13R)`.  Comparing lengths proves the first inequality in
(5); substituting (4) gives the second.

For (6), apply THM-780 to this one vector `P` with `d=11`,
`beta=1/12`, and `alpha=1/13`.  The gap is `1/156`, and the phase-pigeonhole
proof lands inside the **strict** alpha-safe set, giving (6).  Tightness says
`G_P` lies in the union of the `2R`-danger teeth.  Each tooth has measure
`1/(13R)`, so (6) requires at least (7) occupied teeth.  ∎

Every occupied top tooth has nonempty blocker sets at both endpoints.  Choose
one retained-speed owner on its left endpoint and one on its right.  The two
owners must be distinct.  Therefore one of the `11*10=110` ordered flank types
occurs on at least

```text
ceil(N/110)                                               (8)
```

occupied teeth.

Indeed, if all retained speeds were strictly safe at one endpoint, openness
would cross outside the deleted runner's closed tooth, contradicting
tightness.  If one speed `s<2R` blocked both endpoints, its phase change across
the tooth would be `s/(13R)<2/13`.  Two different nearest integers would need
phase separation at least `11/13`; the same nearest integer makes `s`
dangerous throughout the tooth, contradicting its occupied `P`-safe interior.
Thus the endpoint blocker sets are disjoint, and pigeonhole proves (8).

This is a genuine asymptotic incidence constraint.  The scalar floor in (6)
is tiny (`156^11=1331590174384995440787456`), but it is uniform and converts
unbounded height into repeated exact boundary data.

## 3. The forced 13-multiple becomes a moving edge cover

Assume additionally the global deep branch of THM-772, and choose an odd
exception

```text
w=13c.
```

No quotient speed `u in U` is divisible by `13`.  Write the thirteen lifts as

```text
tau_j=(s+j)/13,       j in Z/13Z.
```

Whenever `||cs||>2/13`, the exception `w` is ineligible on every lift, so the
ten quotient runners must cover all thirteen sheets.  Away from an event
`us in Z`, runner `u` kills the two-edge

```text
E_u(s)={j:uj=-floor(us), -ceil(us) mod 13}.               (9)
```

Its endpoints differ by `u^(-1) mod 13`.  Thus the ten labelled Cayley edges
form an edge cover of `Z/13Z` with total degree `20`; the overlap-token excess
is exactly

```text
20-13=7.                                                  (10)
```

At a simple event `us=m`, the labelled edge changes by

```text
{u^(-1)(-m+1),u^(-1)(-m)}
 -> {u^(-1)(-m),u^(-1)(-m-1)}.                           (11)
```

Coverage on both sides forces the departing endpoint to have overlap before
the move and the entering endpoint to have overlap afterward: one of the seven
overlap tokens is transported by the event, through the exact displacement

```text
-2u^(-1) mod 13.                                         (11a)
```

Simultaneous events use the corresponding labelled hyperedge update.

The exact number of `u`-events in the open `w`-ineligible region is, with
`g=gcd(u,c)` and `D=u/g`,

```text
g*(D-2 floor(2D/13)-1).                                  (12)
```

This follows by reducing the phase clock `us` against `cs`: over one common
period the event residues form a complete `D`-grid with multiplicity `g`, and
the two closed radius-`2/13` caps remove `2 floor(2D/13)+1` grid positions.

The event schedule in (11) is a merge of ordinary integer grids.  Its owner
word is an ordinary Christoffel/phase-coset word governed by the same Euclidean
division machinery as THM-778, but it is **not** THM-778's centered midpoint
word.  A faithful finite state is therefore

```text
(phase-coset event word, seven overlap chips, ten labelled Cayley edges). (13)
```

This is a predicate-preserving automaton on each of the `c` ineligible
intervals.  It does not yet prove that the automaton must tear; it replaces an
unbounded continuous collar by an exact finite moving-edge-cover problem.

### 3.1 A transitivity-style collision energy (PROVED, but not a quotient)

Let `d_j>=1` be the thirteen sheet degrees in a covered generic chamber and
put `e_j=d_j-1`.  Since `sum d_j=20`, the seven excess chips have collision
energy

```text
K(d)=sum_j binom(e_j,2)=(sum_j d_j^2-34)/2.              (14)
```

At a simple event write the edge slide as `{a,c}->{c,b}`.  Then only
`d_a'=d_a-1` and `d_b'=d_b+1` change, so

```text
K(d')-K(d)=d_b-d_a+1,        Delta(-K)=d_a-d_b-1.        (15)
```

The second identity is exactly the endpoint-defect form of THM-785's
cyclic-triangle flux.  The rescaled coordinate

```text
X_sheet=8K in 8 Z_(>=0),      Delta X_sheet=8(d_b-d_a+1) (16)
```

has the step-eight quantization of THM-787.  Its zero locus is the maximally
spread pattern: seven sheets of degree two and six of degree one.  This is a
genuine algebraic bridge from the tournament flow work to the collar carrier,
but `K` is incidence energy, not a tournament cyclic-triangle count.

Nor is it a Lyapunov function.  The covered degree pairs `(d_a,d_b)=(3,1)`,
`(2,1)`, and `(2,2)` respectively decrease, preserve, and increase `K`.
More decisively, the exact height-24 atlas contains the two cores

```text
(1,2,3,4,5,6,7,8,20,23),
(1,2,3,4,5,6,8,10,11,14).
```

They have the identical sheet-indexed initial degree vector

```text
(6,1,1,2,1,2,1,1,1,1,1,1,1),       K=10, X_sheet=80,
```

but first tear at `4/23` on sheet `1` and at `2/11` on sheet `7`,
respectively.  Thus the labelled event word and transition instance are
indispensable.  THM-791's Hamiltonian-path coordinate is complementary in the
general metagraph, but here chronology is transitive with one path and the
fixed cyclic tie path has no additional resolving power.  The theorem-bearing
object must retain the labelled event word as well as chip incidence.  The next
subsection shows that, once the initial chip state and event word are fixed,
the full labelled edge multiset in (13) is redundant for the coverage test.

### 3.2 The exact cut-current and quadratic-energy cocycle

There is a sharper state description.  Work in generic chambers, applying all
updates at the same rational event before entering its immediately-right
chamber.  Let `d^0` be the sheet-degree vector immediately right of `2/13` and
put

```text
e^0=d^0-1,                  sum_j e^0_j=7.
```

For the event `k/u`, with `r_u=u^(-1) mod 13`, define its departure, entry,
and signed root current by

```text
a_(u,k)=-(k-1)r_u,          b_(u,k)=-(k+1)r_u,
z_(u,k)=delta_(b_(u,k))-delta_(a_(u,k)).                 (17)
```

For a chamber immediately right of `s`, group equal rational events and set

```text
C(s)=sum_(u in U, 2u<13k, k/u<=s) z_(u,k).               (18)
```

Then the seven-chip state and the coverage predicate are exactly

```text
e(s)=e^0+C(s),
edge cover at s  iff  C_j(s)>=-e^0_j for every sheet j.  (19)
```

Equivalently, for every sheet cut `S subset Z/13Z`, put
`C_S=sum_(j in S) C_j`.  This is net labelled current entering `S`, and

```text
e_S(s)=e^0_S+C_S(s).                                    (20)
```

All cut inequalities `e^0_S+C_S>=0` are equivalent to the thirteen singleton
inequalities in (19).  At a first tear, every newly missing singleton has
`C_j=-e^0_j-1`: its initial chip capacity is overdrawn by exactly one.  Thus
the moving-edge obstruction is a labelled walk by roots
`delta_b-delta_a` in the `A_12` root lattice, killed when it exits the integer
simplex

```text
{e in Z_(>=0)^13 : sum e_j=7},                           (21)
```

whose unlabelled state space has `binom(19,7)=50,388` points.  The runner
labels and phase-coset word generate the roots.  Discarding the resulting root
word destroys the predicate; retaining the roots allows the labels themselves
to be compiled out of the coverage test.

On this affine root walk, `K` has the integrated identity

```text
K(s)=K(2/13+)+<e^0,C(s)>+||C(s)||^2/2.                  (22)
```

Indeed, `K(e)=(sum e_j^2-7)/2` and `sum C_j=0`.  More locally, if a grouped
event has total root increment `z`, then, whenever both adjacent chambers are
covered and `e` denotes the pre-event state,

```text
Delta K=<e,z>+||z||^2/2.                                (23)
```

Formula (15) is the simple-root-step specialization.  This identifies what
the scalar energy retains: the squared displacement and its pairing with the
initial chip allocation.  It forgets the thirteen cut coordinates that decide
exit.

If `supp(e)={j:e_j>0}`, then

```text
7-|supp(e)| <= K(e),                                    (24)
```

with equality exactly when every chip stack has height at most two.  In
particular `K=0` means seven unit chips on seven distinct sheets, but it is not
a safe barrier: the exact height-24 atlas has `14,184` initial covers with
`K=0`, and all of them tear.  The cut-current vector, rather than energy alone,
is an exact predicate-preserving quotient once the event roots are fixed.
More precisely, (13) reduces for predicate testing to

```text
(initial excess vector e^0, grouped root-increment word, current cursor).    (25)
```

The next root in (17) updates `e` directly, so no edge reconstruction is
needed to decide coverage.  Runner labels can therefore be compiled out after
they generate the root word.  They remain essential on the arithmetic side,
where they constrain which root words are realizable.  The two energy-liar
cores prove that the root-event word itself cannot be removed from (25).

This separates the uniform problem into two languages.  The grouped sums of
at most ten `A_12` roots form a finite alphabet, and (21), with one added dead
state, recognizes the regular **safe-current language**.  The ten rational
clocks generate a much thinner arithmetic language of realizable grouped root
words.  Uniform collar closure is therefore an intersection problem: show
that the divisor-complete arithmetic language has no word whose whole relevant
factor remains in the 50,388-state safety automaton.  Unbounded quotient
height enlarges the generating clocks and word length, not the chip state
space.

## 4. The `w=13`, quotient-height-24 automaton tears by a short prefix

There is no tight two-sheet packet

```text
A=2U union {13,y},       |U|=10,       U subset {1,...,24},
```

with distinct positive speeds.  This finite-exact conclusion does not require
primitivity, divisor completeness, the collar hypothesis, or a bound on `y`.

Indeed, THM-769 would make the odd exception `13` eligible at every point of
`G_U`.  For `2/13<s<11/13`, however, all thirteen quotient lifts satisfy

```text
||13(s+j)/13||=||s||>2/13.
```

Thus every lift must lie outside `G_U`, so the ten edges (9) must cover all
thirteen sheets throughout that interval.

The companion exact certificate enumerates the

```text
binom(23,10)=1,144,066
```

ten-subsets of `[1,24]` avoiding `13`.  It initializes immediately to the
right of `2/13` and applies all events `k/u` with `2u<13k<11u`, grouping equal
rational times before testing the next chamber.  There are 117 event groups.
At an event a runner's closed danger set is the union of its left and right
edges, so coverage of all open chambers is equivalent to the required closed
predicate, including simultaneous events.

Exactly `101,850` cores cover the initial chamber, but **every one tears by
the event at `3/8`**, after only the first `38` of the `117` grouped events.
Here a tear at `s` means that the immediately-right open chamber is uncovered.
Exactly three cores reach that last event:

```text
(1,2,3,5,6,7,8,9,10,11),
(1,2,3,5,6,7,8,9,11,20),
(1,2,5,6,7,8,9,10,14,22).
```

In each case the `u=8` slide at `3/8` leaves sheet `3` and enters sheet `6`,
overdrawing the sheet-`3` singleton cut.

The THM-772 subatlas is sharper.  Of `131,183` cores that are
divisor-complete for moduli `2,...,12`, `131,149` are primitive, and
respectively `20,612` and `20,604` are initial static covers.  Every
divisor-complete initial cover tears by `4/11`, after `36` event groups.  The
only two that reach that last event are

```text
(1,2,6,7,8,9,10,14,22,24),
(2,6,8,10,14,16,18,20,22,24).
```

The first is primitive and the second has gcd two.  In both, the `u=22` slide
at `4/11` leaves sheet `5` and enters sheet `12`, overdrawing the sheet-`5`
cut.  Thus the primitive divisor-complete atlas has just the first final
carrier.

The collision-energy stratification of the `101,850` initial covers is

```text
K:      0      1      2      3     4    5    6    7  9  10  11
count: 14184  43927  26518  11036  4838  201  922  193  1  28   2.
```

Immediately before each core's first tear the corresponding histogram is

```text
K:      0      1      2     3     4    5    6    7  9  10  11
count: 15270  45098  26853  9518  4254  204  511  127  1  12   2.
```

So neither small nor large energy separates survivors; the exact separation
is the first violated cut in (19).  The canonical event and decision digests
remain

```text
b7b9f5930d28c2dd7e34464851fe00af941561252585d8889b8737682de01cca
01dafb23a9562fc14f21dbe85bc36d6dcb4e93e4b346dde65a74e82fe2addf0c.
```

Hence every possible bounded core has a chamber containing an uncovered
sheet, which supplies a point of `G_U` where `13` is ineligible.  This
contradicts persistent packet ownership and proves the claim.

## 5. Tournament Analysis and the preserved object

A runner tournament may orient `u->v` by which endpoint moves first, use circle
reflection as the switch/gauge, and the centered event order as the tie
Hamiltonian path.  The ordinary fingerprint is largely transitive and loses
the theorem: tooth multiplicity, left/right blocker sets, sheet degrees,
edge labels, and transported overlap chips disappear.

The challenged vertex sets are:

- top teeth with ordered flank-owner sets for (7)--(8);
- thirteen sheet vertices with labelled Cayley edges for (9)--(13);
- thirteen sheet cuts and `A_12` root currents for (17)--(23);
- proof obligations `(tooth side, owner)` or `(event, overlap chip)`.

The first quotient preserves repeated collar incidence; the second preserves
the exact sheet-cover predicate.  Neither is faithfully replaceable by a raw
runner orientation.  The new target is a finite automaton tear or congruence
obstruction, not another scalar measure or tournament score.

## 6. Scope

This theorem does **not** prove the sporadic collar empty.  It proves that every
hypothetical even-maximum collar row simultaneously carries:

1. a rational maximizer of denominator at most `4R-2` with the center/high-order
   dichotomy (3);
2. a uniform positive safe mass distributed over top teeth, with repeated
   disjoint flank-owner types;
3. in the deep branch, a seven-chip moving edge cover on `Z/13Z`, equivalently
   a cut-constrained `A_12` root-current walk driven by an exact
   Christoffel/phase-coset word;
4. for the forced exception `w=13`, an exact tear for every quotient core
   through height 24, by `3/8` in general and `4/11` in the divisor-complete
   subatlas.

These are the next finite/arithmetic coordinates on which a uniform exclusion
can act.
