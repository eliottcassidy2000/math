---
id: THM-1256
title: COHERENT BLOCKER PHASES ALIGN WITH THE TOOTH WORD, THE CYCLE MINIMUM EXPORTS A THREE-HANDOFF CORRIDOR, AND TOOTHPICK RUNS BREAK
status: PROVED (residual invoice exactly reinterpreted as endpoint-address order; THM-1248 target clearance doubles actual-blocker phase separation; every coherent blocker edge is automatically phase/word aligned because the target's own spoke lies outside its blocking tooth; every binary edge has a target-free covered corridor with five-comb harmonic invoice; at the cycle minimum the next marked tooth is nonadjacent and a canonical three-handoff/rank-two local forest results; two-wall trichotomy reducing the incidence-free/reuse-free branch to a slowest two-cycle with a forced located double-star-plus-bridge spanning tree, protected rank-three forest, and quotient-plus-word nested-carrier tariff; structurally binary six-cycles contain at least two digits of each kind; arbitrary marked inversions remain disjoint adjacent swaps; every immediate backtrack has exact positive detuning and outer speed greater than six times the middle speed; ABAB exclusion; backtrack-DAG height at most four; sharp half-density nonbacktracking-turn floor; mixed-clock beat export; dependency-free exact referee; sorry-free Lean arithmetic core). This is a structural landing theorem, not six-comb noncoverage or LRC(14)
source: codex-2026-07-19 coherent-invoice consumer
depends_on: [THM-1238, THM-1248, THM-1250, THM-1253, THM-1254]
related: [THM-841, THM-848, THM-1196, THM-1252, THM-1255, HYP-7870]
script: 04-computation/lrc14_binary_phase_word_landing_thm1256.py
output: 05-knowledge/results/lrc14_binary_phase_word_landing_thm1256.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCBinaryPhaseWordLanding.lean
script_sha256: 022f2f6693a85088f8f894628c66e101f46c14997f7a343c32d2d6d8f1cab5fb
output_sha256: 46c2f23ecaee4b21f30081f11e79fb5808ca9c650efbc6bb76746530998a9820
formalization_sha256: cfb0fc36ac88162846b3373a02913234ce790b160a2608c7e394ec4f01de9191
---

# THM-1256 — coherent blocker phases align and toothpick runs break

## 1. Setup

Let

```text
G=G_k(c)=[(14k+1)/(14c),(14k+13)/(14c)],   0<=k<c,
c<d_1<...<d_6,
```

and suppose the six strict danger combs cover `G`.  Choose THM-1253's
deletion-minimal cover by individual fast teeth and write its chronological
word as

```text
I_1,...,I_A,       I_a=(alpha_a,beta_a),
owner(I_a)=s_a,    centre(I_a)=n_a/s_a.               (1)
```

Both endpoint sequences are strictly increasing, consecutive teeth overlap,
and minimality gives the decisive separation

```text
beta_a<=alpha_(a+2).                                  (2)
```

Choose every THM-1240 centered-spoke blocker from this same word, as in
THM-1254.  On a simple blocker cycle, let `T_s` be the selected blocker tooth
at the centered phase `t_(v_s)`.  Its owner is `v_(s+1)`.  The `T_s` are
distinct because their owners are distinct.

THM-1254 selects a speed-descent edge

```text
i=v_s -> j=v_(s+1),       d_j<d_i,                    (3)
```

whose relative digit `delta_s` is binary.  Put

```text
T_s       = the n_r-th tooth of s_r=d_j, containing t_i,
T_(s+1)   = the n_0-th tooth of s_0=d_(v_(s+2)), containing t_j.
```

The theorem below identifies exactly what THM-1254's positive residual says,
then lands the phase order on the interval word.

## 2. The residual invoice is endpoint-address order

Suppose first that `T_(s+1)` precedes `T_s`, so the chronological path runs
from `(s_0,n_0)` to `(s_r,n_r)`.  Its reciprocal-centre drift is

```text
Delta=s_0/n_0-s_r/n_r>0.
```

Define the positive endpoint determinant

```text
D_path=n_0 n_r Delta=n_r s_0-n_0 s_r.                 (4)
```

Coherence gives, with `P=P_j`,

```text
P=n_r+k+delta_s,
R=P s_0-n_0 s_r.
```

Consequently

```text
R-D_path=s_0(P-n_r)=s_0(k+delta_s),                   (5)
R>=D_path       iff       P>=n_r.                     (6)
```

The reflected row has the identical interpretation with endpoint-address
surplus `c-k-delta_s`.  Thus THM-1254's unconditional inequality

```text
R>=n_0 n_r Delta
```

is exact and useful, but it is an **address/position invoice**, not a second
piece of overlap mass.  In particular it cannot simply be added to
THM-1253's `H`-debt.  The handoff conservation law

```text
14h_a+14s_as_(a+1)omega_a=s_a+s_(a+1)
```

also points in the opposite direction: centre drift and raw overlap are
complementary parts of the two tooth widths.  Any consumer must retain the
phase segment or a new wall event, not scalarize `R` as measure.

## 3. Every coherent blocker edge is phase/word aligned

For the edge (3), THM-1248 gives

```text
Q_j(t_i-t_j)=theta_i-delta_s,
Q_j=c+d_j>0,             5/28<theta_i<23/28,
|theta_i-delta_s|>5/14.                              (7)
```

Hence

```text
delta_s=0   iff   t_j<t_i,
delta_s=1   iff   t_i<t_j,                            (8)
```

with the quantitative separation

```text
5/(14Q_j)<|t_i-t_j|<23/(28Q_j).                       (9)
```

The old central band alone gave `5/(28Q_j)`; the target-wall leg doubles the
uniform lower separation.  More importantly, mismatch is impossible for an
actual coherent blocker edge.  Write

```text
T_s=J=(a_L,a_R),       t_i in J.
```

At its own centered spoke `t_j`, runner `d_j` has depth greater than `1/4`,
so `t_j` lies strictly outside the closure of `J`.  If `delta_s=0`, (8) puts
`t_j` to the left of `t_i`, hence actually `t_j<a_L`; if `delta_s=1`, then
`a_R<t_j`.  The next marked tooth `T_(s+1)` contains `t_j`.  In the first
case its left endpoint is earlier than `a_L`; in the second its right
endpoint is later than `a_R`.  Both endpoint sequences in (1) increase with
word position.  Therefore `T_(s+1)` occurs on exactly the same side of `T_s`
in the word as `t_j` does in phase order.

This argument uses only that the relative digit is integral: `delta<=0`
selects the left side and `delta>=1` the right.  Thus it applies to **every**
edge of the coherent lasso, not only the binary speed descent.

### 3.1 A literal target-free covered corridor

The chronological subword between the two marked teeth literally covers the
closed segment between `t_i` and `t_j`: consecutive selected intervals
overlap, so their union is one interval containing both phases.  At the first
wall `w` of `J`, the carrier and `d_j` are safe while the strict six-cover
remains active.  The adjacent selected tooth overlaps `J` on the inward side,
giving the actual charged seam

```text
omega>=1/[14 lcm(d_j,h)].                            (10)
```

THM-1248's sharper target-wall calculation gives

```text
|t_j-w|>5/[14(c+d_j)].                              (10a)
```

For the binary minimum-speed edge, (9) also gives
`|t_i-t_j|<23/[28(c+d_j)]<6/(7d_j)`, so the segment from
`w` to `t_j` cannot reach the next `d_j`-tooth.  It is a literal interval
`C_j` on which both `c` and `d_j` are safe and the other five fast combs
cover.  The sharp one-interval discrepancy

```text
mu(C_j intersect D_h)<=|C_j|/7+6/(49h)
```

therefore yields the phase-located harmonic side invoice

```text
sum_(h!=j)1/h >=(7/3)|C_j|>5/[6(c+d_j)].             (10b)
```

### 3.2 The cycle minimum forces a nonadjacent marked tooth

Let `d_j` be the minimum speed on the blocker cycle and write the local
cycle triple as `d_i->d_j->d_h`.  Then `d_h>d_j`.  If the marked
`d_h`-tooth containing `t_j` were adjacent to `J`, that one tooth would have
to span the wall-to-spoke corridor (10a).  Hence

```text
1/(7d_h)>5/[14(c+d_j)],
5d_h<2(c+d_j)<4d_j,                                  (11)
```

contradicting `d_h>d_j`.  The marked teeth are therefore nonadjacent.

Let the four consecutive local owners beginning at the opposite neighbour
of `J` be

```text
h_0 -> d_j -> h_1 -> h_2.                            (11a)
```

These give three distinct chronological handoff occurrences.  If
`h_0!=h_1`, the two walls of `J` already give the rank-two forest
`{h_0,d_j},{d_j,h_1}`.  If `h_0=h_1`, then `h_2!=d_j` by the ABAB exclusion
proved below, so `{d_j,h_1},{h_1,h_2}` is a rank-two forest.  Thus the cycle
minimum canonically exports a phase-located three-handoff window and two
distinct projected seam edges; no arbitrary second-edge choice is needed.

### 3.3 The incidence/reuse-free residual is a canonical six-label path

THM-1248's two-wall owner pigeonhole gives a recursive split before any
metric gain argument.  If a lasso on `v` labels has neither a wall owner on
the lasso nor a repeated outside owner label, its `2v` wall occurrences
inject into only `6-v` outside labels.  Thus `2v<=6-v`; since the cycle has
length at least two, necessarily

```text
v=2,        (L,m)=(0,2).                             (11b)
```

The lasso is a slowest-rooted two-cycle and its four walls have the four
outside labels as distinct owners.  Section 3.2 makes the two marked target
teeth nonadjacent.  If exactly one tooth lay between them, that tooth would
own both facing walls, repeating an outside owner.  Hence there are at least
two intermediate teeth.  Including the two outer wall owners produces a
literal chronological subword

```text
h_L -> d_b -> h_1 -> ... -> h_2 -> d_a -> h_R,       (11c)
```

where `h_L,h_1,h_2,h_R` are the four distinct outside labels.  It contains a
canonical spanning path through all six labels and at least five pairwise
disjoint handoff occurrences, attached to the two centered blocker phases.

Equivalently, the four target-wall seams form the rank-four forest

```text
{d_b,h_L}, {d_b,h_1}, {d_a,h_2}, {d_a,h_R},          (11d)
```

two disjoint three-vertex stars.  The coherent subword between their facing
owners joins the two components; choose its first crossing handoff `{u,v}`.
It is a fifth full, disjoint seam and (11d) plus `{u,v}` is a forced located
spanning tree on all six labels.  In particular

```text
H_fast >= 1/c+(7/12)[1/lcm(d_b,h_L)+1/lcm(d_b,h_1)
          +1/lcm(d_a,h_2)+1/lcm(d_a,h_R)+1/lcm(u,v)]. (11e)
```

The two `d_b`-wall seams and the bridge all lie in THM-1244's protected
`d_a`-safe component: the whole `d_b`-target tooth is protected, and the
subword remains `d_a`-safe until the facing wall of the `d_a`-target tooth.
They form a located rank-three forest on the five blocker labels.  If

```text
H_B=sum_(h!=a)1/h,       g_B=gcd(d_b,h_L,h_1,h_2,h_R),
```

then its three seam lengths total at least `3g_B/(14d_6^2)`.  The protected
Hunter rearrangement improves branch-specifically to

```text
H_B >=(6d_a-c)/[3d_a(c+d_a)]+7g_B/(4d_6^2),          (11e')
```

and only one unlocated edge remains in a five-label spanning-tree extension.

This path also produces the first exact **quotient-plus-word nested-carrier
tariff**.  On the descent `d_b->d_a`, orient so the marked `d_b`-tooth
containing `t_a` precedes the marked `d_a`-tooth.  Write their word positions
as `p<q`; then `q-p>=3`.  Every internal handoff

```text
W_r=I_r intersect I_(r+1),       p+1<=r<=q-2,
```

lies wholly inside the target-free wall-to-spoke corridor: two-step
separation (2) puts its right endpoint before the target wall, while the
earlier marked phase lies before its left endpoint.  These handoffs are
pairwise disjoint and

```text
|W_r|>=1/[14 lcm(s_r,s_(r+1))].                      (11f)
```

Writing `C_ba>5d_a` for THM-1248's cleared target-wall numerator and
`Q_a=c+d_a`, the corridor length is exactly

```text
L_a=C_ba/[14d_aQ_a].
```

Because the other five combs cover the corridor and every `W_r` is counted
twice, the interval discrepancy gives

```text
L_a+sum_r|W_r| <= 5L_a/7+(6/49)sum_(h!=a)1/h,

sum_(h!=a)1/h
 >= C_ba/[6d_aQ_a]+(49/6)sum_r|W_r|
 > 5/(6Q_a)+(7/12)sum_r 1/lcm(s_r,s_(r+1)).          (11g)
```

Unlike the earlier scalar tariffs, (11g) couples one exact centered quotient
clearance to a full internal seam selected by chronological word order.

The residual two-cycle also has a finite exact clock/digit split.  Put

```text
R=(c+d_b)/(c+d_a)>1.
```

If the descent digit is zero and `m=delta_(a->b)>=1`, the two edge equations
give

```text
R=(m-theta_up)/theta_down,
(28m-23)/23 < R < (28m-5)/10.                       (11h)
```

If the descent digit is one and `p=-delta_(a->b)>=0`, then

```text
R=(p+theta_up)/(1-theta_down),
(28p+5)/23 < R < (28p+23)/10.                       (11i)
```

Here the descent magnitude lies in `(5/14,23/28)` while
`theta_up in (5/28,23/28)`.  In particular, if the ascent is binary too,

```text
R<23/10.                                             (11j)
```

In that binary-ascent subbranch the facing gap between the two nonadjacent
marked target teeth is safe for both `d_a` and `d_b`.  Its rational walls and
positivity give length at least `1/[14 lcm(d_a,d_b)]`; the four outside combs
cover it.  The one-interval discrepancy now gives the additional exact
invoice

```text
sum_outside 1/h >=1/[4 lcm(d_a,d_b)].                (11j')
```

The nonbinary residual guardrail has `R=10/3` and ascent digit `2`, showing
that the other strips are genuinely live.

This is a sharper state, not a contradiction.  THM-1248's exact guardrails
realize it on lonely packets, both with binary edge digits and with a
nonbinary ascent digit.  The unconditional nested object is the target-free
five-comb corridor on the descent edge; a four-comb facing-gap invoice needs
the extra hypothesis that the ascent is binary.

### 3.4 Arbitrary marked inversions are still adjacent swaps

The earlier mismatch lemma remains useful for arbitrary marked pairs not
related by a blocker edge.  If an earlier word tooth contains the later
phase and a later word tooth contains the earlier phase, their interiors
overlap, so (2) forces them to be consecutive.  Two such inversions cannot
share a mark, since that would reverse two outer teeth separated by (2).
Hence the phase order of arbitrary marked teeth differs from word order by a
product of disjoint adjacent transpositions.  For the coherent blocker lasso,
the preceding argument proves that none of its own adjacency edges is such a
mismatch.

### 3.5 A structurally binary six-cycle has two digits of each kind

There is one genuinely six-specific finite word exclusion.  Suppose the
cycle has length six and every edge satisfies the structural binary
condition

```text
d_(s+1)<=c+d_s.                                      (11k)
```

Then every digit is `0` or `1`.  It is impossible for either digit to occur
only once.  Indeed, rotate a putative singleton word so that the other five
edges give

```text
t_0>t_1>...>t_5,       Q_s=c+d_s.
```

The sharp phase floor and centered-spoke errors imply

```text
t_0-t_5 >(5/14)sum_(k=1)^5 1/Q_k,
t_0-t_5 <=1/(2Q_0)+1/(2Q_5).                         (11l)
```

Condition (11k) gives `Q_k<=Q_0+kc` for `k<=4`, while the closing edge and
`d_5>c` give `Q_5>=max(2c,Q_0-c)`.  With `x=Q_0/c>2`, (11l) would force

```text
(5/14)sum_(k=1)^4 1/(x+k)
 <1/(2x)+1/[7 max(2,x-1)].                           (11m)
```

The reverse inequality always holds.  Multiplying by `14x`, its strict form
is

```text
5 sum_(k=1)^4 x/(x+k) > 7+2x/max(2,x-1).
```

For `2<x<=5/2`, the left side is greater than `19/2` and the right side is at
most `19/2`.  For `5/2<=x<=3`, the left side is at least
`25(1/7+1/9+1/11+1/13)>10`, while the right side is at most `10`.  For
`x>=3`, the left side is at least `319/28>10` and the right side is again at
most `10`.  Reversing phase order handles the other singleton digit.

Therefore a structurally binary six-cycle has at least two zeroes and at
least two ones.  After the trivial all-one/all-zero monotonic words are
removed, this excludes the twelve weight-one/weight-five words and leaves
exactly `C(6,2)+C(6,3)+C(6,4)=50` binary types.

## 4. Toothpick self-similarity breaks after one backtrack

Now forget the marked cycle and inspect any internal triple in the full word.
Suppose it is an immediate backtrack

```text
h -> j -> h.                                          (12)
```

The two outer `h`-teeth are distinct.  Write their addresses as
`M_-<M_+`, put `r=M_+-M_->=1`, and let `omega_-,omega_+` be their two proper
overlaps with the middle `j`-tooth.  The gap between the inward walls of the
outer teeth is `(7r-1)/(7h)`, whereas the middle tooth has width `1/(7j)`.
Therefore

```text
omega_-+omega_+
 =1/(7j)-(7r-1)/(7h)
 =[h-(7r-1)j]/(7jh)>0.                                (13)
```

In particular

```text
h>(7r-1)j>=6j.                                        (14)
```

Equation (13) is THM-1252's detuned toothpick law, now valid for **every**
immediate backtrack in the minimal word, not only a marked centered target.
It rules out

```text
A -> B -> A -> B,                                     (15)
```

because the first triple would give `A>6B` and the second `B>6A`.
Equivalently, if

```text
e_a={s_a,s_(a+1)},       a=1,...,A-1,                 (16)
```

then every run of one unordered edge symbol has length at most two.  The
ladder can make one toothpick return, but it cannot repeat the return in the
opposite orientation.  Longer star-shaped recurrence remains possible, so
this is a precise break in literal self-similarity rather than a claim of
termination.

There is nevertheless a genuine well-founded coordinate.  Orient every
backtrack pair from its middle owner to its repeated outer owner,

```text
j -> h.                                               (16a)
```

Equation (14) makes this a strict ratio ascent by more than six.  The
backtrack graph is therefore acyclic.  Since THM-1233 places every fast speed
in

```text
c<d_i<2345c,
6^4=1296<2345<7776=6^5,
```

every directed backtrack path has at most four edges.  This does not yet say
that the third blocker emitted at a pair-sum beat follows the next edge of
this graph; it supplies the bounded height that such a recursive transport
would consume.

## 5. The full word has a two-channel internal ledger

Let

```text
B=#{a: e_(a-1)=e_a},
T=#{a: e_(a-1)!=e_a},                                 (17)
```

where `a=2,...,A-1`.  Thus `B` counts immediate backtracks and `T` counts
nonbacktracking turns.  Every internal position has exactly one type:

```text
B+T=A-2.                                              (18)
```

The ABAB exclusion says backtrack positions are never consecutive, whence

```text
B<=floor((A-1)/2),
T>=floor((A-2)/2).                                    (19)
```

The six-label surjectivity argument of THM-1253 independently gives `T>=4`.
THM-1250's private mass gives at least

```text
S=sum_i ceil(d_i/(7c))                                (20)
```

selected owner teeth.  Combining these facts yields the uniform sharpened
turn floor

```text
T>=max(4,floor((S-2)/2)).                             (21)
```

The two types retain different proof data.

* A turn `u->j->v`, `u!=v`, is a rank-two `V` made of two distinct,
  disjoint chronological seam occurrences.
* A backtrack `h->j->h` has two disjoint occurrences of one lcm edge, the
  positive detuning (13), and `h>6j`.  Applying THM-1238 to `(j,h)` supplies
  a safe pair-sum beat in `G` with a literal third blocker and nonzero
  mixed-clock residue.  Its oriented return lies in the height-four
  backtrack DAG.

Thus every internal tooth is either a rank-two wall fork or an exact
detuned toothpick return that exports an off-word beat obligation.  Since
backtracks cannot be adjacent, at least half of a long word consists of the
rank-two branch.

This is an operation-level refinement of THM-1253's full invoice, not an
extra scalar invoice.  Adjacent internal positions can share a handoff, so
one must not count two new seam masses per position.  The valid quantitative
consumer remains

```text
H_fast>=1/c+(7/12)sum_(a=1)^(A-1)1/lcm(s_a,s_(a+1)),
F_6>=(c/16)sum_(a=1)^(A-1)1/lcm(s_a,s_(a+1)).         (22)
```

The new content is that every summand now lies in a word whose local
recurrence cannot stay in a two-label toothpick fibre.

## 6. Farey, `H`-drift, and alternate-gap interpretation

The exact Farey operation in THM-841 retains an ordered neighboring pair;
count-only toothpick recursion fails.  Here the ordered state is

```text
(word position; owner speed; tooth address; centered phase mark;
 phase/position inversion bit; left and right handoff germs).             (23)
```

Formula (5) is the `H`-drift's functional warning: endpoint determinant is a
transport coordinate, while (22) is measure.  Formula (13) is the true
metric failure mass of a toothpick return.  Formula (15) shows why a ladder
recursion must change its owner pair after at most one return.

The third-blocker beat emitted by a backtrack is not yet a smaller safe gap:
the pair-sum clock is integral at its own beat.  A well-founded descent still
has to transport that blocker to a complete safe component or prove that it
creates a genuinely new seam occurrence.  THM-1255 independently closes the
carrier-`41` finite face; it does not turn the present all-scale beat into an
automatic carrier descent.

The precise remaining implication is therefore one of the following.

1. Show that iterating the third-blocker beats from the `B` branch follows
   the height-four backtrack DAG, or strictly decreases a phase-aware
   Farey/address height when it leaves that DAG.
2. Show that one iteration lands outside the already invoiced `A-1`
   handoffs, creating independent excess in (22).
3. Start from the always-aligned, target-free cycle-minimum corridor
   (10a)--(11a) and either descend to a smaller complete gap/address cell or
   convert its `>5/[6(c+d_j)]` five-comb invoice into metric excess beyond
   the three already-counted local handoffs.

No one of these consumers is proved here.

## 7. Tournament and alternate-carrier audit

There are two transitive tournaments on the marked teeth: chronological
tooth order and centered-phase order.  Their edge-disagreement graph on
arbitrary marks is a matching of adjacent chronological teeth, so the
relative permutation is a product of disjoint adjacent swaps.  Section 3
proves the sharper fact that every actual blocker edge lies in the
**agreement** relation.  The blocker relation remains a typed directed cycle
inside that agreement graph.  Collapsing these three relations to one runner
tournament erases the proof.

For the full word, the pair observable is the unordered edge symbol `e_a`.
Equality of consecutive symbols is the switch: equality gives a directed
toothpick return whose large endpoint is forced by (14), while inequality
gives a rank-two `V`.  Chronological position is the tie Hamiltonian path.
The ABAB exclusion says the binary switch word has no adjacent `1`s.

We challenged runners, marked phases, selected teeth, tooth boundaries,
target-free corridors, handoff cells, edge symbols, backtrack positions,
Farey beats, and proof obligations as vertices.  The faithful object is
(23) together with the entire minimal word, the two-order agreement bit, the
target-free wall-to-spoke interval, and the third-blocker beat attached to
each backtrack.  It preserves strict cover truth, phase order, lcm mass,
detuning, and the off-word obligation.  It destroys alternative unselected
blockers and does not yet supply a well-founded global height.

## 8. Verification and scope

The dependency-free exact referee checks `979,836` endpoint residual/address
rows, `14,550` central binary phase rows (`10,500` surviving the sharp
actual-blocker filter and `4,050` rejected), `11,000` exact target-free
corridor and five-comb invoice rows, `4,647,600` cycle-minimum marked-
adjacency exclusions, `713` two-wall spanning-path rows, `81,750` nested
quotient/word tariffs, all `64` structurally binary six-cycle words (excluding
the `12` singleton-digit words and retaining the exact `50` balanced words),
and `1,174` slowest-two-cycle digit strips.  It further checks `3,542`
minimal interval chains, `20,834` marked interval pairs, `27,637` explicit
phase/order mismatches, `395,342` safe-target automatic-alignment witnesses,
the adjacent-swap matching law on `61,336` arbitrary marked phase words,
`340` exact same-provider tooth triples, `63,276` ABAB-free owner words, and
`87,984` same-provider four-windows in which the next edge must turn.  The
maximum unordered-edge run is exactly two and the turn-floor slack attains
zero.  Normal and optimized outputs are byte-identical.

The sorry-free Lean module kernel-checks (4)--(6), both binary phase signs
and their doubled separation, automatic coherent-edge alignment, the sharp
post-target-tooth corridor, the adjacent marked-tooth speed drop and
cycle-minimum contradiction, the five-comb harmonic invoice, the nested
handoff tariff, the four-comb facing-gap invoice, the protected three-seam
coefficient, the two normalized binary-six-cycle ratio obstructions, both
two-cycle digit strips, the three- and five-handoff counts, and the
same-provider next-turn dichotomy.  It also checks the general mismatch-
overlap and nonconsecutive-separation lemmas, adjacent-swap matching,
aligned pair coverage, exact detuning, the `h>6j` implication, ABAB
contradiction, the five-edge compact-box chain obstruction, and the natural-
number turn floor.  Compact minimal-subcover extraction, iteration from a
pair to a full subword, first-wall continuation, and the rational-wall-to-lcm
topology remain explicit paper providers.  There are no proof placeholders
or `native_decide` calls.

Frozen artifact hashes are

```text
source         022f2f6693a85088f8f894628c66e101f46c14997f7a343c32d2d6d8f1cab5fb
output         46c2f23ecaee4b21f30081f11e79fb5808ca9c650efbc6bb76746530998a9820
formalization  cfb0fc36ac88162846b3373a02913234ce790b160a2608c7e394ec4f01de9191
```

THM-1256 does not prove six-comb noncoverage or LRC(14).  It prevents the
new positive residual from being miscounted as `H`-mass, then proves the
strongest honest landing currently available: every coherent blocker edge,
not merely a selected binary edge, is a literal centered-phase corridor; the
cycle minimum supplies a target-free interval, a five-comb harmonic invoice,
and a canonical three-handoff/rank-two window; arbitrary marked order defects
remain disjoint adjacent swaps; and a minimal tooth word cannot repeat a
two-label toothpick return.  The open step is the well-founded descent or
new metric conversion of this aligned typed cell.
