---
id: THM-1256
title: BINARY CENTERED PHASES LAND ON THE TOOTH WORD, AND TOOTHPICK EDGE RUNS BREAK AFTER TWO STEPS
status: PROVED (residual invoice exactly reinterpreted as endpoint-address order; binary digit/phase order; aligned covered subword versus adjacent marked seam; marked phase permutation is a product of disjoint adjacent swaps; every immediate backtrack has exact positive detuning and outer speed greater than six times the middle speed; ABAB exclusion; backtrack-DAG height at most four in the compact box; sharp half-density nonbacktracking-turn floor; mixed-clock beat export; dependency-free exact referee; sorry-free Lean arithmetic core). This is a structural landing theorem, not six-comb noncoverage or LRC(14)
source: codex-2026-07-19 coherent-invoice consumer
depends_on: [THM-1238, THM-1248, THM-1250, THM-1253, THM-1254]
related: [THM-841, THM-848, THM-1196, THM-1252, THM-1255, HYP-7870]
script: 04-computation/lrc14_binary_phase_word_landing_thm1256.py
output: 05-knowledge/results/lrc14_binary_phase_word_landing_thm1256.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCBinaryPhaseWordLanding.lean
script_sha256: 61179e44b4dfb5ff2cef61897010b23d18c15d848a01d465f09b1e131de64013
output_sha256: d2a2ab94ac2b0ae206656f0ea91fd5d649cc1eb032732e0f58d609368f74708e
formalization_sha256: 19838e1e8d9f5e3903dcd8addcec0482bbe90f14aaf5fe5e568e228d3e288d09
---

# THM-1256 — binary phases land on the tooth word and toothpick runs break

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

## 3. The binary digit is the centered-phase order

For the edge (3), THM-1248 gives

```text
Q_j(t_i-t_j)=theta_i-delta_s,
Q_j=c+d_j>0,             5/28<theta_i<23/28.          (7)
```

Hence

```text
delta_s=0   iff   t_j<t_i,
delta_s=1   iff   t_i<t_j,                            (8)
```

with the quantitative separation

```text
5/(28Q_j)<|t_i-t_j|<23/(28Q_j).                       (9)
```

Compare (8) with the positions of `T_s,T_(s+1)` in (1).  There are only two
possibilities.

### 3.1 Alignment: a literal covered phase corridor

If tooth order agrees with phase order, the chronological subword between
the two marked teeth literally covers the closed segment between `t_i` and
`t_j`.  Indeed consecutive selected intervals overlap, so their union is one
interval containing both endpoint phases.

Runner `d_j` is strictly dangerous at `t_i` and has depth greater than
`1/4` at `t_j`.  The segment therefore exits the marked `d_j`-tooth `T_s`.
At its first wall, `c` and `d_j` are safe while the strict six-cover remains
active.  A selected tooth of another label covers that wall and overlaps
`T_s` on the inward side.  By (2) it is the adjacent tooth in the aligned
subword.  This is an actual THM-1253 handoff and has length at least

```text
1/[14 lcm(d_j,h)].                                    (10)
```

Thus alignment turns the abstract binary edge into a phase-located charged
seam on the exact centered-phase corridor.

### 3.2 Mismatch: the two marked teeth are adjacent

Suppose an earlier word tooth contains the later phase and a later word tooth
contains the earlier phase.  Their interiors must overlap: the later left
endpoint lies before the earlier phase, while the earlier right endpoint lies
after the later phase.  Equation (2) then forces the two marked teeth to be
consecutive in the full word.  Their overlap itself is the charged seam

```text
omega>=1/[14 lcm(s_r,s_0)].                           (11)
```

This proves the promised landing dichotomy:

> Every binary speed-descent blocker edge either aligns with a literal
> centered-phase corridor and its first-exit seam, or is a one-step marked
> inversion whose two blocker teeth directly overlap.

There is a useful global strengthening.  Any inversion between two marked
phases forces their teeth to be consecutive by the same argument.  Two
inversions cannot share a marked tooth: three consecutive reversed marks
would reverse the two outer marks, whose teeth are separated by (2).
Therefore the phase order of all marked teeth differs from their word order
by a product of **disjoint adjacent transpositions**.  The inversion graph is
a matching.  For a simple blocker cycle of length at least three, its
mismatched cycle adjacencies are consequently matching edges; the two-cycle
counts the same unordered inversion twice and is the only duplicate-edge
exception.

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
3. Use the matching law for marked inversions to force enough aligned binary
   corridors that their quantitative phase lengths (9), together with the
   exact union-length identity, exceed the available tooth widths.

No one of these consumers is proved here.

## 7. Tournament and alternate-carrier audit

There are two transitive tournaments on the marked teeth: chronological
tooth order and centered-phase order.  Their edge-disagreement graph is not
an arbitrary tournament switch.  Section 3 proves it is a matching of
adjacent chronological teeth, so the relative permutation is a product of
disjoint adjacent swaps.  The blocker relation remains a typed directed
cycle.  Collapsing these three relations to one runner tournament erases the
proof.

For the full word, the pair observable is the unordered edge symbol `e_a`.
Equality of consecutive symbols is the switch: equality gives a directed
toothpick return whose large endpoint is forced by (14), while inequality
gives a rank-two `V`.  Chronological position is the tie Hamiltonian path.
The ABAB exclusion says the binary switch word has no adjacent `1`s.

We challenged runners, marked phases, selected teeth, tooth boundaries,
handoff cells, edge symbols, backtrack positions, Farey beats, and proof
obligations as vertices.  The faithful object is (23) together with the
entire minimal word and the third-blocker beat attached to each backtrack.
It preserves strict cover truth, phase order, lcm mass, detuning, and the
off-word obligation.  It destroys alternative unselected blockers and does
not yet supply a well-founded global height.

## 8. Verification and scope

The dependency-free exact referee checks `979,836` endpoint residual/address
rows, `14,550` binary phase rows, `3,542` minimal interval chains, `20,834`
marked interval pairs, `27,637` explicit phase/order mismatch witnesses,
the adjacent-swap matching law on `61,336` marked phase words, `340`
exact same-provider tooth triples, and `63,276` ABAB-free owner words.  The
maximum unordered-edge run is exactly two and the turn-floor slack attains
zero.  Normal and optimized outputs are byte-identical.

The sorry-free Lean module kernel-checks (4)--(6), both binary phase signs,
the mismatch-overlap and nonconsecutive-separation lemmas, the adjacent-swap
matching obstruction, aligned pair coverage, the exact detuning identity,
the `h>6j` implication, the ABAB contradiction, the five-edge compact-box
chain obstruction, and the natural-number turn floor.  Compact
minimal-subcover extraction, iteration from a pair to a full
subword, first-wall continuation, and the lcm endpoint quantum remain
explicit paper topology providers.  There are no proof placeholders or
`native_decide` calls.

Frozen artifact hashes are

```text
source         61179e44b4dfb5ff2cef61897010b23d18c15d848a01d465f09b1e131de64013
output         d2a2ab94ac2b0ae206656f0ea91fd5d649cc1eb032732e0f58d609368f74708e
formalization  19838e1e8d9f5e3903dcd8addcec0482bbe90f14aaf5fe5e568e228d3e288d09
```

THM-1256 does not prove six-comb noncoverage or LRC(14).  It prevents the
new positive residual from being miscounted as `H`-mass, then proves the
strongest honest landing currently available: binary blocker order is either
a literal centered-phase corridor or an adjacent charged seam; all marked
order defects are disjoint adjacent swaps; and a minimal tooth word cannot
repeat a two-label toothpick return.
