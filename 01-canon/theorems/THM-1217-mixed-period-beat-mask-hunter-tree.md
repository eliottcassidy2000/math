---
id: THM-1217
title: THE MIXED-PERIOD BEAT-MASK HUNTER TREE — quotient-clock lifting, intersection credit, and exact cyclic-run escape
status: PROVED (analytic finite-mask theorem; exact referee; Lean cardinality/tree and corrected-packet kernel)
source: codex-2026-07-19-S82
depends_on: [THM-1179, THM-1192, THM-1216]
related: [THM-1194, THM-1197, THM-1198, HYP-7715, HYP-7845]
script: 04-computation/lrc14_mixed_period_beat_mask_tree_referee_codex_S82.py
output: 05-knowledge/results/lrc14_mixed_period_beat_mask_tree_referee_codex_S82.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCMixedPeriodBeatMaskTree.lean
script_sha256: 3f30fd41b68b524f258e6df3fa3a040534bdc0c6bb80dae1e2106a3cfacbc492
output_sha256: 1dcd2f0ca1a9057199e76081d3fdb2bc2422cc61893128732a3e92ef49cfa349
formalization_sha256: d9c074f3340925c623a2e2b5914eef491831af10596ad08589d61f660ce10fdf
---

# THM-1217 — the mixed-period beat-mask Hunter tree

## Statement

Let `q>0` be a sum or positive-difference beat denominator for six fast
integer speeds `b1,...,b6`, with the defining pair labelled `b5,b6`:

```text
q=b6-b5  or  q=b6+b5.                                  (1)
```

Put

```text
gi=gcd(bi,q),        Qi=q/gi,
d0=gcd(q,b1,...,b6), L=q/d0=lcm_i Qi,
A(Q)=2ceil(Q/14)-1.                                    (2)
```

The strict danger mask of `bi` lifts from `Z/Qi Z` to a mask `Mi` on the
master clock `Z/L Z`, and

```text
|Mi|=Ci=(L/Qi) A(Qi).                                  (3)
```

The defining masks agree, `M6=M5`, and all masks contain zero.  Hence only
the five obligations `M1,...,M5` are relevant.  Write

```text
U=M1 union ... union M5,
C0=1+sum_i(Ci-1).                                      (4)
```

Then

```text
|U| <= C0.                                             (5)
```

In particular, if `C0<L`, every block of `C0+1` consecutive beat numerators
contains a residue outside `U`.

There is a potentially sharper overlap-aware form.  For any spanning tree `T`
on the five mask obligations, put

```text
C_T=sum_i |Mi|-sum_{ij in E(T)} |Mi intersect Mj|.      (6)
```

Then

```text
|U| <= C_T.                                            (7)
```

Thus, if `C_T<L`, every `C_T+1` consecutive beat numerators contain a safe
residue.  Optimizing (6) is the maximum-spanning-tree problem whose edge
weight is mask-intersection cardinality.

Finally, when `U` is proper, let `ell(U)` be the longest cyclic consecutive
run contained in `U`.  The exact threshold is

```text
B_run=ell(U)+1:                                        (8)
```

every `B_run` consecutive numerators escape.  This statement has no analogue
when `U=Z/L Z`; a full clock has no finite escape run.  Likewise the
cardinality suppliers in (5) and (7) require `C0<L` and `C_T<L`, respectively.

If a closed slow safe gap

```text
G_k(a)=[(14k+1)/(14a),(14k+13)/(14a)]                  (9)
```

contains at least one of these threshold numbers of `q`-beat points, the
point outside `U` is safe for `a,b1,...,b6`.  This is a mixed-gcd theorem; it
does not assert that some useful beat denominator or sufficiently long beat
block must occur in every putative LRC(14) counterexample.

## Proof of the master-clock formula

For each prime `r`, write `v=v_r(q)` and `w_i=v_r(b_i)`.  Then

```text
v_r(Qi)=v-min(v,w_i)=max(v-w_i,0).
```

The valuation of the least common multiple is therefore

```text
max_i v_r(Qi)
 =v-min(v,min_i w_i)
 =v_r(q/d0).
```

This proves `L=lcm_i Qi`.  In particular every `Qi` divides `L`, so the
reduced mask has a well-defined lift to the master clock.

Write `bi=gi ui`.  Then `gcd(ui,Qi)=1`, and at a beat point `t=p/q`,

```text
||bi p/q||<1/14
iff 14 min(ui p mod Qi,-ui p mod Qi)<Qi.                (10)
```

Multiplication by the unit `ui` permutes `Z/Qi Z`.  The strict symmetric
window in (10) has `A(Qi)` residues.  Each reduced residue has `L/Qi` lifts,
which proves (3).  If (1) is a difference, `b6` is congruent to `b5` modulo
`q`; if it is a sum, `b6` is congruent to `-b5`.  Since (10) is sign
invariant, `M6=M5` in either case.

## Common-zero and Hunter-tree proofs

Adjoin the five masks one at a time.  Because every row contains residue
zero, each row after the first has intersection of cardinality at least one
with the preceding union.  Four units of overlap are therefore forced:

```text
|U|+4 <= sum_i |Mi|,
```

which is exactly (5).

For (7), root `T` and order its vertices so every parent precedes its child.
When a child mask `M_j` is adjoined, the union exposed so far contains its
parent `M_i`.  Therefore

```text
|M_j intersect (earlier union)| >= |M_i intersect M_j|.
```

Summing the four insertion identities gives (7).  Since any consecutive block
of at most `L` integers has distinct residues modulo `L`, either strict
cardinality inequality forces an escape.  The run statement (8) is immediate
from the definition of the longest cyclic run, with properness essential.

## Corrected `a=79` residual packet

The theorem closes a residual packet that survives the current coarse gates:

```text
a=79,
(b1,...,b6)=(140,210,350,420,490,770),
q=b6-b5=280,
d0=70, L=4,
(Q1,...,Q6)=(2,4,4,2,4,4).                             (11)
```

On `Z/4Z`, the five relevant lifted masks are

```text
({0,2},{0},{0},{0,2},{0}),
U={0,2}.                                                (12)
```

Their sizes sum to seven.  The tree

```text
(M1,M4), (M1,M2), (M1,M3), (M1,M5)
```

has intersection credit `2+1+1+1=5`.  Consequently

```text
C0=3,   C_T=2,   ell(U)=1,
B_common=4, B_tree=3, B_run=2.                          (13)
```

For `k=11`, the slow gap supplies exactly

```text
P_280(11)={40,41,42}.                                  (14)
```

The tree bound already forces escape, and the unique escaping numerator in
this block is `p=41`.  At `t=41/280`, the complete thirteen-speed packet

```text
(1,2,3,4,5,6,79,140,210,350,420,490,770)               (15)
```

has least residues

```text
(41,82,123,116,75,34,121,140,70,70,140,70,70),         (16)
```

all at least `280/14=20`.  Thus (15) is lonely at this exact time.

The row is genuinely beyond the coarse filters used to locate it:

```text
140/79 < 13/6,
70/79 < 397/432,
79 sum_i 1/bi = 21804/13475 > 1.                        (17)
```

It therefore survives the universal first-tooth gate, the older common-gcd
ratio gate, and the harmonic-crowding test, while the mixed mask tree closes
it directly.  The common-period premise of THM-1216 fails because the gcds
alternate between `140` and `70`.

## Exact census and guardrails

The dependency-free referee exhausts all unordered five-mask packets on every
master clock `L<=48`, retaining every distinct unit mask at each divisor
period.  It checks `734765` packets:

```text
proper unions                         582852
full unions                           151913
common-zero cap active (C0<L)         538523
tree cap active (C_T<L)               582796
tree strictly sharper than C0         572063
exact run strictly sharper than C_T   582733
```

Kruskal's maximum-tree result is independently compared with all `125`
labelled five-vertex spanning trees on `382` deterministic packets.  The
first full union is the trivial `Q_i=1` clock, but full unions do **not**
arise only from period one: the first without `Q_i=1` occurs at

```text
L=16, periods=(2,16,16,16,16), sizes=(8,3,3,3,3).       (18)
```

This guardrail is realized, rather than merely an abstract choice of masks.
Take `q=16` and

```text
(b1,...,b6)=(17,35,53,71,88,104),
```

with defining difference `104-88=16`.  The five obligations are

```text
{0,1,15}, {0,5,11}, {0,3,13}, {0,7,9},
{0,2,4,6,8,10,12,14};                                  (19)
```

their union is all of `Z/16Z`, and both `C0` and the optimal tree cap equal
`16`.  This is a full **beat-mask** clock, not a construction covering an
entire slow gap and not an LRC counterexample.

Thus neither positivity of every reduced period nor the presence of a common
zero implies a hole.  Strict cap inequalities or direct properness must be
proved, not inferred.

## Structural and tournament audit

The faithful vertices are not runners.  They are master residues, quotient
clocks, and lifted mask obligations.  The quotient-clock incidence diagram
preserves the exact beat predicate, while its nerve records overlaps and the
cyclic word records block phase.  Passing only to mask sizes destroys quotient
placement and run length.

On mask obligations the pairwise observable is

```text
w_ij=|Mi intersect Mj|.
```

It is symmetric, so any tournament orientation is merely a gauge.  Ordering
by `(period,label)` gives the transitive score histogram `(0,1,2,3,4)`, no
directed cycles, five singleton SCCs, one Hamiltonian path, and ten flips under
gauge reversal.  The useful tournament-related object is the undirected
maximum-weight spanning tree: its four edges are exactly the overlap credits
that survive root-first insertion.  The challenged assumption is therefore
that tournament vertices or arcs should be runners; quotient masks and proof
obligations are the information-preserving choice here.

## Formalization boundary

`LRCMixedPeriodBeatMaskTree.lean` kernel-checks the fibre-cardinality model,
the five-mask common-zero ledger, the arbitrary rooted-tree intersection
ledger, subtraction-free block consumers, and the properness guard for exact
run suppliers.  It also proves the exact master masks, union, witness hole,
period data, beat-gap membership, thirteen least residues, and frontier
arithmetic for (11)--(17).  Prime-valuation identification of the master
clock, the analytic strict-window count, and the slow-gap/block conversion
remain paper-level provider lemmas rather than being overstated as formalized.

## Scope

THM-1217 replaces a single common reduced period by the exact master-clock
object and exposes two layers of usable overlap beyond density: a Hunter tree
and the cyclic run.  It closes (11), but it does **not** prove that every open
LRC(14) packet has `C_T<L`, a proper union, or a beat block long enough to use
the resulting threshold.  Full mixed-period unions such as (18) are the
precise obstruction the next structural lemma must exclude or bypass.
