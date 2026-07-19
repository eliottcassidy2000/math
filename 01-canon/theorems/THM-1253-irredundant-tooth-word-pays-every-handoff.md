---
id: THM-1253
title: IRREDUNDANT TOOTH WORD PAYS EVERY HANDOFF — the full chronological seam invoice
status: PROVED (minimal interval-chain separation; all raw consecutive handoffs pairwise disjoint; full unweighted and twelve-piece functional multiplicity-excess charges; factor-three upgrade over the Cayley-averaged owner-reuse debt; sharp four-turn word skeleton; dependency-free exact referee; sorry-free Lean interval/arithmetic core)
source: codex-2026-07-19-S78 continuation
depends_on: [THM-1166, THM-1198, THM-1250]
related: [THM-1178, THM-1199, THM-1248, THM-1252, THM-1254]
script: 04-computation/lrc14_full_chronological_seam_invoice_thm1253.py
output: 05-knowledge/results/lrc14_full_chronological_seam_invoice_thm1253.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCFullChronologicalSeamInvoice.lean
script_sha256: 6e2a1866562d3f7097e8c730185fd709ebe12ec5f10aca0372f1ae5a8a3d51e4
output_sha256: df3f3b6b28725ba0c350854e1561699be537f5dd4f724ff305e2069247714e90
formalization_sha256: 36d53a66ca2f4c5ea9fe8e8437a6685f45ee42512ea509e791255cd3baf1b279
---

# THM-1253 — the full chronological seam invoice

## 1. Setup

Let

```text
G=[(14k+1)/(14c),(14k+13)/(14c)]
```

be a complete closed safe gap of the positive integer speed `c`.  Suppose
the strict danger combs of six distinct speeds

```text
c<d1<...<d6
```

cover `G`, and put `H=sum_i 1/d_i`.  As in THM-1250, choose a
deletion-minimal finite subcover by individual open danger teeth and order it
by left endpoint:

```text
I_a=(alpha_a,beta_a),       a=1,...,N.                (1)
```

No selected interval contains another, so both endpoint sequences are
strictly increasing.  Connected coverage and minimality also give

```text
alpha_(a+1)<beta_a,                                  (2)
```

and every raw consecutive overlap

```text
W_a=(alpha_(a+1),beta_a)                              (3)
```

lies in `int(G)`.  THM-1250 retained a spanning tree of these handoffs, then
used Cayley averaging to charge repeated owner pairs.  The missing interval
order observation below charges every occurrence directly.

## 2. Minimality separates all handoffs

For every `a<=N-2`, one has

```text
beta_a<=alpha_(a+2).                                  (4)
```

Indeed, if `alpha_(a+2)<beta_a`, then `I_a` and `I_(a+2)` overlap.  Because
their left and right endpoints are ordered, their union contains all of
`I_(a+1)`:

```text
I_(a+1) subset I_a union I_(a+2).
```

This makes `I_(a+1)` deletable, contrary to minimality.  Equality in (4) is
allowed and is exactly why the strict open-endpoint convention causes no
problem: the adjacent handoff intervals then approach the common wall from
opposite sides but do not meet.

It follows from (4) and endpoint monotonicity that the entire family

```text
W_1,...,W_(N-1)                                       (5)
```

is pairwise disjoint.  This is stronger than pairwise disjointness only for
repetitions of one unordered owner pair.

If `s_a,n_a` are the speed and centre address of `I_a`, the same raw
tooth-endpoint calculation as THM-1250 gives

```text
|W_a|
 =[s_(a+1)(14n_a+1)-s_a(14n_(a+1)-1)]
    /[14s_as_(a+1)]
 >=gcd(s_a,s_(a+1))/[14s_as_(a+1)]
 =1/[14 lcm(s_a,s_(a+1))].                           (6)
```

Adjacent selected teeth have different speeds, since distinct teeth of one
comb are disjoint.

## 3. Coverage excess pays the full word

Let `C(t)` be the number of active fast danger combs at `t`.  Since the six
combs cover `G`,

```text
integral_G(C-1)
 =sum_i mu(G intersect D_(d_i))-|G|.                  (7)
```

Every `W_a` lies where `C>=2`, and the `W_a` are pairwise disjoint.  Hence

```text
sum_a |W_a|
 <=sum_i mu(G intersect D_(d_i))-|G|.                 (8)
```

THM-1166's exact singleton discrepancy is

```text
mu(G intersect D_d)<=|G|/7+6/(49d),
|G|=6/(7c).                                          (9)
```

Summing (9) and substituting in (8) yields

```text
sum_a |W_a| <=(6/49)(H-1/c).                         (10)
```

Combining (6) and (10) proves the full occurrence invoice

```text
H >=1/c+(49/6)sum_a |W_a|
  >=1/c+(7/12)sum_(a=1)^(N-1)
                    1/lcm(s_a,s_(a+1)).              (11)
```

Equation (11) has the same scale-covariant coefficient as THM-1250's
five-edge located tree, but now applies to every transition occurrence.
It improves THM-1250's multiplicity-averaged coefficient `7/36` by the exact
factor three.  No forest inequality is needed: spatial disjointness replaces
graphic acyclicity.

The same separation also sharpens the functional `H`-drift of THM-1199.
In the normalized slow-gap coordinate, its six-bin density satisfies
`f>=3/4`, while physical-to-normalized measure contributes `7c/6`.
Therefore every raw handoff has weighted mass at least

```text
(7c/6)(3/4)|W_a| >=c/[16 lcm(s_a,s_(a+1))].          (12)
```

Because the `W_a` are disjoint, their weighted masses lie directly in the
weighted coverage excess.  With THM-1198's twelve-piece one-comb envelope
`Pbar`, this upgrades the old tree functional to the full word:

```text
F_6:=sum_i Pbar(6d_i/(7c))-1
 >=(c/16)sum_(a=1)^(N-1)1/lcm(s_a,s_(a+1)).          (13)
```

Thus the literal functional form requested by the near-tiling/H-drift route
and the harmonic invoice (11) now charge exactly the same chronological
clock word.  They remain different upper envelopes and must not be added as
independent mass.

## 4. Scalar and word-skeleton consumers

THM-1198's private mass still forces

```text
n_i>=ceil(d_i/(7c)),
N=sum_i n_i.                                         (14)
```

Let `g0=gcd(d1,...,d6)`.  Since every transition has

```text
1/lcm(s_a,s_(a+1))>=g0/d6^2,
```

(11) gives the factor-three strengthened scalar recurrence

```text
H >=1/c+[7g0/(12d6^2)]
           [sum_i ceil(d_i/(7c))-1].                 (15)
```

There is also a sharp structural count.  Write `e_a={s_a,s_(a+1)}` for the
unordered transition edge.  An internal position is a nonbacktracking turn
exactly when `e_(a-1)!=e_a`, equivalently
`s_(a-1)!=s_(a+1)`.  The first edge sees at most two labels, and every edge
change introduces at most one new label.  Because the word visits all six
labels, it contains at least

```text
6-2=4                                                  (16)
```

nonbacktracking turns.  The star walk
`1,0,2,0,3,0,4,0,5` shows sharpness.  This is the smallest combinatorial
skeleton on which the two-sided wall germs of THM-1252 can attach.  Equation
(16) alone is not a metric credit—point-only private witnesses at coincident
walls are possible—but all adjacent handoffs on both sides are already paid
by (11).

## 5. Tournament and carrier audit

The runner-speed tournament is still transitive, with score histogram
`(0,1,2,3,4,5)`, no directed cycles, singleton SCCs, and one Hamiltonian
path.  It forgets multiplicity and every wall location.  Orienting every
missing pair after projecting the tooth word likewise destroys the fact used
in (8).

The faithful binary observable here is instead the chronological transition
edge `e_a`, with the wall order as tie Hamiltonian path.  Consecutive equal
edge symbols are immediate backtracks; changes are the four forced turns in
(14).  The proof-bearing carrier is the ordered interval family with both
endpoint sequences, not a tournament:

```text
(G; (alpha_a,beta_a,s_a,n_a)_(a=1)^N;
 W_a; exact lcm clock of every W_a).                  (17)
```

We challenged runners, individual teeth, tooth boundaries, unordered owner
pairs, overlap cells, edge runs, private components, and proof obligations as
vertices.  Projecting to owner pairs preserves the multiplicities but loses
the separation (4); projecting to a spanning tree preserves Hunter legality
but loses two thirds of the available word weight.  The ordered wall cells
preserve exactly the LRC predicate consumed by (7)--(11).

## 6. Verification and scope

The dependency-free referee checks `1,716` exact three-interval redundancy
rows, `3,997` minimal endpoint chains and `15,341` pairs of separated
handoffs, `460,607` positive gcd/lcm tooth numerators, all `112,320`
surjective adjacent-distinct owner words of lengths six through eight, the
sharp four-turn floor, `110,160` private occurrence rows, `357,201` scalar
debt rows, and every unweighted and functional coefficient above.  Normal and optimized outputs are
byte-identical.

The sorry-free Lean module kernel-checks the three-interval redundancy and
handoff-separation lemmas, the exact coverage-excess rearrangement, the full
lcm occurrence debt, its factor-three comparison and scalar consumer, and
common-dilate covariance.  Compact subcover extraction, deletion minimality,
and measure integration remain explicit paper topology providers.

Frozen artifact hashes are

```text
source         6e2a1866562d3f7097e8c730185fd709ebe12ec5f10aca0372f1ae5a8a3d51e4
output         df3f3b6b28725ba0c350854e1561699be537f5dd4f724ff305e2069247714e90
formalization  36d53a66ca2f4c5ea9fe8e8437a6685f45ee42512ea509e791255cd3baf1b279
```

THM-1253 does not exclude six-comb covers or prove LRC(14).  It removes the
Cayley one-third loss from the positioned owner-reuse route and shows that a
failed oriented-germ transport is automatically another fully charged seam
occurrence.  The remaining task is to lower-bound the total lcm word weight
sharply enough, or combine the full invoice with THM-1252/1254's two-sided
and coherent centered circuits to obtain a well-founded address descent.
