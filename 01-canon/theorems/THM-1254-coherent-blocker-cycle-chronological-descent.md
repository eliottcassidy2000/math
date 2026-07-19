---
id: THM-1254
title: A COHERENT BLOCKER CYCLE UNCONDITIONALLY PAYS A CHRONOLOGICAL DRIFT INVOICE — binary speed descent selects the original or reflected tooth-word descent
status: PROVED (common-irredundant-subcover blocker selection; marked-cycle-tooth injection; forced speed-descent edge and binary relative digit; same-edge original/reflected position dichotomy; unconditional full chronological-drift invoice; carrier-at-most-1171 general-digit guardrail; dependency-free exact referee; sorry-free Lean arithmetic core). This is a coupling theorem, not six-comb noncoverage or LRC(14)
source: codex-2026-07-19 coherent-selection frontier audit
depends_on: [THM-1198, THM-1240, THM-1248, THM-1250]
related: [THM-841, THM-848, THM-850, THM-1127, THM-1156, HYP-7870]
script: 04-computation/lrc14_coherent_blocker_chronology_thm1254.py
output: 05-knowledge/results/lrc14_coherent_blocker_chronology_thm1254.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCCoherentBlockerChronology.lean
script_sha256: 6b970ba8da33a877ce37b1cf26bbfba84c5536b3a68f098c7daa7ea1eeef2ee3
output_sha256: 7986c5ae7828a063beee4a7b56077a9ef15e9e42632f43c24fb11e392d0d5a90
formalization_sha256: 7d824b37dca5d451d4e5e426a760fb40a5b580c19ea6e2ec1d72370e5d073f87
---

# THM-1254 — a coherent blocker cycle unconditionally pays one tooth-word drift invoice

## 1. Setup

Let

```text
G=G_k(c)=[(14k+1)/(14c),(14k+13)/(14c)],   0<=k<c,       (1)
c<d_1<...<d_6,
```

and suppose the six strict danger combs of the `d_i` cover `G`.  For each
label `i`, let

```text
Q_i=c+d_i,
t_i=P_i/Q_i in int(G)                                  (2)
```

be THM-1240's centered carrier spoke.  Thus `c` and `d_i` are safe at
`t_i`, while at least one of the other five fast labels is dangerous.

THM-1250 chooses a deletion-minimal finite cover of `G` by individual open
fast teeth.  Order its teeth by their left endpoints:

```text
I_1,...,I_N.                                           (3)
```

Their right endpoints are strictly increasing as well, and every
consecutive pair has a positive raw overlap in `int(G)`.  Every `t_i` lies
in at least one selected tooth.  Choose one such tooth and write

```text
A_i in {1,...,N},
t_i in I_(A_i),
b(i)=owner(I_(A_i)).                                   (4)
```

Since `d_i` is safe at `t_i`, one has `b(i)!=i`.  Therefore `b` is a
loopless functional graph and has a directed cycle.  The important new
point is that all blocker choices in this cycle have been made from **one
common irredundant tooth word**.  THM-1248 applies without alteration,
because it permits any strict blocker at each centered spoke.

## 2. The cycle marks distinct teeth in one linear order

Write a simple selected cycle of runner labels as

```text
v_0 -> v_1 -> ... -> v_(r-1) -> v_0,       2<=r<=6.    (5)
```

Let `T_s=I_(a_s)` be the tooth chosen at `t_(v_s)`.  Its owner is
`v_(s+1)`, with indices modulo `r`.  Let `N_s` be its integer tooth address,
so

```text
T_s is the N_s-th tooth of speed d_(v_(s+1)).          (6)
```

The marked teeth `T_s` are distinct: their owners `v_(s+1)` are distinct.
Consequently the cyclic word of integer positions

```text
a_0,a_1,...,a_(r-1),a_0                               (7)
```

has at least one strict descent and at least one strict ascent.

This elementary observation is the missing coupling.  The blocker cycle is
not an abstract relation superimposed on THM-1250's tree.  It is marked on
the same linearly ordered tooth carrier that produced every located Hunter
credit.

## 3. A cyclic position descent is an exact typed mixed circuit

Choose `s` with

```text
a_(s+1)<a_s.                                           (8)
```

The forward chronological subword from `T_(s+1)` to `T_s` runs from

```text
(s_0,n_0)=(d_(v_(s+2)),N_(s+1))
```

to

```text
(s_r,n_r)=(d_(v_(s+1)),N_s).                          (9)
```

Every consecutive step is an actual positive tooth overlap.  In the
canonical lift its reciprocal-centre drift is

```text
Delta=s_0/n_0-s_r/n_r>0.                              (10)
```

The next blocker edge of the same selected cycle is

```text
v_(s+1) -> v_(s+2).                                   (11)
```

At the centered spoke of `v_(s+1)`, its blocker is literally the initial
tooth `T_(s+1)` of (9).  Hence the closing centered edge in THM-1250's mixed
identity has

```text
P=P_(v_(s+1)),          N=N_(s+1)=n_0.               (12)
```

The preceding cycle edge `v_s->v_(s+1)` has THM-1248 relative digit

```text
delta_s=P_(v_(s+1))-(k+N_s).                         (13)
```

Therefore the formerly uncontrolled address-product term collapses exactly:

```text
P n_0-N n_r
 =N_(s+1)(P_(v_(s+1))-N_s)
 =N_(s+1)(k+delta_s).                                 (14)
```

Substituting (14) into THM-1250's mixed-circuit identity gives the sharper
common-stalk form

```text
R=P s_0-N s_r
 =N_(s+1)N_s Delta+s_0(k+delta_s).                   (15)
```

All addresses in (14)--(15) are positive.  Hence

```text
k+delta_s>=0
  ==> R>=N_(s+1)N_s Delta>0,                          (16)
```

while the other branch is the strict addressed-product descent

```text
P n_0<N n_r  iff  k+delta_s<0.                        (17)
```

THM-1248 gives `delta_s>=-586`, so (17) can occur only when

```text
k<=585.                                               (18)
```

Thus the affine sign-reversal guardrail in THM-1250 is not generic once the
blockers are selected coherently from the chronological word.  On every
marked position descent it is either a full positive drift invoice or a
bounded lower-gap address.

## 4. A speed descent makes the invoice unconditional

The vertices of (5) carry distinct numerical speeds.  Choose `s` so that

```text
d_(v_(s+1))<d_(v_s).                                  (19)
```

Such an edge exists: take a maximum-speed vertex of the cycle and follow its
outgoing edge.  In particular

```text
d_(v_(s+1))<=c+d_(v_s)=Q_(v_s),
```

so THM-1248's target-at-most-source-clock lemma applies to this **same**
edge and gives

```text
delta_s in {0,1}.                                     (20)
```

The two teeth `T_s,T_(s+1)` are distinct.  There are exactly two cases.

* If `a_(s+1)<a_s`, Section 3 applies in the original word.  Since `k>=0`
  and `delta_s>=0`, its factor satisfies

  ```text
  k+delta_s>=0.                                       (21)
  ```

* If `a_(s+1)>a_s`, reflect the circle by `t |-> 1-t`.  The reflected word
  reverses the marked positions, so the **same edge `s`** is now a position
  descent and Section 3 applies there.  The carrier-gap and address laws are

  ```text
  k'=c-k-1,
  P_i'=Q_i-P_i,
  N_s'=d_(v_(s+1))-N_s.                               (22)
  ```

  If `M_s=k+N_s`, then

  ```text
  M_s'=Q_(v_(s+1))-1-M_s,
  delta_s'=1-delta_s,                                 (23)
  k'+delta_s'=c-k-delta_s>=0,                         (24)
  ```

  because `c-k>=1` and `delta_s<=1`.

Thus one of the two located versions of the same typed mixed circuit always
pays the full chronological drift:

```text
R>=n_0 n_r Delta>0.                                   (25)
```

There is no bounded-carrier fallback in this conclusion.

> **Unconditional coherent-cycle invoice theorem.**  Every hypothetical
> six-comb cover of a complete `c`-safe gap admits a blocker cycle coherently
> marked on one irredundant tooth word.  A speed-descent edge of that cycle
> yields, according only to the order of its two consecutive marked teeth,
> either an original or reflected located mixed circuit satisfying (25).

The conclusion still does not say that this invoice exceeds the available
Hunter mass; it supplies the previously missing unconditional sign and
location coupling.

## 5. The `1171` arbitrary-digit guardrail remains valid

For completeness, the older scale split remains useful when one deliberately
chooses arbitrary position-descent and position-ascent edges rather than the
binary speed-descent edge above.  Reflecting as in (22)--(23), an original
descent `s` and an original ascent `t` have respective factors

```text
k+delta_s,                    c-k-delta_t.             (26)
```

If either factor is nonnegative, the corresponding circuit pays its full
drift.  If both fail, the full THM-1248 bank

```text
-586<=delta_s,delta_t<=587                              (27)
```

and integrality imply

```text
k<=585,               c-k<=586,               c<=1171. (28)
```

This is a secondary general-digit lemma, not an alternative in the primary
theorem.

## 6. What the older viewpoints contribute

The theorem is the common operation-level lesson of several older probes.

* **Hamiltonian-path `H`-drift.**  THM-848 shows that averaged scalar drift
  is not a transport state; one needs the target-coupled directional orbit.
  Here the blocker cycle alone is the averaged shadow, while
  `(T_s,a_s,N_s)` couples each direction to its actual target tooth.  Formula
  (15) is the resulting exact functional form.
* **Toothpick and Farey ladders.**  THM-841 refutes count-level toothpick
  self-similarity but proves that an ordered Farey-neighbour pair is an
  operation-congruent address.  THM-1127 likewise keeps the torus strip and
  phase.  The common tooth word here is the analogous recursive carrier:
  order and addresses update, while a total tooth count does not.
* **Fano/`chi_7`.**  THM-1156 and THM-1247 organize seam obligations but do
  not place them.  No Fano symmetry is needed in (14).  Exact owner and tooth
  placement repairs the liar fibres left by the contracted `q=15` stalk.
* **The `j=4` flood tail.**  THM-850 shows that Heawood currents recover a
  root-edge field, while rerooting destroys every nonconstant edge-only
  scalar.  The repair is again a common operation carrier: the legal
  extension word together with cyclic position.  Equations (4)--(14) are the
  slow-gap version of that repair.
* **Kakeya needles.**  Direction alone is insufficient.  The marked tooth
  retains direction, affine offset, owner, strict side, and chronology.  The
  proof uses no planar Kakeya dimension theorem.

These are not independent analogies added after the proof.  They explain why
the one new choice in (4)—select every blocker from the same minimal
subcover—preserves exactly the coordinate all older quotients discarded.

## 7. Tournament and alternate-carrier audit

On the cycle sources, take the pair observable

```text
C(i,j)=A_j-A_i,                                        (29)
```

orient `i->j` when `C(i,j)>0`, and break impossible marked-cycle ties by
speed.  The tie Hamiltonian path is increasing marked-tooth position.  This
position tournament is transitive, with score histogram
`(0,1,...,r-1)`, no directed triangles, singleton SCCs, and one Hamiltonian
path.  The blocker relation (5) is a directed cycle.  Overlaying the two
typed relations forces a backward blocker edge and hence the mixed circuit
of Section 3.  Collapsing the types to one tournament would erase the proof.

We challenged runners, centered spokes, selected teeth, tooth boundaries,
overlap intervals, relative digits, cycle obligations, and marked positions
as vertices.  The faithful carrier is

```text
(one irredundant tooth word;
 six marked centered spokes and their selected teeth;
 blocker-cycle incidence;
 owner, absolute tooth address, relative digit, and reflected address).   (30)
```

It preserves literal cover truth at every sampled phase, chronological
subpaths, the exact `d_i`/`c+d_i` interface, and reflection.  It destroys
unselected alternative blockers, global safe-set measure away from the word,
and any claim that the resulting positive invoice is already large enough.

## 8. Verification and scope

The dependency-free exact referee checks every cyclic position word of
length `2..6`; every pair of speed-order and marked-position-order words of
length `2..6`; both binary digits on the maximum-speed vertex's forced
descent edge; the reflection formulas; the product-gap collapse; and all
`1,378,276` pairs of allowed relative digits in the secondary two-sided
carrier bound.  The Lean module kernel-checks (14), the simplified mixed
identity (15), the reflection laws, the same-edge original/reflected binary
invoice dichotomy, the one-sided `585` cutoff, and the secondary two-sided
`1171` guardrail.  Compact minimal-subcover extraction, selection of a tooth
covering each centered phase, and selection of the maximum-speed cycle vertex
remain explicit paper topology providers.  There are no proof placeholders
or `native_decide` calls.

Frozen hashes are

```text
script         6b970ba8da33a877ce37b1cf26bbfba84c5536b3a68f098c7daa7ea1eeef2ee3
output         7986c5ae7828a063beee4a7b56077a9ef15e9e42632f43c24fb11e392d0d5a90
formalization  7d824b37dca5d451d4e5e426a760fb40a5b580c19ea6e2ec1d72370e5d073f87
```

THM-1254 does not prove six-comb noncoverage or LRC(14).  It makes the
THM-1248 blocker cycle and THM-1250 chronological tree coherent by choice,
then uses one forced binary speed-descent edge to eliminate the free affine
sign term unconditionally: the marked-tooth order selects either the original
or reflected full drift invoice.  The bound `c<=1171` survives only as a
general-digit guardrail for noncanonical edge choices.
