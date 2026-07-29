---
id: THM-2883
title: Ranked-apex hitting closure of all THM-741 roots
status: >
  PROVED + FINITE-EXACT + VERIFIED + INDEPENDENTLY PROOF-AUDITED.
  The exact ranked-apex recursion closes all 602 rank-impossible THM-741
  roots.  Together with the disjoint 584 top-four-positive and 816
  rank-feasible roots, this proves all 2002/2002 THM-741 roots.  THM-738
  supplies the small-speed chamber.  General LRC(14) remains OPEN.
source: codex-2026-07-29; independent audit root/lrc-rank-impossible-overlap-2026-07-29
depends_on:
  - THM-731-covering-middle-order-x-integral-autocorrelation-discrepancy
  - THM-732-disc-v-bernoulli-edge-pair-dedekind-form-exact-certificates-far-element-tail
  - THM-735-bonferroni-simultaneous-multi-peel-defeats-the-clustered-non-isolated-wall
  - THM-738-near-AP-three-slot-closure-all-1001-bodies-in-1-14
related:
  - THM-741-near-AP-four-slot-closure-all-2002-bodies-in-1-14
  - THM-2881-rank-impossible-pair-residual-closure-of-thirty-eight-thm741-roots
verification:
  - 04-computation/lrc14_thm741_all_root_pure_tail_top_four_atlas_codex_20260728.py
  - 05-knowledge/results/lrc14_thm741_all_root_pure_tail_top_four_atlas_codex_20260728.out
  - 04-computation/lrc14_thm741_all_rank_feasible_head_atlas_codex_20260728.py
  - 05-knowledge/results/lrc14_thm741_all_rank_feasible_head_atlas_codex_20260728.out
  - 04-computation/lrc14_thm741_residual_apex_hitting_closure_codex_20260729.py
  - 05-knowledge/results/lrc14_thm741_residual_apex_hitting_closure_codex_20260729.out
---

# THM-2883 — ranked-apex hitting closure of all THM-741 roots

## 1. Statement

For every nine-element set

```text
E subset {1,...,14}
```

and every four distinct positive integers `a,b,c,d` outside `E`, the
thirteen-speed family

```text
E union {a,b,c,d}
```

has positive lonely-time measure at threshold `1/14`.  Equivalently, all
`C(14,9)=2002` roots in the statement of THM-741 are uniformly closed.

This proves global **THM-741**.  It does **not** prove general LRC(14):
THM-741 is the near-AP chamber in which at least nine speeds already lie in
`{1,...,14}`.

## 2. Exact root partition and the tail estimate

For a root `E`, let

```text
G_E={t in R/Z : ||et||>1/14 for every e in E},
m_E=|G_E|,                    r_E=#components of G_E,
c_E(w)=|G_E intersect D_w|,   w>=15,
```

where `D_w={t:||wt||<=1/14}`.  Boundary conventions do not affect measure.
Put `s=99/70>sqrt(2)`.  The THM-731/732/735 chain gives the strict
one-comb estimate

```text
c_E(w)<u_E(w):=m_E/7+s r_E/(7w).                     (1)
```

The exact scan reconstructs all `2002` roots and seals their global top four
coverages.  It reproduces the disjoint partition

```text
584  top-four-positive roots,
816  rank-feasible roots,
602  rank-impossible roots.                          (2)
```

The first `584` roots satisfy

```text
m_E-q_1-q_2-q_3-q_4>0
```

and close by the union bound.  The second `816` are closed by the finite
ranked-head addendum to THM-741.  Its corrected semantic census is `63,265`
union-bound-dangerous quadruples; the stored computation checks `79,930`
literal endpoint controls because the exhaustive legacy `K<=10` layer keeps
`16,665` harmless extra controls.

It remains to close exactly the `602` rank-impossible roots.  Their
newline-joined lexicographic body digest is

```text
a4ad9f9b5a8ce16103450ac05684cd13ec33637e8e3737218a9186086cd639d4. (3)
```

## 3. The global top-ten apex gate

For each of the `602` roots, scan `15<=w<=600`, let `q_14` be the fourteenth
ranked coverage there, and put

```text
T_14(E)=s r_E/[7(q_14-m_E/7)].
```

Every row has `q_14>m_E/7`.  The verifier scans every additional integer
through `ceil(T_14)-1`; for all later integers, `(1)` is strictly below the
finite `q_14`.  Hence the first fourteen stored values are the global ranks
over every `w>=15`.  The largest tail-first value is

```text
989
```

at `E={1,2,8,9,10,11,12,13,14}`.

Write the global ranks as `q_1>=...>=q_14`, with speed as the deterministic
tie-breaker.  Exact arithmetic gives

```text
m_E-(q_11+q_12+q_13+q_14)>0                          (4)
```

on every root.  The minimum margin is

```text
67759/5045040
```

at `E={1,2,3,5,7,8,9,11,13}`.  Thus every quadruple whose individual
coverage sum can reach `m_E` contains at least one of the ten highest-ranked
speeds.  These ten speeds are the **apices**.

## 4. Recursion after fixing an apex

Fix an apex `a` and form the literal carrier

```text
C_(E,a)=G_E minus D_a,        m_(E,a)=|C_(E,a)|,
c_(E,a)(w)=|C_(E,a) intersect D_w|,  w!=a.            (5)
```

The verifier again seals the global top three carrier coverages by `(1)`.
Across the `6020` root-apex pairs the largest tail-first value is `1166`.
Let these global values be `p_1>=p_2>=p_3`.

There are three exact classes:

```text
direct:   m_(E,a)-p_1-p_2-p_3>0,                     4272;
rank2:    direct fails and
          R_2:=m_(E,a)-p_1-p_2>m_(E,a)/7,             1657;
failed:   neither inequality holds,                    91. (6)
```

Every direct apex closes all three remaining slots by the union bound.

For a rank2 apex, define the finite third-rank head

```text
H_(E,a)={w>=15, w!=a : c_(E,a)(w)>=R_2}.              (7)
```

If a triple is not contained in `(7)`, its two other coverages are at most
`p_1,p_2`, so its coverage sum is strictly below `m_(E,a)`.  It therefore
suffices to enumerate the triples in `(7)` whose exact individual sum is at
least `m_(E,a)`.

The largest head has `5711` speeds and the largest tail-first value is
`168888`, both at

```text
E={1,4,5,6,7,8,9,10,13},   apex a=23 (global apex rank 6).
```

Output-sensitive enumeration leaves exactly `29,622` literal triples.  Every
literal carrier after their three combs has positive measure.  The minimum is

```text
183/43792
```

at

```text
E={1,3,4,7,9,10,11,12,13},  a=16,  triple={17,19,23},
```

with six components; the minimum component count over all obligations is
four.

## 5. Weighted transversal recursion

The final step is the following reusable combinatorial lemma.

> **Weighted transversal lemma.**  Let `S` be a ground set with nonnegative
> weights `c:S->R`, let `m>0`, and consider `k`-element targets.  Let `A` be
> a finite apex set and let `C subset A` be the recursively closed apices:
> every target containing an element of `C` has already been closed on its
> literal carrier.  Put `F=A minus C`.  If the sum of the `k` largest
> weights on
>
> ```text
> F union (S minus A)=S minus C
> ```
>
> is strictly below `m`, then every target whose weight sum can reach `m`
> meets `C`; consequently every target is closed.

Indeed, a target meeting `C` uses its apex certificate.  A target avoiding
`C` lies in `S minus C`, so its weight sum is below `m`; the union bound
closes it.

Apply the lemma with `k=4`, `A` the global top ten, and `C` the `4272+1657`
direct or rank2-closed apices.  The `91` nominally failed apex branches form
`F`.

For each root, filter its global top fourteen to keep the failed apices among
ranks `1,...,10` and every rank `11,...,14`.  The first four entries of this
filtered list are exactly the four largest weights on `S minus C`.  The
reason ranks beyond fourteen cannot enter is important: ranks `11,12,13,14`
always remain in `S minus C`, so they already supply four allowed sentinels,
and every rank `>=15` has no larger coverage than rank fourteen.

The resulting complement margin is positive on all `602` roots.  Its minimum
is

```text
57684167/7467740280
```

at

```text
E={1,3,4,5,7,9,10,11,13},
allowed top four speeds=(23,17,53,72).                 (8)
```

Thus every individually dangerous quadruple hits a recursively closed apex.
The `91` nominal failures require no new carrier branch, and all `602`
rank-impossible roots are closed.

The lemma is stated independently of this rung because the same weighted
transversal recursion can organize higher-slot finite-rank arguments.  No
higher rung is claimed here.

## 6. Exact implementation and validity controls

The verifier is rational-exact.  For the unit-periodic tooth indicator it
uses the primitive

```text
P(x)=floor(x)/7+min({x},1/14)+max(0,{x}-13/14),
```

so for a carrier with components `[l_i,r_i]`,

```text
c(w)=sum_i [P(w r_i)-P(w l_i)]/w.                     (9)
```

The scalar `Fraction` implementation and an integer-vector implementation
of `(9)` agree on two values for every one of the `6020` apex carriers
(`12040` controls).  Every local apex subtraction agrees with the full
endpoint subtraction (`6020` controls).  For every rank2 carrier, the first
literal triple agrees between cached-pair and simultaneous subtraction, and
also with direct thirteen-comb reconstruction (`1657+1657` controls).

The vector arithmetic is overflow-audited.  The exhaustive root census gives

```text
max_E m_E=353/1176.
```

Every chosen apex has coverage above `m_E/7`, hence its carrier mass is

```text
m'<6m_E/7<=353/1372<2/7.                              (10)
```

With common endpoint denominator `D`, the nonnegative vector numerator is
`N=14Dw c(w)<4Dw`.  The runtime guard `4Dw<2^63` therefore protects the
axis sum.  It also gives `Dw<2^61`; endpoint products, the scaled primitive,
and the `14`-fold remainder terms all remain in signed `int64`.

An independent audit replayed `144` vector/scalar values, `12`
sequential/simultaneous/full-family carriers, and `17,600` synthetic
rank2-pruning instances containing `44,554` equality witnesses.  All passed.
The auditor also checked every strict inequality and the weighted-transversal
orientation.

Canonical SHA-256 ledgers are:

```text
rank-impossible apex gate
  af026722b290a2e578ded6f105722a104836a2cab31dee2272d836a7a7c14ecf
apex classification
  24a4b9b07424b23f6fe334f84c45ef45efa0e2c5ea695d402d1be88687bd87e0
failed-apex subatlas
  2f97c76ccbf5ae6f53c9c32b7a8366a3bd62a4b83957c63fb96283c9cbec3743
weighted-transversal gate
  f6c5f21c5a16651b69a681d2709a498b83cabf6bc845cce27dfd37ff0b8cc000
rank2 dangerous triples
  9df78b9841ba5fbac29b49965dd45498656d5093d377d16c8929699bcb308243
rank2 literal endpoint ledger
  dbef39aa9b7d681706321835dfd34143c1303392c58ea5788a41c24aa05780a4.
```

Normal and optimized full replays are byte-identical to the stored output:

```bash
python3 04-computation/lrc14_thm741_residual_apex_hitting_closure_codex_20260729.py \
  --workers 6 --literal-rank2
python3 -O 04-computation/lrc14_thm741_residual_apex_hitting_closure_codex_20260729.py \
  --workers 6 --literal-rank2
```

The LF-normalized script and output hashes are inserted at the promotion
checkpoint:

```text
script a5f3dcc1a23defea4b3dc067675d83141f1866022d6d01946617a97de69e5b0e
output 3d15be210ca1637a5a63942248758af48b90d2f469c7a4244320d8bcaec5b24c.
```

## 7. Consequence

Equations `(2)` and `(8)` give the exact pure-tail census

```text
584+816+602=2002/2002.                               (11)
```

If one of the four added speeds is at most `14`, the completed family has at
least ten speeds in `{1,...,14}` and is closed by THM-738.  Hence every
whole THM-741 root is closed, proving THM-741 as stated.

General LRC(14) remains **OPEN**. ∎
