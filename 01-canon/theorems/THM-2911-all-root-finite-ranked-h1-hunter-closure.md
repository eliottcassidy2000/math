---
id: THM-2911
title: "All-root finite ranked H1 Hunter closure"
status: >
  PROVED + FINITE-EXACT + VERIFIED.  Of the 6,180 scalar-hard marked
  suffixes satisfying the ranked-r4/H1 entry inequality, all 5,999 with
  exact discrepancy cutoff at most 15,000 close: 5,975 by an ordered
  literal-child pivot and 24 by maximum-spanning-tree Hunter repair.
  Branchwise composition with THM-2905 and THM-2904 closes 138 route
  roots.  Exact set difference against the proved 88-root baseline adds
  93 roots, taking the proved union to 181 and the residual to 3,251.
source: root-2026-07-29
depends_on:
  - THM-735-bonferroni-simultaneous-multi-peel-defeats-the-clustered-non-isolated-wall
  - THM-2892-eight-body-five-slot-heavy-triangle-closure
  - THM-2893-complement-cap-finite-core-flag-lemma
  - THM-2897-partition-cap-tropical-convolution-and-alternating-pair-ladder
  - THM-2899-all-root-ranked-suffix-scalar-census
  - THM-2901-all-hard-exact-global-pair-cap-and-route-partition
  - THM-2903-one-hard-actual-h3-link-and-bad-triangle-closure
  - THM-2904-all-root-ranked-h1-hunter-closure
  - THM-2905-all-hard-hunter-star-envelope-census
related:
  - THM-2902-omission-six-ranked-h1-depth-one-closure
  - THM-2907-paircap-exception-h4-heavy-link-child-census
verification:
  - 04-computation/lrc14_j6_ranked_h1_thm2911/verify.py
  - 04-computation/lrc14_j6_ranked_h1_thm2911/scout.py
  - 04-computation/lrc14_j6_ranked_h1_thm2911/hunter_two_star_exceptions.py
  - 04-computation/lrc14_j6_ranked_h1_thm2911/hunter_all_24_star_survivors.py
  - 05-knowledge/results/lrc14_j6_ranked_h1_thm2911/locked.out
  - 05-knowledge/results/lrc14_j6_ranked_h1_thm2911/ordinary_optimized_replay.out
---

# THM-2911 -- all-root finite ranked H1 Hunter closure

**PROVED + FINITE-EXACT + VERIFIED.**

## 1. Statement

Consider the `14,806` scalar-hard marked suffixes in the exact all-root
`j=6` ledger of THM-2899.  Exactly `6,180` satisfy the strict ranked
four-singleton inequality

```text
R_4=q_1+q_2+q_3+q_4<6h/7.                                (1)
```

For each such row, let `N` be the exact discrepancy cutoff in `(6)` below.
Every one of the `5,999` rows with

```text
N<=15,000                                                   (2)
```

admits no allowed five-cover of its literal carrier.  The remaining `181`
rows have larger `H_1` cutoffs and are not claimed by this finite-`H_1`
certificate.

The finite certificate closes every hard branch on `127` bodies; adjoining
THM-2899's five scalar-terminal bodies gives `132` finite-`H_1` route
roots.  On hard branch keys, its exact composition with THM-2905's `G_5` certificate and
THM-2904's hostile-centre pivot certificate is

```text
finite r4/H1 branches                                  5,999
G_5 branches                                           2,964
hostile-centre pivot branches                            279
three-route branch union                               6,118

finite r4/H1 route roots (127 hard + 5 scalar)           132
three-route terminal roots                               138. (3)
```

Set difference against the proved `88`-root union through THM-2904 gives

```text
three-route roots intersecting the proved baseline        45
new roots                                                  93
proved union                                              181
official residual                           3,432-181=3,251. (4)
```

This is a strict finite-exact improvement, not a proof of LRC(14).

## 2. Why the rank-four gate makes the problem finite

Fix one marked suffix.  Its theorem-bearing data are

```text
C       literal post-apex carrier,
h=mu(C) carrier mass,
P       actual excluded gate prefix,
r       number of interval components of C,
c(w)    mu(C intersect D_w), for allowed w not in P,
q_i     globally sealed singleton order statistics.       (5)
```

Every allowed four-set has union at most `R_4`, by subadditivity and rank
selection.  Apply the `s=1`, `k=5` row of THM-2893/2897.  Under `(1)`,
every label of a hypothetical five-cover belongs to

```text
H_1={w not in P:c(w)>=h-R_4}.
```

This conclusion is stronger than merely forcing one distinguished centre:
all five labels lie in the same finite core.

Put

```text
epsilon=h-R_4-h/7=6h/7-R_4>0,
gamma=(99/70)r/7.
```

The strict interval discrepancy estimate

```text
c(w)<h/7+gamma/w
```

therefore seals the whole core at

```text
N=ceil(gamma/epsilon)-1
 =ceil((99/70)r/[7(6h/7-R_4)])-1.                        (6)
```

Every allowed integer through `N` is scanned, with the actual prefix `P`
excluded.  Equality at the `H_1` threshold is retained.  Thus `(2)` is a
complete finite universe, not a sampling horizon.

## 3. Ordered literal-child certificate

Order `H_1` by decreasing parent coverage `c(w)`, breaking ties by the
integer label.  Every hypothetical five-cover has a unique earliest core
label `x`; its other four labels lie in the strict later suffix

```text
H_1^(>x)={y in H_1:x precedes y}.                         (7)
```

Let

```text
C_x=C minus D_x,             h_x=mu(C_x),
d_x(y)=mu(C_x intersect D_y).
```

If fewer than four later labels remain, the pivot branch is vacuous.
Otherwise, if the four largest exact values `d_x(y)` on the later suffix
have sum strictly below `h_x`, then no four of those labels cover `C_x`.
This closes every five-set assigned to `x`.

The verifier first uses the inherited parent coverages as safe upper
bounds.  Only when that bound is inconclusive does it compute the literal
child values `d_x(y)`.  This changes cost, not logic: every pivot is tested
on its actual suffix, and no earlier label is silently reintroduced.

This ordered proof allocation closes `5,975` of the `5,999` finite rows.
Its smallest positive strict margin is

```text
1/3,166,020.                                               (8)
```

The order in `(7)` is an acyclic gauge on proof branches.  It does not
orient the symmetric overlap graph and is not a tournament on runners.

## 4. The 24 nonstar Hunter repairs

The ordered scalar test leaves exactly `24` rows.  Each has one unresolved
pivot, at index zero.  For a candidate five-set `Q` containing that pivot,
put

```text
A_w=C intersect D_w,                  w in Q,
i(x,y)=mu(A_x intersect A_y).
```

For every spanning tree `T` on `Q`, Hunter's tree inequality gives

```text
mu(union_(w in Q) A_w)
 <=sum_(w in Q)c(w)-sum_(xy in T)i(x,y).                  (9)
```

Indeed, root `T` and add its vertices parent-first.  When a child set is
added, its intersection with the preceding union contains its intersection
with its parent.  Summing the resulting one-step union bounds proves `(9)`.
Consequently the strongest pairwise tree invoice uses a
maximum-weight spanning tree:

```text
Psi_5(Q)=sum_(w in Q)c(w)
         -max_T sum_(xy in T)i(x,y).                      (10)
```

Exact star pruning reduces the enormous nominal five-set universe to just
`54` star-hostile sets.  All pair intersections in those sets are computed
from literal residual subtraction.  Deterministic Kruskal maximizes the
tree credit in `(10)`, and every one of the `54` sets satisfies

```text
Psi_5(Q)<h.
```

The smallest Hunter margin is

```text
1,799,771/75,716,368.                                     (11)
```

Thus all `24` remaining finite rows close.  A separate exact
branch-and-bound in the shard computation finds no witness and serves as a
hostile control; the theorem uses the explicit Hunter repair above.

## 5. Exact finite census

The complete census is

```text
seven-body roots                                        3,432
marked gate branches                                   41,415
scalar-hard branches                                   14,806
rank-four H1 eligible                                   6,180
finite cutoff N<=15,000                                 5,999
ordered literal-child closures                          5,975
maximum-tree Hunter closures                               24
witnesses                                                    0
deferred N>15,000 rows                                     181. (12)
```

The exact `H_1` core-size quantiles on the finite rows are

```text
minimum 5, p25 8, median 16, p75 36,
p90 83, p95 141, p99 283, maximum 444.                  (13)
```

The maximum core would have

```text
binom(444,5)=140,578,526,088
```

raw five-subsets.  The ordered pivot and tree sidecar are therefore the
mechanism that makes the finite theorem practical; the cutoff alone is not.

## 6. Three-route recomposition

All branch keys retain

```text
(body, gate size, rank, apex, excluded prefix).           (14)
```

On the union of the three branch certificates, the exact nonempty atomic
membership classes are

```text
G_5 only                                      55
finite H_1 only                            2,875
finite H_1 and G_5                         2,909
hostile-centre pivot only                     64
finite H_1 and hostile-centre pivot          215
                                               ---
                                             6,118.       (15)
```

There is no `G_5`--pivot overlap because THM-2904 evaluates the
`G_5`-surviving rows.  The pair-cap exceptions occur in none of the three
classes.

THM-2899 closes the scalar-soft marked suffixes and THM-2892 supplies the
eight-body chamber.  Requiring every scalar-hard suffix of a body to lie in
the union `(15)` gives `133` hard-terminal bodies; adding THM-2899's five
vacuous scalar-terminal bodies gives the `138` route roots in `(3)`.

The finite `H_1` certificate alone, unioned with the proved `88`-root
baseline, already gives `180` roots.  The full three-route join contributes
one further root that no route closes alone:

```text
E=(1,2,3,5,6,8,13).                                      (16)
```

Its exact two-branch anatomy is

```text
rank 1, apex 22, prefix (22):
  finite H_1 cutoff 6,930;
  G_5 margin -174401/12612600;

rank 2, apex 21, prefix (22,21):
  hostile-centre pivot closure;
  no finite-r4/H1 entry;
  G_5 margin -383/254800.                                (17)
```

Thus `(16)` is a genuine compositional gain, rather than a root counted
twice under two names.  Exact artifact-derived set difference gives `(4)`.

## 7. Relation to THM-2904

The two `H_1` constructions are complementary and must not be conflated.

- THM-2904 starts after `G_5` failure.  Concavity determines a hostile
  singleton threshold `lambda>h/7`; only the maximum singleton of a cover
  is forced into its small core, of size at most `13`.  The inherited pair
  cap controls four leaves.
- This theorem starts from the stronger rank-four inequality `(1)`.
  It forces all five labels into a possibly much larger core, then uses
  exact literal children and nonstar overlap.  Its cutoff is deliberately
  limited to `15,000`.

The source, target, and destroyed information differ: THM-2904 preserves a
small possible-centre set but forgets nonstar child geometry; this theorem
preserves the complete all-label finite core and literal suffixes but
defers long discrepancy tails.  Their branchwise join in `(16)` is exactly
the connection that neither terminal-root list reveals.

## 8. Verification

The locked verifier checks all source and output hashes, the complete
32-shard hash battery and aggregate `3,432`-body count, core histograms,
exception anatomy, strict margins, Hunter repairs, and explicit root lists.
It invokes the hardened THM-2905 parser and independently reparses the
hash-pinned THM-2904 ledger.  It retains the full key `(14)`, reconstructs
the proved `88`-root baseline from canonical artifacts, and recomputes every
set union and difference in `(15)`--`(17)`.

All 32 shards were freshly replayed under ordinary Python with their
recorded worker schedule and were byte-identical to the locked optimized
artifacts.  The final verifier was run under both ordinary and optimized
Python with byte-identical output.  All guards are explicit; no
theorem-bearing check disappears under `python3 -O`.  Repository-text
hashes normalize CRLF to LF and reject lone carriage returns, including
the transitive Hunter/G5 audit imports.

Canonical directory layout:

```text
04-computation/lrc14_j6_ranked_h1_thm2911/
  scout.py
  hunter_two_star_exceptions.py
  hunter_all_24_star_survivors.py
  verify.py

05-knowledge/results/lrc14_j6_ranked_h1_thm2911/
  all_s00_of32.out ... all_s31_of32.out
  hunter_all_24_star_survivors.out
  ordinary_optimized_replay.out
  thm2904_hostile_centre.ledger.out
  locked.out
```

The explicit `132`-root finite-`H_1` list and `93`-root additive list are
printed in `locked.out`.  The source, output, semantic, branch-key, and
root-set hashes are pinned by `verify.py` and printed in `locked.out`.

Reproduction of the compact theorem check:

```bash
python3 04-computation/lrc14_j6_ranked_h1_thm2911/verify.py
python3 -O 04-computation/lrc14_j6_ranked_h1_thm2911/verify.py
```

Replaying the 32 interval shards is intentionally a separate, expensive
operation; the exact commands and worker schedule are recorded in the
package README and replay attestation.

## 9. Equality and scope boundary

All closure inequalities are strict.  Equality at the finite-core
membership threshold is included, not discarded.  The `181` rows with
`N>15,000` are deferred from the finite-`H_1` route solely by the stated
resource cutoff.  They are not cover witnesses, no extrapolation from the
finite rows is claimed, and other branch routes may close some of them.

The theorem does not close THM-2901's `52` pair-cap exceptions, does not
settle every branch failing the three sufficient certificates, and does
not prove LRC(14). ∎
