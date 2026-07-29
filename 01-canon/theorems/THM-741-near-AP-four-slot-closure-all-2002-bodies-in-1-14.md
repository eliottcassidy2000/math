---
id: THM-741
title: NEAR-AP FOUR-SLOT CLOSURE — every 13-speed family with AT LEAST 9 speeds in {1,…,14} satisfies LRC(14). Equivalently, for EVERY 9-element body E ⊆ {1,…,14} (all C(14,9)=2002) and all v₁<v₂<v₃<v₄ not in E, {E,v₁..v₄} is lonely. Proof = the THM-735 Bonferroni tree at j=4: legs J4 (one inequality, all four ≥ V₁(E)) / J3 (per-v₁ exact bodies) / J2 (per-(v₁,v₂)) / J1 (per-(v₁,v₂,v₃) tail) / bottom (exact-ℚ sweeps of covering quadruples via lcm-multiples) — with PROVED P1/P2 LEMMA-SKIPS at every level (subtrees where the next Bonferroni threshold already fires from the parent's exact data close without computing the child body; sound because P1/P2 are one-level bounds off exact data)
status: CLAIMED globally; 595/2002 root bodies now proved uniformly.  The live direct resume ledger reached 290/2002 clean bodies at the last pull but is not harvested.  Exact addenda prove all 21/21 whole flood bodies.  An all-root coverage atlas proves every pure-tail extension 15<=a<b<c<d for exactly 584/2002 roots; THM-738 closes their complementary small-speed extensions, making 584 whole-root closures.  Exactly ten flood roots lie in that atlas set, so the eleven separately repaired flood exceptions raise the proved union to 595 roots.  These proofs use literal root carriers plus the THM-735(ii) covariance cap, not Fano transport.  The remaining 1407-root discharge and global LRC(14) stay open.
source: kind-pasteur-2026-07-13-S128 (cont.5); exact flood/completed-family addenda codex-2026-07-15/18; all-root atlas codex-2026-07-28
depends_on:
  - THM-735   # the simultaneous multi-peel lemma (j=4,3,2,1 legs) + P1/P2 peel lemmas (THM-733)
  - THM-731 / THM-732 / THM-366
related:
  - THM-734 (j=2, 364 bodies), THM-738 (j=3, 1001 bodies) — this is rung j=4 of the same ladder
  - THM-737/739/740 (opus: pack-clock / cluster-clock / two-cluster — tiling the j≥7 seam from the coherent side)
  - klein-S295 (HYP-6580): the LRC(13) localization — the residual is PLACEMENT (the good set reaching
    the middle [1/14,13/14]), the AP-cluster the unique end-trapped case, and the uniform version a
    stability theorem around the AP. In that frame the ladder is the QUANTITATIVE middle-reach from
    the bounded side: rung j proves middle-reach (indeed L>0) for every covering family with ≥(13−j)
    in-window speeds — the "practical closure" klein's consolidation names, grown one rung per run
  - MISTAKE-122 (j≤6), MISTAKE-141 (exact thresholds), HYP-6540 (calibration)
verification:
  - 04-computation/lrc14_j4_flood_portable_pruner_codex_S14.py
  - 05-knowledge/results/lrc14_j4_flood_56_exact_codex_S16.out
  - 05-knowledge/results/lrc14_j4_flood_57_exact_codex_S14.out
  - 05-knowledge/results/lrc14_j4_flood_67_exact_codex_S15.out
  - 04-computation/lrc14_j4_flood_reroot_shadow_codex_S16.py
  - 05-knowledge/results/lrc14_j4_flood_reroot_shadow_codex_S16.out
  - 04-computation/lrc14_j4_next_anchor_shadow_frontier_codex_S16.py
  - 05-knowledge/results/lrc14_j4_next_anchor_shadow_frontier_codex_S16.out
  - 04-computation/lrc14_j4_34_ancestral_carrier_frontier_probe_codex_S16.py
  - 05-knowledge/results/lrc14_j4_34_ancestral_carrier_frontier_probe_codex_S16.out
  - 04-computation/lrc14_j4_three_small_frontier_workload_codex_S16.py
  - 05-knowledge/results/lrc14_j4_three_small_frontier_workload_codex_S16.out
  - 05-knowledge/results/lrc14_j4_three_small_frontier_workload_codex_S16.jsonl
  - 04-computation/lrc14_j4_34_a15_pure_tail_exact_codex_20260717.py
  - 05-knowledge/results/lrc14_j4_34_a15_pure_tail_exact_codex_20260717.out
  - 04-computation/lrc14_j4_34_a16_pure_tail_exact_codex_20260717.py
  - 05-knowledge/results/lrc14_j4_34_a16_pure_tail_exact_codex_20260717.out
  - 04-computation/lrc14_j4_34_a17_pure_tail_exact_codex_20260717.py
  - 05-knowledge/results/lrc14_j4_34_a17_pure_tail_exact_codex_20260717.out
  - 04-computation/lrc14_j4_34_a18_pure_tail_exact_codex_20260717.py
  - 05-knowledge/results/lrc14_j4_34_a18_pure_tail_exact_codex_20260717.out
  - 04-computation/lrc14_j4_34_a19_pure_tail_exact_codex_20260718.py
  - 05-knowledge/results/lrc14_j4_34_a19_pure_tail_exact_codex_20260718.out
  - 04-computation/lrc14_j4_34_exact_top_four_coverage_codex_20260728.py
  - 05-knowledge/results/lrc14_j4_34_exact_top_four_coverage_codex_20260728.out
  - 04-computation/lrc14_j4_37_top_five_exception_codex_20260728.py
  - 05-knowledge/results/lrc14_j4_37_top_five_exception_codex_20260728.out
  - 04-computation/lrc14_j4_six_more_top_four_coverage_codex_20260728.py
  - 05-knowledge/results/lrc14_j4_six_more_top_four_coverage_codex_20260728.out
  - 04-computation/lrc14_j4_four_ranked_exception_coverage_codex_20260728.py
  - 05-knowledge/results/lrc14_j4_four_ranked_exception_coverage_codex_20260728.out
  - 04-computation/lrc14_j4_label1_head_tail_partition_codex_20260728.py
  - 05-knowledge/results/lrc14_j4_label1_head_tail_partition_codex_20260728.out
  - 04-computation/lrc14_thm741_all_root_pure_tail_top_four_atlas_codex_20260728.py
  - 05-knowledge/results/lrc14_thm741_all_root_pure_tail_top_four_atlas_codex_20260728.out
  - 04-computation/lrc14_thm741_sharded_resume_runner_codex_20260717.py
  - 05-knowledge/results/lrc14_thm741_sharded_resume_runner_codex_20260717.out
---

# THM-741 — the near-AP four-slot closure (2002 bodies, overnight run)

Statement and method as in the title. New at this rung: (a) the general-Q bottom now has THREE
slots available to cover Q(E) before v₄; (b) the **v₀-upper-bound node closure** (the design that
survived testing — see below).

## Design findings (cont.5, recorded per MISTAKE discipline)

1. **The naive per-level "lemma-skips" are algebraically VACUOUS.** The intended skip — "at index v,
   does the NEXT leg's Bonferroni threshold already fire for all remaining slots ≥ v+1, using P1/P2
   bounds?" — can never fire inside the loop ranges: the skip point sits at ~k·√2·r/m-type slopes
   ABOVE the very thresholds (V₂, V₃, v₀) that bound the loops (e.g. level 1: fires only when
   3·√2·v·m < (24/7)·m·v — impossible since 3√2 > 24/7). First smoke test confirmed 0 firings across
   758,921 nodes. Lesson: a subtree-skip must compare the child threshold against the CURRENT index,
   and the child threshold's slope (≈0.28–0.57 per unit) only catches the index beyond the parent
   threshold — i.e. outside the loop. Removed.
2. **The optimization that works: bound v₀ from ABOVE, not the loop from below.** With E₂ exact,
   P1/P2 give v₀(E₃) ≤ v₀ᵘ = √2·(v₃m₂+(15/7)r₂)/(6·((6/7)m₂−8r₂/(49v₃))); "v₄ > v₀ᵘ ⟹ tail" is
   sound, and the bottom candidate set (covering v₄ ∈ (v₃, v₀ᵘ], via lcm(Qb)-multiples) is computable
   WITHOUT the exact E₃ body. The exact level-3 subtract happens only at nodes with candidates;
   sweeping covering v₄ ∈ (v₀, v₀ᵘ] is redundant-but-sound. Smoke body: 758,921 E₃ subtracts → ~10³
   (the candidate nodes), the dominant cost replaced by ~4 Fraction ops per node.

## Evidence log

- [x] design validated on smoke body (pre-fix 278.6 s, post-fix probe below); regression subtree
      (E={1..9}, v₁=10) checked as (V₂, E₂-count, v₃-nodes, sweeps) = (154, 143, 7537, ≥27) — the
      ≥ is because v₀ᵘ ≥ v₀ sweeps a superset of THM-738's 27
- [ ] probe calibration (4 sample bodies incl. a true flood) + extrapolation
- [x] live resume run launched; `290/2002` bodies clean at the last pull, ledger not yet harvested here
- [ ] all 2002 bodies clean; tight census; verdict

## Portable fixed-`E2` pruning addendum (codex-S14/S16): three flood bodies closed

The original driver hard-codes a Windows scratch path.  Its earlier reported
`171/2002` body ledger is not present in this checkout; the live runner later
advanced to `290/2002`, likewise not yet harvested here.  The portable
companion above hash-guards the original exact interval kernel, accepts a
platform-neutral optional JSONL state path, and adds the following exact
screen before constructing the fragmented `E3` interval list.

Let `G=G(E2)` have `r` components and measure `m`, put `a=v3`, and write
`m_a=|G\D_a|`.  For every `b>a`, THM-732's discrepancy bound gives

```text
|G intersect D_b| <= m/7 + sqrt(2) r/(7b).
```

Therefore

```text
|G \ (D_a union D_b)|
 >= m_a - m/7 - sqrt(2) r/(7b).                         (A1)
```

Using the existing rational majorant `S2=99/70>sqrt(2)`, all integer
`b>S2*r/(7(m_a-m/7))` close whenever the denominator is positive.  Before
computing `m_a`, P2 supplies the monotone lower bound

```text
m_a-m/7 >= 5m/7 - 8r/(49a),                             (A2)
```

so one binary search locates the first `a` after which every `b>a` is
certified.  This retains the fixed `E2` Kakeya-comb carrier: the expensive
`G(E2 union {a})` endpoint list is built only for the finite survivors of
(A1)--(A2).  The proof is a union bound plus THM-732/P2; no Fano symmetry or
body transport is inferred.

For the flood core

```text
E_(5,7)={5,7,8,9,10,11,12,13,14},
```

the exact tree has

```text
root: r=16, m=581453/2522520, V1=131;
E1/E2/v3 nodes:                         121 / 10,929 / 525,362;
P2-preclosed v3 nodes:                                    181,445;
exact-m3 nodes / closed without E3:              343,917 / 339,348;
residual E3 nodes / exact bottom sweeps:             4,569 / 28,847;
fallback nodes:                                                   0.
```

All `28,847` bottom measures are positive.  The smallest is

```text
7/858
```

at `{1,2,3,4,5,7,8,9,10,11,12,13,14}`.  Hence **the `(5,7)` flood body is
closed uniformly over all four added speeds**.  This is one of the 21 flood
bodies, not a quotient representative: the remaining 20 require their own
exact runs.  THM-741 therefore remains `CLAIMED`; this addendum establishes
one rigorous flood row and a portable pruning invariant for the rest.

The same unmodified exact driver closes a second, genuinely different row,

```text
E_(6,7)={6,7,8,9,10,11,12,13,14}.
```

Here the exact tree has

```text
root: r=20, m=24289/105105, V1=164;
E1/E2/v3 nodes:                       154 / 17,296 / 1,028,036;
P2-preclosed v3 nodes:                                      354,361;
exact-m3 nodes / closed without E3:                673,675 / 664,115;
residual E3 nodes / exact bottom sweeps:               9,560 / 73,323;
fallback nodes:                                                   0.
```

All `73,323` bottom measures are positive.  The smallest is `97/4004`, at
`{1,2,3,5,6,7,8,9,10,11,12,13,14}`.  Hence the `(6,7)` flood body is also
closed uniformly over all four added speeds.  This is a literal second edge,
not a Fano or tournament transport of `(5,7)`: the eight-dimensional Heawood
cycle sector and the metric interval data forbid that quotient.  Thus two of
the 21 flood bodies were exact at that checkpoint.

The same driver now closes the third edge of the small-label triangle,

```text
E_(5,6)={5,6,8,9,10,11,12,13,14}.
```

Its exact tree has

```text
root: r=24, m=563009/2522520, V1=203;
E1/E2/v3 nodes:                       193 / 26,764 / 1,952,406;
P2-preclosed v3 nodes:                                      669,288;
exact-m3 nodes / closed without E3:              1,283,118 / 1,263,298;
residual E3 nodes / exact bottom sweeps:              19,820 / 191,000;
fallback nodes:                                                    0.
```

All `191,000` bottom measures are positive.  The smallest is
`57191/2522520`, at
`{1,2,3,4,5,6,8,9,10,11,12,13,14}`.  Hence `(5,6)` is uniformly closed as
well.  Three of the `21` literal flood bodies are now exact; eighteen remain,
and THM-741 remains `CLAIMED`.

## Re-root descent and the five-small shadow (codex-S16)

There is nevertheless a rigorous transport, but it acts on **completed
families**, not on Fano-equivalent root bodies.  Put

```text
H={8,9,10,11,12,13,14},        E_e=H union e
```

for an edge `e` of `K_7`.  If the uniform four-extension theorem has been
proved for `E_e`, then every 13-speed family `S` containing `E_e` is closed:
literally `S=E_e union X` with `|X|=4`, so it is one of the already certified
extensions.  In particular, the exact `(5,6)`, `(5,7)`, and `(6,7)` rows cast upward
**anchor shadows** across every other choice of root edge.

This also exposes a descent obstruction for the Fano carrier.  Given any
small-label triple `{a,b,c}`, the single actual family

```text
H union {a,b,c,15,16,17}
```

is a four-extension of each of `E_ab,E_ac,E_bc`.  Hence a scalar attached only
to the root edge and claimed to descend to completed-family data must agree on
the three edges of every triangle of `K_7`.  The seven Fano triangles alone
give seven disjoint 3-edge components: their equality system has rank `14`
and quotient dimension `7`.  The other `28` triangles are connected on all
`21` edges (degree `8`, `84` transition adjacencies), so their equality system
has rank `20` and quotient dimension `1`.  Therefore

> every root-edge-only invariant of completed flood families is constant.

The same `28` noncollinear triples index THM-850's Heawood hexagons, but the
maps are different: alternating cycle currents recover an edge **carrier**,
whereas re-root equalities coequalize presentation charts.  This proves why
the Fano/Heawood coordinates can organize inputs yet cannot be a nonconstant
symmetry quotient of the final LRC predicate.

The anchor shadows have a concrete metric payoff.  Consider a completed
13-speed family containing `H` and exactly five labels `K subset {1,...,7}`;
its last speed is then one integer `w>=15`.  Of the `C(7,5)=21` sets `K`,
exactly `18` contain one of `(5,6),(5,7),(6,7)` and are already in an anchor
shadow.  The three residual sets are

```text
12345, 12346, 12347.                                      (B1)
```

For `P=H union K`, let `G(P)` have `r` components and measure `m`.  THM-732
and `S2=99/70>sqrt(2)` give, uniformly for the remaining speed,

```text
|G(P) minus D_w| >= 6m/7-S2*r/(7w)>0
whenever w>S2*r/(6m).                                     (B2)
```

The exact verifier sweeps every integer from `15` through the floor of this
cap.  Across (B1) the three caps are respectively

```text
122, 97, 83,
```

for `260` exact interval sweeps in total.  Every measure is positive; the
smallest swept value is `16607/840840`, at `K=12345,w=23`.  At the first
integer beyond each cap, the rational lower bounds in (B2) are also checked
strictly positive, so the conclusion is genuinely uniform in unbounded `w`,
not a bounded experiment.

Every six-small-label set contains at least two of `{5,6,7}`, hence contains
one of the three anchors.  Notably, the former exceptional family
`H union {1,2,3,4,5,6}` is exactly where the new `(5,6)` full-body run attains
its minimum `57191/2522520`.  Consequently:

> **Five-small shadow.** Every 13-speed family containing `H` and at least
> five labels from `{1,...,7}` is strictly lonely.

The unresolved flood tail is therefore confined to completed families with
at most four small labels, equivalently at least two speeds above `14`.  This
does not close any one of the eighteen remaining root bodies uniformly, and THM-741
remains `CLAIMED`.

Tournament Analysis uses the three residual bases (B1), not runners and not
root edges, as vertices.  The pair observable is the exact cap difference
`C(K')-C(K)`, `C=S2*r/(6m)`; orient from lower to higher cap, using
lexicographic order for ties.  All caps are distinct, so the tournament is
transitive: score histogram `{0,...,2}:1`, no directed 3-cycles, three
singleton SCCs, and one Hamiltonian path, namely

```text
12347 -> 12346 -> 12345.
```

This quotient retains the one-speed proof horizon but destroys the interval
comb and exact clearance.  Challenging the old vertex choice is essential:
root edges are presentation charts glued by re-root triangles, not intrinsic
vertices of a completed-family transport.

## Next-anchor set system and the four-small exact closure (codex-S16)

The three certified anchors form the triangle

```text
A={56,57,67}
```

on `T={5,6,7}`; put `O={1,2,3,4}`.  A small-label set `K` is outside the
containment shadow exactly when `|K intersect T|<=1`.  Therefore the exact
shadow census by `k=|K|` is

| `k` | all `K` | shadowed | unshadowed |
|---:|---:|---:|---:|
| 2 | 21 | 3 | 18 |
| 3 | 35 | 13 | 22 |
| 4 | 35 | 22 | 13 |
| 5 | 21 | 18 | 3 |
| 6 | 7 | 7 | 0 |

Indeed the unshadowed count is
`C(4,k)+3 C(4,k-1)`.  This makes the marginal value of a new edge completely
explicit.  For an `OO` edge inside `{1,2,3,4}`, its gains in the five rows
`k=2,...,6` are

```text
C(2,k-2)+3 C(2,k-3) = (1,5,7,3,0),                       (B3)
```

whereas for an `OT` edge the gains are

```text
C(3,k-2) = (1,3,3,1,0).                                  (B4)
```

Thus every `OO` edge has the maximum possible four-small gain `7`, versus
`3` for every `OT` edge.

The low-CPU verifier also computes exact root, `E1`, and `E2` data without
entering any `v3` tree.  Ordering first by four-small gain and then by the
exact `E2` node count gives

```text
OO: 34:62282, 24:78826, 23:105089, 14:120308,
    13:128042, 12:150592;

OT: 47:21428, 37:25232, 45:32118, 35:37872,
    27:41163, 46:43067, 36:45861, 17:48818,
    25:60258, 15:71554, 26:78388, 16:100967.               (B5)
```

For the bicriterion `(four-small gain high,E2 low)`, the single-edge Pareto
frontier is exactly `{34,47}`.  If the next computation must close a whole
flood edge, `(3,4)` is therefore the coverage-first choice: it has gain `7`
and the smallest `E2` bank, `62,282`, among all gain-seven edges.  This is a
proof-work recommendation, not a proof that its deeper tree is smallest.

There is a stronger no-go against running that whole body merely to clear the
four-small stratum.  For each of the `13` presently unshadowed four-sets, put
`P=H union K`; two external speeds `15<=a<b` remain.  If `G(P)` has `r`
components and measure `m`, choose the least `V` with

```text
2/V < 5m/(S2*r).                                          (B6)
```

For `a,b>=V`, two applications of THM-732 leave at least

```text
5m/7-(S2*r/7)(1/a+1/b)>0.
```

For each `15<=a<V`, compute `G(P+a)` exactly, with component count `r_a` and
measure `m_a`.  Then THM-732 closes every

```text
b > S2*r_a/(6m_a),                                        (B7)
```

leaving only the integer pairs below that cap for exact evaluation.  Every
one of the `1,788` intermediate `a` nodes has `m_a>0`.  The finite obligations
are

```text
K:     1234 1235 1236 1237 1245 1246 1247 1345 1346 1347 2345 2346 2347
pairs: 5793 4137 3295 2263 1832 2752  971 2476 2209 1223  884  938  410,
```

for exactly `29,183` candidate pairs.  Exact rational interval subtraction
finds all `29,183` strictly positive, with zero covering failures.  The global
minimum is

```text
482219/29008980
```

at `K=1235`, `a=15`, `b=46`, i.e. for the family
`{1,2,3,5,8,9,10,11,12,13,14,15,46}`.  As an independent kernel check, five
flat-index quantiles in each of the thirteen `K`-banks (`65` samples total)
were replayed using full interval subtraction instead of the sparse measure
path.  All `65` agree exactly; the canonical sample manifest has SHA-256

```text
ce5122bd547452c916da92da507ea449753d870883c004508bbc67394ef78a89.
```

Consequently:

> **Four-small shadow.** Every 13-speed family containing
> `H={8,9,10,11,12,13,14}` and at least four labels from `{1,...,7}` is
> strictly lonely.

Together with the containment triangle, this closes all `35` four-small
completed families uniformly over their two unbounded external speeds.  At
this checkpoint the unresolved flood tail was confined to at most three small labels,
equivalently at least three speeds above `14`.  This closes no additional
whole flood body: exactly three of the `21` literal edges are still certified
uniformly, and eighteen remain.

If two more `OO` anchors were desired for the smaller-label frontier, the
three perfect-matching pairs have combined exact `E2` counts

```text
12+34:212874,       13+24:206868,       14+23:225397;
```

the matching `(13,24)` is smallest by this first-layer proxy.

Tournament Analysis uses the `18` unresolved root edges only as proof-job
vertices.  Since the four-small stratum is now closed directly, the pair
observable is the lexicographic vector
`(Delta gain_3,Delta gain_4,-Delta E2)`; orient toward larger gain and then
smaller `E2`, with lexicographic edge order for a true tie.  It is transitive, with score
histogram `{0,...,17}:1`, no directed 3-cycles, eighteen singleton SCCs,
`81` flips against the raw lex gauge, and one Hamiltonian path, beginning

```text
34 -> 24 -> 23 -> 14 -> 13 -> 12 -> 47 -> ... -> 16.
```

This scheduler retains containment gain and first-layer work but destroys
`v3` geometry, bottom margins, and final-family identity.  The then-next direct
frontier has the `22` unshadowed three-small `K`-bases as vertices; runners,
Fano flags, and root charts discard that final-family quotient.  If a whole
edge must be run next, `(3,4)` remains the coverage-first choice: it gains five
three-small sets and has the least `E2` count among all gain-five edges.

## Three-small exact closure (codex-S16)

The anchor triangle shadows exactly the thirteen three-subsets containing at
least two of `{5,6,7}`.  The other `22` bases have the form

```text
P=H union K,       |K|=3,       |K intersect {5,6,7}|<=1.
```

For three ordered external speeds `15<=a<b<c`, let the exact good set of `P`
have measure `m` and `r` components.  Choose the least `V_a` with

```text
3/V_a < 4m/(S2 r).                                      (B8)
```

Three applications of THM-732 show that every `a,b,c>=V_a` leaves positive
measure.  For each `15<=a<V_a`, subtract `D_a` exactly, obtaining `(m_1,r_1)`,
and choose the least `V_b` with

```text
2/V_b < 5m_1/(S2 r_1).                                  (B9)
```

Every `b,c>=V_b` is then positive.  Finally, after exact subtraction at each
`a<b<V_b`, with survivor `(m_2,r_2)`, THM-732 closes every

```text
c>S2 r_2/(6m_2).                                        (B10)
```

The exact bank therefore includes every integer
`b<c<=floor(S2 r_2/(6m_2))`; in particular an integral equal-cap endpoint is
not skipped.  These are the `4,5,6` residual-density coefficients from the
THM-735/732 ladder.  P2 is not used in this three-small computation.

All `4,408` one-external and `359,236` two-external nodes have positive exact
measure.  The terminal bank contains exactly `1,357,920` triples.  Every
terminal sparse subtraction is strictly positive, with global minimum

```text
28470499/1520554035
```

at `K=134,(a,b,c)=(17,23,37)`.  Five deterministic flat-index samples in each
base bank, `110` total, agree with full subtraction through the same pinned
core.  These are sparse/full cross-path checks, not an independent kernel.
The row, cross-path, and global certificate manifests are respectively

```text
7632b1ef6cfaadc11f8091a7ac3d6fbebe49b353944d60f2c19561256c97658f,
3a4e288500174602ea83cd48db388f7160de6256c668b282a301d94db52fe886,
b05c56708781287bdb958c7937d5e70a781803a4ccb9a3226fc6c83ed957b886.
```

The canonical 22-row JSONL is hash-linked to the exact source and premise;
normal and optimized render-only runs reproduce the stored report byte for
byte.  The stored JSONL is complete and clean.  Operationally, its malformed-
last-line reader does not truncate a partial fragment before a later append,
so such a crash residue should be removed before resuming; this caveat does
not affect the frozen 22-row artifact.  The zero ledgers are empty.  Had a zero occurred, it would have been a
zero-measure/tight-certificate alarm rather than, by itself, proof of a
covering family.

Consequently:

> **Three-small shadow.** Every 13-speed family containing
> `H={8,9,10,11,12,13,14}` and at least three labels from `{1,...,7}` is
> strictly lonely.

For each of the eighteen root-edge bodies not already certified, the only
unresolved completions now add four speeds all greater than `14`; equivalently,
the fixed root edge supplies the completed family's only two small labels.
This common shadow still closes no fourth whole body: the whole-body count
remains `3/21`, and THM-741 remains `CLAIMED`.  The live direct runner
independently reached `290/2002` clean bodies at the last pull; its unharvested
resume ledger is not a substitute for a completed theorem.

Tournament Analysis uses the 22 residual completed-family bases as proof-job
vertices.  Orient by fewer terminal triples, then fewer `b` and `a` nodes,
with lexicographic tie-breaking.  The tournament is transitive, has score
histogram `{0,...,21}:1`, no directed triangles, 22 singleton SCCs, 196 flips
from lexicographic order, and one Hamiltonian path beginning

```text
347 -> 247 -> 345 -> 346 -> ... -> 134 -> 123.
```

It retains exact finite workload and final-family base but destroys interval
geometry and margins.  Runners, Fano flags, and root-edge charts do not even
preserve the completed-family bank being decided.

## Ancestral-carrier envelope for the pure `(3,4)` tail (codex-S16)

Fix `E=H union {3,4}` and write the remaining ordered speeds as
`15<=a<b<c<d`.  Let `A` be an exact current survivor and let `C` be any
ancestor carrier with `A subset C`, measure `m_C`, and `r_C` components.  If
`s` ordered needles `w_i` remain, THM-732 gives the reusable one-sided bound

```text
|A minus union_i D_(w_i)|
 >= |A|-s m_C/7-(S2 r_C/7)sum_i 1/w_i.                  (B11)
```

If the first possible speed is `v`, the reciprocal sum is largest at
`v,v+1,...,v+s-1`.  Applying P2 to that first needle and charging the other
`s-1` needles to `C` also gives

```text
|A minus union_i D_(w_i)|
 >=6|A|/7-8r_A/(49v)-(s-1)m_C/7
   -(S2 r_C/7)sum_(i=1)^(s-1)1/(v+i).                   (B12)
```

Taking the maximum of (B11)--(B12) over all ancestors is sound.  It is a proof
certificate, not an equality quotient: it retains ancestor containment and
the `(m,r)` tradeoff while discarding later-needle overlap geometry.

The exact root is

```text
r=28,       m=433607/2522520.
```

The common-threshold tree starts at `V1=308`.  Merely using distinct ordered
speeds improves the union cutoff to `306`; the P2/ancestor envelope improves
it to `291`.  Thus 17 of the 293 baseline `a` branches close algebraically.
Exact sparse first-child measures close 63 more, leaving 213 branches that
need a literal `G1` carrier; the baseline exact `E2` bank has 61,379 nodes.

A deliberately bounded `15<=a<=40` audit performs no `c` subtraction and no
`d` sweep.  It reduces 4,796 standard `b` nodes to 4,391 sparse measures and
2,853 literal `G2` constructions.  On those survivors the standard `c`
frontier has 244,652 nodes, whereas the ancestor envelope retains 144,695,
closing 99,957 without constructing `G3`.  Normal and optimized runs
byte-match the frozen exact output.

Tournament vertices are the ancestor levels `C0,C1,C2`, not runners or Fano
flags.  The pair observable is the nodewise smaller certified `c` horizon,
then residual workload; the tournament is transitive with Hamiltonian path
`C1->C0->C2`.  Different ancestors win in real nodes, so keeping the entire
ancestor antichain is useful.  This bounded prepass proves only the screened
subtrees.  It evaluates no final margins and does not close the whole `(3,4)`
body or THM-741. ∎

## Exact first branch of the pure `(3,4)` tail (codex-2026-07-17)

The bounded prepass above did not evaluate any final margins.  The first
literal-carrier branch is now complete.  Fix

```text
E={3,4,8,9,10,11,12,13,14},       a=15,
```

and let `15<b<c<d`.  Exact subtraction at `a=15` gives

```text
r_1=26,       m_1=184909/1261260,       V2=189.
```

Thus `b>=189` is closed by THM-735's three-needle common-threshold leg.  The
remaining `173` values of `b` generate `11,177` `c` obligations.  The monotone
P2/fixed-`E2` inequality closes `2,849` of them without an exact `m3`.  Of the
other `8,328`, exact sparse `m3` closes `7,166` without constructing a `d`
bank.  The final `1,162` `c` carriers generate exactly `17,198` integer `d`
sweeps below their rational caps.  Every sweep is strictly positive, with
minimum swept clearance

```text
32953849/624660036
```

at `(b,c,d)=(17,19,23)`.  Values of `d` above each cap are covered by the
strict fixed-`E2` discrepancy inequality, including the correct integral-cap
endpoint convention.  Hence

> `E union {15,b,c,d}` is lonely for every `15<b<c<d`.

The terminal and certificate ledgers have SHA-256 hashes

```text
3594fb4a07d9ee79780f7c99cf4cf2427b0a921282f8c2c19249c46c339602b2,
2ce412ad92743b74400e3e4cba57dfb1c1c75a3dc4b919a998d762df050f180a.
```

For a second evaluation path, the first, middle, and last terminal of every
active `b` bank were rebuilt as a full thirteen-comb union.  All `177` samples
agree with nested sparse subtraction; their manifest is
`623f110f3f6a0f26b275cfe2097e80d01cb11705008228e70ec2734d9c02f4cc`.
Normal and `-O` runs are byte-identical.  Source/output hashes are
`0d63edea45a523ddd71494e02d4500fe9725cf4e424d11a3a70d9dd6e2a9f0fc`
and `d50adea8ccabae192c95f235067e240e0530d80e52164291b07f632326d948fc`.

The proof-bearing vertices are `b`-indexed obligations with their nested
survivor components, rational cap, certificate type, and exact terminal
margins.  Their workload tournament is transitive on the `61` active banks,
but discards interval geometry.  A Fano/`chi_7` flag identifies `(3,4)` only
as one of the 21 edge addresses and supplies no metric transport.  In Kakeya
language each `D_w` is a translated periodic one-dimensional needle comb;
shared adaptive component incidence, not an isolated needle or a dimension
estimate, is the faithful statistic.

Combined with the earlier root/first-measure screen, the `(3,4)` body now has
17 root-closed branches, 63 exact-measure-closed branches, this first complete
literal-carrier branch, and `212` literal-`G1` branches still open.  This does
not close a fourth whole flood body or global THM-741.

## Exact second branch and the nontransported carrier (codex-2026-07-17)

The next first-speed branch also closes, but not by monotone transport from
`a=15`.  For `a=16`, exact first subtraction gives

```text
r_1=26,       m_1=29921/210210,       V2=194.
```

The two first-child survivor sets are genuinely incomparable.  Their common
measure and directional differences are

```text
|G15 intersect G16| = 105857/840840,
|G15 minus G16|     =   4019/194040,
|G16 minus G15|     =    419/25480.
```

Thus neither carrier contains the other.  What transports is only the
certificate schema: the THM-735 common-threshold leg, then P2 against each
fresh exact `E2`, then the fixed-`E2` rational `d` cap.

For `a=16`, the `177` values `17<=b<194` generate `11,786` `c` obligations.
Every P2 predicate was evaluated linearly—not imported from the `a=15`
cutoffs—and its truth set was checked to be a terminal interval.  This closes
`2,986` nodes.  Of the `8,800` exact sparse-`m3` nodes, `7,579` close without a
`d` bank.  The remaining `1,221` carriers generate exactly `18,182` finite
terminal sweeps.  There are no fallback nodes and no integral cap endpoints.
Every sweep is positive; the minimum swept clearance is

```text
6999703/133617120
```

at `(b,c,d)=(19,23,32)`.  Consequently

> `E union {16,b,c,d}` is lonely for every `16<b<c<d`.

The terminal and certificate ledgers have SHA-256 hashes

```text
20daaa1ec85df38843115c5577f4eb2a48080f51e60dcdeec13cb272e8f3b1b4,
001c43a3d7e830473451c2eee33732484e3b00dfd6b97933f7715869bcf9996a.
```

First/middle/last terminals in every active `b` bank give `187` full-union
crosschecks, all agreeing with nested sparse subtraction; their manifest is
`1652d5afe562c4222bec2e500c15d838de1c2f8102deaae8d409a5f84c5b4ab3`.
Normal and `-O` runs are byte-identical.  Source/output hashes are
`0445c144120bec1d48ea16c7001ff723ebc6bbd6431ad13b6f057b7c49e44020`
and `0d90a88451b868f65f139fc3ee38c17ad2d0dfeb93cf1d187cb4846a47d2f104`.

This is an instructive Kakeya/Fano guardrail.  `D_16` is a new periodic comb,
not a monotone extension or translate of `D_15`; in fact the later branch has
more `b`, `c`, and terminal obligations.  Fano/`chi_7` still identifies only
the common root edge `(3,4)`.  The predicate-preserving carrier must retain
the freshly recomputed nested components, phase, certificate type, rational
caps, and margins.  The scheduler tournament on the `63` active `b` banks is
transitive telemetry and loses that geometry.

After the root and exact-measure screens, two literal-carrier branches are now
complete and `211` remain.  No fourth whole flood body and no global THM-741
claim follows.

## Exact third branch and the nonmonotone workload (codex-2026-07-17)

The same schema closes `a=17`, again from a freshly constructed carrier.  Exact
first subtraction gives

```text
r_1=26,       m_1=758281/6126120,       V2=223.
```

An independent bidirectional comparison with `a=16` gives

```text
|G16 intersect G17| =   992991/9529520,
|G16 minus G17|     = 1090283/28588560,
|G17 minus G16|     =     25831/1319472.
```

Both directional differences are positive.  Thus the proof schema transports,
but neither the carrier nor an assumed monotone horizon does: `G17`, every P2
cutoff, and every fixed-`E2` cap were recomputed literally.

The `205` values `18<=b<223` generate `15,795` `c` obligations.  Linear audits
show each P2 truth set is a terminal interval and close `4,001` nodes.  Of the
remaining `11,794` exact sparse-`m3` nodes, `10,065` close before constructing
a `d` bank.  The final `1,729` carriers generate exactly `30,507` finite
terminal sweeps.  There are again no fallback nodes and no integral cap
endpoints.  Every sweep is positive; the minimum is

```text
2503360059/62849566780
```

at `(b,c,d)=(23,31,37)`.  Therefore

> `E union {17,b,c,d}` is lonely for every `17<b<c<d`.

The terminal and certificate ledgers have SHA-256 hashes

```text
54fd735e7e441d50b3595fa69332a4b752902511a0077e61eac205b065101e7f,
fb8edfdddc36a69969b6799b116a8eeffbcdfde2960aabbb7163791bb4bbae39.
```

First/middle/last terminals from every active `b` bank give `209` independently
rebuilt full-union crosschecks; their manifest is
`b58c644435d860c1fa9a0911524d00b3cb0c031a4f5b99f32628dfebefee42a2`.
Normal and `-O` runs are byte-identical.  Source/output hashes are
`120074e9617d50f353427d94f3babfb42a77cebff8dab3af74d1c445522fbc53`
and `1bc560811bf531f41ac3f0f4e0763d9a98ffa5f5dc605903c6e68130c8ce2f79`.

The tournament on the `70` active `b` obligations orients by minimum swept
clearance, then workload and label.  It is transitive, with one Hamiltonian
path, and is scheduler telemetry rather than a symmetry quotient.  The exact
carrier retains the nested component stalk, phase, certificate type, caps, and
margins; runner labels, residues, isolated Kakeya needles, wall events, and
Fano flags do not preserve it.  In particular `D_17` changes the stalk in both
directions relative to `D_16`, while `chi_7` supplies only their common root-edge
address `(3,4)`.

The workload has increased sharply from `a=16` despite the adjacent first
speed, so no monotone-workload inference is available.  Three literal-carrier
branches are now complete and `210` remain.  This is still one pure tail: it
does not close a fourth whole flood body or global THM-741.

## Exact fourth branch and adjacent-branch oscillation (codex-2026-07-17)

The next branch also closes from fresh exact data.  For `a=18`, first
subtraction gives

```text
r_1=26,       m_1=17747/120120,       V2=187.
```

The carrier again fails to nest with its predecessor:

```text
|G17 intersect G18| = 1564249/14294280,
|G17 minus G18|     =     30761/2144142,
|G18 minus G17|     =      45637/1191190.
```

Both directional differences are positive.  Every `G18` interval, P2 cutoff,
and fixed-`E2` cap was therefore recomputed; no `a=17` result row was reused.

The `168` values `19<=b<187` generate `10,762` `c` obligations.  Linear P2
audits close `2,750`, leaving `8,012` exact sparse-`m3` nodes.  Of these,
`6,987` close before a `d` bank.  The final `1,025` carriers generate `14,399`
finite sweeps, with no fallback node and no integral cap endpoint.  All are
strictly positive.  The minimum is

```text
518788/9144135
```

at `(b,c,d)=(21,29,32)`.  Hence

> `E union {18,b,c,d}` is lonely for every `18<b<c<d`.

The terminal and certificate ledger hashes are

```text
44bd2e2dea062e5090243177f9c8871b0ce820767f972df7278cb5e56d670b08,
f25008da5e5009ae1d3ee15063c90086b599af689ba842e4119e9e03ef594f06.
```

There are `165` independently rebuilt first/middle/last full-union samples,
with manifest
`86046df7baa90ce91798f30df2a0bd8c6babac70d6e830651c65a97703a895a3`.
Normal and `-O` runs are byte-identical.  Source/output hashes are
`e9c0f84f6eb0d0ac81fb4f5241db9163843c3345228ad518cc1702ada460dc5b`
and `827e17e0fa3a3828dbe45de01931fe9285a4a0ce13a29562fd6099d8661b07c1`.

The scheduler tournament on the `56` active `b` obligations is transitive
under minimum clearance, workload, and label, with one Hamiltonian path.  It
retains proof-job telemetry but not the component geometry.  The actual
predicate-bearing object remains the nested component stalk with phase,
certificate type, caps, and margins.  Isolated Kakeya needles, runner labels,
residues, wall events, and Fano flags destroy that object; `chi_7` again gives
only the common root-edge address `(3,4)`.

The sequence itself is structural evidence: from `a=16` through `a=18`, the
fresh cutoffs are `194,223,187` and terminal workloads are
`18,182,30,507,14,399`.  Both rise and then fall while consecutive carriers
remain incomparable.  Thus neither a monotone horizon nor monotone workload
can substitute for exact branch geometry.  Four literal-carrier branches are
now complete and `209` remain.  No fourth whole flood body or global THM-741
claim follows.

## Exact fifth branch and cutoff-level nonmonotonicity (codex-2026-07-18)

The next fresh carrier tree closes `a=19`.  Exact first subtraction gives

```text
r_1=26,       m_1=6388289/47927880,       V2=207.
```

It cannot be obtained by nesting the `a=18` carrier.  The literal interval-set
comparison is

```text
|G18 intersect G19| = 5691269/47927880,
|G18 minus G19|     =       15793/544635,
|G19 minus G18|     =       11617/798798.
```

Both directional differences are positive.  More strongly, both adjacent
second-level horizons were rebuilt on all `167` common `b` labels.  Relative
to `a=18`, the fresh `V3` horizon moves down/equal/up on `25/3/139` labels and
the P2 cutoff moves down/equal/up on `23/6/138`.  The full paired-bank digest is

```text
d7bc67012c8d5e161dfcefff0840ecef186e4341b86b295653fb6544921addae.
```

Thus even cutoff direction is nonmonotone; no adjacent cutoff or result row is
transported.

The `187` values `20<=b<207` generate `13,332` `c` obligations.  Linear P2
audits close `3,361`, leaving `9,971` exact sparse-`m3` nodes.  Of those,
`8,597` close before a `d` bank.  The final `1,374` carriers generate exactly
`21,302` finite terminal sweeps, with no fallback and no integral cap endpoint.
Every sweep is positive.  The minimum is

```text
368611/8351070
```

at `(b,c,d)=(23,32,52)`.  Consequently

> `E union {19,b,c,d}` is lonely for every `19<b<c<d`.

The terminal and certificate ledgers have hashes

```text
91a878c51d96cebeaa015598b4c78f14c8a6b0fa181dc87827108a22c7e7aaeb,
86533640f52a442c138fcea45068fe06fd91fb9e8de353ee187cfcca6166a567.
```

First/middle/last terminals from every active `b` bank give `187` independently
rebuilt full-union crosschecks, all exact, with manifest
`1c9d0dec8d4f6a8065021930182cce6763506885d09465bbd57252c5fd3dedd3`.
Normal and `-O` runs are byte-identical.  Source/output hashes are
`a63afe6693107e1fc4fdc10c84f7ecdb62673c226bc5817e1ec9e804c5b0b416`
and `9b841ab368e676278450b69581c699f4487e62cc57f0c7af82195b63f0523e7d`.

Tournament Analysis uses the `63` active `b` obligations, with observable
`(minimum swept clearance,-sweeps,-b)`.  Its lexicographic switch gives a
transitive tournament, one Hamiltonian path, and path digest
`fe91efac9d31909e0b3317cadb8ef41130ae40d17ff8c8ebb09c44235ddc9e52`.
This is scheduling telemetry only.  The proof carrier is the nested component
stalk with phase, certificate type, cap, and exact margin.  Runners, residues,
isolated Kakeya needles, wall events, and Fano flags erase necessary data;
`chi_7` supplies only the common root-edge address `(3,4)`.

At that branch-by-branch checkpoint, five literal-carrier branches were
complete and `208` remained.  That certificate by itself did not close a
fourth whole flood body or global THM-741.

## Exact top-four root-coverage envelope closes the whole `(3,4)` tail (codex-2026-07-28)

The branch tree above retains the full nested component stalk, which is needed
when one estimates each later comb against the carrier left by the earlier
ones.  For this root, however, a coarser object is unexpectedly decisive:
the **individual coverage profile of every future comb against the unchanged
root carrier**.

Let

```text
E={3,4,8,9,10,11,12,13,14},       G=G(E),
r=28,                              m=433607/2522520,
c(w)=|G intersect D_w|.
```

For four distinct future speeds `15<=a<b<c<d`, the ordinary union bound gives

```text
|G minus (D_a union D_b union D_c union D_d)|
 >= m-c(a)-c(b)-c(c)-c(d).                              (B13)
```

The exact verifier computes `c(w)` for every integer `15<=w<=472`.  Sparse
subtraction, full subtraction, a direct ten-comb union, and an independent
endpoint-incidence sum over the disjoint teeth of `D_w` agree on all `458`
values.  The four largest are

```text
rank:       1                    2                     3                4
w:         17                   23                    19               32
c(w):      85973/1786785        240307/5801796        14017/363090     11521/315315.
```

No uncomputed speed can enter this list.  Indeed THM-735(ii), through the
THM-731/732 covariance-discrepancy chain, gives

```text
c(w) <= m/7+sqrt(2)r/(7w)
     <  m/7+(99/70)r/(7w) =: u(w).                       (B14)
```

The rational crossing against the fourth value is exact:

```text
(99/70)r / (7(c(32)-m/7)) = 33297264/70523,
472 < 33297264/70523 < 473,
u(472)-c(32) = 1301/347266920 > 0,
c(32)-u(473) = 1093/50618568 > 0.                        (B15)
```

Since `u(w)` decreases, every `w>=473` has `c(w)<c(32)`.  The displayed four
values are therefore the global top four over every integer `w>=15`.  The
distinctness of the four added speeds is the only ordering information used.
Their total possible root coverage and the residual in (B13) are

```text
c(17)+c(23)+c(19)+c(32) = 308603791/1873980108,
m-that sum                  = 135228493/18739801080 > 0. (B16)
```

Consequently

> **Fourth whole flood body.** For every four integers
> `15<=a<b<c<d`, the family
> `E union {a,b,c,d}` is strictly lonely.

The literal family formed from the four maximal individual coverages,

```text
{3,4,8,9,10,11,12,13,14,17,19,23,32},
```

has exact lonely measure `805109/16702140` and `16` components by both direct
thirteen-comb union and nested subtraction.  This is a hostile cross-check,
not the reason (B16) is uniform.  Conversely, the older ancestral root bound
(B11) at the four first possible speeds `15,16,17,18` is
`-81555377/62537475<0`; the new gain comes from ranking exact individual comb
coverages rather than charging all combs by the same worst-case discrepancy.

At that checkpoint this closed one literal root edge, raising the exact flood
count from `3/21` to `4/21`.  It uses no Fano/`chi_7` transport.  The earlier
five full first-speed computations remain independent nested-carrier controls
but are strictly subsumed for this body.  The other seventeen pure tails
remained; the next section applies the same envelope across six more roots.
Global THM-741 and LRC(14) remain open. ∎

## The top-four envelope closes six more literal roots (codex-2026-07-28)

The preceding proof is an instance of a reusable finite-envelope lemma.  For
any fixed root good set `G` with measure `m` and `r` components, put
`c(w)=|G intersect D_w|` and

```text
u(w)=m/7+(99/70)r/(7w).
```

Suppose an exact scan of `15<=w<=W` has four largest values
`q1>=q2>=q3>=q4`, and `u(W+1)<q4`.  By THM-735(ii),
`c(w)<u(w)<=u(W+1)<q4` for every `w>W`.  Thus `q1,...,q4` are the global top
four coverages over every future integer speed.  Consequently

```text
m-(q1+q2+q3+q4)>0                                      (B17)
```

is a uniform proof for all four distinct future speeds.  This lemma retains
the complete one-comb coverage profile against the root but discards all
later nested-carrier geometry.  It is therefore a sufficient envelope, not
an equivalence.

The generic verifier uses the common exact bank `15<=w<=482` for six further
literal roots.  On every one of the `6*468` root/speed pairs, sparse
subtraction, full subtraction, direct ten-comb union, and an independent
two-pointer tooth-incidence sum agree exactly.  Write

```text
T=(99/70)r/(7(q4-m/7));
```

then `u(floor(T))>q4>u(floor(T)+1)` is checked in exact rational arithmetic.
The complete certificates are:

| root edge | `r` | `m` | top-four speeds | top-four sum | margin in (B17) | exact `T` |
|---|---:|---:|---|---:|---:|---:|
| `24` | 30 | `31789/194040` | `23,22,17,16` | `312668333/1972610640` | `3499547/657536880` | `214053840/484061` |
| `35` | 24 | `652/3465` | `17,19,23,46` | `3436749929/18739801080` | `1988211/416440024` | `10820304/31291` |
| `36` | 28 | `504167/2522520` | `19,17,23,21` | `1786984753/9369900540` | `57162379/6246600360` | `49945896/153311` |
| `45` | 24 | `514471/2522520` | `17,23,19,32` | `1274636801/6814473120` | `6497515/384406176` | `342486144/981901` |
| `46` | 28 | `518957/2522520` | `17,19,23,32` | `349312099/1873980108` | `2532941/131047560` | `33297264/69023` |
| `47` | 20 | `261193/1261260` | `17,19,23,32` | `2395349819/12493200720` | `575561731/37479602160` | `28540512/96835` |

The largest crossing is the `46` row:

```text
482 < 33297264/69023 < 483,
u(482)-q4 = 14089/1418497080 > 0,
q4-u(483) = 389/27075048 > 0.
```

Hence the common endpoint `482` is not a rounded numerical convenience; it
is the exact minimal suite-wide finite cutoff.  Direct thirteen-comb union
and nested subtraction also agree on the family formed from each row's four
maximal individual coverages.  Their exact lonely measures are respectively

```text
98391511/1972610640, 87938369/1249320072,
7150467/109589480,   2647649/40085136,
177789/2783690,      1969187/30834720.
```

Roots `13` and `16` are complete negative controls through the same
four-path bank and exact tail crossing.  Their global top-four margins are

```text
-25012943/1249320072,       -6867281/312330018,
```

while their literal top-four families still have positive exact lonely
measure `906053/21917896` and `12226321/312330018`.  Thus the verifier does
not confuse positivity of one extremal family with positivity of the
root-coverage union envelope.  Root `37`, the closer exceptional profile, is
reserved for a separate higher-rank audit and is not claimed here.

Therefore six additional whole flood bodies are proved:

> for each `e` in `{24,35,36,45,46,47}` and every
> `15<=a<b<c<d`, the family `H union e union {a,b,c,d}` is
> strictly lonely.

Together with `34,56,57,67`, the exact whole-body census at that checkpoint was

```text
{24,34,35,36,45,46,47,56,57,67}: 10/21.
```

At this checkpoint the eleven root edges
`12,13,14,15,16,17,23,25,26,27,37` remain outside this certificate; the
next section closes `37` by retaining one more coverage rank.  No
Fano/`chi_7` transport is used.  Global THM-741 and LRC(14) remain open. ∎

## Top-five envelope plus one exact exception closes `(3,7)` (codex-2026-07-28)

The top-four method has a sharp repair when its margin fails but the fifth
coverage is separated.  Put

```text
E={3,7,8,9,10,11,12,13,14},       G=G(E),
r=20,                              m=53619/280280,
c(w)=|G intersect D_w|.
```

Four independent exact paths compute `c(w)` for `15<=w<=293`.  The first
five ranks are

```text
rank:       1                    2                  3                  4                  5
w:         19                   17                 21                 23                 46
c(w):      134663/2662660       2578/51051         38609/840840       25331/552552       79435/1933932.
```

The THM-735(ii) cap from (B14) crosses the fifth value at

```text
(99/70)r / (7(c(46)-m/7)) = 547026480/1860739,
293 < 547026480/1860739 < 294,
u(293)-c(46) = 1829953/39664945320 > 0,
c(46)-u(294) = 733/947626680 > 0.                        (B18)
```

Thus these are the global first five ranks for all integer `w>=15`.  The
plain top-four margin is genuinely negative:

```text
m-c(19)-c(17)-c(21)-c(23)
 = -3183347/2082200120 < 0.                              (B19)
```

But distinctness now makes the residual obligation a singleton.  Any
four-speed set other than the literal top-four set omits at least one of
those four ranks, so its individual-coverage sum is at most
`c_1+c_2+c_3+c_5`.  Exact simplification gives

```text
m-c_1-c_2-c_3-c_5
 = 843411/260275015 > 0.                                 (B20)
```

The union bound therefore closes every nonexceptional quadruple.  For the
sole exception, nested subtraction and a direct thirteen-comb union agree:

```text
E union {17,19,21,23}
has 8 good components and measure 137897/2299080 > 0.    (B21)
```

The large value in (B21) exhibits the missing coordinate in the failed
top-four envelope: its four high individual coverages overlap heavily.
There is no need to bound those overlaps uniformly because the fifth-rank
gap isolates exactly one quadruple.

Consequently

> **Eleventh whole flood body.** For every four integers
> `15<=a<b<c<d`, the family
> `{3,7,8,9,10,11,12,13,14,a,b,c,d}` is strictly lonely.

The companion verifies all `279` finite coverages by sparse subtraction,
full subtraction, direct ten-comb union, and independent tooth incidence;
it also verifies both constructions of (B21), the exact tail crossing, and
the failed plain-envelope control (B19).  Combined with the preceding
six-root theorem, this raises the exact flood count from `10/21` to `11/21`.
The ten root edges `12,13,14,15,16,17,23,25,26,27` remain.  This proves no
Fano/`chi_7` transport and not global THM-741 or LRC(14). ∎

## Ranked finite heads close `(2,3),(2,5),(2,6),(2,7)` (codex-2026-07-28)

The isolated top-five repair extends to a finite ranked-exception lemma.  Let
the globally ranked individual root coverages be

```text
q1>=q2>=... .
```

If, for some `K`,

```text
m-q1-q2-q3-q_(K+1)>0,                                  (B22)
```

then every four-speed set not wholly contained in the first `K` ranks is
positive by the union bound: at most three of its speeds can occupy the first
three ranks, and its fourth coverage is at most `q_(K+1)`.  Only the finite
head of `C(K,4)` literal quadruples remains.  This is a recursive sharpening
of the top-four profile: preserve progressively more ranked atoms only until
the tail becomes uniformly harmless, then restore full interval geometry on
that finite head.

For the four remaining non-`1` root edges, exact coverage banks give:

| root | `r,m` | `K` | finite `W` | `q_(K+1)` | exact cap crossing `T` | gate (B22) | preceding gate | `C(K,4)` |
|---|---|---:|---:|---|---|---:|---:|---:|
| `23` | `30,358427/2522520` | 36 | 1796 | `c(350)=29/1225` | `535134600/297953` | `30353/410960550` | `-667243/3328780455` | 58,905 |
| `25` | `26,409261/2522520` | 22 | 1175 | `c(117)=317/11466` | `92756664/78919` | `37841/551170620` | `-2069/91861770` | 7,315 |
| `26` | `30,413747/2522520` | 10 | 840 | `c(25)=3303/107800` | `267567300/318211` | `69547/612411800` | `-81687/306205900` | 210 |
| `27` | `22,52147/315315` | 23 | 924 | `c(92)=206/7245` | `225648423/244061` | `15721/367447080` | `-1269581/58424085720` | 8,855 |

Here

```text
T=(99/70)r/(7(q_(K+1)-m/7)),
floor(T)=W,
u(W)>q_(K+1)>u(W+1).                                   (B23)
```

Thus the listed head and outside rank are global over every integer `w>=15`,
not artifacts of a truncated sort.  The negative preceding gates prove that
each `K` is the first possible rank for this particular three-largest-plus-
outside certificate.  The plain top-four margins are also all negative:

```text
root 23: -2106724/156165009,     root 25: -5458325/468495027,
root 26: -1226821/122482360,     root 27: -1002853/76488984.
```

So none of these four bodies was already closed by the preceding top-four
lemma.

The verifier evaluates every one of the `75,285` head quadruples.  It caches
full exact carriers through the first three speeds, subtracts the fourth
sparsely, replays five deterministic terminal quantiles per root by full
subtraction, and rebuilds each global minimum both as a direct thirteen-comb
union and as four full nested subtractions.  Exactly `79,31,21,36` head
quadruples, respectively, have nonpositive individual-coverage union margins;
all of those and all remaining head quadruples are strictly positive.  The
global head minima are

| root | minimum lonely measure | four speeds | components |
|---|---:|---|---:|
| `23` | `1068173/39939900` | `16,19,22,25` | 20 |
| `25` | `11739671/340723656` | `17,19,23,25` | 14 |
| `26` | `287374099/9369900540` | `17,19,22,46` | 20 |
| `27` | `39526373/1338557220` | `17,19,26,46` | 18 |

For roots `23` and `27` the actual head minimum even lies outside the
nonpositive-union-margin subbank.  This is a useful guardrail: the ranking
selects the finite proof universe, but it does not predict the metric minimum
after comb overlaps are restored.

Consequently

> for each `e` in `{23,25,26,27}` and every
> `15<=a<b<c<d`, the family `H union e union {a,b,c,d}` is
> strictly lonely.

The ten earlier bodies, the independent root-`37` top-five certificate, and
these four ranked-head bodies give

```text
{23,24,25,26,27,34,35,36,37,45,46,47,56,57,67}: 15/21.
```

The exact residual is the six-edge star at label `1`:

```text
12,13,14,15,16,17.
```

This is a structural change in the frontier, not an asserted symmetry:
literal root carriers were computed throughout, and no Fano/`chi_7`
transport was used.  Global THM-741 and LRC(14) remain open. ∎

## A finite-head/tail partition closes the label-1 star (codex-2026-07-28)

It remains to close the pure tails of the six roots

```text
E_b={1,b,8,9,10,11,12,13,14},       2<=b<=7.
```

Write `G_b=G(E_b)`, let `m_b=|G_b|` and `r_b` be its number of interval
components, and put

```text
c_b(w)=|G_b intersect D_w|,
u_b(w)=m_b/7+(99/70)r_b/(7w).                           (B24)
```

THM-735(ii), through the THM-731/732 covariance-discrepancy chain, gives
`c_b(w)<u_b(w)`.  Fix the exact horizon

```text
W=2500,                    u_*=u_b(2501).
```

For every `15<=w<=2500`, the companion computes `c_b(w)` in four exact ways:
sparse subtraction from `G_b`, full subtraction, direct ten-comb union, and
an independent two-pointer sum over the teeth of `D_w`.  All `6*2486`
comparisons agree.  If `q_1>=q_2` are the top two values in this finite bank,
the exact data are:

| root | `r_b,m_b` | `q_1,q_2` (speed:value) | `u_*` | `m_b-q_1-q_2-2u_*` |
|---|---|---|---:|---:|
| `12` | `32,319927/2522520` | `23:453587/11603592`, `46:212279/5801796` | `182859895/8832351528` | `4947963661/507860212860` |
| `13` | `30,3319/25740` | `23:453587/11603592`, `17:67493/1786785` | `65750513/3154411260` | `25633591579/2466749605320` |
| `14` | `30,335047/2522520` | `23:240307/5801796`, `17:135481/3573570` | `944979467/44161757640` | `61555969441/5755749079080` |
| `15` | `26,6716/45045` | `23:41075/892584`, `19:36649/798798` | `527231/22531509` | `9555282629/918985147080` |
| `16` | `30,365567/2522520` | `19:377149/7987980`, `23:47539/1054872` | `1021309987/44161757640` | `15408825262/2412336011085` |
| `17` | `22,384011/2522520` | `19:128503/2662660`, `23:25331/552552` | `1038897919/44161757640` | `106901554667/9649344044340` |

Every displayed margin is positive and every row also has `q_2>u_*`.
Consequently, if at least two of the four added speeds exceed `W`, their
total root coverage is at most `q_1+q_2+2u_*`.  For three tails this follows
from `q_1+3u_*<q_1+q_2+2u_*`; for four tails it follows from
`4u_*<q_1+q_2+2u_*`.  Thus every branch with at least two tail speeds is
closed.

Suppose exactly one speed `t` exceeds `W`, and let `S` be the other three.
The root union bound closes the branch unless

```text
sum_{w in S} c_b(w)+u_* >= m_b.                          (B25)
```

Equality is deliberately included in the dangerous bank.  There are only
`1,6,1,4,6,3` such triples for roots `12,...,17`.  For each one, construct
the literal residual carrier

```text
G_{b,S}=G_b minus union_{w in S}D_w
```

with measure `m_{b,S}` and `r_{b,S}` components.  Applying (B24) to this
residual carrier gives, for every `t>W`,

```text
|G_{b,S} minus D_t|
 > 6m_{b,S}/7-(99/70)r_{b,S}/(7*2501).                  (B26)
```

The minimum of the right side over each complete dangerous bank is:

| root | dangerous triples | minimum in (B26) | attaining `S` | `m_{b,S},r_{b,S}` |
|---|---:|---:|---|---|
| `12` | 1 | `71621894273/1608224007390` | `19,23,46` | `29780131/551170620,22` |
| `13` | 6 | `34858373223/959291513180` | `17,19,23` | `958105/21917896,14` |
| `14` | 1 | `1865084760919/54679616251260` | `17,19,23` | `154798927/3747960216,16` |
| `15` | 4 | `15348916166/365505456225` | `19,23,25` | `6349621/125266050,18` |
| `16` | 6 | `156275963542/4556634687605` | `17,19,23` | `6484031/156165009,16` |
| `17` | 3 | `12865251479/300437451930` | `17,19,23` | `1052099/20593188,12` |

All are strictly positive, so every exactly-one-tail branch closes.

It remains only to take all four speeds in the finite head.  The root union
bound closes every quadruple except those satisfying

```text
sum_{w in S}c_b(w) >= m_b.                               (B27)
```

Again equality is included.  The complete candidate counts and exact minimum
survivors are:

| root | candidates in (B27) | minimum survivor | attaining speeds | components |
|---|---:|---:|---|---:|
| `12` | 2,589 | `17069726431/693372639960` | `17,19,23,37` | 16 |
| `13` | 7,609 | `48080163/2027405380` | `17,19,23,37` | 16 |
| `14` | 2,581 | `162578519/7415750160` | `19,23,32,37` | 16 |
| `15` | 7,290 | `544970477/18539375400` | `19,23,25,37` | 16 |
| `16` | 7,522 | `1707671939/77041404440` | `17,19,23,37` | 16 |
| `17` | 5,087 | `13420153/455992020` | `17,19,23,31` | 16 |

Every candidate is rebuilt both by nested exact subtraction and by direct
thirteen-comb union; the full carriers, component counts, and measures agree.
The candidate banks themselves are independently generated by specialized
descending-rank loops and by an arity-generic recursion whose exact suffix
sum is an upper envelope for every pruned branch.  Both enumerators return
the same labelled tuples.  The respective dangerous/head recursion-node
counts are

```text
12: 4/2614,  13: 12/7642, 14: 4/2604,
15: 7/7320,  16: 10/7556, 17: 6/5112.                   (B28)
```

As hostile controls, the ordinary top-four margins for all six roots are
negative:

```text
-3653017/275585310, -25012943/1249320072,
-50908597/3747960216, -68144053/3747960216,
-6867281/312330018,  -228625/11639628.                  (B29)
```

Thus this is a genuine residual-carrier repair, not a rediscovery of the
earlier top-four envelope.

The script hash is

```text
40c09ae5383d5243e22595dbba1304a5209e00a613730f55f7382cd6e25f8811.
```

Normal and optimized full replays produce the same stored report, whose hash
is

```text
1aff874ff7cfbe23c7fed5394f13d407d4143f94a2292b366ed7d0c8267432bc.
```

The combined coverage/triple/quadruple manifest is

```text
556e997f0fb7b270943611ab560ce72a32f077c97ed3e6354da8d6176d2c76e4. (B30)
```

An independent from-scratch endpoint/bisect audit reproduced all six
candidate counts and all twelve minima; its canonical consequence-ledger
digest is

```text
967ccc76a736256f074c7435046ff0a9d1812da9d369e19a108d9e7e82529901.
```

Therefore every pure tail of roots `12,13,14,15,16,17` is strictly lonely.
The earlier **Three-small exact closure** handles every non-pure completion:
if an added speed is at most `14`, the final family has at least three labels
from `{1,...,7}`.  Combining the six rows here with the preceding fifteen
whole-root certificates gives:

> **Complete flood/H-core stratum.** Every thirteen-speed family containing
> `H={8,...,14}` and at least two labels from `{1,...,7}` is strictly lonely.
> Equivalently, all `21/21` literal flood roots are closed.

Families containing `H` with only zero or one such label lie outside
THM-741's nine-in-window hypothesis.  This addendum closes one 21-body
sub-stratum, not all `2002` bodies.  Global THM-741 and LRC(14) remain open. ∎

## All-root pure-tail top-four atlas: 584 exact strata (codex-2026-07-28)

Let

```text
E range over all C(14,9)=2002 nine-subsets of {1,...,14},
G_E be its good set, m_E=|G_E|, r_E=#components(G_E),
c_E(w)=|G_E intersect D_w|,                    w>=15. (B31)
```

THM-735(ii), via the proved THM-731/732 discrepancy bound and the rational
majorant `sqrt(2)<99/70`, gives the strict all-speed estimate

```text
c_E(w)<u_E(w):=m_E/7+(99/70)r_E/(7w).                (B32)
```

This gives a general **finite-rank compactness lemma** for the arithmetic
tooth combs.  If `q_k(N)` is the `k`th largest value among
`c_E(15),...,c_E(N)`, and

```text
q_k(N)>m_E/7+(99/70)r_E/[7(N+1)],
```

then every global top-`k` coverage occurs at a speed at most `N`.  Indeed,
every unscanned tooth comb has strictly smaller coverage by `(B32)`.  The
infinite pure-tail obstruction is therefore a finite ranked head whenever
the chosen rank stays above the limiting density `m_E/7`.

For every root `E`, the companion evaluates all `586` exact coverages
`c_E(w)`, `15<=w<=600`, and sorts them as
`q_1(E)>=...>=q_586(E)`.  It then checks

```text
q_4(E)>m_E/7,
T_E:=(99/70)r_E/[7(q_4(E)-m_E/7)]<601.               (B33)
```

Thus `u_E(w)<q_4(E)` for every unscanned `w>=601`: the four highest values
in the finite bank are the four highest coverages over every integer
`w>=15`.  The largest crossing threshold over all `2002` roots is

```text
max_E T_E
 =39783744/67829
```

and is attained at
`E={1,2,4,8,10,11,12,13,14}`.

Define the exact top-four margin

```text
M_E=m_E-q_1(E)-q_2(E)-q_3(E)-q_4(E).                 (B34)
```

The exhaustive sign census is

```text
M_E>0: 584,             M_E=0: 0,             M_E<0: 1418.
```

For any one of the `584` positive roots and any four distinct added speeds
`15<=a<b<c<d`, the union bound and the global ranking give

```text
|G_E minus (D_a union D_b union D_c union D_d)|
 >=m_E-c_E(a)-c_E(b)-c_E(c)-c_E(d)
 >=M_E>0.
```

Hence every pure-tail extension of each of these `584` roots is lonely.
This already closes each root uniformly: if an arbitrary four-speed
extension contains a speed at most `14`, then its final thirteen-speed
family has at least ten speeds in `{1,...,14}` and is lonely by the proved
THM-738 three-slot theorem.  Therefore

> **The atlas plus THM-738 close 584/2002 whole THM-741 root bodies.**

Exactly ten of the `21` flood roots belong to this atlas set:

```text
(2,4),(3,4),(3,5),(3,6),(4,5),(4,6),(4,7),(5,6),(5,7),(6,7).
```

The other eleven flood roots are precisely

```text
(1,2),(1,3),(1,4),(1,5),(1,6),(1,7),
(2,3),(2,5),(2,6),(2,7),(3,7),
```

and were closed by the top-five, ranked-head, and residual-carrier addenda
above.  Consequently the current proved union contains exactly

```text
584+(21-10)=595
```

distinct whole THM-741 roots, leaving `1407` roots outside these
certificates.

The least positive pure-tail margin is

```text
47/1513512
```

at `E={1,2,4,5,6,8,11,12,14}`.  As a hostile boundary, the most negative
top-four margin is

```text
-9158777/174053880
```

at `E={1,3,4,5,7,8,10,11,13}`.  No conclusion is drawn for any of the
`1418` nonpositive rows.

The computation covers exactly `2002*586=1,173,172` rational entries.  Its
canonical coverage-manifest hash is

```text
63fd4a08965c9e5f2665dde27ae8db792cd446bc4578290ab4aa14b01fc469f7,
```

and the lexicographically ordered positive-body digest is

```text
93ed30b15748a90ea78aa1f392e2335cb284fa683f98e8d068b8a4f0f6af7a54.
```

The script and stored-output hashes are

```text
script 52448dce08d8e71149d2b19e4c9dd933c6c6b487ce5ec4e00d668c538ac3b7ce
output 6b28577862808f13782c18206cf0124783fbad7cf328a7e386ea514ce1e5e02e.
```

For all `8008` decisive root/speed pairs, sparse subtraction, full
subtraction, direct ten-comb union, and independent tooth-incidence summation
agree.  Each of the `584` extremal thirteen-speed families is additionally
rebuilt by both direct union and full nested subtraction.

An independent implementation classified endpoint breakpoints by midpoint
cells and used a closed periodic tooth primitive rather than the companion's
interval subtraction.  It reproduced the entire sign census, both digests,
all three extremal rows, and the strict global minimum

```text
min_E(q_4(E)-u_E(601))=196297/739641084>0.
```

The exact atlas itself proves only the stated pure-tail strata; the small-speed
chamber in the whole-root conclusion is imported from THM-738.  No conclusion
is drawn from the top-four test for the `1418` nonpositive roots, although
separate addenda close the entire flood sub-stratum.  This does not promote
the stale direct-run ledger or prove global THM-741 or LRC(14). ∎

## Sharded-runner integrity repair (codex-2026-07-17)

The live-pull inline sharding patch in `a00592c74` accidentally de-indented the
full-run block out of `main`; the module failed immediately with `NameError`
and the nominal per-shard writer still targeted the shared state file.  The
mathematical engine has therefore been restored byte-for-byte to its validated
hash `5aa81d9d...dac96`, which also restores every historical companion's hash
guard.  Sharding now lives in a separate orchestration wrapper.

That wrapper uses lexicographic-index sharding rather than process hashes,
writes one file per shard, rejects partial or malformed JSONL rows, coalesces
only identical duplicates, and emits a canonical merge only after all `2002`
bodies are present.  Its self-test proves the partitions for `1<=s<=16` are
disjoint, exhaustive, and balanced; six shards have sizes
`334,334,334,334,333,333`.  This repairs the computation path but does not
promote the unharvested `290/2002` ledger or add any unrun body claim.
