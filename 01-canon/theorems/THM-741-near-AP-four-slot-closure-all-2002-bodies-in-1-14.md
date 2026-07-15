---
id: THM-741
title: NEAR-AP FOUR-SLOT CLOSURE — every 13-speed family with AT LEAST 9 speeds in {1,…,14} satisfies LRC(14). Equivalently, for EVERY 9-element body E ⊆ {1,…,14} (all C(14,9)=2002) and all v₁<v₂<v₃<v₄ not in E, {E,v₁..v₄} is lonely. Proof = the THM-735 Bonferroni tree at j=4: legs J4 (one inequality, all four ≥ V₁(E)) / J3 (per-v₁ exact bodies) / J2 (per-(v₁,v₂)) / J1 (per-(v₁,v₂,v₃) tail) / bottom (exact-ℚ sweeps of covering quadruples via lcm-multiples) — with PROVED P1/P2 LEMMA-SKIPS at every level (subtrees where the next Bonferroni threshold already fires from the parent's exact data close without computing the child body; sound because P1/P2 are one-level bounds off exact data)
status: CLAIMED (kind-pasteur-2026-07-13-S128 cont.5) — the OVERNIGHT run (owner-directed): multiprocessing over bodies, resume-safe per-body JSONL checkpoint (scratchpad; harvested to 05-knowledge/results on completion), probe-calibrated. Regression built in: the subtree (E={1..9}, v₁=10) must reproduce THM-738's body-{1..10} numbers (V=154, 143, 7537, 27) exactly. Upgrades to PROVED when all 2002 bodies close clean
source: kind-pasteur-2026-07-13-S128 (cont.5)
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
- [ ] overnight run launched (heavy-first, resume-safe JSONL, cpu−2 workers)
- [ ] all 2002 bodies clean; tight census; verdict

## Portable fixed-`E2` pruning addendum (codex-S14/S16): three flood bodies closed

The original driver hard-codes a Windows scratch path, and its reported
`171/2002` body ledger is not present in this checkout.  The portable
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
completed families uniformly over their two unbounded external speeds.  The
unresolved flood tail is now confined to at most three small labels,
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
`v3` geometry, bottom margins, and final-family identity.  The next direct
frontier has the `22` unshadowed three-small `K`-bases as vertices; runners,
Fano flags, and root charts discard that final-family quotient.  If a whole
edge must be run next, `(3,4)` remains the coverage-first choice: it gains five
three-small sets and has the least `E2` count among all gain-five edges.

## Ancestral-carrier envelope for the pure `(3,4)` tail (codex-S16)

Fix `E=H union {3,4}` and write the remaining ordered speeds as
`15<=a<b<c<d`.  Let `A` be an exact current survivor and let `C` be any
ancestor carrier with `A subset C`, measure `m_C`, and `r_C` components.  If
`s` ordered needles `w_i` remain, THM-732 gives the reusable one-sided bound

```text
|A minus union_i D_(w_i)|
 >= |A|-s m_C/7-(S2 r_C/7)sum_i 1/w_i.                  (B8)
```

If the first possible speed is `v`, the reciprocal sum is largest at
`v,v+1,...,v+s-1`.  Applying P2 to that first needle and charging the other
`s-1` needles to `C` also gives

```text
|A minus union_i D_(w_i)|
 >=6|A|/7-8r_A/(49v)-(s-1)m_C/7
   -(S2 r_C/7)sum_(i=1)^(s-1)1/(v+i).                   (B9)
```

Taking the maximum of (B8)--(B9) over all ancestors is sound.  It is a proof
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
