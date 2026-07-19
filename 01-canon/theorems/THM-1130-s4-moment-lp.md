---
id: THM-1130
title: The S4 moment LP is an exact fixed-sextuple union bound; the supplied 792-core run is only an 18-sextuple-per-core sample, and its hardest family exhausts the q<=40 atlas
status: CORRECTED.  The S1--S4 identities and rational basic-solution LP are exact for each fixed sextuple.  The full script does traverse all 792 cores, but tests only 18 heuristic sextuples per core; it finds 70 cores with a positive LP3 sample and 5 with a positive LP4 sample, while negative sampled cores are not certified.  The former core-reduction, ~97%, trivial-enumeration, and likely-S5-closure claims are withdrawn.  One exact positive sample actually covers every q=15..40 core-safe obligation, so no deeper moment hierarchy on that same atlas can certify it; an enlarged chart supplies a witness immediately.
source: kind-pasteur-2026-07-18-S128 (cont.67; owner: try S4 first, then run r=6 on the 70 cores)
depends_on:
  - THM-1122    # the S₃ moment LP this extends
  - THM-1111    # the MST prune the moment ladder supersedes
script: 04-computation/lp4_r6_kps_S128c67.py, lp4_r6_chunk_kps_S128c67.py (+ .out)
---

# THM-1130 — the S₄ moment LP

> **All-scale scope audit (accepting codex-S73; MISTAKE-164).**  Everything below
> is stated inside THM-1102's *candidate* box `KB=333`, whose width-16 max-T scan
> is telemetry and not a proved uniform tail.  So the survivor counts describe a
> bounded computation, not a finite reduction of uniform `r=6`.  The moment-LP
> inequality itself (Section I) is unconditional and scale-free; only the
> per-core survivor census is box-relative.  Uniform `r=5` and `r=6` are both
> open, and per codex-S73 only `r<=4` is uniformly closed in this hierarchy.
> I accept this audit: my scaling-ratio arguments sampled and inferred all-scale,
> which is the same error I had been naming in others' samples without applying
> it to my own ratios.

> **Family-quantifier audit (codex-S67; THM-1122 correction propagated).**
> Each core receives only one top-cardinality, three greedy, and fourteen seeded
> random sextuples.  LP4 is exact for each of those 18 families, but a negative
> result says nothing about the untested sextuples.  The frozen main output did
> complete all 792 cores; the earlier statement that the core scan was partial
> was also wrong.  Counts below are positive-*sample* telemetry, never survivor
> counts or an exact reduction.


## (I) The extension

THM-1122 bounded the union by a three-moment LP. Adding the quadruple overlaps

> S₄ = Σ_{i<j<k<l} |Aᵢ ∩ Aⱼ ∩ A_k ∩ A_l|

gives a fourth equation, so the bound becomes

> **max Σ_d n_d subject to Σ_d n_d·C(d,j) = S_j for j = 1,2,3,4, n_d ≥ 0, d = 1..6.**

Four equality constraints put the optimum at a basic solution with at most four nonzero
n_d, so it is computed exactly by enumerating the fifteen 4-subsets of {1,…,6} and solving a
4×4 rational system. The whole bound remains solver-free and certifiable — S₄ costs fifteen
extra bit-ANDs per sextuple.

## (II) What the complete seeded sample actually says

The main run traverses all 792 cores and tests exactly 18 deterministic heuristic
sextuples at each core.  Its honest summary is

```text
positive LP3 sample in 70 cores,
positive LP4 sample in  5 cores,
negative sampled cores certified: 0.
```

The five LP4-positive sampled cores are

```text
[1,2,3,5,7,8,9], [1,2,3,5,8,9,11], [1,2,4,5,6,7,11],
[1,4,5,6,7,9,11], [2,3,5,6,7,9,11].
```

Thus S4 is empirically much sharper on this candidate generator.  It does not
remove 65 exact core branches, because no exact core branches were enumerated.

## (III) The ladder of bounds

| bound | worst sampled margin at r=6 | sampled scope |
|---|---|---|
| MST (spanning tree) | +22 | — |
| LP with S₁,S₂ only | **+77.3** (worse than MST) | — |
| LP with S₁,S₂,S₃ | +14.8 | positive sample in 70 cores |
| **LP with S₁,…,S₄** | **+20/9** | positive sample in 5 cores |

The shape is informative: the *pairwise* moment term alone is worse than the spanning tree,
because moments without higher terms discard the combinatorial structure the tree exploits.
On the seeded bank, S4 sharply improves S3.  That observation ranks fixed-family
certificates; it does not establish that moment depth is the global proof lever.

## (IV) The exact atlas wall

For the hardest sampled core and sextuple

```text
P = [1,2,3,5,7,8,9],
K = [234,300,253,168,216,211],
```

the `q=15..40` core-safe atlas has `n=190` obligations and exact moments

```text
(S1,S2,S3,S4) = (262,110,46,8).
```

The LP4 optimum is `1730/9`, giving margin `20/9`.  More importantly, the
actual union of the six masks has size exactly `190`: the family covers the
whole atlas.  Therefore S5, S6, or the complete multiplicity distribution on
this unchanged atlas cannot prove a safe obligation; deeper moments converge
to equality, not to a negative bound.  The full thirteen-speed row is primitive
and Covering, but `t=4/41` is already `1/14`-lonely (and THM-1121's universal
atlas also supplies `t=4/53`).  The missing coordinate is **atlas enlargement
or adaptation**, not another moment.

The hardened scripts keep `Fraction` margins exact for classification, seed
random candidates by global core index so chunk/main runs agree, and label the
five rows as positive samples.  Their outputs no longer use “survivor” or
“closes” terminology.

## Named next
- Treat the five positive-sample rows as heuristic stress cases, not a residual
  search space.
- Enlarge or adapt the rational-obligation atlas; the first exact wall already
  disappears at denominator 41.
- For a true finite reduction use THM-1135's harmonic ratio box, then prove a
  universal weighted/adaptive certificate on its bounded strata.
- `r=6` remains open.  Neither a sampled core traversal nor S5 on the old atlas
  can close it.
