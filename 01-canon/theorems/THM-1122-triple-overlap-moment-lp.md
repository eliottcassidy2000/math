---
id: THM-1122
title: THE TRIPLE-OVERLAP CORRECTION AS A MOMENT LP — an exact per-family union bound plus a seeded candidate audit, not a core closure. Writing n_d for the number of points lying in exactly d kill-sets gives exact moments S_j = Σ_d n_d C(d,j), and maximizing Σ_d n_d subject to S₁,S₂,S₃ gives a rigorous upper bound for each fixed family. The supplied all-core script nevertheless tests only 18 heuristic sextuple candidates per core. Its 722 negative sampled cores are not certified, while its positive sample at [1,4,6,7,8,9,12] proves that the LP₃ bound alone cannot uniformly close r=6.
status: (I) PROVED — the moment identities, LP relaxation, and exact basic-solution computation. (II) MEASURED only on seeded core/family samples. (III) CORRECTED 2026-07-18 — the original "722 of 792 cores certified" and projected residual-enumeration claims were false because the candidate generator is not exhaustive. The audit covers all cores but not all sextuples within a core.
source: kind-pasteur-2026-07-18-S128 (cont.66; owner: work the triple-overlap correction and close the gap)
depends_on:
  - THM-1111    # the MST prune this tightens, and the gap it was asked to close
  - THM-1102    # the r=6 enumeration wall
  - THM-1121    # exact finite-horn closure; supersedes the proposed residual enumeration
script: 04-computation/moment_lp_kps_S128c66.py, lp3_r6_allcores_kps_S128c66.py (+ .out)
---

# THM-1122 — the triple-overlap correction

## (I) The formulation

The gap in THM-1111 was that the MST bound charges each set only its best *pairwise*
overlap with a predecessor, so a point lying in four or five kill-sets is under-counted. The
fix is to stop reasoning about the sets and reason about the **multiplicity distribution**.

Let n_d = #{x ∈ bits(P) : x lies in exactly d of A₁,…,A_r}. Then

> |⋃Aᵢ| = Σ_d n_d,  S₁ = Σ_d n_d·d,  S₂ = Σ_d n_d·C(d,2),  S₃ = Σ_d n_d·C(d,3),

where S₁ = Σ|Aᵢ|, S₂ = Σ_{i<j}|Aᵢ∩Aⱼ|, S₃ = Σ_{i<j<k}|Aᵢ∩Aⱼ∩A_k| are all directly
computable from the masks. So an upper bound on the union is the LP

> **maximise Σ_d n_d subject to the three moment equations, n_d ≥ 0, d = 1..r**,

and **coverage requires that optimum ≥ n**. Three equality constraints means the optimum is
attained at a basic solution with at most three nonzero n_d, so it is computed exactly by
enumerating triples of d-values and solving a 3×3 rational system — no LP solver, no
floating point.

## (II) The triple term is what does the work

Displayed margins (bound − n) over the seeded heuristic candidates in the first script:

| r | MST bound | LP with S₁,S₂ only | LP with S₁,S₂,S₃ |
|---|---|---|---|
| 4 | −4 | −4.0 | **−10.0** |
| 5 | −2 | +4.8 | **−10.0** |
| 6 | +22 | +77.3 | **−2.0** |

These are observations on the displayed candidates, not extrema over all cores and all
families. Two things are worth noting. First, the **pairwise-only** LP is *worse* than the spanning
tree at r = 5 and 6 — moments without the triple term throw away the combinatorial
information the MST was exploiting. Second, with S₃ included the LP became **exact** on the
r=4 and r=5 cases examined: LP₃ = 106 against an actual union of 106, and 136 against 136.

## (III) It does not close r=6

On a 14-core sample the r=6 margin was −2 and the correction looked like it had closed the
gap. The second script then visited all **792 seven-speed cores**, but tested only one
top-cardinality, three greedy, and fourteen seeded-random sextuple candidates per core:

> **worst margin +14.8**, at core [1,4,6,7,8,9,12] with n = 248 and LP₃ = 262.8
> **70 of 792 cores** have a sampled candidate with LP₃ ≥ n

Thus the LP₃ bound by itself cannot uniformly reject every sextuple: the displayed positive
candidate is already a witness to that methodological limitation. Conversely, LP₃ < n for
all 18 sampled candidates at one of the other 722 cores says nothing about the untested
sextuples, so it does **not** certify that core. This is the
second time in two sessions that a promising sample was overturned by the wider run — the
pattern is consistent enough to be worth naming: **on this problem, a bound that looks
conclusive on a dozen cores has a real chance of failing on eight hundred.**

## (IV) What remains useful

The exact moment formulation remains a valid certificate for any fixed sextuple and a useful
way to rank candidate families. The sampled distribution — 70 cores with at least one
nonclosing candidate and 722 without one among the 18 tested — is honest reconnaissance.
It is not a reduction in the exact search space. The finite horn that motivated that proposed
search is instead closed uniformly by THM-1121's exact weighted atlas.

## Named next
- If this lens is pursued, either enumerate all sextuples in a stated finite domain or derive
  a universal maximization theorem; sampled negative margins cannot discharge a core.
- S₄ or a hybrid moment/MST/weighted-dual bound may sharpen the fixed-family certificate,
  but it needs the same quantifier audit.
- Use the 70 positive-sample cores as a heuristic structure bank, not as an exact survivor set.
