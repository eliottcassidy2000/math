---
id: THM-1122
title: THE TRIPLE-OVERLAP CORRECTION AS A MOMENT LP — a real tightening that closes 722 of 792 r=6 cores but NOT all of them. (I) THE RIGHT FORMULATION: writing n_d for the number of points of bits(P) lying in exactly d of the kill-sets, the union is Σ_d n_d while the overlap sums are its moments — S₁ = Σ n_d·d, S₂ = Σ n_d·C(d,2), S₃ = Σ n_d·C(d,3). So an upper bound on the union is the linear program **max Σ_d n_d subject to those three equations and n_d ≥ 0**, and coverage requires that optimum ≥ n. With three equality constraints the optimum sits at a basic solution with ≤ 3 nonzero n_d, so it is solved EXACTLY by enumerating triples of d-values — no solver needed, exact rational arithmetic. (II) IT IS A GENUINE TIGHTENING, and the triple term is what does the work: on adversarial sets the worst margins are MST **−4 / −2 / +22** at r = 4/5/6, pairwise-only LP **−4 / +4.8 / +77.3** (i.e. *worse* than MST — moments without triples are weaker than a spanning tree), and the full triple LP **−10 / −10 / −2**. At r=4 and r=5 the LP is **exact** on the cases examined: LP₃ = 106 = the actual union, and 136 = the actual union. (III) BUT IT DOES NOT CLOSE r=6: run across ALL 792 seven-speed cores with an adversarial sextuple search, the worst margin is **+14.8** at core [1,4,6,7,8,9,12] (n = 248, LP₃ = 262.8), and **70 of 792 cores** still admit LP₃ ≥ n. My initial 14-core sample showed −2 and looked conclusive; the full run refutes it. (IV) WHAT IT IS WORTH ANYWAY: 722 of 792 cores are now certified outright with no enumeration, an ~11× reduction in cores, which composes with the MST prune's ≥2000× reduction on the tail — together these should bring r=6 from ≈140 days to well under an hour, making the enumeration finally practical
status: (I) PROVED — the moment identities are exact and the LP is a valid relaxation, so its optimum is a genuine upper bound on the union; solving at basic solutions is exact. (II) MEASURED, with the r=4/r=5 exactness observed on the cases examined rather than proved. (III) MEASURED over all 792 cores with a heuristic adversarial search — the 70 failing cores are witnessed, so "does not close r=6" is established; the exact maximum of LP₃ per core is not. (IV) is a projection from the two measured reduction factors, not a completed run
source: kind-pasteur-2026-07-18-S128 (cont.66; owner: work the triple-overlap correction and close the gap)
depends_on:
  - THM-1111    # the MST prune this tightens, and the gap it was asked to close
  - THM-1102    # the r=6 enumeration wall
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

Worst margins (bound − n) over adversarial searches:

| r | MST bound | LP with S₁,S₂ only | LP with S₁,S₂,S₃ |
|---|---|---|---|
| 4 | −4 | −4.0 | **−10.0** |
| 5 | −2 | +4.8 | **−10.0** |
| 6 | +22 | +77.3 | **−2.0** |

Two things are worth noting. First, the **pairwise-only** LP is *worse* than the spanning
tree at r = 5 and 6 — moments without the triple term throw away the combinatorial
information the MST was exploiting. Second, with S₃ included the LP became **exact** on the
r=4 and r=5 cases examined: LP₃ = 106 against an actual union of 106, and 136 against 136.

## (III) It does not close r=6

On a 14-core sample the r=6 margin was −2 and the correction looked like it had closed the
gap. Running the full **792 seven-speed cores** with an adversarial sextuple search:

> **worst margin +14.8**, at core [1,4,6,7,8,9,12] with n = 248 and LP₃ = 262.8
> **70 of 792 cores** still admit LP₃ ≥ n

So coverage is not ruled out everywhere, and r=6 is not closed by this bound. This is the
second time in two sessions that a promising sample was overturned by the wider run — the
pattern is consistent enough to be worth naming: **on this problem, a bound that looks
conclusive on a dozen cores has a real chance of failing on eight hundred.**

## (IV) What it is worth anyway

The correction is not a proof, but it is a large practical gain:

- **722 of 792 cores certified outright**, with no enumeration at all — an ~11× reduction in
  the number of cores that need any sextuple search.
- This composes with THM-1111's MST prune, which cut the per-core tail by ≥2000×.

Together these take r=6 from ≈3.6 × 10¹² sextuples (≈140 days) to a projected ~10⁸ over 70
cores — well under an hour. The enumeration that was infeasible two sessions ago is now
routine; it simply has not been run yet.

## Named next
- **Run r=6** on the 70 surviving cores with the MST prune. That is now a short job and it
  would close r=6 outright.
- **Push the LP further**: adding S₄ = Σ_{i<j<k<l}|Aᵢ∩Aⱼ∩A_k∩A_l| costs one more moment and
  one more constraint, and the jump from S₂ to S₃ was worth 79 units of margin at r=6. If S₄
  buys another 15 it closes the remaining 70 cores and r=6 needs no enumeration at all.
  This is the cheapest remaining lever and should be tried before the run.
- The 70 failing cores are worth inspecting as a set — if they share structure (they all
  contain 1, or all have a particular max), that structure is the real obstruction.
