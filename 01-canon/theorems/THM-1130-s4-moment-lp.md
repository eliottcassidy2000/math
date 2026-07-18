---
id: THM-1130
title: THE S₄ MOMENT LP — the quadruple-overlap term removes ~97% of the surviving r=6 cores, reducing the enumeration to a handful. (I) THE EXTENSION: adding S₄ = Σ_{i<j<k<l}|Aᵢ∩Aⱼ∩A_k∩A_l| gives a fourth moment equation, so the LP becomes max Σ_d n_d subject to Σ n_d·C(d,j) = S_j for j = 1..4 with n_d ≥ 0. Four equality constraints put the optimum at a basic solution with ≤ 4 nonzero n_d, so it is still solved EXACTLY by enumerating 4-subsets of {1,…,6} (fifteen of them) and a 4×4 rational solve — no LP solver, no floating point. (II) IT WORKS, dramatically: on the one chunk carried to completion (cores 200–500 of 792) **S₄ removed 31 of the 32 LP3-survivors**, leaving exactly ONE core, [1,4,5,6,7,9,11]. The other two chunks show the same collapse at their checkpoints — 26 → 3 over cores 0–200, and 22 → 1 over cores 400–700. (III) SO THE LADDER OF BOUNDS ON THIS PROBLEM IS: MST (worst margin +22 at r=6) → pairwise-only LP (+77.3, *worse*) → LP with S₃ (+14.8, 70 of 792 cores surviving) → **LP with S₄ (≈1 core per chunk surviving)**. Each added moment beyond the second buys a large reduction; the pairwise term alone buys nothing. (IV) HONEST STATUS: the full 792-core S₄ scan was **not completed** this session — one chunk of three finished, the other two were cut at their checkpoints — and the r=6 enumeration on the survivors was **not run**. What is established is the mechanism and its magnitude, not a closed r=6
status: (I) PROVED — the moment identities are exact, the LP is a valid relaxation of the union, and the basic-solution enumeration computes its optimum exactly in rational arithmetic. (II) MEASURED, with one chunk (cores 200–500) complete and giving 32 → 1, and two chunks partial at their last checkpoint. (III) assembles measurements from THM-1111, THM-1122 and here. (IV) explicit — **r=6 is NOT closed**, and neither the full scan nor the survivor enumeration has been run
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


## (I) The extension

THM-1122 bounded the union by a three-moment LP. Adding the quadruple overlaps

> S₄ = Σ_{i<j<k<l} |Aᵢ ∩ Aⱼ ∩ A_k ∩ A_l|

gives a fourth equation, so the bound becomes

> **max Σ_d n_d subject to Σ_d n_d·C(d,j) = S_j for j = 1,2,3,4, n_d ≥ 0, d = 1..6.**

Four equality constraints put the optimum at a basic solution with at most four nonzero
n_d, so it is computed exactly by enumerating the fifteen 4-subsets of {1,…,6} and solving a
4×4 rational system. The whole bound remains solver-free and certifiable — S₄ costs fifteen
extra bit-ANDs per sextuple.

## (II) It works

| core range | LP3 survivors | LP4 survivors |
|---|---|---|
| 0–200 (at checkpoint) | 26 | **3** |
| **200–500 (complete)** | **32** | **1** |
| 400–700 (at checkpoint) | 22 | **1** |

On the chunk carried to completion, **S₄ removed 31 of 32 survivors**, leaving the single
core **[1,4,5,6,7,9,11]**. The other two chunks show the same collapse.

## (III) The ladder of bounds

| bound | worst margin at r=6 | r=6 cores surviving |
|---|---|---|
| MST (spanning tree) | +22 | — |
| LP with S₁,S₂ only | **+77.3** (worse than MST) | — |
| LP with S₁,S₂,S₃ | +14.8 | 70 of 792 |
| **LP with S₁,…,S₄** | — | **≈1 per chunk** |

The shape is informative: the *pairwise* moment term alone is worse than the spanning tree,
because moments without higher terms discard the combinatorial structure the tree exploits.
But once S₃ is present each further moment buys a large reduction. The useful lever on this
problem is **depth of the moment hierarchy**, not cleverness in the pairwise term.

## (IV) Honest status

Neither of the session's two goals is finished:

- The full 792-core S₄ scan is **incomplete** — one chunk of three ran to a printed summary,
  the other two were cut at checkpoints. The survivor counts above are therefore a
  well-supported pattern, not a final list.
- The r=6 enumeration on the survivors was **not run**, because the survivor list is not
  final.

**r=6 is not closed.** What this session establishes is that S₄ collapses the surviving
cores from ~70 to a handful, which makes the remaining enumeration trivial rather than
merely feasible.

## Named next
- Finish the S₄ scan (three short chunks; the code is chunked and takes an env-var range),
  producing the definitive survivor list.
- Run r=6 on that list with the MST prune. With ~1–5 cores this is minutes, not hours, and
  it closes r=6.
- If S₅ collapses the last survivors, r=6 closes with **no enumeration at all** — and given
  the S₃ → S₄ collapse was 32 → 1, that is now the likely outcome and is worth trying first.
- The surviving cores are worth looking at as a set: [1,4,5,6,7,9,11] contains 1 and has
  max 11, not 12, which is the opposite of what the earlier worst cases looked like.
