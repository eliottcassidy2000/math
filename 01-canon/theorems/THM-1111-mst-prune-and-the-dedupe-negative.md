---
id: THM-1111
title: THE PAIRWISE-OVERLAP (MST) PRUNE — valid, strong, and NOT sufficient; plus the residue/mask dedupe REFUTED. Both were my own named-next proposals from THM-1102, and they came out opposite ways. (I) THE DEDUPE IS DEAD: the certificate sees a killer only through its kill-set, so killers with equal kill-masks are interchangeable — but across ten sampled seven-speed cores the number of distinct masks equals the number of killers **exactly** (dedupe factor 1.000, 176/176 and 189/189). No two killers below KB are certificate-identical, and the proposed orders-of-magnitude saving is zero. (II) THE MST PRUNE IS VALID: for any sets, |⋃Aᵢ| ≤ Σ|Aᵢ| − Σ_{i≥2} max_{j<i}|Aᵢ ∩ Aⱼ|, and maximising the subtracted term over orderings is exactly the **maximum spanning tree** of the intersection graph, so coverage requires **Σ|Aᵢ| − MST_max ≥ n** — strictly stronger than the Σ|Aᵢ| ≥ n condition used at r=4 and r=5. (III) AND IT IS STRONG: of 1983 random sextuples that passed Σ ≥ n, **zero** passed Σ − MST ≥ n, i.e. at least a ~2000× further reduction on the random tail — which would bring r=6 from ~140 days toward hours. (IV) BUT IT IS NOT SUFFICIENT, and this is the honest verdict: against *adversarial* sextuples (top-by-size, greedy max-marginal, local search) the margin Σ − MST − n is **+2 at r=4, −2 at r=5, +36 at r=6**. Positive margins mean the bound does not rule out coverage, so the enumeration is reduced but not eliminated. My first read of (III) — that coverage might be provably impossible — was premature and is withdrawn here. (V) WHAT THE SEARCH DOES SHOW: the best sextuples found reach only **0.957, 0.932, 0.971** of n in actual union size at r = 4, 5, 6 — close to covering but never covering, which is consistent with the exhaustive r=4 and r=5 runs having found zero uncertified families
status: (I) REFUTED — measured, factor exactly 1.000 on every sampled core. (II) PROVED — the inequality is elementary (each Aᵢ beyond the first loses at least its best overlap with a predecessor; the ordering-optimal subtracted sum is the max spanning tree). (III) MEASURED on random sextuples, giving a LOWER bound on the reduction only, since a zero count bounds the rate from one side. (IV) MEASURED on a 40-core sample per r with a heuristic adversarial search — the positive margins are witnessed, so "not sufficient" is established; the exact maximum of Σ − MST is not. (V) MEASURED, and it is a search result, not a proof that no sextuple covers
source: kind-pasteur-2026-07-18-S128 (cont.65; owner: work the pairwise-overlap prune and the residue-vector dedupe)
depends_on:
  - THM-1102    # the r=6 enumeration wall, and the two proposals this tests
related:
  - THM-1071    # (III) the positive correlation that makes overlaps large is what powers the prune
script: 04-computation/prune_dedupe_kps_S128c65.py, mst_adversarial_kps_S128c65.py, mst_across_r_kps_S128c65.py (+ .out)
---

# THM-1111 — the MST prune, and the dedupe negative

## (I) The residue/mask dedupe is refuted

The idea was sound in principle: the certificate sees a killer only through its kill-set, so
killers with identical masks are interchangeable and could be enumerated once with
multiplicity. Measured across ten seven-speed cores:

| core | #killers | #distinct masks | dedupe factor |
|---|---|---|---|
| [1,4,6,8,9,10,12] | 176 | 176 | **1.000** |
| [1,3,4,6,7,10,11] | 189 | 189 | **1.000** |
| … (all ten) | | | **1.000** |

Every killer below KB has a distinct kill-mask. The saving is exactly zero, and the
sextuple count scales by 1.000⁶ = 1. Proposal dead.

## (II) The MST prune is valid

For any sets A₁,…,A_r and any ordering,

> |⋃ Aᵢ| = Σᵢ |Aᵢ \ ⋃_{j<i} Aⱼ| ≤ Σᵢ ( |Aᵢ| − max_{j<i} |Aᵢ ∩ Aⱼ| ).

The subtracted sum is the weight of a spanning tree (each i > 1 attaches to one
predecessor), and maximising over orderings gives exactly the **maximum spanning tree** of
the complete graph weighted by |Aᵢ ∩ Aⱼ| — Prim's order attains it. Hence

> **coverage of bits(P) requires Σ|Aᵢ| − MST_max ≥ n.**

This is strictly stronger than the Σ|Aᵢ| ≥ n condition that made the r=4 and r=5 runs
feasible, and it is exactly the place where THM-1071(III)'s positive correlation pays: the
overlaps are large, so MST_max is large.

## (III) It is very strong on the random tail

Of **1983** random sextuples that passed Σ ≥ n across six cores, **0** passed
Σ − MST ≥ n. That is at least a ~2000× further reduction of the tail, which on its face
would take r=6 from ≈140 days toward hours.

## (IV) But it is not sufficient — my first read was wrong

A zero count on random samples says nothing about the adversarial cases, and those are the
ones that matter. Searching deliberately (top-6 by kill-set size, greedy max-marginal
selection, and local search on the score):

| r | worst margin Σ − MST − n |
|---|---|
| 4 | **+2** |
| 5 | **−2** |
| 6 | **+36** |

Positive margins mean the MST bound does **not** rule out coverage, at r=4 and especially
at r=6. So the enumeration is reduced, not eliminated. I had briefly read (III) as evidence
that coverage might be provably impossible; that reading is withdrawn.

## (V) What the adversarial search does show

The best sextuples found reach, in actual union size,

| r | max |⋃Aᵢ| / n |
|---|---|
| 4 | 0.957 |
| 5 | 0.932 |
| 6 | 0.971 |

Close to covering, never covering. That is consistent with the exhaustive r=4 and r=5 runs
finding zero uncertified families, and it suggests the true obstruction is real but sits in
the gap between the MST bound and the actual union — the bound is loose by roughly the
margin above.

## Named next
- The gap in (V) is the object worth attacking: actual unions top out near 0.95n while the
  MST bound allows n + 36. A second-order correction — subtracting triple overlaps, or a
  fractional-relaxation bound rather than a spanning tree — would close much of that gap and
  could turn the prune into a proof.
- Practically, r=6 should now be attempted **with** the MST prune: the reduction is at least
  ~2000× on the tail, and the surviving adversarial sextuples are few enough to enumerate
  directly if they can be characterised (they are the ones with large, weakly-overlapping
  kill-sets, i.e. killers divisible by *different* small moduli).
- Do not spend further effort on mask dedupe.
