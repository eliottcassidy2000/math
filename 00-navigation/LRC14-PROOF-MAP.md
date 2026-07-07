# LRC(14) proof map — the two routes, their obligations, and the tractable path

**opus-2026-07-06-S121.** Assembly-owner map reconciling the two LRC(14) proof threads, with the
critical correction that the "J-K reduction" is a **citation**, not an unbuilt bridge. Supersedes
the S120 "the gap thread is unwired" framing: the gap thread reaches the top level through a
*citable* published reduction, so the fleet's `(C)` work is on the LRC(14) path after all — it
just needs to be *wired* (as a cited hypothesis), exactly like LRC(≤13).

**Target.** `LRC14Statement := ∀ v : Fin 13 → ℤ, (∀ i, v i ≠ 0) → ∃ t, Lonely 14 v t`
(13 nonzero speeds, threshold `1/14`). LRC(≤13) is settled (owner directive) and enters as a
citation.

---

## Route 1 — witness / good-period density (the Lean skeleton's DAG)

`LRCFourteenSkeleton`, the "official" conditional DAG.

```
LRC14Statement
  ⟸ Mreach(v) ≥ 1/14  for all v            [GREEN: lonely_of_Mreach_ge, compactness handoff]
  ⟸ positive good-period density ρ*(shape) > 0   [OPEN: thm527_partA_density_pos_implies_reach
                                                   — the slow/fast change of variables + criterion C]
  ⟸ ρ* > 0 over all k=8..13 cluster shapes  [OPEN: the 1/7 global-witness floor; the old 2/7
                                             uniform-ρ* floor is REFUTED (kps zeros)]
  + sieve reductions                        [GREEN: no-mult-of-14, all-odd, saturated]
```

**Open critical path:** the `k = 8..13` witness/density floor — hard analysis (Riesz-product /
Bedert-style). Named Lean obligations: `thm527_partA_density_pos_implies_reach`,
`gK8_concentration_extremality`, `doublet_Rtail_uniform_bound`, the witness-floor cases.

---

## Route 2 — J-K reduction → rank-2 torus → 1-D gap (the recent-work thread; owner's intent)

This is the "LRC(14) → n=12 rigidity" route. Its top link is a **citation**.

```
LRC14Statement (near-extremal / rigidity behaviour)
  ⟸ [J-K REDUCTION — CITE Jain–Kravitz / Giri–Kravitz 2024:                    ← CITATION, not a proof obligation
       "Relative Lonely Runner Spectra" / "The structure of Lonely Runner spectra";
       the accumulation points of the LRC spectrum are governed by the relative
       spectra of 2-dimensional (rank-2) subtori]
  ⟸ (A)  no coupled proper rank-2 subtorus U has M(U) ∈ (1/13, 2/25)
  ⟸ (A) ⟸ (C):                                                                 [draft opus-S101]
       · projection floor                    [GREEN: LRCTorusProjection.torus_point_of_projection, S99]
       · pigeonhole rigidity lemma           [OPEN in Lean: "all full-support directions tight ⟹ rank≤1";
          2×2 core GREEN (LRCRankRigidity.dep_of_two_proportional, S102); the infinite-pigeonhole
          wrapper (finitely many Sym-12 orderings, infinitely many directions ⟹ two share) is unformalized]
  ⟸ (C)  the 1-D 12-speed Farey gap (THE CRUX): no integer 12-family has M ∈ (1/13, 2/25); only the
         dilated AP attains 1/13.
```

**Status of (C) (the crux), from the fleet:**
- canonical mediant family excluded — **PROVED** (mac-mini THM-632, parity; N=12 machine-checked).
- base+outlier, dilated-AP, order-2/3 species — **swept empty** at N=12 (mac-mini, kps ~150k families).
- structural signature (opus-S120): every gap member = `(N−2)`-AP + exactly 2 defects; the 2-defect
  world is closed.
- **RESIDUAL:** rule out gap members with `≥ 3` defects — the Freiman-stability step
  `M ∈ (1/13,2/25) ⟹ 10 of 12 form a dilated AP` (quantitative inverse-sumset; Fan–Sun/Freiman-3k−4).

**Open critical path:** (i) the `≥3`-defect Freiman-stability residual (the crux); (ii) the
infinite-pigeonhole rigidity wrapper; (iii) *wire* the J-K citation + (A)⟸(C) + (C) into a
top-level conditional theorem (parallel to Route 1's `lrc14_from_witness_floor`).

---

## Which route, and what to do

- **Route 2 is the owner's intent and where the momentum is**, and its top link (J-K) is a
  *citation*, so it is closer than S120 suggested. Its remaining *mathematical* obligation is one
  clean statement — the `(C)` Freiman-stability `≥3`-defect exclusion — plus two *formal* wiring
  tasks (the pigeonhole wrapper; the citation-conditional top-level theorem).
- **Route 1's** remaining obligation is the analytic density floor — genuinely hard, and the
  refuted 2/7 uniform floor shows it needs the corrected 1/7 global-witness estimate.
- **Recommendation.** Drive Route 2: (a) close (C) via the 2-defect signature + the ≥3-defect
  Freiman step [math, the crux]; (b) formalize the pigeonhole rigidity wrapper [S102 + Finite
  pigeonhole]; (c) build the citation-conditional top-level theorem `lrc14_gap_route` that reads
  `[J-K citation] + (A) ⟹ LRC14`, with (A)⟸(C) already reducible — so the fleet's (C) proofs
  finally register against `LRC14Statement`.

The single sentence: **LRC(14) closes when (C) closes**, because everything above (C) — the
projection floor, the rigidity lemma, and the J-K reduction — is either GREEN, provably clean, or
a citation. The crux is the `≥3`-defect Freiman-stability exclusion.

**Prior-work anchors:** Jain–Kravitz / Giri–Kravitz 2024 (the rank-2 accumulation reduction —
Route 2's top link); Fan–Sun arXiv:2306.10417 (the spectrum-gap gcd template for (C)); Bedert
arXiv:2511.16636 (Riesz products — Route 1's floor technique).

**Citation caveat (honest).** The J-K top link is used here at the *structural* level established
by web search: "the accumulation points of the LRC spectrum are governed by the relative spectra
of rank-2 subtori" (Jain, *Relative Lonely Runner Spectra*; Giri–Kravitz, *The structure of Lonely
Runner spectra*, arXiv:2304.01462). Before wiring it as a Lean citation-hypothesis, the exact
statement and the dimension bookkeeping (13-speed LRC(14) → rank-2 in the draft's `(ℝ/ℤ)¹²` → the
12-speed 1-D gap `(C)`) must be pinned against the paper — the PDF did not extract cleanly this
session. Treat "J-K reduces the LRC(14) gap to rank-2 / to (C)" as a well-supported lead to
confirm, not yet a certified citation, and keep the owner's LRC(≤13) citation policy in mind: a
formal caveat belongs only in the final write-up, not in the working DAG.
