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
- ~~structural signature (opus-S120): every gap member = `(N−2)`-AP + 2 defects~~ **REFUTED
  (opus-S122):** N=7 has a 3-defect gap member `{1,3,4,5,7,13,18}` (`M=3/23`, longest-AP
  `{1,3,5,7}`=4=N−3). Both N=7 members share `M=3/23` (order 2) but have 2 vs 3 defects — so the
  *defect count is not the governing parameter; the order is* (two families, same order, different
  defect counts). The Freiman/defect framing was a mis-read of 3 examples.
- **RESIDUAL — the DEFECT-STRATIFICATION whittle (opus-S123, synthesizing S120 + mac-mini-S31 +
  kps-S41).** Stratify 12-families by `d = 12 − (longest sub-AP)`. The strata partition *all*
  families, and each is excluded from the open gap `(1/13, 2/25)`:
  - `d = 0`: the whole set is a dilated AP ⟹ `M = 1/13` (the boundary; the AP).
  - `d = 1`: dilated 11-AP + 1 outlier ⟹ mac-mini's ladder ⟹ `M ∈ {rungs} ∪ {plateau}`, min
    non-AP `= 2/25` (at `{1..11,24}`) ⟹ `M ≥ 2/25` (verified S123). **[open: prove the ladder bound]**
  - `d = 2`: dilated 10-AP + 2 outliers ⟹ `M ≥ 2/25` (verified S123; mac-mini's 2-outlier / +0.007
    margin). **[open: prove the 2-outlier bound]**
  - `d ≥ 3`: **GREEN — kps `LRCMod25Floor`** (`mod25_covering_floor`/`loose_of_mod25_covering`):
    a `≥3`-defect family is loose (`M ≥ 2/25`) via a `(ℤ/25)*` rotation (no-mult-of-25 branch) or a
    small-denominator clearance (mult-of-25 branch). The `25 = 5²` is `2N+1` at `N=12` and is the
    n-specific reason this closes at `12` but not at `7` (`2·7+1 = 15 = 3·5`, not a prime power —
    which is exactly why the `n=7` 3-defect member `{1,3,4,5,7,13,18}` slips through).
  So **(C) at `N=12` reduces to the `d=1` and `d=2` bounds** — a *finite-defect* residual, not the
  infinite-order gauntlet. The order view (kps HYP-4557) is a different, infinite slicing of the
  same set; the defect view is finite and has `d≥3` already GREEN.
- **TWO-MODULUS DICHOTOMY / sharpening of the mod-25 leg (mac-mini S32, HYP-4622).** The mod-25
  clearing rotation used above (the `d≥3` GREEN and kps-S41) works *exactly* when the family is
  **not** a full transversal mod 25: a clearing `c∈(ℤ/25)*` exists ⟺ the unit-speeds **miss** one
  of the ten `±`-classes; then the explicit witness is `t = a^{-1}/25` for a missed class `[a]`
  (⟹ `M ≥ 2/25`). So the mod-25 rotation branch is not "all no-mult-of-25 families" (that
  over-claims — the AP itself is a full transversal and is *not* cleared) but precisely the
  **non-transversals**, at *any* defect count. The genuine residual across all of `d=0,1,2` is the
  **full-transversal** families (= mac-mini S7 "gap member ⟹ full transversal mod 25"), for which
  the top wall is `M<2/25 ⟹ dilated AP` — a mod-13 rigidity (bottom wall). **0 counterexamples in
  ~15k structured + adversarial full-transversal families**; the AP is the unique full transversal
  with `M<2/25`. Frame: `2N+1=25` closes the top (miss a class), `N+1=13` pins the bottom (AP);
  a gap member must thread both and cannot. This makes the mod-25 hypothesis a *decidable residue
  condition* (Lean-ready against `LRCMod25Floor`) and isolates the `d=1,2` residual as the
  full-transversal-⟹-AP rigidity.

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

The single sentence: **LRC(14) closes when (C) closes**, and `(C)` now reduces (opus-S123) to the
**`d=1` and `d=2` defect bounds** — `d≥3` is GREEN (kps `LRCMod25Floor`), `d=0` is the AP. So the
finite residual is: *a dilated 11-AP + 1 outlier, and a dilated 10-AP + 2 outliers, are never in
the open gap at `N=12`* — both are ladder-law statements with `M ≥ 2/25`, verified and awaiting
proof. Everything above `(C)` — projection floor, rigidity lemma, J-K reduction — is GREEN,
provably clean, or a citation.

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
