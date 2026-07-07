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

**Status of (C) — the current synthesis (opus-S126, superseding the layered history below).**

`(C)` has converged, across all agents, to a single statement and a single *finite* residual:

> **`(C)` ⟺ the AP `{1,…,12}` is the unique 12-integer-family (up to dilation) with `M < 2/25`.**
> Equivalently (kps-S42, via the settled `LRC(≤13)` floor `M ≥ 1/13`): no family attains an
> `order-k ≥ 2` value in the gap — both gap edges `1/13, 2/25` are `k=1` Kravitz rungs, and every
> interior value `3/38, 4/51, …` has `k ≥ 2`.

**The three-case split (mac-mini-S32 pair-blocking, from kps-S41's mod-25 core; = opus-S124's
mod-25 dichotomy).** Split 12-families by their residues mod 25:

1. **NON-blocker** — the unit-speeds miss one of the ten antipodal `±`-pairs mod 25 (`⟺` a
   clearing rotation `c ∈ (ℤ/25)*` exists, `c = a⁻¹` for a missed pair `{a,−a}`). Then
   `M ≥ 2/25` at `t = c/25`. **GREEN** — kps `LRCMod25Floor` (`loose_of_mod25_covering`) + mac-mini
   `LRCMod25Transversal` (THM-634, the explicit miss-a-pair witness).
2. **BLOCKER (full transversal)** — no mult of 25, and the `±`-residues cover all ten pairs
   (`= (ℤ/25)*`). **THE RESIDUAL.**
3. **mult of 25** — a speed `≡ 0 mod 25` sits at residue 0 for every rotation; clears at a small
   denominator (`M ≥ 2/11, 2/17, …`). **EASY.**

**Case 2 is a FINITE COVERING SYSTEM (kps-S43) — the key reduction.** The blockers are
*defect-agnostic* (they span every defect count `d ≥ 1`, not just `d = 1,2` — correcting opus-S123's
"`d≥3` GREEN via mod-25": there are `d≥3` blockers, e.g. `{1,2,3,4,6,7,8,9,10,11,13,55}`, `d=5`,
that mod-25 does *not* clear). Every non-AP blocker has `M ≥ 1/12`, and:

> **a finite set of moduli `q ∈ {6,…,39}` clears every non-AP blocker** — verified 0 uncleared of
> 27 218 (sample, height ≤ ~110). Each clearance is a `rational_point_margin` certificate at
> `t = c/q` (the same Lean atom as `LRCMod25Floor`, just at `q` instead of 25).

So **`(C) = case 1 (GREEN) + case 3 (easy) + case 2 [finite `q ≤ 39` covering + the AP exception]`**,
and *every branch is a margin certificate.* The AP is the unique uncovered blocker because it is
the global `M`-minimizer (`1/13`, unique since `13` is prime = the tight locus), so it has no slack
at any modulus; every other blocker (`M ≥ 1/12`) has a clearing `q ≤ 39`. **The crux is now a
finite, Lean-ready covering system — not an analytic rigidity.** (opus-S125's two-modulus factoring
— `q=13` collision-clears the bottom, `q=25` clears the top — is the special case `q ∈ {13, 25}`;
kps-S43's full `q ∈ {6..39}` achieves 0 residual.)

**Formalized so far (all GREEN, kernel-pure):** `LRCMod25Floor` (kps, case-1 core),
`LRCMod25Transversal`/THM-634 (mac-mini, case-1 miss-a-pair witness), `LRCLadderD1`/THM-633
(mac-mini, the `{1..11}+x` ladder — a sub-family of case 2), `LRCBinderInfeasible` (opus, the
mediant `k=2` parity gate), `LRCSubfamilyCap` (opus, the plateau `M ≤ M(subfamily)`).

**Open critical path for `(C)`:** prove the covering is **uniform over all heights** — every
non-AP blocker clears at some `q ≤ Q₀` (`39` on the sample). This is a *finite* residue condition
(clearing at `q` depends only on `v_i mod q`), an **Erdős-covering-system-flavored** statement, not
analysis — plus the AP-exception (immediate from `M`-minimality) and the easy case 3. Then wire the
`q≤Q₀` covering + case 1 + case 3 into a Lean theorem `M(V) < 2/25 → V = dilated AP`.

**Open critical path for the top level:** (i) the pigeonhole rigidity wrapper (finishes (A)⟸(C));
(ii) *wire* `[J-K citation] + (A)⟸(C) + (C)` into a top-level conditional theorem (parallel to
Route 1's `lrc14_from_witness_floor`); (iii) pin the exact Jain–Kravitz statement.

<details><summary>Superseded layered history (S120–S125)</summary>

- canonical mediant family excluded — PROVED (mac-mini THM-632, parity).
- opus-S120 "gap member = (N−2)-AP + 2 defects" — REFUTED (opus-S122: 3-defect member
  `{1,3,4,5,7,13,18}` at N=7). Defect count does not govern; order (or the blocker/non-blocker
  split) does.
- opus-S123 defect stratification (`d=0/1/2/≥3`) — subsumed by the blocker split; its "`d≥3` GREEN
  via mod-25" attribution was imprecise (kps-S43: `d≥3` blockers exist, clear at small denom, not
  mod-25).
- opus-S124 mod-25 dichotomy = the case-1/case-2 split; opus-S125 two-modulus = `q ∈ {13,25}` slice
  of kps-S43's covering.
</details>

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

The single sentence: **LRC(14) closes when (C) closes**, and `(C)` now reduces (opus-S124) via
the **mod-25 dichotomy** to a *single* residual — the near-AP moat:
- **cleared** (some `c ∈ (ℤ/25)*` has all `v_i·c mod 25 ∈ [2,23]`) ⟹ `M ≥ 2/25`, LOOSE — **GREEN,
  now with the existence half too**: kps `loose_of_mod25_covering` handed the rotation `c`;
  mac-mini `LRCMod25Transversal.loose_of_misses_pair` (THM-634, S33b) *produces* it — from the
  decidable hypothesis "misses an antipodal pair `{a,−a}` mod 25", the explicit `c = a⁻¹`, `t =
  a⁻¹/25`. So branch (a) — *every non-transversal family is loose, at any defect count* — is fully
  machine-checked. (This also corrects the "d≥3 GREEN" filing: ~2% of `d≥3` families are
  transversals, not cleared by the rotation; the clean line is transversal / non-transversal, i.e.
  kps-S43 "defect-agnostic".)
- **non-cleared** (⟺ `±{v_i mod 25}` covers all 20 units of `(ℤ/25)*`) ⟹ `M = 1/13` (AP) or
  `M ≥ 1/12` (plateau) — **the residual (b)**, verified over 50k families, `0` in the gap. The
  `d=0` (AP boundary) and `d=1` (`{1..11,x}`, THM-633) slices are done; the open piece is the
  saturated `d ≥ 2` plateau `M ≥ 1/12`.

Neither branch meets `(1/13, 2/25)`. So `(C)` = **(b): a mod-25-saturated 12-family is the AP or
the plateau** — the AP-rigidity heart, with the "spread" families peeled off GREEN, and the target
pinned to a finite residue-covering condition. (This subsumes the S123 defect stratification: one
dividing line instead of four strata.) Everything above `(C)` — projection floor, rigidity lemma,
J-K reduction — is GREEN, provably clean, or a citation. Formalize (b) via mac-mini's ladder +
opus-S115's subfamily-cap plateau, restricted to mod-25-saturated families, and `(C)` closes.

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
