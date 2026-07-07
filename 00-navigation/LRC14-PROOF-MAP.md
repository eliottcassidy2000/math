# LRC(14) proof map — the two routes, their obligations, and the tractable path

**opus-2026-07-06-S121.** Assembly-owner map reconciling the two LRC(14) proof threads, with the
critical correction that the "J-K reduction" is a **citation**, not an unbuilt bridge. Supersedes
the S120 "the gap thread is unwired" framing: the gap thread reaches the top level through a
*citable* published reduction, so the fleet's `(C)` work is on the LRC(14) path after all — it
just needs to be *wired* (as a cited hypothesis), exactly like LRC(≤13).

> ## ⛔ CORRECTION BANNER (opus-2026-07-06-S130) — READ FIRST: ROUTE 2 HAS TWO BROKEN LINKS
>
> A careful correctness audit found that **Route 2 does not prove LRC(14)**, at BOTH ends:
>
> 1. **TOP LINK INVALID (MISTAKE-117).** The "J-K reduction" is NOT a valid published reduction.
>    Giri–Kravitz (arXiv:2304.01462) study the **accumulation points** of the spectrum
>    (`acc(S(n)) = S(n-1)`), NOT the **supremum** that the LRC bounds; the abstract says verbatim
>    *"Rather than attack this conjecture, we study the structure of the sets S(n)."* Controlling
>    rank-2 subtori (accumulation-point data) does NOT bound the sup. So `(A) ⟹ LRC(14)` is
>    UNWARRANTED — Route 2 is **disconnected from LRC(14) at the top**; even a full proof of
>    `(C)`/`(A)` would not prove LRC(14) via this citation. The S121 claim below ("the J-K reduction
>    is a citation, not an unbuilt bridge… the gap thread reaches the top level through a citable
>    published reduction") is **RETRACTED**.
> 2. **BOTTOM LINK NOT FINITE (MISTAKE-116).** `(C)` does NOT reduce to a finite covering `{2..Q0}`
>    (mac-mini-S36, verified opus-S130): compressed varying-k families `≡ AP mod lcm(2..Q0)` escape
>    every `q ≤ Q0` and clear only at `nextprime(Q0)`; the covering modulus is UNBOUNDED. "Every
>    non-AP clears at some `q`" is EQUIVALENT to the analytic `(G)`, not a finite reduction.
>
> **What still holds.** `(C)` is TRUE (verified opus-S130: unique AP at `1/13`, empty gap). The Lean
> is sound (all valid conditional implications / honest `Prop` obligations). The rank-2 rigidity +
> torus projection + covering-endpoint theorems are correct math. What is RETRACTED is the FRAMING
> that LRC(14) was "nearly closed via 3 obligations." **Reroute:** Route 1 (bound `Mreach ≥ 1/14`
> DIRECTLY — the sup) is the correctly-aimed project route; the recognized external route is Tao's
> finite reduction (2018) + Malikiosis–Santos–Schymura (2025) + Rosenfeld/Trakulthongchai sieving.
> The Route-2 material below is retained as correct conditional math / a spectrum-structure study,
> NOT as a proof route. See MISTAKE-116, MISTAKE-117.

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

> **Structured-branch tool (kps-S52/S53, HYP-4697/4707) — a *sup* fact, so it lives on Route 1. FORMALIZED GREEN.**
> The **coarse / scale reduction** `M(v) ≥ M(K) − A/L` (for `vᵢ = aᵢ + L·kᵢ`, all `kᵢ ≥ 1`,
> `|aᵢ| ≤ A`, `K = {distinct kᵢ}`) grounds the **multi-scale** families directly in the
> **settled LRC(≤13)**: if the 13 speeds cluster into `≤ 12` groups at any scale `L`
> (equivalently, carry a close pair at scale `L`), then `M(K) ≥ 1/13` ⟹ `M(v) > 1/14` for
> `L > 182·A`; `13`-distinct-cluster families **descend** (smaller height). **Now machine-checked
> (kps-S53): `LRCCoarseReduction.lean` — `reach_transfer_coarse` + `lonely14_of_coarse_le12`,
> kernel-pure, in the manifest.** So the witness/density floor is only needed for the
> **single-scale / compressed** families — the same bounded-ratio domain as **Tao (2018)**'s
> finite reduction. Over that residue the floor is a **dichotomy** (kps-S53
> `lrc_singlescale_density_floor`): the AP `{1..13}` is the unique tight family (`ρ = 0`,
> isolated witness, `M = 1/14`) closed by **rigidity/three-gap**, and the spread bulk (`ρ > 0`,
> perturbations jump to `M ≥ 1/13`) closed by **decorrelation** — two problems, not one uniform
> estimate. Survives opus-S130 because it bounds the **sup**, not accumulation points.
> Reflections `the-coarse-reduction-is-a-sup-fact-…-kps-S52`,
> `…-formalized-and-the-density-node-is-reduced-to-single-scale-kps-S53`.

> **↳ ESCAPE FAMILIES COLLAPSE TO THE MOAT (klein-S152, HYP-4711) — 2.5 open pieces → 1.**
> The coarse bound fails ONLY when the coarse part `K` is the AP (`M(K)=1/14`, no slack): the
> perturbed dilated AP `vᵢ = aᵢ + L·(d·i)` (mac-mini's S36/S37 escape / L-lift family; the r=13
> crux kps's "descend" can't close — descent loses `A/L` at the slackless AP). **NEW:** these carry
> an explicit **AP-inherited conjugate witness** — the dilated AP's `φ(14)=6` witnesses `t_c=c/(14dL)`
> each bind one antipodal pair; a shift `δ=O(A/L²)` keeps the pair `≥1/14` iff
> `a_{i₊}/v_{i₊} ≥ a_{i₋}/v_{i₋}`, and the conjugate `c↦14−c` flips that test, so one always works ⟹
> `M(v) ≥ 1/14` uniformly (verified 200/200 base; slope predicts winner 100/100; permuted 200/200;
> the families are in fact LOOSE, true `M≈0.1–0.25` — `a≠0` decorrelates — but the witness gives the
> `1/14` floor without needing that lift). So the escape families are **not** a separate decorrelation
> obstruction; answers mac-mini-S38's open item (a) (the "sharper base-structure bound"). **The whole
> LRC(14) residual is now ONE object — the MOAT:** `{1..13}` is the unique single-scale 13-family
> (up to dilation) with `M<1/13` (the 13-family `(C)`-analog, used directly as a sup bound = Route 1).
> Not a proof (moat open; witness verified, not yet formal — needs `L₀(A)≈200A`). Reflection
> `the-escape-families-collapse-to-the-moat-the-ap-carries-its-own-witness-klein-S152`.

> **✅ ROUTE 1 IS THE LIVE ROUTE (owner directive, S130) — and its density floor is ROBUST (opus-S130).**
> After the Route-2 audit, the owner directed effort here. The density floor was worked
> correctness-first and found genuinely sound (NOT a 2/7-style artifact). The load-bearing quantity is
> `μ_{1/7}(E) = meas{x : max-gap of {frac(e·x) : e∈E} > 1/7}`, and it reduces to:
> (i) **AP-minimality** `μ_{1/7}(E) ≥ μ_{1/7}({1..k})` — strongly verified (40 aggressive adversarial
>     descents + structured adversaries, NONE below the AP; `μ_{1/7}` is dilation/translation
>     invariant). Proof = three-gap equidistribution (the AP orbit is maximally spread). OPEN lemma.
> (ii) **Exact AP constants** via three-gap piecewise-linear breakpoints — RIGOROUS (opus-S130,
>     `lrc_mu17_exact_threegap`): `μ_{1/7}({1..k})` = 1 (k≤7, pigeonhole), 691/735, 247/294, 38/49,
>     1381/2205, 13823/24255, **477/1078** (k=8..13). The `k=13` value **exactly matches the canon's
>     `rhoGlobFloorRat(13)=477/1078`**. Min over the relevant range `≈ 0.44 ≫ m_P=0.0565`.
> WHY 1/7 WORKS (where 2/7 failed): `1/7≈0.143` sits well below the typical max-gap of `k` well-spread
> points (`~H_k/k≈0.34`), so the good-set stays a large majority — a structural reason, not luck.
> REMAINING (honest, both bounded analysis on a correctly-aimed route — no wrong-object flaw):
> (A) prove AP-minimality (three-gap lemma); (B) the finite-`Vmax` error budget `O(#arcs/Vmax)` for
> Part A (`LRCWitnessPartA` has the arithmetic glue). Scripts: `lrc_witness17_floor_probe`,
> `lrc_ap_minimizes_mu17`, `lrc_mu17_exact_threegap` (all `_opus_S130`).
>
> **Convergence (kps-S53 + opus-S130):** kps's coarse reduction sends multi-scale families to
> LRC(≤13), leaving the **single-scale** residue where the floor is a *dichotomy* — near-AP
> (rigidity/**three-gap**) + spread (decorrelation). My `μ_{1/7}` three-gap work IS that near-AP
> branch (AP-minimality + exact constants); mac-mini-S38's `reach_decorr` handles the spread/escape
> branch. The pieces fit.

---

## Route 2 — J-K reduction → rank-2 torus → 1-D gap (RETIRED as a proof route, opus-S130; correct conditional spectrum math)

This is the "LRC(14) → n=12 rigidity" route. Its top link is a **citation**.

```
LRC14Statement (near-extremal / rigidity behaviour)
  ⟸ [J-K REDUCTION — CITE Jain–Kravitz / Giri–Kravitz 2024:                    ← CITATION, not a proof obligation
       "Relative Lonely Runner Spectra" / "The structure of Lonely Runner spectra";
       the accumulation points of the LRC spectrum are governed by the relative
       spectra of 2-dimensional (rank-2) subtori]
  ⟸ (A)  no coupled proper rank-2 subtorus U has M(U) ∈ (1/13, 2/25)
  ⟸ (A) ⟸ (C):                                                                 [WIRED opus-S129: LRCRoute2Assembly.torus_loose_of_rank2, GREEN]
       · projection floor                    [GREEN: LRCTorusProjection.torus_point_of_projection, S99; consumed via LRCTorusReduction.torus_loose_of_loose_direction]
       · pigeonhole rigidity lemma           [GREEN: 2×2 core LRCRankRigidity.dep_of_two_proportional (S102) AND
          the infinite-pigeonhole wrapper LRCRankRigidity.dep_of_infinite_common_proportional (S102, finitely many
          Sym-12 orderings, infinitely many (1,N) directions ⟹ two share) — BOTH formalized]
       · the C-bridge                         [the sole OPEN input: "a not-loose (1,N) direction is (after centering)
          proportional to its finite dilated-AP ordering vector" = (C) in wrapper-ready form; carried as `CBridge` Prop.
          Centering is what makes the classifier finite (Sym 12, not Sym 12 × ℚ).  S129 composes floor+wrapper: rank-2 + CBridge ⟹ torus loose]
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

> **Refinement (klein-S144, HYP-4611, validated ~140k families).** The node is *not* uniform over
> *all* heights, and it should not be: a family `≡ AP mod lcm(2..39)` is obstructed at every `q ≤ 39`
> (matches the AP's residues) — a genuine covering gap — but its entries differ from `{1,…,12}` by
> `~10¹⁶`, so it is **non-compressed, carries a far element, and peels** (THM-620 / THM-608) before
> reaching `(C)`. The correct statement is **uniform over COMPRESSED families** (`max ≤ 13·min`, the
> post-peel `(C)` domain): those clear at **`q ≤ 31`, 0 gaps to height 650,000** (a compressed family
> cannot be `≡ AP mod` a large `L` without a far entry). So the node = *compressed non-AP ⟹ cleared
> at `q ≤ 39`* (covering) **⊕** *non-compressed ⟹ peels* (composition). NB kps-S44's circulating
> "min-clear-mod ≤ 14" is too low — the covering **set** needs `q` up to `38` (CRT-lifts `≡ AP mod
> lcm(2..12,25)`); `31` on compressed. Scripts: `lrc14_covering_uniformity_klein_S144.py`,
> `lrc14_covering_compressed_uniformity_klein_S144.py`.

> **⚠ CORRECTION (mac-mini-S36/S37, HYP-4667, verified).** klein-S144's key step —
> *"a compressed family cannot be `≡ AP mod` a large `L` without a far entry"* — is **FALSE**.
> Counterexample (verified at `L = lcm(2..39)`): `V = {i + L·kᵢ : i=1..12}` with **all** `kᵢ ≥ 1`
> and varying (e.g. `k = (1,2,1,2,…)`). Then `V ≡ AP mod L` (obstructed at every `q ≤ 39`), yet
> `V` is **compressed** (`max/min = 2 ≤ 13`) with **no far entry** — every entry is lifted to `~L`,
> so they cluster at one scale. klein assumed *some* entries stay `~1` (a far element that peels);
> lifting *all* entries uniformly keeps it compressed. As a 13-lift `vᵢ = i + 13·Kᵢ`, the
> `Kᵢ ~ L/13 ~ 10¹⁴` are huge yet the family is compressed — *"compressed"* bounds the lift **range**,
> not its **values**. This family is non-DilatedAP, compressed (does **not** peel by kps-S49's own
> `max > 13·min` criterion), fails all `q ≤ 39`, and clears **only at `q = 41 = nextprime(39)`**
> (loose, `M ≥ 3/41 > 2/25`, so (G) holds). **General:** for every `Q₀`, the `≡ AP mod lcm(2..Q₀)`
> varying-`k` families escape `{2..Q₀}` and clear at `nextprime(Q₀)` — the covering modulus is
> **unbounded**. So `CoveringComplete` as wired (`∃ q`, no bound) is **`= (C)` exactly**, an honest
> open obligation, **not** a finite `q ≤ 39` residue check — the latter is provably impossible.
> The covering is a **reformulation** of the crux, not a reduction; `(C)` = *"every non-AP 12-family
> has `M ≥ 2/25`"* remains the open **analytic** core (the escape families approach `2/25⁺`, so any
> proof must be tight in the margin). Scripts: `lrc_escape_verify` / `lrc_covering_regress` /
> `lrc_escape_at_Q39` `_macmini_S36/S37`. Reroute: analytic (density floor / theta-sum Cohn-Elkies /
> Fan-Sun's *actual* n=4 method if it extends), not a finite covering pile.

**Open critical path for the top level (updated opus-S129).** (i) DONE — the pigeonhole rigidity
wrapper was already GREEN (`dep_of_infinite_common_proportional`, S102), and (A)⟸(C) is now WIRED
(`LRCRoute2Assembly.torus_loose_of_rank2`, S129: rank-2 + `CBridge` ⟹ torus loose, composing the
GREEN projection floor + GREEN rigidity wrapper). (ii) DONE — `[J-K] + (A)⟸(C) + (C)` is wired into
the top-level conditional `LRCRoute2Assembly.lrc14_via_route2` (`JKReduction` + `CBridge` ⟹
`LRC14Target`), parallel to Route 1's `lrc14_from_witness_floor`. The remaining OPEN inputs are now
just two `Prop` obligations: the **`CBridge`** (= (C) in centered wrapper-ready form — discharge via
the covering system + a centering lemma) and the **`JKReduction`** citation (pin the exact
Jain–Kravitz dimension bookkeeping against the paper). (iii) still pending: pin the exact J-K
statement.

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
