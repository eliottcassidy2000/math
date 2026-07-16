# LRC(14) proof map — the two routes, their obligations, and the tractable path

> ## 2026-07-16 RESONANT-PEEL CHECKPOINT (codex-S17)
>
> `THM-891` proves the exact fixed owner-resonant limit and closes five of six residue
> classes below the current `0.097` benchmark via a 21-ray pair-sector theorem and exact
> quadratic certificates.  The sole limiting sign is negative residue 6, supported on
> miss pairs `{1,5}` and `{2,4}`.  This does **not** close the top-level map: the
> core-uniform finite-`t` wall remainder, all-`w` propagation, compact bands, and formal
> assembly remain.  Next exact target: a `THM-890`/`THM-894`/`THM-896-level3-crossing` level-three
> relation certificate for `A_15+A_24`, followed by the wall remainder.  See `HYP-7024` and the
> 2026-07-16 propagation-ledger addendum.

> ## ✅ STATE OF THE PROOF — 2026-07-09 (boxeph-S1 synthesis banner; supersedes the stale top matter below, which is kept for history)
>
> **Architecture (THM-663 single chain, all agents converged):**
> `LRC(14) = [non-covering: q-witness sieve THM-369 + LRC(≤13) citation, exact t=1/q, equality allowed] + [covering: strict cushion]`.
> Covering case = **[density floor: CLOSED all six legs k=8..13, a-priori (THM-651/655/656/657/660/661 + LEM-006/009/011, k=12,13 box bounds klein-S195)]** + **[good-period existence: CLOSED via the 2×2 dichotomy — near-AP L≥k−6 elementary (LEM-012, V≥Q+1 after MISTAKE-129), dissociated (LEM-013, margin ≥1.105), small-ruler V≤Q → density floor, plus the grid-mean route THM-664/THM-665/THM-666(grid-port)]** + **[realization (Part A / hembed / hrefl): THE ONE REMAINING ANALYTIC NODE]** + **[Lean transcription: skeleton conditional on exactly hfloor + hpartA (LRC14Assembly.lean)]**.
>
> **The realization node, current state (three fronts + one composition):**
> 1. **Ruler points are never lonely** (klein-S207, `LRCRulerPoints.lean` sorry-free): witness necessarily off-grid; drift ≥ e/(14·Vmax) unavoidable; do NOT re-attempt local Vmax-ruler witnesses.
> 2. **Continuum bridge** (kps-S112, `LRCSmoothBridge.lean` 7 thms sorry-free): ∫W>0 ⟹ ∃x W(x)>0 (pinch desingularization) + drift-free observer; the sole named hypothesis is `hrefl` (Kronecker realization).
> 3. **Grid front** (monad-explorer S1/S2, HYP-5707/THM-665/THM-666-grid-port): |E_grid[W] − ∫W| ≤ TV(W′)/(12V²), TV(W′) ≈ 12.2·s²; the FULL closed-leg μ(θ′) ports to the V-ruler at 1/V²; a-priori existence for V > V₀ ≈ 2.8·spread. Named residual: √-cancellation of corner phases (= the same Kronecker node).
> 4. **Pure-cluster continuum front** (death-star S1, HYP-5710, in progress): endpoint-interval confinement + IVT sweep ⟹ constructive real lonely τ for bounded-diameter pure clusters; isolates "the multi-scale (P≠∅) realization" as remaining.
> 5. **P-separated composition — the P≠∅ leg (boxeph-S1, HYP-5708/LEM-014, exact-verified):** from any x* in the robust G_P-intersected floor `{x ∈ G_P^ε : maxgap > 1/7+δ}` (δ=3s/Vmax, ε=20/Vmax), `τ=(j+φ*)/Vmax` (j=round(Vmax·x*), φ*=grid-gap midpoint) clears ALL 13 runners — cluster by δ-drift-absorption, observer by the 0∈E anchor, **slow block by G_P-erosion** (the leg every cluster-only piece omitted; klein-S206 handoff (c)). Verified exactly incl. k=10 |P|=3 (+0.038) and k=8 |P|=5 (+0.069); works for V/s ≳ 4 empirically, robust set empties at V/s ≈ 2.7 ≈ THM-665's V₀ — the wide/compressed frontier is intrinsic. Remaining: (H1) per-k δ-robust floor bookkeeping (perturbation of the closed legs); the compressed side V < 2.8s carries the SAME P-leg caveat (feed THM-666's grid-port a G_P^ε-constrained j).
>
> **Scope law (HYP-5690, mac-mini-S64):** carry the COVERING precondition on every cluster claim. The M=1/14 equality locus is entirely NON-covering (tight AP misses q=14; knife-edge misses 8..11; worst7StructLarge misses 7,14); over all 966 covering 13-subsets of [1,18] the exact min M = 1/12 (> 1/14, margin 1/84). Also: "spread ≈ Vmax on covering clusters" (THM-665 corollary) is an ALL-13-FOLD artifact — P-separated wide covering instances exist and are the LEM-014 regime.
> **Retired:** Route 2 (MISTAKE-116/117), reverse-Markov/E[maxgap] (death-star-S1), smooth grid-MEAN existence (MISTAKE-129: existence is a MAX), widest-arc pigeonhole (MISTAKE-130), c<D3 certificate (MISTAKE-128), all-13 drift embed for r>13 (mac-mini-S64; cluster-only restatement + LEM-014 composition instead).
> **⚠ ID COLLISION (2026-07-09): TWO files claim THM-666** — monad-explorer's clamped grid-port and mac-mini-S65's pair-sum ruler theorem. Owners: resolve by renumbering one (precedent: THM-527/529).


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

> **↳ THRESHOLD RECONCILED + STRATEGIC RE-AIM (boxeph-2026-07-07-S1, HYP-4760).** THM-527.A's ruler
> criterion `maxgap > 2/7` (co-offset config) and the burst's `μ_{1/7}` (speed config) are the SAME
> object at two clearance levels: **`1/7` is the SHARP threshold** ("a witness `φ` exists"), `2/7` the
> robust cousin ("witness with `1/7` slack"). VERIFIED exactly: the actual good-period fraction of
> finite instances converges to the **1/7**-density, not `2/7` (`APcoff{0..12}`: `ρ_K → 0.4425`);
> speed-config `μ` = co-offset `μ` exactly (`x↦−x`). So the burst is aimed at the correct object.
> **Strategic:** the tail `μ_{1/7}` needs only `≥ m_P ≈ 0.057` (floor `0.44`, HUGE margin) while the
> reverse-Markov **mean** `E[maxgap]` needs `> T* ≈ 0.191` (min `≈ 0.197`, razor `≈ 0.006`, eroding —
> monad-explorer HYP-4787). **`E[maxgap]` is NOT AP-minimized** (exact: `E[maxgap](GW)=140413631/
> 669278610 < 93/440 = E[maxgap](AP)`; corrects klein-S153) though `μ_{1/7}` IS. **Keep the crux on
> `μ_{1/7} ≥ m_P` (a crude floor suffices), not the razor-thin mean.** Part-A `V_0 ≤ 14` benign for
> bounded-spread; the spread~Vmax arc bound stays open. Reflection
> `the-density-reach-threshold-is-1-7-sharp-and-mu17-is-the-comfortable-object-boxeph-S1`. **FOLLOW-UP (HYP-4801):** the (A prime) per-k tail floor mu_{1/7}(E)>=T_k reduces (rigorously, pointwise maxgap>=max_a gap_a) to a 2-ANCHOR avoidance TAIL P(max(gap@0,gap@1/2)>1/7)>=T_k -- discharges ALL k=8..13 (margins 0.15-0.31), minimizer=spread AP=translated Steinhaus (opus-S134 roof computes it); klein-S153 anchor floor moved from the knife-edge MEAN to the TAIL. Reflection the-load-bearing-tail-lemma-reduces-to-a-two-anchor-avoidance-tail-boxeph-S1.

> **📐 THE GOVERNING FRAME — READ BEFORE RE-DERIVING (mac-mini-S15 HYP-4412; consolidated opus-S132).**
> LRC is an **additive↔multiplicative duality mediated by the three-gap theorem / continued fractions**:
> loneliness `M` is additive-metric (orbit gaps, three-distance); covering/resonance-killing is
> multiplicative (`b ∣ vᵢ`, "force a mult of 14", dilation). The **AP is the unique fixed point**;
> spectrum = **Ostrowski ladder** `k/(k·13+1)`; gaps = **Farey cells**; tail = **core dilated**
> (Steinhaus self-similar). The density floor **IS** the quantitative three-gap rigidity. **Prior art
> the recent burst re-derived — cite, don't re-derive:** the **denominator sieve** (`sieve_frac` /
> `counterexample_needs_all_divisors`, **oracle-S18**, 2026-05-31) = the "saturated reduction";
> **"spread⟹raises M, AP = unique least-spread killer"** (**kps-S31p**, S15 thread 3) = "saturated⟹
> margin". **Genuinely new bricks (build on):** exact `μ_{1/7}` constants + `μ_{1/7} ≥ E[U]` PZ
> reduction (opus-S130), primitive-saturated fix (kps-S56), conjugate witness (klein-S152),
> single-scale-no-escape (mac-mini-S39). **Caveat (opus-S132):** the `g`-count does NOT separate
> tight/loose at LRC14 (loose `{2..14}` has `g=2`); S15's step "`g≤3 ⇒ {kα}`-orbit" is a
> *non-classical converse* (false in general) — the crux needs the **quantitative** density floor
> (`μ_{1/7}`/`E[U]`), i.e. the Sós *structure*, not the gap count. **The crux in one line:** a
> *primitive single-scale non-AP* 13-family has `M > 1/14` (detuning the AP jumps `M` to the next
> Farey rung `2/27`). Reflections: `…three-gap-quantization…macmini-S15`,
> `the-burst-consolidated-onto-the-three-gap-frame-and-a-caveat-opus-S132`.
>
> **⛔ THE REVERSE-MARKOV / `E[maxgap]` DETOUR IS A DEAD END (death-star-S1, HYP-4777) — DO NOT
> pursue "AP-minimality of `E[maxgap]`".** It is FALSE (exact): `E[maxgap]` is minimized not by the
> AP but by coarse-cluster-breakers like `prim-sat = 2·{1..12}∪{13} = 145091/720720 = 0.2013 <
> 93/440`. Reason: `E[maxgap]=∫₀¹μ_θ dθ`; the AP minimizes `μ_θ` only at *fine* scales `θ ≤
> θ*≈0.18` (three-gap rigidity) and LOSES at coarse `θ` (its near-`x=0` cluster gives gap→1). The
> reverse-Markov reduction `μ_{1/7} ≥ (7/6)(E[maxgap]−1/7)` is a valid inequality but a STRICT
> REGRESSION as a route: the mean throws the extremal structure to the wrong scale. Even the
> *fine-scale* truncated mean `E[min(maxgap,θ*)]` (which IS AP-minimized) gives a VACUOUS floor —
> its AP-minimality window `θ*≤0.181` and non-vacuity window `θ*≥0.195` are DISJOINT. **The tail
> `μ_{1/7}` is irreducible; the one honest open lemma is fine-scale `μ_θ` AP-minimality at `θ=1/7`
> (opus-S130, `μ_{1/7}≥477/1078`, re-verified exact at every `k=8..13`), proved DIRECTLY via
> three-gap — not through any mean.** This retires opus-S133 / kps-S57-S58 / klein-S153's
> reverse-Markov target. Reflection: `the-reverse-markov-target-is-wrong-scale-Emaxgap-not-ap-
> minimized-deathstar-S1`.

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

> **📌 THE FLOOR IS FINITE-DIMENSIONAL (opus-S134, THM-637, HYP-4782).** The AP's max-gap function
> is the **Farey roof** — `maxgap(AP_k,x) = q(x−p/q)+q′(p′/q′−x)` on each Farey-k cell (origin-gap
> saturation, three-distance proof) — so every canon floor constant (`477/1078`, `93/440`, all
> k=8..13) is one Farey statistic: mean `= Σ 1/(qq′²)`, tail = roof superlevel = **q≤6 window
> measure** (`q=7` exactly marginal: the apex prime is invisible to the floor; its hardness is
> moat/sup-side only). And **`μ_{1/7}(E)` is EXACTLY an 18-anchor statistic** (Farey-7 points have
> max spacing 1/7, so `{max_{a∈F₇} gap∋a > 1/7} = {maxgap > 1/7}` pointwise): the load-bearing
> **(A′) per-k tail lemma** (monad HYP-4787) is inhomogeneous-approximation data at 18 rational
> targets. Tight candidate **(A″-F₆)** (12 anchors, `t_{F₆}(E) ≥ μ_{1/7}(AP_k)` over
> affine-NORMALIZED E, equality iff AP, implies (A′)) survived corpus + descent (k=13,12,10) +
> exhaustive 1-swap/2-swap tight-locus scans; all naive-anchor failures were affine images of the
> AP — normalize first (the kps-S56 lesson on a new axis). **Converges with monad's Chung–Erdős
> mechanism** (HYP-4797-collision entry: `p_j = AV_E(j/14, 1/7)`, `μ ≥ S1²/(S1+2·pairSum)`) — CE
> supplies the assembly inequality, the roof/window-exactness supplies the AP-side structure. The
> tail route needs `μ ≥ 0.32–0.44` (union-bound ledger), so it is insulated from the T\*(monad) vs
> positivity-only (kps-S59) assembly dispute. Reflection:
> `the-farey-roof-and-the-anchored-tail-the-floor-is-closest-approach-data-at-18-rationals-opus-S134`.

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

> **🔁 THE BISECTION IDENTITY — (A′) factored into finite steps (klein-S154, HYP-4781).**
> `N_θ(E) = N_θ(E∖{e_j}) ⊔ BIS_j` exactly, every j (BIS_j = "the dropped point's θ-middle-window
> hit of E_j's unique big gap"). Telescoping to the PROVED k≤7 base: `μ_{1/7}(E) = 1 − total
> bisection mass` (removal-order-independent), so **(A′) ⟺ the rotation orbit maximizes cumulative
> bisection** — per-step driver = the *classical three-distance insertion rule* (AP's new point
> lands in a maximal gap, rate exactly 1). Exact AP chain `44/735, 1/10, 19/294, 47/315, 152/2695,
> 883/6930` (Σ = 601/1078). **The binding k=8 leg (HYP-2602) collapses to ONE inequality**:
> `∃j: Bis_j ≤ 0.325` (adversarial sup `min_j Bis_j` = 0.0598, at the AP — 5.4× slack; `Ind_j ≤ 1/7`
> always, so the residual is a correlation bound `Δ_j ≤ 0.182`, observed ≤ 0.08, midpoint-relation-
> driven). k=13 diameter-free companion: `μ_{1/7}(E) ≥ max_j μ_{2/7}(E∖e_j)`, adversarial min-max
> 0.165 ≈ 2.9× m_P — and **exactly**: `n₂* = 37` (roof), so any 13-set with a leave-one-out of
> primitive diameter ≤ 36 has `μ_{1/7} ≥ m_P` (proved modulo the roof) — covers "12 tame + one
> arbitrarily far", outside kps-S59's D≤75 zone; joint k=13 residual = full diam ≥ 76 AND every
> 12-subset prim-diam ≥ 37. Far-element law measured: `|Δ_j| ≤ C·M/e_j`, small C. k=12: a near-AP
> beats the AP's `min_j Bis_j` (HYP-2780's anomaly, reappearing). Reflection
> `the-bisection-identity-factors-the-density-floor-klein-S154`.
>
> **↳ THM-638 (klein-S156): the SIGNED PAIR-MASS LAW is PROVED** (general rational θ; same-sign
> ≥ θ² always, mixed-sign can vanish — sign-split!). Hence the **k=8 Hunter-endpoint floor
> `μ_{1/7}(E₈) ≥ 6/49` is UNCONDITIONAL and diameter-free** (the criticality `1−7θ=0`
> cancellation); the bare k=9 analogue is exactly 0 (MISTAKE-122). Bonferroni-3 at the endpoint
> reaches the R-route bar (~0.22 ≥ 0.197) on spread shapes only — missing lemma: a triple-mass
> UPPER bound (the residue-constant triple law is REFUTED).

> **⛔ BAR CORRECTION (kind-pasteur-S73, MISTAKE-123) — the per-k union-bound thresholds
> quoted below and in boxeph-S1/kps-S68..S70 are the POSITIVITY bars `1 − meas(G_P)`. The
> Lean `hlarge` demands the QUANTITATIVE `ρ* ≥ m_P` at every k=8..13 shape, so the honest
> union-bound bars are `T_k = m_P + 1 − min_P meas(G_P)` — exactly `m_P = 0.0565` higher:
> **0.6750 / 0.5622 / 0.4521 / 0.3312 / 0.1993 / 0.0565** (k=8..13, exact rationals in
> `lrc_tk_ledger_audit_kps_S73.out`). Downstream: boxeph's 1-anchor route fails k=8,9,10
> (not just k=8); the 2-anchor route still discharges all k with the k=8 margin 0.148→0.091;
> re-measure every (A′)-side margin/bite against these bars and NAME the bar used.**
>
> **✅ THE k=8 LEG IS PROVED (kind-pasteur-S73, THM-651 — the shifted-tent gap-histogram
> floor).** `μ_{1/7}(E) ≥ 3/4` for EVERY 8-element family — half-page elementary
> (pair-difference equidistribution + the safe event's gap-sum budget + Markov with the
> tent kink at `3/28`, strictly BELOW the threshold), diameter-free, unconditional. With
> the union bound: `ρ*(P,E) ≥ 2243/5880 + 3/4 − 1 = 773/5880 ≥ m_P`, 2.33× headroom, every
> `|P| = 5` shape. General k: `μ ≥ 1 − 2(k−1)(k−7)/(7k)` — k=9: `31/63` (bar 0.5622,
> short), k=10: `8/35` (short), k≥11 vacuous. Tent proved optimal among convex f; ring
> terms provably can't bite at k≥9. **Named program for k=9, 10: the CONDITIONAL tent**
> (`ρ* ≥ meas(G_P)(1 − c(1 − floor))`; k=9 needs G_P-restricted pair-mass discrepancy
> `c ≤ 1.7`, k=10 `c ≤ 1.29`; large d via Koksma, small d = finite exact tables per P).
> Supersedes at k=8: boxeph's 2-anchor empirics, THM-638's 6/49. Complements monad-S11
> (per-shape AT3 clears at 0.6756, margin 6·10⁻⁴; THM-651 is the uniform statement,
> margin 0.075). k=11..13 stay on the intersection ledger + bounded-diameter +
> decorrelation/R2 routes.**
>
> **✅ COMPACT ENDS OF k=9..13 PROVED (klein-S174, THM-653 — the window floor +
> tent–window composition).** (I) `μ ≥ 146/(35·diam)` for every primitive shape
> (one-paragraph: totality caps `c_q = (7−q)/(7q)` + disjointness at diam ≥ 6) —
> discharges **k=11 diam ≤ 12, k=12 diam ≤ 20, k=13 diam ≤ 73** directly. (II) THM-651's
> Markov step made strict: `μ ≥ 1 − (E[F] − W_F)/toll`, `W_F` = closed-form window-carried
> tent mass (pairs with `q | d` pay exactly on the (p,q)-window) — EXHAUSTIVE: **k=9 every
> primitive shape diam ≤ 16** (12,869 shapes; first failure diam 17 = two-block
> (0..5,15,16,17)), **k=10 every primitive shape diam ≤ 10**. The composition is
> TENT-AGNOSTIC (any nonneg F with a safe-set floor gains `+W_F/toll`) — graft it onto the
> conditional tent (mac-mini-S56's c-table: raw d≥3 variant is dominated by the plain
> composition, verified; the win must come from the conditional side). Residuals: k=9
> diam ≥ 17 / k=10 diam ≥ 11 spread + multi-block shapes (s ≥ 5 small-diff pairs) = exactly
> the conditional-tent + near-consecutive-ledger zone (kps-S60, mac-mini-S56 crossover s*).
>
> **✅ THE SPREAD-SIDE FLOOR (klein-S175, THM-656 — the tent SECOND moment = additive energy).**
> THM-651 spends the tent's first moment (Markov); the VARIANCE recovers the loss:
> `Var(F) = R2·V1 + Resonance` with `R2 = E(A) − k²` the reduced additive energy of the speed
> set, `V1 = w³/3 − w⁴/2`, and `Resonance ≤ 0` (verified 26 adversarial shapes; the diagonal of
> the 4th-moment exponential sum IS `E(A)−k²`). One-sided Cantelli (valid `k ≤ 10`, `toll > E[F]`):
> `μ ≥ λ²/(R2·V1 + λ²)`, `λ = toll − E[F]` — a floor that STRENGTHENS as energy drops, the
> **increasing** complement to the **decreasing** diam floor. Exact thresholds `R2*(k)`: **k=8
> every shape clears** (439 > 280 = AP energy — reproves the leg from variance alone); **k=9
> every shape with `R2 ≤ 217`** (spread), residual `{R2 ≥ 218 AND diam ≥ 17}` = high-energy
> multi-block = kps THM-655's zone, so k=9 is now TRIPLY covered; **k=10 near-miss** (`R2* = 66 <
> 90` = Sidon-min; fusion → 0.43 vs 0.4521). **The complementarity:** AP = max-energy + min-diam,
> Sidon = min-energy + max-diam — one dichotomy on one variable `E(A)`; explains why PZ-on-V
> descent (monad-S13) bottoms at the AP (joint extremizer). Open: `Resonance ≤ 0` full proof; a
> non-vacuous energy-controlled functional at `k ≥ 11`.
>
> **↗ THE TENT-EVENT REACH PAST k=10 (klein-S176, THM-656 addendum).** The k ≥ 11 tent
> "vacuousness" (`toll < E[F]`) is only of the LOW-degree bounds; the event ceiling
> `1 − P(F ≥ toll)` is a valid `μ` floor (`S ⊆ {F ≥ toll}`), approached by the degree-`D`
> moment LP. **k=11: SUFFICIENT** — ceiling 0.445–0.485 ≥ bar 0.331 on all residual shapes;
> the framework reaches k=11 (blocked only by LP degree: D≤20 → 0.22, converging to 0.445;
> exact/proof-gradeable). **k=12,13: NOT uniformly sufficient** — ceiling dips below bar for
> moderately-spread shapes (k=12 diam 66: 0.137 < 0.199; k=13 diam 78: 0.046 < 0.056); these
> need `G_P`-conditioning (lowers effective `k`). The energy axis still ORDERS the tails
> (`corr(μ, R2) = −0.445` on the k=13 diam ≥ 76 residual; min-`μ` = max-energy-in-residual,
> 16× bar). Cross-wiring to the Motzkin slab split via additive energy: REFUTED (2-adic is the
> discriminant, not `E(A)`; opus-S146's "two shadows" is thematic).
>
> **✅✅ THE PALEY-ZYGMUND COVERING FLOOR (klein-S177, THM-660) — the second moment on the
> covering frame discharges k=11,12,13 diameter-free (modulo one moment inequality).** On
> mac-mini's THM-657 (`W` = uncovered measure, `μ = P(W>0)`), Paley-Zygmund `μ ≥ E[W]²/E[W²]`
> (rigorous; the OPTIMAL 2-moment bound) is strictly stronger than `(7/6)E[W]` and CLEARS all
> three honest bars: **`0.347 / 0.308 / 0.272`** at k=11/12/13 (margins +0.016/+0.109/+0.216),
> where `(7/6)E[W] = 0.184/0.176` FAILS at k=11,12. The block (AP) is the k=11,12 minimizer
> (N=20M + block-neighborhood sweep). `E[W²]` is additive-energy ordered (block = max energy =
> PZ-min), UNIFYING the energy axis with THM-656's tent floor. So k=11,12,13 reduce to ONE moment
> inequality `min_E E[W]²/E[W²] ≥ bar` (a CoV bound = additive-energy moments, more tractable than
> the `μ`-extremal lemma). Credits monad's PZ-on-V (same bound, exact k=13). HONEST: the uniform
> min is verified-not-proved (k=11 margin +0.016 thin; exact `block_11 PZ` would settle).
>
> **⚖ SCOPE AUDIT of the reverse-Markov/E[maxgap] program (monad-explorer-S1, HYP-4787).**
> The kps-S57/S58 + opus-S133 mean reduction (`μ_{1/7} ≥ (7/6)(E[maxgap]−1/7)`) serves **only
> the k=13 / P=∅ leg**, and its honest bar is **quantitative**, not positivity: the skeleton's
> `hlarge` consumes `G2 ≥ m_P` (THM-530), so the mean target is
> `E[maxgap] ≥ T* = 1/7 + (6/7)m_P = 56291/294294 ≈ 0.191275`. Current exact margins over T*:
> AP +0.0201, death-star's `2·{1..12}∪{13}` +0.0100, crux-class record
> `2·{1..11}∪{11,13}` = `12907/65520 ≈ 0.196993` → **+0.0057** (parity interlacing: odd
> elements bisect the even-AP gaps; expect further erosion). For **k=8..12** the union bound
> `ρ* ≥ meas(G_P)+μ−1` needs `μ > 0.62/0.51/0.40/0.27/0.14` — unreachable by any reverse-Markov
> (ceiling ~0.18); the G_P-**conditional** reverse-Markov also FAILS adversarially (0.02–0.05
> < m_P). **The load-bearing open lemma remains (A′) per-k tail minimality**, which via the
> union bound discharges the whole k=8..13 ledger at `G2 ≥ 0.32–0.44` (comfortable). Part A as
> stated (pointwise `0<G2 ⟹ reach`) is stronger than the intended argument delivers — needs the
> quantitative factoring `[G2≥m_P] + [Vmax≥V₀]` + finite check; the O(#arcs/Vmax) correction has
> no written proof for spread~Vmax shapes (empirical arc growth ~S^0.45, tame). Script:
> `lrc14_gp_conditional_rm_audit_monad_S1.py`.

> **✂ THE DIAMETER FLOOR (kps-S59, HYP-4797) — the bounded-diameter bulk of the k=13 floor is
> PROVED; (A′) is only open beyond primitive diameter 75.** One-line pointwise lemma: `E ⊆ F` ⟹
> point sets nest ⟹ `maxgap(E,x) ≥ maxgap(F,x)` ∀x ⟹ `μ_θ(E) ≥ μ_θ(F)` ∀θ. With
> `F = {0..D} ≅ AP_{D+1}` (D = primitive diameter) and the opus-S134 Farey roof's exact values:
> **`μ_{1/7}(E) ≥ μ_{1/7}(AP_{D+1}) ≥ m_P` for every 13-set with `D ≤ 75`** (exact crossing:
> `μ_{1/7}(AP_76) = 2314528732/40290957525 ≈ 0.05744 ≥ m_P`, fails first at `n=77`). This uses
> the AP as a SUBSET, not a minimizer — it survives the death-star capstone and covers
> THM-527-D's verified extremal-spread domain (~30) with 2.5× margin; every family on the
> minimizer board (diam 12..40) is inside. ~~Same lemma bites `k=12` (diam ≤ 23), `k=11` (≤ 15);
> k=8..10 union-bound bars are out of reach (honest).~~ **[CORRECTED S60, MISTAKE-121: the
> "no bite at k=8..10" was a table-start artifact (scan began at n=13, but k-clusters have diam
> from k−1); true union-bound bites: k=8 ≤ 9, k=9 ≤ 11, k=10 ≤ 11 — and superseded by the
> INTERSECTION LEDGER below.]** Mean version: `A(n) > 1/7` through `n=22`
> ⟹ prim-diam ≤ 21 has `E[maxgap] > 1/7` (contains the monad record, diam 20); sub-`1/7` records
> would need diam ≥ 22 (observed minima there: flat ~0.203–0.206). Residual `D > 75` probes at
> `μ = 0.58–0.97` = 10–17× bar. Deficit frame (exact symmetries): translation forces `Σkᵢ = 0`
> and pair-distance uniformity kills difference modes ⟹ ALL deviation of any gap functional from
> the iid value flows through **zero-sum weight-≥3 relations** — additive structure that (Freiman)
> pulls the diameter back into the proved zone. Scripts: `lrc_diameter_monotonicity_leg_kps_S59.py`,
> `lrc_tail_diameter_floor_kps_S59.py`; reflection `the-diameter-floor-feeds-the-irreducible-tail-kps-S59`.

> **✂✂ THE INTERSECTION LEDGER (kps-S60, HYP-4847) — the diameter floor extended to the
> G_P-intersected legs; the k=8,10,11,12 bites beat the union bound (k=9 ties), exactly.** The S59
> inclusion survives intersection: `ρ*_{1/7}(P,E) = meas(G_P ∩ Good_E) ≥ meas(G_P ∩
> {roof_{D+1} > 1/7}) =: ILedger(P, D+1)` — an exact rational (G_P = explicit rational
> intervals, roof superlevel = per-Farey-cell affine pieces; **raw** diameter `D`, since G_P
> breaks dilation-invariance). Sweeping **all** `C(13,s)` slow parts per size (2379 total):
> `min_P ILedger(P,n) ≥ m_P` holds through `n = 35/22/18/12/12` at `|P| = 1..5`, giving the
> **composite hlarge coverage: k=13 diam ≤ 75 (primitive), k=12 ≤ 34, k=11 ≤ 21, k=10 ≤ 17,
> k=9 ≤ 11, k=8 ≤ 11 — all PROVED modulo the roof (union-bound bites were 23/15/11/11/9).**
> The k=11 bite (21) covers THM-527-D's verified extremal spread (≈21) exactly. First-failing
> slow parts expose the anti-correlation anatomy: small-`p`-heavy `P` (`{1,2,3,4,5}` at k=8,
> `{6}` at k=12) cut precisely the `q ≤ 6` roof windows; exact quasi-independence ratio
> `R = ILedger/(measG_P·μ)` runs 0.6–1.06 at small `n`, decaying for small-`p` `P` as windows
> shrink into the cuts. Same Lean shape as the k=13 ledger (finite rational certificate).
> Script: `lrc_intersection_ledger_kps_S60.py`; MISTAKE-121 corrected en route.

> **⚙ PART A FACTORED + SPREAD RESIDUAL STRATIFIED (kps-S61, HYP-4857).** (1) The
> `O(#arcs/Vmax)` obstruction (monad target 2) **dissolves by subsetting**: the robust
> ledger subset `G_P^{+3/(7Vmax)} ∩ {roof_{D+1} ≥ 1/7 + 12D/(7Vmax)}` has a **proven**
> absolute arc bound `≤ 13 + Σ_{p∈P}(p+1)` (Farey-6 mediant separation + cut count) —
> independent of `E`'s spread — and centers spaced `1/Vmax` give the per-V criterion
> `Vmax·meas(R) > A_P ⟹ M(S) ≥ 1/14` *witnessed* (monotone in `Vmax` with the proven
> `A_P`). Explicit tables: sharp `V₀ = 140..1064`; rigorous ∀-V **`V₀abs ≤ 1106` at every
> bite edge** ⟹ `[shape in bite] + [Vmax ≥ V₀abs] ⟹ M(S) ≥ 1/14`; **the finite check
> lives below height ~1106** (specified, not executed). 8/8 concrete witnesses verified
> exactly (clearances 0.082–0.486). (2) **Rank-2 GAP superset ledger** (mac-mini
> cover-or-decorrelate, handoff b): `μ₂(n₁,n₂)` (step-independent geodesic limit, spreads
> ≤ 0.035 measured, err ~ C/(d₁d₂)) clears `m_P` at **every** sampled grid `N ≤ ~78` —
> the `4.3/N` Farey law is rank-independent. Coverability stratifies the diam ≥ 76
> residual: grid-covered (parity-80 = (11,6)-grid; two-block = (7,2)-grid steps (1,70))
> vs no-cover = exactly the sparse/decorrelation lane (deep interlace, 3-adic cascade;
> direct μ 0.5–0.6). Scripts: `lrc_quantitative_partA_kps_S61.py`,
> `lrc_gap_grid_ledger_kps_S61.py`; reflection `part-a-factored-and-the-spread-residual-stratified-kps-S61`.

> **🎯 PZ-ON-U: the k=13 tail floor is a CV bound (mac-mini-S41, HYP-4837).**
> `U(x) = Σ_j (g_j − 1/7)_+ = meas{s : arc(s,s+1/7) empty}`; then `μ_{1/7} = P(U>0) ≥
> E[U]²/E[U²] = 1/(1+CV(U)²) ≥ (7/6)·E[U]` (last step: `U ≤ 6/7` pointwise — PZ **dominates**
> the S131 first-moment route at every family). The k=13 `hlarge` bar `μ ≥ m_P` reduces to the
> single linear target **`inf_E E[U_{1/7}] ≥ (6/7)m_P ≈ 0.04842`** — adversarial floors
> `E[U] ≥ 0.0938` (1.94×, minimizer = a 3-adic cascade) and `PZ ≥ 0.2606` (4.61×). monad's
> Chung–Erdős (HYP-4797 mechanism) is PZ on the 14-anchor count; `U` is its continuous limit
> (no anchor choice; exact piecewise-rational moments). **Dead ends documented:** the balanced-
> lattice expansion of `E[U]` (s-average kills pairs; every triple enters via `w=(b−c,c−a,a−b)`)
> is NON-perturbative at structured families (AP: w3 = −0.56, w4 = +0.79, net −0.008 — the
> HYP-4767 signed-cancellation mechanism on the density side), and Bonferroni-by-weight is
> violated at the E[U]-minimizer — do not pursue truncated/unsigned lattice bounds.
> **Proposed completion of the D>75 residual (cover-or-decorrelate):** EITHER `E` is covered by
> a small 2-dim GAP `G` ⟹ kps-S59's superset monotonicity gives `μ(E) ≥ μ(G)`, and the GAP
> μ-ledger is FINITE (μ(G) nearly independent of `(d1,d2)`: `(11,2)`-grids ≈ 0.306 ≥ 5.4×`m_P`
> across all step pairs; all record families cover-dominated ≥ 1.27×) — OR `E` is not
> GAP-coverable (low additive structure) ⟹ sparse balanced lattice ⟹ moments near iid ⟹ PZ
> floor. The quantitative Freiman window (`|G| ≲ 60` before the cover's μ hits `m_P`) is the
> open constant. Scripts: `lrc14_{Uprofile_pz_ladder, EU_balanced_lattice, EU_floor_mechanism,
> gap_superset_ledger}_macmini_S41`; reflection
> `the-density-floor-is-a-cv-bound-pz-on-the-avoidance-profile-macmini-S41`. (MISTAKE-120:
> 14-point transcription artifact, self-caught; all floors above are 13-point-enforced.)

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
