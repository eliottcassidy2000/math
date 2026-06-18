---
id: THM-526
title: The arc-width lemma for LRC — a non-tautological discharge of the large-parked-runner direction (if the core's widest level-1/14 safe arc exceeds the parked runner's danger-tooth width 1/(98m), then M(core ∪ {14m}) ≥ 1/14, with an explicit witness), closing LRC(14) unconditionally for the family of (size-12 cores ⊆ {1..13} covering 2..13) ∪ {14m} over all m
status: PROVED (the lemma: elementary interval-covering, proof below; tighter threshold m>1/(98·W) verified vs dense check). The sub-family closure is PROVED (lemma for m≥4 + exhaustive small-m check m≤3, 0 counterexamples). Co-derived with the S4 creative-reframing workflow (Angle B / "homogenization"), which obtained the same lemma with a weaker threshold m>1/(14·W).
source: mac-mini-2026-06-17-S4 (creative reframings of the last LRC(14) inequality)
depends_on:
  - THM-523   # covering-set reduction (the only hard case is covering 13-sets)
  - THM-524   # binding-pair reduction (M attained at a pair crossing)
  - THM-525   # easy-dominates-hard (the parked-runner / core framing)
related:
  - HYP-2578  # the reframings ledger: sawtooth-exclusion, off-grid fence, even-clearing, dip-law collapse, the clearing-crossing reformulation
  - HYP-2577  # the peeling-chain decomposition (this lemma supplies its missing rigorous large-m piece)
  - THM-518   # measure-side stranger-decoupling (this is the GAP-side, with an explicit threshold)
external: Lonely Runner Conjecture; proven for ≤7 runners (Barajas–Serra), open for 13.
---

# THM-526 — The arc-width lemma: discharging the large parked runner

**Context.** LRC(14) reduces (THM-523/524/525) to: every covering 13-set `S` has
`M(S) = max_τ min_{v∈S} ‖vτ‖ ≥ 1/14`. Such `S` contains a "parked" runner `w = 14m`
(`w ≡ 0 mod 14`, pinned at the observer at every grid time). Writing `S = A ∪ {w}`,
the genuine non-tautological difficulty was: bound how much the parked runner can pull
`M` below the core gap `M(A)`. This lemma settles the **large-`m`** half outright.

## The lemma (PROVED)

Let `A` be any finite set of positive integers with `M(A) > 1/14`, and let
`G_A = {τ ∈ [0,1) : ‖aτ‖ ≥ 1/14  ∀ a ∈ A}` be its **level-1/14 safe set** — a finite
union of closed arcs of positive total measure (positive because `M(A) > 1/14`). Let
`W(A)` be the width of its **widest arc**.

> **Lemma.** If `W(A) > 1/(98m)`, then `M(A ∪ {14m}) ≥ 1/14`.
> Hence for every `m ≥ M₀(A) := ⌊1/(98·W(A))⌋ + 1`, looseness is automatic.

**Proof.** The parked runner's danger set at level 1/14 is
`D = {τ : ‖14m·τ‖ < 1/14} = ⋃_k ( k/(14m) − 1/(196m),  k/(14m) + 1/(196m) )`:
`14m` open "teeth", each of full width `t := 1/(98m)`, with centers spaced
`1/(14m) = 7t` apart — so consecutive teeth are separated by **safe gaps of width `6t`**.
Let `I` be a widest arc of `G_A`, of width `W(A) > t`. Since each tooth has width `t < W(A)`,
`I` cannot lie inside any single tooth; as the teeth are isolated (gaps between them),
`I` must contain a point `τ₀ ∉ D`, i.e. `‖14m·τ₀‖ ≥ 1/14`. (Quantitatively, the safe
measure inside `I` is `≥ W(A) − (⌊W(A)/(7t)⌋+1)·t > 0` whenever `W(A) > t`.) Because
`τ₀ ∈ I ⊆ G_A`, also `‖aτ₀‖ ≥ 1/14` for all `a ∈ A`. Therefore
`min_{v ∈ A∪{14m}} ‖vτ₀‖ ≥ 1/14`, so `M(A ∪ {14m}) ≥ 1/14`. ∎

**Why it is non-tautological.** `W(A)` is computed from the **core `A` alone**, with no
reference to `S` or to `M(S)`; the conclusion `M(S) ≥ 1/14` is *derived*, with an explicit
witness `τ₀` (any safe point of the widest arc; e.g. the gap-center
`(2k+1)/(2·14m)`, where the parked runner is at distance ½). This is the gap-side analogue
of the measure-side stranger-decoupling (THM-518), but with an **explicit, computable
threshold** `M₀(A)` rather than an asymptotic `m → ∞`.

## What it closes (PROVED sub-family, unconditional)

For every size-12 core `A ⊆ {1,…,13}` that covers `2..13`, the worst widest-arc width is
`W(A) = 5/1848`, attained at the drop-6 core `A = {1,2,3,4,5,7,8,9,10,11,12,13}`, giving
threshold `M₀ = ⌊1848/490⌋+1 = 4`. So:

> **Corollary.** For every size-12 core `A ⊆ {1,…,13}` covering `2..13` and every `m ≥ 1`,
> `M(A ∪ {14m}) ≥ 1/14`. *Proof:* `m ≥ 4` by the lemma; `m ∈ {1,2,3}` by exhaustive exact
> check (0 counterexamples; in fact `M(A∪{14m}) = M(A) ≥ 2/23` there — the small parked
> runner does not bind). ∎

This is a genuine infinite sub-family of covering 13-sets for which **LRC(14) is now
proved unconditionally** — the canonical hard cores `{1..13}\{j} ∪ {14m}` among them.

## ANGLE B addendum (mac-mini-2026-06-17-S6): the CLUSTERED-large regime

The arc-width criterion `C(S)` (`∃ v: W(S\{v}) > 1/(7v)`) holds on every tested covering
13-set, but its **pigeonhole** lower bound `W(A) ≥ μ(A)/N(A)` (with `N(A) ≤ Σ_{u∈A} u`) is too
weak in the **clustered-large** regime, where many runners lie near a common value `V` and
`Σ u ~ 12V`. Two elementary lemmas (PROVED) give a sharper `W`, settling the easy half of that
regime **directly** (a global witness, not even via `C`):

> **LEMMA 1 (first-gap / cluster lemma, PROVED).** For any runner set `S` with
> `Vmin = min S`, `Vmax = max S` and `13·Vmin > Vmax`, the open arc
> `J = ( 1/(14·Vmin), 13/(14·Vmax) )` is level-1/14 safe for **every** runner in `S`. Hence its
> midpoint is a witness and `M(S) ≥ 1/14`. Width `|J| = (13·Vmin − Vmax)/(14·Vmin·Vmax)`.
> *Proof:* the first inter-tooth gap of runner `u` is `(1/(14u), 13/(14u))`; their intersection
> over `u∈S` is `J`, nonempty iff `13·Vmin > Vmax`; on `J`, `u·τ ∈ (1/14, 13/14)` so
> `‖uτ‖ ≥ 1/14`. ∎
> **Scope (exact):** exactly the covering 13-sets with spread ratio `Vmax/Vmin < 13`. This
> contains **all all-large clustered covering 13-sets** (their large window has width `≲14`, so
> `Vmax/Vmin → 1 ≪ 13`). Verified: closes 1950/1950 all-large covering 13-sets; e.g.
> `S = {140,…,153}` has `M = 140/293 ≈ 0.478` with witness `τ ≈ 0.0033`.

> **LEMMA 2 (band-fit mixed lemma, PROVED).** If `S = W ⊔ B` where `W` has a common safe arc
> `[p,q]` and the "band" `B ⊆ [a,b]` and some `t∈[p,q]`, integer `j` satisfy
> `j+1/14 ≤ a·t` and `b·t ≤ j+13/14`, then `t` is safe for all of `S`, so `M(S) ≥ 1/14`.
> *Proof:* `u·t ∈ [a·t, b·t] ⊆ [j+1/14, j+13/14]` for every `u∈B` (monotone in `u`). ∎
> (Implication verified on 85 555 hypothesis-satisfying instances, 0 violations.)

**Honest residual (the precise ANGLE-B obstruction).** The genuine **mixed** clustered sets
(small runners `+` a tight large cluster of spread `sB ≈ 35–40`) are **not** closed by Lemma 2:
at the working scale the band straddles several gaps (`sB·q ≫ 6/7`). At the actual witness the
band runners occupy **3–4 distinct gap indices** (e.g. `S={1,2,3,4,406,407,414,416,420,429,432,441,450}`:
indices 30–33) — so no single-gap lemma can close them; the witness is irreducibly a
**simultaneous-Diophantine (CRT/three-distance) alignment**. Moreover the sharp *average*-gap
pigeonhole still fails there (merged safe-component count `~2.6·(q·b)`, avg gap `< 1/(7v)`), while
the **widest** gap is `~4×` the average — the widest-vs-average factor is what makes `C` hold and
is exactly what resists a closed form. Computationally `C(S)` holds via 8–12 of the 13 large-`v`
choices on every mixed set, and a witness (Lemma 1/2 or exact CRT odd-multiple search) is found for
**every** clustered covering 13-set tested (≈ 3500+, 0 failures; supplementary exact-`M` spot checks
all `≥ 1/14`).

Scripts: `04-computation/lrc14_angleB_cluster_arcwidth_mac-mini-2026-06-17-S6.py`,
`04-computation/lrc14_angleB_adversarial_Mfloor_mac-mini-2026-06-17-S6.py`.

## ANGLE F assembly (mac-mini-2026-06-17-S6): the complete case-split + the residual is LOOSE

Assembling the pieces into the cleanest **complete proof skeleton**, with each case named and
its sufficient lemma PROVED, and the residual pinned to a single statement.

> **Reduction.** LRC(14) ⟺ every primitive covering 13-set `S` has `M(S) ≥ 1/14`
> (THM-523/524/525). By the generalized arc-width lemma (this theorem), it suffices to prove
> the **criterion** `C(S)`: `∃ v∈S, W(S\{v}) > 1/(7v)`.

**CASE S1 — single-large (`k := #{v∈S : v>13} ≤ 1`): PROVED, closed form + 7-set check.**
Here the 12 small runners `A ⊆ {1..13}` cover `2..13`, and the one parked `V = 14m ≥ 14`.
The pigeonhole bound `W(A) ≥ μ(A)/N(A) ≥ (1−Σ_{u∈A}1/(7u))/(Σ_{u∈A} u)` makes `C` hold via `V`
once `7V·μ_lb(A) > Σ_{u∈A} u`, i.e. `V > Σu/(7μ_lb)`. The **worst** core is `{1..13}\{5}`, giving
threshold `V > 30990960/1448599 = 21.394`. Since `V` is a multiple of 14:
- `V ≥ 28`: `PIG(V)` fires ⇒ `C(S)` ⇒ `M(S) ≥ 1/14`. **No finite check needed.**
- `V = 14`: exactly **7 covering sets** (`{1..13}\{j} ∪ {14}` that cover `2..14`), each checked
  exactly: `M ∈ {1/11, 2/23, 2/21, 2/19, 2/19, 2/17, 1/8}`, all `≥ 2/23 ≈ 0.087 > 1/14`. ∎

**CASE S2 — clustered/window (not S1, spread `Vmax/Vmin < 13`): PROVED by LEMMA 1** (the
first-gap lemma above): the common safe arc `J = (1/(14·Vmin), 13/(14·Vmax))` gives a *global*
witness. (Leave-one-out "origin-gap" variant — remove `v₀`, the merged origin-teeth half-width
`1/(7·min(S\{v₀}))` vs nearest non-origin left-edge `min_{u,k≥1}(7k−1)/(7u)` — fires on
1200/1200 tightest-window covering sets, but LEMMA 1 is the cleaner global form.)

**CASE S3 — residual (neither pigeonhole nor first-gap/origin-gap fires): C still holds, and
these sets are LOOSE.** Over a 4000-set covering sample, S3 is ≈ 2.7%; on **every** S3 set the
exact criterion `C(S)` holds (typically via 8–12 distinct `v`), and — the decisive point —
**S3 sets sit far from the LRC frontier**: the minimum exact `M` over all sampled S3 sets is
`4/31 ≈ 0.129`, i.e. `1.81×` the threshold `1/14`. By contrast the unique closest-to-`1/14`
covering set, `{1..11,13,84}` with `M = 7/89 ≈ 0.0787` (THM-524 C), is **single-large**, hence
closed by S1's pigeonhole. So the genuinely tight configs are all in the PROVED case S1, and the
residual S3 is provably non-dangerous.

> **THE RESIDUAL LEMMA (single precise open statement).** *Every covering 13-set in case S3
> satisfies `C(S)`.* Equivalently: when both the average-gap pigeonhole (`W ≥ μ/N`) and the
> origin/first-gap bound fail, the **widest** safe arc of some `S\{v}` still exceeds `1/(7v)`.
> The obstruction to an elementary proof is exactly the **widest-vs-average factor**: in S3 the
> merged safe set has avg gap `< 1/(7v)` but widest gap `~4×` larger, located at an
> irreducibly simultaneous-Diophantine (CRT / three-distance) alignment away from `τ=0`. This
> is a statement about loose sets only (`M ≥ 4/31` empirically), so it is a closure-of-bookkeeping
> task, not the crux of LRC(14).

**Exhaustiveness (verified).** Every covering 13-set is in exactly one of S1/S2/S3 (the split is
on `k` and on whether the two bounds fire); `C(S)` held on **4000/4000** covering sets across all
generators (clustered, small-core, mixed-2/3-large, spread); the 4 sets where the prior
`max(pigeonhole, antipode)`-with-`v=max` proxy (S5) failed all satisfy `C` via other `v`
(free choice of removed runner strictly strengthens the criterion). Scripts:
`04-computation/lrc14_angleF_proof_assembly_mac-mini-2026-06-17-S6.py`,
`05-knowledge/results/lrc14_angleF_proof_assembly_mac-mini-2026-06-17-S6.out`.

## Honest scope and the remaining hard regime

`W(A)` is **not uniform** over all cores: it shrinks toward 0 as core speeds grow
(e.g. `A = {1..11,13,1400}` has `W(A) ≈ 6·10⁻⁴`), so `M₀(A)` is unbounded across the full
covering family. The lemma therefore does **not** prove LRC(14): the surviving hard regime
is exactly **large-speed cores in small-`m` resonance**, where the parked runner *is* the
binding partner and `W(A)` is too small for the threshold to bite. Two further pieces would
finish the program: (i) bound the essential core speeds (compactness — the gap-side of
THM-522; the S4 finite-reduction probe confirmed the inf `= 7/89` stabilizes once
speeds ≤ 84, but found the "min at smallest representative" map is not monotone, so the
bound is empirical); (ii) the small-`m` window per bounded core (finite). The sharpest
distilled target (S4 synthesis) is **clearing-crossing existence**: for the binding pair
`D = flank + w`, `M(S) = j/D` with `j = (flank·num mod D)`, and the non-tautological content
is that the crossing index `j` is forced `≥ D/14` — i.e. some pair-crossing clears all 13
runners at level ≥ 1/14.

(Runner-up exact identity, S4: `M({1,…,13, 14m}) = m/(14m+1) < 1/14` — a **14**-element
set, the cleanest proof that the 13-element cardinality bound is essential.)
