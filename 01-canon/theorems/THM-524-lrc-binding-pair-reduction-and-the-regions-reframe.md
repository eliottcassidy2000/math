---
id: THM-524
title: The binding-pair (sawtooth-envelope) reduction of the LRC gap, and the regions/sections reframe — M(S)=max_τ min_v‖vτ‖ is attained at a binding-PAIR crossing τ=k/(v_a±v_b) or an all-odd peak τ=½, so LRC(N) reduces from an n-runner condition to ~O(n²) pairwise SWITCHES (find one pair-crossing with gap≥1/N and others clear); the regions/SDR view is the on-grid special case, the complement=reversal Z/2 is the one exact tournament bridge, and the covering hard core condenses to the closed form M=7m/(84m+5)
status: PROVED (the binding-pair reduction: elementary sawtooth lower-envelope lemma, proof below; the candidate set is exactly the M-tool's, stress-validated 0/300 grid-beats-candidate + dense-grid validation). VERIFIED (the reduction identity: 0 counterexamples over ~11,000 random configs n=2..13 + named cases). The covering closed form M=7m/(84m+5) and the tournament ledger are VERIFIED/honest-analogy as marked.
source: mac-mini-2026-06-17-S2 (user's "regions of the loop, not runners" reframe; the binding-pair switch view + tournament analogy)
depends_on:
  - THM-523   # the q-witness / covering-set reduction (gap side); this supplies the off-grid engine the section lens lacks
related:
  - THM-522   # kind-pasteur measure-side compactness
  - HYP-2571  # covering closed form M=7m/(84m+5); SDR≠loneliness; overtaking-tournament trivial (the honest ledger)
  - OPEN-Q-104
external: Lonely Runner Conjecture; the optimal τ for integer speeds is rational at a bounded denominator (classical); view-obstruction.
---

# THM-524 — The binding-pair reduction and the regions/sections reframe of LRC

**Setup.** `M(S) = max_{τ∈[0,1)} min_{v∈S} ‖vτ‖` is the LRC gap; LRC(N) is `M(S) ≥ 1/N`
for all primitive `S` (here N=14, |S|=13). The user's reframe: track the **N sections** of
the loop, not the runners. This note records what that buys (and what it doesn't), with the
real reduction proved.

## A. The binding-pair reduction (PROVED)

> **Lemma (sawtooth lower-envelope).** For positive integers `v_1,…,v_m`, the maximum of
> `f(τ)=min_i ‖v_iτ‖` over `[0,1)` is attained at a τ\* of one of two types:
> **(crossing)** `‖v_aτ\*‖=‖v_bτ\*‖=f(τ\*)` for two indices `a≠b` (a *binding pair*), forcing
> `τ\* = k/(v_a+v_b)` or `k/(v_a−v_b)`; or **(peak)** `‖v_iτ\*‖=½` (so `τ\*=(2k+1)/(2v_i)`),
> which occurs iff all speeds are odd (then `M=½`).

*Proof.* Each `‖v_iτ‖` is a triangle wave, continuous and piecewise-linear with breakpoints
at `τ=j/(2v_i)` (peaks `½` at odd j, troughs `0` at even j). On any maximal interval `I`
between consecutive breakpoints of *all* the waves, each `‖v_iτ‖` is affine, so `f|_I = min`
of affine functions is **concave** piecewise-linear, with interior kinks exactly where two
affine pieces cross. A concave PL function on `I` attains its max at an interior kink or at
an endpoint of `I`. An interior kink is a crossing of two pieces that are both the running
minimum: `‖v_aτ\*‖=‖v_bτ\*‖=f(τ\*)`, where `v_aτ−p = ±(v_bτ−q)` for integers `p,q`, i.e.
`(v_a∓v_b)τ\* ∈ ℤ`, giving `τ\*=k/(v_a±v_b)`. An endpoint is a breakpoint `j/(2v_i)`: a trough
(`f=0`, never the max unless `f≡0`) or a peak (`f=½` iff that wave is the min there, i.e. all
others are `≥½` — only possible if every speed is odd and `τ\*=½`). Taking the global max over
all `I` gives the claim. ∎

**Consequence (the reduction).** `M(S) = max( g(S,½), max over pairs (a,b)∈S and integers k
of g(S, k/(v_a±v_b)) )`, where `g(S,t)=min_v‖vt‖`. So **LRC(N) at S ⟺ this max ≥ 1/N**:
check the `≤ C(n,2)·2` candidate **pair crossings** (plus the all-odd peak), each with an
`O(n)` *others-clear* scan. A 13-runner condition becomes **~78 pairwise SWITCHES** + clearance
— **polynomial**, `O(n²·d·n)`. This is the precise "regions/switches, a few conditions"
the reframe asked for, and it is the rigorous foundation of the exact-M tool
(`lrc14_gap_M_exact`, validated vs dense grid; envelope candidate-set stress 0/300).

Two facts make it sharp: **(i)** the optimum always has `≥2` binding runners (min over 2000
random = 2), so a binding pair always exists; **(ii)** *others-clear* is independent and
necessary — e.g. `{1..13}` has a pair-only gap `½` at pair (1,3) but true `M=1/14`, so the
pair value alone overstates `M`.

## B. The regions/sections layer is the ON-GRID special case (VERIFIED)

At a grid time `τ=a/14` (a a unit mod 14), runner i sits in **section** `r_i = v_i·a mod 14`;
the observer is lonely iff no `r_i=0`. This is exactly the THM-523 q=14 witness:
> some grid time is lonely ⟺ every `v_i ≢ 0 (mod 14)` (verified 4000/4000),
and `(ℤ/14)^* = {1,3,5,9,11,13}` permutes the grid witnesses (`gridM(a·S)=gridM(S)`).

The user's **"each runner its own section" = perfect SDR** = residues distinct *and* nonzero
mod 14. **Honest caveat (HYP-2571): SDR is strictly stronger than loneliness**, not a
characterization — 133 thirteen-runner configs are lonely-on-grid but not SDR (a section
doubly occupied, e.g. residues `{1,2,2,4,…,13}`). And the section lens is **blind off-grid**:
the covering hard core `{1..11,13,84}` has `gridM=0` (runner 84≡0 sits in section 0 forever)
yet true `M=7/89` off-grid. So the regions view *describes the easy on-grid case beautifully
and is exactly the q-witness, but the binding-pair reduction (A) is what decides the hard
cases*. `{1..13}` is the maximally-spread perfect SDR at every unit a.

## C. The covering hard core condenses to a closed form (VERIFIED)

For the covering family `{1..11,13} ∪ {84m}` (84 = lcm(12,14), the minimal covering big
element): **`M = 7m/(84m+5)`**, strictly increasing in m, limit `1/12`, minimum at m=1 =
**`7/89`**. `{1,…,11,13,84}` is the unique closest-to-1/14-from-above covering set (452-set
scan). The whole 1/14 margin is one cross-multiplication: `7/89 > 1/14 ⟺ 98 > 89` (i.e.
`9=3²>0`). The exotic near-counterexamples thus reduce to a **finite checklist**: the small
binding flank lives in `{2,4,5,13}`, so "beats 1/14?" = check (small flank)×(big element)
sum-crossings. This is the gap-side companion to THM-523's covering reduction.

## D. The tournament ledger (honest: one exact bridge, the rest analogy)

- **EXACT — complement = reversal.** In `(ℤ/14)^* ≅ ℤ/2 × ℤ/3`, the element `a=13≡−1` sends
  section `r ↦ 14−r` — literally the *complement* of the section assignment, the order-2
  reversal, matching the tournament complement `T↦T^op` (fixing the antipode section 7, the
  self-complementary locus). Moreover the binding pairs of the SDR case `{1..13}` —
  `(1,13),(5,9),(3,11)` — all **sum to 14**, i.e. are complement pairs `v↔N−v`: the optimum's
  switch *is* the `v↦−v` reversal. And "M is a pairwise switch" mirrors the project's
  single-arc-flip model. These are genuine.
- **ANALOGY ONLY.** (i) The `ℤ/3` part `{1,9,11}` (cube roots of 1 mod 14, `9=3²`) echoes the
  project's Φ₃ / `H=3^α` world, but "multiply sections by 9" is not a proven functor to the
  3-cycle relabeling — same abstract `ℤ/2×ℤ/3`, no functorial Φ₃ bridge. (ii) **The overtaking
  "position tournament"** (`i→j` iff `frac(v_iτ\*) > frac(v_jτ\*)`) is **always transitive**
  (H=1, Ω empty) — `frac(·)` is a total order — so the naive "loneliness = #Hamiltonian paths
  (Rédei)" hope is **dead**. (iii) The section conflict graph (`i∼j` iff same section) gives
  "perfect SDR ⟺ edgeless ⟺ χ=1", suggestively like "Ω edgeless ⟺ ideal gas", but the
  M-optimum is a 2-body pair, not an independent-set property.

## E. Status & the next step

- A (binding-pair reduction): **PROVED** (elementary; generalizes to LRC(N) verbatim) +
  identity verified 0/11,000.
- B (regions = on-grid q-witness; SDR strictly stronger; off-grid blindness): **VERIFIED**.
- C (covering closed form `7m/(84m+5)`, finite checklist): **VERIFIED** (small-flank-palette
  conjectural, HYP-2571).
- D (tournament ledger): complement=reversal **EXACT**; Φ₃ and conflict-graph **analogy**;
  overtaking tournament **trivial** (proved).
- **Next:** tie *grid-sharpness* (`gridM=M`) to the binding *complement* pair `(v,N−v)`, and
  prove grid-sharpness fails *exactly* for covering sets — turning the section lens into a
  genuine off-grid criterion and merging it with THM-523. Then the residual `inf M` over
  covering sets (`≥7/89`?) is the same compactness frontier as THM-522.
