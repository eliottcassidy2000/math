---
source: klein-2026-07-23-S407
status: LRC(14) working session (owner: "can the snippet be the key to LRC(14)?"). Computational, honest.
  Confirms/extends kps-S131 (1/25 trivially loose for a fixed 13-set ⇒ snippet is a GENERAL bound) and
  answers mac-mini-S169 (Σv² does NOT predict the SOS crossing degree). Integrates opus-S4b owner clue
  (1/25 = 24-speed floor) and the FORK. No new proof; a map of where the snippet does/doesn't fit.
tags: [lrc14, snippet, riesz, sos, curvature, second-moment, free-energy, fork, session]
---

# LRC(14) session — the snippet is a GENERAL bound, not a config; curvature is the weight, not the difficulty

**klein-2026-07-23-S407.** Owner reframe (relayed via opus-S4b): the snippet is an outside LLM's CRITICAL STEP
in a Lonely Runner proof — LRC confirmed. I worked the exact object `L(S)=meas{τ:‖vτ‖>1/14 ∀v∈S}` and the
Riesz/SOS certificate. Findings, robust ones first.

## 1. "> 1/25" is TRIVIALLY loose for a fixed 13-speed set (confirms kps-S131, direct data)
Threshold sweep of the Riesz certificate `∫MR/∫R < 1 ⇒ L>0` at the danger threshold `H`:
| core | `H=1/14` (conj) | `H=1/25` (snippet) | `H=1/26` (classical) |
|---|---|---|---|
| `{1..13}` tight | ratio 1.06 **stalls**, L=0 | ratio **0.12** fires, L=0.329 | 0.11 fires |
| `{1..13}\6∪56` loose | 1.09 stalls, L=0.0056 | 0.16 fires, L=0.320 | 0.14 fires |
| `{1..13}\6∪98` reson | 1.07 stalls, L=0.0052 | 0.22 fires, L=0.312 | 0.20 fires |
So at `H=1/25` the method fires with **huge margin** (L≈0.32, ratio≈0.15) on every 13-speed core; the crux is
the **knife-edge at `H=1/14`** (tight cores have L=0 exactly there). ⇒ `{1..13}` clears 1/25 trivially, so it is
the **second-moment REFERENCE/normalization** (`2457=3·Σ₁¹³k²=3·819`), **not the certified hard instance**
(independent of kps-S131, same conclusion). The snippet's real content is a **general/uniform** statement over
arbitrary 13-speed `V` — the union-bound-beating (Bedert) regime, not a fixed-config check.

## 2. Σv² does NOT predict the SOS crossing degree (answers mac-mini-S169)
Using mac-mini's sound SOS certificate `λ_min(T_M^{(N)})<1 ⇒ L>0` on the 13 drop-j cores `{1..13}\{j}∪{56}`
(`Σv²=3955−j²`): `corr(Σv², λ@N=80)=−0.53`, but `corr(L, λ@80)=−0.75`. The **hardest** core (`j=6`, crossing
degree ≈60, `λ@80=0.909`) has **middling** `Σv²=3919` and the **minimum** `L=0.0056`. So the SOS difficulty
tracks the **loneliness/resonance structure `L`**, not the scalar curvature. `Σv²` is the *weight* (kps-S131);
it is **decoupled** from the difficulty — matching klein-S405 ("curvature ⊥ additive/resonance structure") and
the standing guardrail "target `L=p₀`, the genuine order is RESONANCE, not the moment."

## 3. The FORK (kps-S131) under my data + the owner's 24-speed reframe (opus-S4b)
`X=(2457/6592)log_B−log_A=0.045725 > 1/25`. Two readings:
- **(a) wider-gap bound** (`gap≥1/25>1/26`, an unconditional union-bound improvement, INCOMPLETE for the full
  conjecture which needs `1/14`). My §1 data supports (a) **if the speed count is 13**: for 13 speeds the method
  reaches ≈`1/14`, so a general `>1/25` result is a weak sub-bound (Bedert regime), not the conjecture.
- **(b) surplus clearing a tight floor** (could be full): opus-S4b notes `1/25 = ML({1..24}) = 1/(24+1)` EXACTLY
  — the tight floor for **24 speeds**. Then `>1/25` clears the 24-speed floor by 0.0057, a genuine proof step.
- **The tension (unresolved):** the WEIGHT says 13 (`2457=3·Σ₁¹³k²`), the FLOOR says 24 (`1/25=ML({1..24})`).
  Candidate reconciliation: a REDUCTION where a 25-runner/24-speed target floor `1/25` is cleared while a
  13-speed sub-config's second moment (819) is the curvature weight — but that is exactly the **reduction** whose
  soundness opus-S4b flags as the load-bearing risk. Cannot settle without the clause stating what ">1/25" concludes.

## 4. Soundness (kps-S130/S131): the free-energy reading is the only sound one, and my curvature supports it
`∫M log R` is signed ⇒ not a direct certificate (kps-S130). The sound formulation is a **free energy / soft-max**
`X=(1/β)log∫e^{β g(τ)}dτ`, `g=min_v‖vτ‖` — which satisfies `X ≤ ML` (since `g≤ML ⇒ ∫e^{βg}≤e^{βML}`), so
**`X>1/25 ⇒ ML>1/25` soundly**, and `X→ML` as `β→∞`. My `Σv²=½E''(0)` is the **Laplace saddle-width** of such a
free energy (curvature of the tent energy at the resonance) — so the second-moment weight is *native* to the
free-energy reading (R2-valid), corroborating kps. Caveat (klein-S406): the *spectral* free energy `log∫R` is
freq-blind at 2nd order, so the `Σv²` must ride the **spatial** saddle (g's local curvature), not the amplitude
free energy — the sound functional is a saddle-point Laplace bound on the loneliness integral, not `∫M log R`.

## Verdict on the owner's question
The snippet is almost certainly **not the key to LRC(14)** in the sense of a new certificate the repo lacks: for
13 speeds `>1/25` is trivially loose, and the sound content is a *general wider-gap* bound (Bedert regime), which
is weaker than the `1/14` the repo already reaches per-config (mac-mini's SOS). Its genuine value to the repo is
the **T1 certified-artanh rigor layer** and the **free-energy/soft-max soundness frame** (a Lean-portable way to
bound a log-rate above a floor). If the external LLM claimed a FULL proof and its step is reading (a) at 13 speeds,
the proof is **incomplete at this step** (reaches `1/25`, not `1/14`) — the classic universality/insufficient-floor
failure mode opus-S4b warned of. The one clause around eq(27) settles it.
→ kps-S130/S131, opus-S4b, mac-mini-S169 (SOS), klein-S405/S406 (curvature, tension), THM-515/518.
