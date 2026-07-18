# The self-similar resonance is a scaling law to a fixed frame base, not a genuine recursion — `Q_s(d) = C·d²` with `C` frame-fixed

*boxeph-2026-07-18-S97. Testing whether the S96 "w plays a Lonely Runner against the mode lattice —
self-similarity one level up" is a GENUINE RECURSION (LRC(14)-resonance reducible to a lower-count LRC,
inductable) or a WALL. Verdict, from three verified computations: **NOT a genuine recursion.** It is a
one-step SCALING LAW terminating at a FIXED, settled LRC(≤6) base — single-generator (Dirichlet), no
tower. But this pins the resonant obstruction to a FIXED frame-local constant, which is a real advance.
LRC(14) not closed. Verified S97 computation (`lrc14_resonance_recursion_test_boxeph_S97.py`).*

## The two ways it could have been a recursion — both fail

For the density-route family `E_t = {1,…,6, t}` (≤6 bounded + one far element `t`), section `s`, with the
endpoint sum `S(n) = Σ_k ε_k e(n p_k)`:

**(1) Avoidance recursion — DEAD by forcing.** THM-729 / klein-S279 (line 57): the peel scale is
`w = d ≥ diam` — it **is** the removed far element / the scale-separation parameter, not a free choice.
So "choose a non-resonant `w`" (which *would* be a lower LRC) cannot apply: `w = d` is forced, and `d`
generates its own lattice `dℤ`, so the peel is maximally resonant. Verified in S96.

**(2) Structural recursion — DEAD by measurement.** Does the resonant structure at the forced peel
descend into a shrinking tower of LRC instances, or terminate? Three tests:

- **TEST A — scaling law, not a tower.** `|S(ta)|/t → |ν̂(a)|` for a **fixed** 7-vector `ν̂`
  (`|ν̂(1,2,3)| = 0.164, 0.235, 0.316`, matched to 3 digits by the measured `|S(ta)|/t` from `t=240` to
  `1920`). So `S(ta) = t·ν̂(a) + o(t)`: the family reproduces the **same fixed comb** `ν̂` at every scale,
  scaled by `t`. No shrinking sequence — a scaling law that terminates in **one step**.

- **TEST C — single scaling generator (Dirichlet, not a multi-runner LRC).** The mode comb lives on
  `tℤ + 7·{−k,…,k}` — the far-element lattice `t`, decorated by the **fixed** section modulus `7`
  (from `1/14 = 1/(2·7)`). Tested at `t = 241` (coprime to both `60` and `7`): the big teeth are still at
  `n ≡ 0, ±7, ±14 (mod t)`, with `n mod 60` showing **no pattern** — the frame `{1,…,6}` contributes
  **no independent scaling generator**. One scaling generator (`t`) + a fixed finite decoration ⟹
  "`w` vs the mode lattice" is a **1-runner Diophantine (Dirichlet) problem**, not a ≥2-runner LRC. A
  single runner has no LRC content (Dirichlet is unconditional); there is no smaller LRC to induct on.

- **TEST B — the base is fixed and settled.** `ν̂ = DFT(ν)`, `ν` = the far element's signed section
  measure, `ν̂ = [0, 0.164, 0.235, 0.316, 0.316, 0.235, 0.164]` (symmetric; `ν̂(0)=Σν=0`, a signed-count
  conservation law). `ν` is determined by the bounded frame — an **LRC(≤6)** object (the frame's
  loneliness is the empty section `s`). The "one level up" is literally one level, to the bounded frame,
  and it is **settled and fixed** — it does not iterate `LRC(13)→LRC(12)→⋯`.

## The exact consequence: `Q_s(d) = C·d²`, `C` fixed and frame-derived

`Q_s(t) = Σ_{ℓ≠0}|S(ℓt)|²/ℓ²` samples **only** exact multiples of `t`, where `|S(ℓt)| = t·|ν̂(ℓ)|`. Hence

> **`Q_s(t) = C·t²`, `C = 2·Σ_{ℓ≥1}|ν̂(ℓ mod 7)|²/ℓ²` — a FIXED constant depending only on the ≤6-frame.**

Verified both ways: `Q_s(t)/t² → 0.1365` (direct bilinear form, converged by `t=240`), and the
reconstruction `2·Σ_ℓ|ν̂(ℓ)|²/ℓ² = 2·0.0682 = 0.1364` matches to 3 digits. So `Error = √Q_s/t → √C ≈
0.369`, a **fixed constant** — not shrinking, no descent. This is the quantitative form of "not a
recursion": the resonant residual is a constant, computed from a fixed 7-vector, reached in one step.

## Why this is still an advance (the reframe)

The density route wanted `Q_s = o(r²)` (any power-saving). At the forced resonant peel this is now proved
**FALSE with the exact constant**: `Q_s(d) = C·d²`, `C ≈ 0.1365 > 0`. But the route never actually needed
`Error → 0` — it needs `Error < Φ_∞(frame)` (the reduced family's good-set floor). Both sides are now
**fixed constants of the bounded ≤6-frame**:

> The density row at the resonant peel closes ⟺ **`√C < Φ_∞(frame)`** — a fixed, finite, frame-local
> inequality between two LRC(≤6)-derived constants. Not an open recursion, not an open analytic bound.

This converts the S96 wall into a **decidable comparison**. (The boundary `√C = Φ_∞` is exactly tightness
`M = 1/14`, where `μ_0 = 0` — consistent with the deep-well family being both resonant and extremal.)
Computing `Φ_∞(frame)` in klein's normalization and checking the inequality is the concrete next step;
it requires the frame's main-term definition, not a new idea.

## Net (honest)

- **Answer to the question:** the self-similar resonance is **NOT a genuine recursion**. It is a scaling
  law `S(ta) = t·ν̂(a)`, single-generator (Dirichlet against `tℤ`, decorated by the fixed `7`),
  terminating in one step at a fixed, settled LRC(≤6) base `ν̂`. There is no LRC tower to induct on; the
  "Lonely Runner one level up" is a metaphor for a 1-runner Diophantine approximation, not a sub-LRC.
- **Verified new fact:** `Q_s(d) = C·d²` with `C = 2Σ|ν̂(ℓ)|²/ℓ² ≈ 0.1365` a fixed frame constant;
  `Error → √C ≈ 0.369`. So `Q_s = o(r²)` is false at the peel *with the exact constant*.
- **The advance:** the obstruction is now a **fixed frame-local inequality `√C < Φ_∞(frame)`**, a finite
  decidable comparison — strictly more tractable than the (false) `o(r²)` target.

LRC(14) not closed. The recursion hope is closed off cleanly, and the resonant wall is replaced by a
concrete two-constant comparison over the bounded frame.

Cross-links:
[[soft-weyl-closes-density-off-resonance-but-the-far-peel-is-a-lonely-runner-one-level-up-boxeph-S96]],
THM-886 (resonance law of Q_s), THM-729 (density = autocorrelation discrepancy),
[[both-lrc14-routes-bottom-on-the-same-multilinear-cancellation-not-one-mollification-klein-S279]].
