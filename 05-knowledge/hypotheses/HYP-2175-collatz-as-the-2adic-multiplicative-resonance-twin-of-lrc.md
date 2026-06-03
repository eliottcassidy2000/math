---
id: HYP-2175
status: SYNTHESIS + FORMALIZATION — Collatz is the same 2-adic-vs-multiplicative resonance question
  as the repo's LRC/doubling work; a precise dictionary, with the cycle resonance equation and the
  2-adic shift formalized in Lean. Dictionary verified small-scale; the conjectures are the open part.
source: claudebox-2026-06-03-S614
related: [HYP-2117, HYP-2155, HYP-2140, HYP-2145, THM-004, HYP-2063]
---

# HYP-2175: Collatz is the 2-adic/multiplicative resonance twin of the Lonely Runner

The Collatz conjecture and the repo's LRC/doubling work are the **same question**: does a
multiplicative orbit RESONATE with the 2-adic structure, and is the resonance trivial except at the
base? The repo's machinery — the doubling map / 2-adic seam (HYP-2117), the resonance energy
(HYP-2155), rigidity = orbit-type (HYP-2140), the divisor blocks (HYP-2145) — is the shared toolkit.

## The dictionary

| structure | Lonely Runner | Collatz (shortcut map `T`) |
|---|---|---|
| dynamics | runner phases `v_i·t mod 1` | `T(n)=n/2` (even) / `(3n+1)/2` (odd) |
| 2-adic side | doubling `x↦2x mod n`; seam at even `n=2q` (HYP-2117) | the `÷2` steps; `v₂(3n+1)` per step = the per-step rigidity-height (S596 H6) |
| multiplicative side | the speeds `{v_i}` | the `×3` factor |
| **RESONANCE = obstruction** | additive `Σ m_i v_i = 0`; the resonance energy `Σ∏\|ĝ\|` (HYP-2155) | multiplicative `2^K = 3^L·∏(1+1/3n_i)`, i.e. `2^K ≈ 3^L` (a cycle) |
| the conjecture | no resonant config dips below `δ` (sidestep the resonance) | no nontrivial resonance (no cycle); every orbit → 1 |
| randomness (**Lemma A**) | circuit-free ⇒ equidistribution ⇒ `G ≥ δ` | balanced parity signature (odd-density `< log₃2`) ⇒ contraction (Tao's almost-all) |
| structure (**Lemma B**) | 3-term/AP, high resonance energy, tight | resonant cycles `2^K≈3^L`, the structured core |
| **binary signature** | tiling signatures (THM-004 `inshat = 2·desc+1`) | the parity vector (Lagarias): bijection `[0,2^K) ↔ {0,1}^K` |
| the `+1` / translation | the AP translation (additive face); IFS `T`-branch `+1` (HYP-2117) | the `+1` in `3n+1` (the same `T`-branch `+1`) |
| orbits | every config reaches a lonely time | every `n` reaches the cycle `{1}` |

## Verified this session

- `v₂(3n+1)` over odd `n` is geometric `P(k)=2^{−k}`, `E[k]=2` ⇒ drift `3/2^2 = 3/4 < 1` (contraction)
  — the per-step 2-adic rigidity-height, the Lemma-A drift.
- A Syracuse `L`-cycle satisfies `2^K·∏n_i = ∏(3n_i+1)`, hence `2^K ≥ 3^L` (a resonance from above);
  only the trivial `L=1` cycle is feasible with small elements — the multiplicative `Σ m_i v_i=0`.
- The parity-vector bijection `[0,2^K)↔{0,1}^K` (verified `K=3..10`): the parity vector is a 2-adic
  binary signature, `T` a 2-adic shift = HYP-2117 doubling dynamics + the `3x+1` `+1`-twist.
- Trajectory odd-density `≈ 0.48 < log₃2 = 0.631` ⇒ contraction (the Lemma-A randomness side).

## Formalized (math-lean, sorry-free)

- `Math/Collatz/Resonance.lean`: **`cycle_resonance`** (`2^K·∏n = ∏(3n+1)`), `three_pow_le_two_pow`
  (`2^K ≥ 3^L`), `even_three_mul_add_one`. The cycle = the multiplicative resonance.
- `Math/Collatz/Parity.lean`: the shortcut map + branch lemmas, `half_modEq`, **`shortcut_mod_pow`**
  (`n ≡ m mod 2^{K+1} ⇒ T(n) ≡ T(m) mod 2^K`) — the 2-adic shift / parity-bijection inductive heart.

First Collatz content in the repo's Lean; built on the same `Nat`/`Fin` machinery as the LRC modules.

## The unifying thesis

Both problems pit a **multiplicative structure** (the speeds; the `×3`) against the **2-adic
structure** (the doubling/clock; the `÷2`), and ask whether they RESONATE. LRC: the resonance is the
additive relation `Σ m_i v_i=0`, large on the AP (the high-energy core), sidestepped by construction.
Collatz: the resonance is `2^K≈3^L` (a cycle), and the conjecture is that no nontrivial one is
realizable. In both, the 2-adic seam (even `n`; the `÷2` contraction) is where the action is, the
"randomness" side (circuit-free; balanced parity) gives the generic bulk (Lemma A / Tao), and the
"structured" resonance is the hard core (Lemma B / cycles). The repo's doubling/resonance/orbit
toolkit transfers directly.

## Open / next

- Formalize the parity-vector bijection in full (iterate `shortcut_mod_pow`), and the
  no-small-nontrivial-cycle bound (the resonance `2^K≈3^L` Diophantine, Steiner/Simons–de Weger).
- The "sidestep the resonance" (S605) analog for Collatz: is there a construction (the inverse
  doubling-IFS tree, HYP-2117) that covers ℕ, dual to the LRC sieve?
- Whether the repo's divisor-block / CRT-split (HYP-2145) gives a `mod 2^a · 3^b` block structure
  for Collatz trajectories.

**Artifacts:** computation inline; math-lean `Math/Collatz/Resonance.lean`, `Parity.lean`.
Builds on HYP-2117 (doubling/IFS/seam), HYP-2155 (resonance energy), HYP-2140 (orbit-type),
HYP-2145 (divisor blocks), THM-004 (binary signatures).
