---
id: HYP-2120
status: EXPLORATION + PARTIAL RESULT — the good-time measure is a relation-lattice theta /
  Molien constant term; corrections are graded by relation LENGTH (each extra term ~ 1/(k+1)
  smaller), so the 3-term count N₃ is the leading-order obstruction. Verified: the AP-translation
  flip (energy-invariant, 3-term-destroying) drives tight↔safe. Full discrepancy bound OPEN.
source: claudebox-2026-06-03-S585
related:
  - HYP-2064
  - HYP-2065
  - HYP-2115
---

# HYP-2120: the LRC good-time measure is a relation-lattice theta function, graded by character length

A two-lens (abstract-functional ⟷ representation-theory) account of the Lemma A (circuit-free /
randomness) vs Lemma B (3-term / structure) dichotomy, with a quantitative backbone for Lemma A.

## The dictionary

| object | abstract-function (Haskell) | representation theory |
|---|---|---|
| runner `v` at time `t` | `ev v t = {v·t}`, curry on `v` ⇒ character | `χ_v(t)=e^{2πivt}`, irrep of `ℝ/ℤ` |
| speed set `S` | `v : Fin k → ℤ` | Laurent poly `p_S=Σ_{v}x^v`, a virtual character |
| **3-term `v_c=v_a+v_b`** | a fold `t·v_c = t·v_a+t·v_b` | **fusion** `χ_a·χ_b=χ_c` |
| circuit-free | no folds | multiplicatively-independent characters |
| tournament arc `a→b` | curry `ev` at witness `t*` ⇒ `sign Im χ_{v_a−v_b}(t*)` | arcs indexed by the **difference set** `S−S` |

## The good-time measure is a theta / Molien constant term (PROVED identity, verified)

Let `δ=1/(k+1)`, `B_i={t:‖v_i t‖<δ}`, and `good = {t : all runners lonely}`. Fourier-expanding
the single-runner indicator `1_{‖vt‖≥δ}=Σ_m ĝ(m)e^{2πi m v t}` (`ĝ(0)=1−2δ`, `ĝ(m)=−sin(2πmδ)/(πm)`):

> **Identity.** `meas(good) = ∫₀¹ Π_i 1_{‖v_i t‖≥δ} dt = Σ_{m∈Λ} Π_i ĝ(m_i) = [z⁰] Π_i G(z^{v_i})`,
> where `Λ = ker(v : ℤ^k→ℤ) = {m : Σ m_i v_i = 0}` is the **relation lattice** and `G(z)=Σ_m ĝ(m)z^m`.

So `meas(good)` is a **theta sum over the relation lattice** — equivalently the multiplicity of the
trivial character (a Molien-type constant term) in `Π_i G(χ_{v_i})`. Validated numerically: the
relation-sum reproduces the direct grid overlap (`{2,3,5}`: 0.061 vs 0.060; `{1,2,3}`: 0.100 vs
0.100), and the **m=0 term is exactly `(1−2δ)^k`** — the independence/Weyl baseline `→ e^{−2}>0`.

## The character-length grading (the Lemma A backbone — verified)

Every nonzero `m∈Λ` contributes `Π ĝ(m_i)`. Since distinct speeds forbid 2-term relations
(`v_a=v_b` impossible), the shortest relation has `ℓ¹`-length **3** (a 3-term sum `a+b=c` or a
doubling `2a=b`). Each coordinate `m_i≠0` carries a factor `ĝ(m_i)=O(sin 2πδ)=O(δ)=O(1/(k+1))`:

> **Grading.** a relation of length `L` contributes `O(δ^L)`. Hence
> `meas(good) = (1−2δ)^k + O(N₃·δ³) + O(δ⁴)`, where `N₃ = #{3-term relations}` (the `ℓ¹=3`
> lattice vectors). **Each extra relation-term costs a factor `~1/b̂(1) ~ (k+1)/2`** (verified:
> 3-term/4-term correction ratio 1.98→2.87→3.80→4.75 at k=4,6,8,10 — grows with k).

**Consequences (all matching the de-risking the user reported):**
- **Circuit-free (`N₃=0`) ⇒ Lemma A:** corrections are `O(δ⁴)`, suppressed relative to the `O(δ³)`
  hard term, so `meas(good)=(1−2δ)^k(1+O(δ))>0` — the `~+0.08` stable margin is the `(1−2δ)^k`
  baseline minus a vanishing tail. Verified: circuit-free good-measure hugs `(1−2δ)^k≈0.133`,
  margin `G−δ>0` through k=8.
- **4-term-rich is safe:** 4-term relations (`a+b=c+d`, the additive energy) are length-4, `O(δ⁴)`,
  suppressed by `~(k+1)/2` vs 3-term — so high energy alone never makes a config hard. Verified.

## The AP-translation flip (verified — isolates 3-term as THE order parameter)

`{m+1,…,m+k}` has **translation-invariant additive energy** but its 3-term count collapses once
`m≥k` (sums exceed the max). At k=6: energy fixed `146`, `N₃ = 9→4→1→0→0` for `m=0,2,4,6,8`, and
the margin flips `G−δ = −0.000→+0.130→+0.190→+0.225→+0.248`. **`{1,…,6}` is tight not because it is
an AP or high-energy, but because it is the maximal 3-term-closed set; its translate `{7,…,12}`
(same energy) is safe.** In the theta picture: translation keeps the energy (length-4 lattice
vectors) but lengthens the minimal relation `ℓ¹`-distance `λ₁(Λ)` from 3 to ≥4.

## Open / next (the two unproven halves)

- **Lemma A discrepancy bound (analytic):** sum the `O(δ⁴)` tail rigorously — bound
  `Σ_{m∈Λ, m≠0} |Π ĝ(m_i)| < (1−2δ)^k` for `λ₁(Λ)≥4`. The grading says each term is small; the
  count of length-`L` relations in `Λ` is the missing ingredient (a lattice point-count / the
  number of length-`L` additive configurations in `S`).
- **Lemma B (the fold delivers 1/(k+1)):** when `N₃>0`, the fusion `χ_a χ_b=χ_c` confines the orbit
  to a subtorus; the theta sum is a delicate resonance (the `{1..6}` corrections cancel the constant
  to exactly 0). The 2025 prime-divisibility machinery plugs in here.

## Formalized crux (math-lean `982473f`)

`Math/LonelyRunner/SumFree.lean` (sorry-free): `sumFree_of_lo_add_lo_gt_hi` and
`consecutive_block_sumFree` — a speed set above its own diameter (`hi < lo+lo`) has **no** 3-term
relation `a+b=c`. This is the machine-checked mechanism of the AP-translation flip: translating a
set above its length destroys all 3-term relations while preserving additive energy, so it lands
in the circuit-free / Lemma A regime regardless of energy.

**Artifacts:** formal — `claude-monad/math-lean` `Math/LonelyRunner/SumFree.lean` (`982473f`);
informal — `04-computation/lrc_character_sieve_exploration_s585.py` (+`.out`),
`07-reflections/lrc-two-lens-character-grading-s585.md`. Builds on HYP-2064/2065 (sieve / thin
core = the same circuit-free/structured split), HYP-2115 (the CRT/sieve thread).
