# The function-field route, tested: the packing is unconditional, but it does NOT force INV

*boxeph-2026-07-18-S91. Owner asked to prove the function-field INV via unconditional packing (the
S90 direction). I built the `𝔽_p[t]` runner problem and tested it. **Result: the packing IS
unconditional and exact — but it does NOT force the AP core, so the function-field route does not
close INV.** This corrects the S90 speculation. Scripts: `lrc_function_field_boxeph_S91.py`,
`lrc_ff_inv_test_boxeph_S91.py`.*

## The function-field runner problem (clean, and it works)

Over `A = 𝔽_p[t]`, `K_∞ = 𝔽_p((1/t))`, `‖x‖ = p^{deg{x}}`. For a rational time `a/Q`
(`deg a < deg Q`), `‖v·a/Q‖ = p^{deg((v·a) mod Q) − deg Q}`. The dictionary:

- **Threshold `1/p`** ⟺ the residue `(va) mod Q` has top coefficient (degree `deg Q − 1`) nonzero.
- **"Covering"** ⟺ the speeds' roots cover all of `𝔽_p`. **Sieve (automatic):** a non-root point `c`
  gives `‖v·1/(t−c)‖ = 1/p` for all `v` (evaluation `v(c) ≠ 0`) — non-covering ⟹ lonely at level 1.
- **Level** `= deg Q`. Covering ⟹ no level-1 loneliness ⟹ must use `deg Q ≥ 2`.
- **Deep well** `= 𝔽_p^* ∪ {t^p − t}`: the constants `𝔽_p^*` (the AP core `{1,…,p−1}`) plus the
  vanishing polynomial `t^p−t = ∏_{c}(t−c)` (the far killer, covering every point). Verified covering,
  no level-1, lonely at level 2, core top-coeffs `= 𝔽_p^*` exactly, killer collides in top-coeff.

## The packing IS unconditional and exact (the good news)

At a level-2 lonely time (`deg Q = 2`), each residue is linear `α_v t + β_v` with `α_v ∈ 𝔽_p^*`.
Since `α_{u−v} = α_u − α_v`, **a difference-closed set has pairwise-distinct top coefficients**, and
`p−1` distinct nonzero values in `𝔽_p^*` **must be all of `𝔽_p^*`** — with **no slack**. The
archimedean "`≤ 1` excess / sub-case B / `q=13val+1` ambiguity" genuinely **vanishes**: the ultrametric
makes the packing exact, exactly as S90 hoped.

> **FF Packing Lemma (proved).** A difference-closed `(p−1)`-family that is level-2 lonely has residue
> top-coefficients `= 𝔽_p^*`; hence it is `γ·𝔽_p^*` (a dilated `𝔽_p^*`) — the AP core.

## But it does NOT force the AP core (the decisive finding)

The Lemma **assumes** difference-closure; it does not derive it — and covering + level-2-loneliness do
**not** derive it. Counterexamples over `𝔽_5` (verified): the families

`{1, 2, 3, t} ∪ {t^5−t}`, `{1,2,3, t+1} ∪ {t^5−t}`, `{1,2,t,t+1} ∪ {t^5−t}`, …

are all **covering, level-2-lonely (`M = 1/5`), yet their 4-core is NOT a dilated `𝔽_5^*`** (top-coeffs
`[1,1,2,3]`, `[1,2,2,3]`, … — repeats). So

> **the naive FF-INV ("covering + level-2 lonely ⟹ AP core") is FALSE.**

**Why (and why this refutes the S90 shortcut).** The archimedean rigidity lives on the *fine*
condition `M ∈ (1/14, 1/13)` — a strict inequality producing the razor-thin band `q = 13val+1`. The FF
value group is **discrete** (`p^ℤ`): there is no value strictly between `1/p` and `1/p²`, so the deep
well sits **at** `1/p` (level 2), not below it, and "level 2" lumps the near-tight deep well together
with genuinely loose covering families. The discreteness that made the *packing* exact also **erases
the tightness gradient** that the *inverse theorem* needs. The unconditional packing is real but
insufficient: it forces `𝔽_p^*` **only** for a difference-closed core, and deriving difference-closure
is the same residual as in `ℤ` (THM-1017's open half).

## Corrected direction

- **S90 correction:** "prove FF-INV via unconditional packing, then lift" does **not** work as a
  shortcut. The packing transfers and becomes exact; the *hard half* (difference-closure of the core)
  does **not** transfer — it is destroyed, not solved, by the ultrametric.
- What the FF picture *does* give: a clean, exact model of THM-1017's **elementary half** (difference-
  closed core ⟹ `𝔽_p^*`/AP + far killer `t^p−t = lcm`-analog covering all points), confirming that the
  lcm-forcing and the packing are the "easy" transferable structure, and isolating **difference-closure
  derivation** (= Tao's n=12 additive rigidity) as the irreducibly archimedean/additive core.
- A finer FF invariant than `deg Q` (to recover the tightness gradient) would be needed for a genuine
  transfer — e.g. counting level-2 lonely times, or the Frobenius action of `t^p−t` on `𝔽_p[t]/(Q)`.
  Open whether any such refinement restores the inverse theorem.

Honest: LRC(14) is not closed, and the function-field shortcut is now **ruled out** — a useful negative.
The irreducible core remains the additive rigidity (difference-closure derivation / Tao n=12), which
does not simplify over `𝔽_p[t]`.

Cross-links: [[freiman-via-resonance-packing-function-fields-and-taos-optimistic-conjecture-boxeph-S90]]
(corrected here), [[THM-1017-ap-core-bridge-reduction]], [[HYP-7310]],
[[the-lrc14-crux-is-sharp-freiman-additive-energy-and-a-discrete-markoff-spectrum-boxeph-S89]].
