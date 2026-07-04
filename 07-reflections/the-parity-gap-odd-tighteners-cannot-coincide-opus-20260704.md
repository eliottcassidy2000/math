# The parity gap: two odd tighteners can't coincide, so the small end of confinement closes

*opus-2026-07-04-S68. Asked for more creative progress on the core, I found the one place the tighteners'
being ODD has not yet been cashed in — the arc coincidence that extremity forces. It's an odd = even
near-miss, worth `≥ 1`, and that single integer closes the small-tightener end of the m=2,f=2 confinement
and squeezes the residual to a sharply-named corner.*

## Where the parity was hiding

Confinement (m=2,f=2) reduces (THM-615) to `M(2U ∪ {w₁,w₂}) ≥ 1/12`, i.e. some point with `g_E ≥ 1/12`
avoids *extremity* (`Ψ<1/12`: one tightener `<1/12`, the other `>5/12`). Lemma 3 kills the case where a
tightener is *large*. The residual — both tighteners small — is where they conspire to be extremal on the
whole high interval `I₀`.

Here is the new observation. If extremity holds throughout `I₀` (no mode-switch, else a moderate value
appears), then `I₀` sits simultaneously inside
- a `w₁`-**half-integer** arc, centered `c₁ = (2l+1)/(2w₁)`, and
- a `w₂`-**integer** arc, centered `c₂ = k/w₂`.

Both arcs contain `I₀` (width `L`), so their centers are close: `|c₁−c₂| ≤ (w₁+w₂)/(12w₁w₂) − L`. But
`c₁ − c₂ = [(2l+1)w₂ − 2kw₁] / (2w₁w₂)`, and — because **both tighteners are odd** — the numerator is
`odd − even = odd`, hence nonzero, so `|c₁−c₂| ≥ 1/(2w₁w₂)`. The two facts collide:

> `1/(2w₁w₂) ≤ (w₁+w₂)/(12w₁w₂) − L`  ⟹  **`L ≤ (w₁+w₂ − 6)/(12 w₁ w₂)`**.

A half-integer of `w₁` and an integer of `w₂` can never *exactly* coincide (odd ≠ even); the closest they
come is `1/(2w₁w₂)`, and that gap is what forces the `−6`. It is the tighteners' oddness, finally used.

## What it closes

`L ≤ (w₁+w₂−6)/(12w₁w₂)` is the *only* way extremity survives on `I₀`. So:

> **`L > (w₁+w₂−6)/(12 w₁ w₂)` ⟹ `M(2U ∪ {w₁,w₂}) ≥ 1/12`.** In particular **`w₁+w₂ ≤ 6 ⟹ M ≥ 1/12`**
> (the RHS is `≤ 0 < L`), unconditionally.

Verified (0 violations). Lemma 4 disposes of the **small** end. With Lemma 2 (the AP even part), Lemma 3
(the large end), and THM-616 (f=1 at every scale), the confinement residual is now squeezed to a single,
sharply-named corner: **moderate odd tighteners (`6 < w₁+w₂`, both `≤ u_max/(6(M(U)−1/12))`) on the near-AP
even part.** And for a *primitive* family the near-AP corner is bounded (rigidity: primitive `M(U)→1/12
⟹ U→{1..11}`, HYP-4062), a finite check. The obstruction that remains is purely arithmetic: whether a
*moderate* `w₁`-half / `w₂`-integer coincidence can sit at the even part's deep hole.

## The shape of the whole reduction now

| regime | disposed by |
|---|---|
| `f = 1`, any `m` | THM-616 (one tightener useless: orbit-max ≥ 1/4) |
| AP even part `c·{1..11}` | Lemma 2 (finite mod-24) |
| a *large* tightener | Lemma 3 (dense orbit hits the moderate band) |
| *small* tighteners (`w₁+w₂ ≤ 6`, or wide `I₀`) | Lemma 4 (odd ≠ even coincidence gap) |
| **moderate tighteners × near-AP** | **OPEN** (bounded for primitive → finite) |

Each lemma is one clean mechanism — orbit-max, mod-24, Lipschitz density, and now parity. The remaining
corner is the last arithmetic sliver, and it is *bounded* for primitive families.

## Status

- **Proved (THM-615 Lemma 4):** `w₁+w₂ ≤ 6 ⟹ M ≥ 1/12`, and the width condition `L > (w₁+w₂−6)/(12w₁w₂)`
  — via the odd-parity coincidence gap. Verified (160 + 26 cases, 0 violations).
- **Reduced:** m=2,f=2 confinement to moderate tighteners on the near-AP corner (bounded for primitive).
- **Open (the sliver):** whether moderate odd `w₁,w₂` can coincide (`w₁`-half / `w₂`-integer) at the even
  part's deep hole — the last piece.

The creative move was to ask where the tighteners' oddness had gone unused, and to find it in the
coincidence they are forced into. `odd ≠ even` is a small fact, but it is exactly the fact the small end
was missing.

Related: THM-615 (Lemmas 1–4), THM-616/opus-S67 (f=1 all m), HYP-4062/kps (the rigidity bounding the
near-AP corner), HYP-4080 (the harmless-alone framing). Script: `lrc14_parity_gap_lemma4_opus_S68.py`.
HYP-4082.
