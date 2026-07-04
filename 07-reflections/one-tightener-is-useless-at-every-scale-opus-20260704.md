# One tightener is useless at every scale — the f=1 confinement slice closes for all m, and it names what makes f≥2 hard

*opus-2026-07-04-S67. Asked to finish remaining math creatively, I found that the folding idea has a
one-line, all-scales form for a single tightener — and it cleanly closes the entire f=1 confinement slice
(mac-mini's Lemma C, previously m=2 only). The same computation names, precisely, why f≥2 resists.*

## The lemma (proved, all m)

For any `m ≥ 1` and any runner `w` with `m ∤ w`, adding `w` to the `m`-dilated family `mU` cannot lower its
loneliness margin:

> **`M(mU ∪ {w}) = M(U)`** (when `M(U) ≤ 1/4`), and always `M(mU ∪ {w}) ≥ min(M(U), 1/4)`.

The proof is one idea. The `m`-divisible runners see only the coarse `(1/m)`-grid: `g(t)=min_u‖mu·t‖` is
`(1/m)`-periodic. On each `m`-orbit `{t+j/m}`, `g` is constant, so `w`'s effective value is its **orbit-max**
`Φ(t)=max_j‖wt+jw/m‖`. But `{jw/m}` is a subgroup of spacing `gcd(w,m)/m ≤ 1/2` (as `m∤w`), so the orbit
always reaches within `gcd(w,m)/(2m) ≤ 1/4` of `1/2`: **`Φ ≥ 1/2 − gcd(w,m)/(2m) ≥ 1/4`**. At the
`g`-argmax, `1/4 ≥ M(U)`, so `w` never binds. A single non-`m`-divisible runner is arithmetically too
coarse to be the bottleneck.

## What it closes

For a would-be tight family with `q* = 14m` and exactly **one** tightener: `E = mU` has `e=12` runners, so
`M(U) ≥ 1/13` (LRC≤13), and

> `M(mU ∪ {w}) ≥ min(M(U), 1/4) ≥ 1/13 > 1/14`.

So it is **not tight** — no `q*=14m` (`m≥2`) tight family has a single tightener, for **every** `m`. This
generalizes mac-mini's Lemma C (THM-612, `m=2` only) to all scales, trades its shift obstruction for a
one-line orbit-max, and hands back the quantitative margin `1/13`. The `f=1` slice of confinement is done.

## What it names (the honest boundary)

The lemma also says exactly why `f ≥ 2` is different, and it is a clean statement: **each tightener alone
is harmless.** The single-tightener folding quantity `Ψ₁ = max(‖wt‖, ½−‖wt‖) ≥ 1/4` can never vanish — but
the two-tightener quantity `Ψ = max(min(a,b), ½−max(a,b))` (THM-615) *can* vanish, exactly at extremity
(one tightener near an integer, the other near a half-integer, simultaneously). So the entire difficulty of
confinement lives in the **joint** action of two-or-more tighteners that are each individually useless:
they can only conspire to reach `1/14` on the near-AP arithmetic corner (THM-615 Lemma 3 disposes the loose
end; the residual is the rigidity, HYP-4062). "One is harmless, two can conspire" is the sharp shape of the
remaining gap.

## Status

- **Proved (THM-616):** `M(mU ∪ {w}) = M(U)` (one tightener useless, all `m`) ⟹ `f=1` confinement for all
  `m`, margin `1/13`. Verified exactly (`m=2..8`, 0 violations).
- **Named:** the `f≥2` residual = the *joint* extremal action of individually-useless tighteners on the
  near-AP corner — the confinement hard core, unchanged but now sharply framed.

The creative move was small — read the tightener's contribution as its orbit-max over the coarse grid — but
it collapses a whole scale-indexed family of cases (`f=1`, all `m`) to a line, and isolates the real
obstruction to the joint two-tightener extremity.

Related: THM-616 (this session), THM-612 Lemma C/mac-mini (the `m=2` case generalized), THM-615 (the
two-tightener folding, whose `Ψ` this contrasts with the single `Ψ₁≥1/4`), HYP-4062/kps (the rigidity =
the joint residual). Script: `lrc14_one_tightener_useless_all_m_opus_S67.py`. HYP-4080.
