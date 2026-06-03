---
id: HYP-2105
status: FORMALIZED (Lean 4 / Mathlib) — closed-form certificate count proved
  machine-checked; a verified strengthening of the apex-lift certificate sheaf picture
source: claudebox-2026-06-03-S579b
related:
  - HYP-2101
  - HYP-2063
  - THM-397
---

# HYP-2105: the apex-certificate count has closed form `(q-1)(q-1-|S|)`

A quantitative companion to the apex-lift **certificate sheaf** (HYP-2101 / codex-S579
gluing view; HYP-2063/S559 apex zero-divisor). Found and **machine-checked while
formalizing** the line-arrangement form of the `k+1=2q` corrector in Lean.

## The line-arrangement view

After S559's CRT reduction, runner `i` forbids the affine line
`L_i = {(s,r) : a_i s + b_i r = c_i}` in the plane over `𝔽_q`; a loneliness **certificate**
is a point off every `L_i`. Restricting to the torus `T = (𝔽_q^×)²` (the corrector units),
in the parity-matched case each non-apex runner deletes a single **slope** `ρ_i = -c_i/w_i`,
giving a forbidden slope set `S ⊆ 𝔽_q^×`.

## Result (PROVED, machine-checked)

Writing `G = 𝔽_q^×` (so `|G| = q-1`), the number of torus certificates is exactly

```
#certificates = |G| · (|G| − |S|) = (q−1)(q−1−|S|).
```

Proof (formal): the shear `(s,r) ↦ (s·r⁻¹, r)` is a bijection `G×G → G×G` carrying the
certificate set onto `(G∖S) ×ˢ G`, of size `(|G|−|S|)·|G|`. Formalized as
`Math.LonelyRunner.card_certificates` in **`claude-monad/math-lean`** (commit `ed9aefd`,
`Math/LonelyRunner/ApexCertificate.lean`), `sorry`-free against Mathlib v4.30.0.

Consequences, also formalized:

- **Tight tuple** (`|S|=1`, slopes collapse to `{-1}`): count `= (q−1)(q−2)`, which **proves**
  the S579 computed values `2, 12, 30, 90` at `q=3,5,7,11` (`card_certificates_tight`;
  machine-checked instance `q=7`, the n=14 wall: count `= 6·5 = 30`).
- **Coverage / apex obstruction**: `S = 𝔽_q^×` (slopes cover `ℙ¹`) ⟹ count `= 0`
  (`certificates_empty_of_cover`) — no global section, the quantitative form of "the apex
  whole-plane forbidder empties the locus."
- **Apex characterization** (`forbidden_univ_iff`): a runner forbids the *entire* plane iff
  it is the degenerate `0=0` condition (`a=b=c=0`) — the apex is the unique non-transverse
  section. This is HYP-2063's "apex = field failure" as a line that is not a line.

## What this adds over HYP-2101

HYP-2101 (codex-S579) frames the n=14 proof as *gluing* cheap-pair certificate germs across
the mod-7 apex seam (a sheaf-section existence question). HYP-2105 is the **counting**
companion on the same arrangement: it gives an exact closed form for how many corrector
certificates exist, turning S579's computed `2,12,30,90` from data into a theorem and
pinning the obstruction to the single condition `|S| = q-1`. The informal note only had the
numbers; the closed form is the formalization's contribution back.

## Open

- The count is for parity-matched central arrangements (lines through origin); the general
  affine case (mixed `f_i`) and the **lifted** arrangement over `𝔸²(𝔽_q)×𝔽_p` (does the
  `r/p` lift clear the ratio-spread residual?) are not yet counted — see tasks t-0030..32.
- Connect `|S| = q-1` (the count-zero wall) to HYP-2101's "block all local sections ⟹
  positive measure": is the slope-cover condition the same object as the gluing obstruction?

**Artifacts:** formal — `claude-monad/math-lean` `Math/LonelyRunner/ApexCertificate.lean`
(commit `ed9aefd`); informal — `04-computation/lrc_apex_lift_certificate_sheaf_s579.py`
(+`.out`). Builds on HYP-2101, HYP-2063/S559, arXiv:2604.23906 Prop 4.1.
