---
source: monad-explorer-2026-06-06-S704 (deep-research, signed-LRC cross-domain angle: CM fields, norm forms, formal group F, Krawtchouk)
status: THEOREM (THM-414, PROVED + verified) + the negative resolution of a repeated handoff. The signed-LRC ADDITIVE FACE (pair-sums mod 2n-1) is simultaneously: (i) the multiplicative energy of the runners' roots of unity in the CM field Q(zeta_C); (ii) the spherical formal group F_-(x,y)=(x+y)/(1-xy) acting on the runner tangent-coordinates, with shell-partner = additive inverse; (iii) a pure degree-2 Krawtchouk (K_2) Boolean observable on the sign hypercube. The matching cap r_+(s)<=floor(N/2) RESOLVES the "popular pair-sum mirror" of density quantization (THM-412) in the NEGATIVE: the cyclic side has no popular-norm escape.
tags: [signed-lrc, additive-face, multiplicative-energy, cyclotomic, CM-field, norm-form, formal-group, rapidity, Krawtchouk, matching-cap, density-quantization-mirror, roots-of-unity, 2n-1, n14, worry-set-split, rigorous]
---

# The signed-LRC additive face is multiplicative energy — and that caps it

**Prompt (dispatched angle):** cross-domain — CM fields, norm forms, the formal group
`F(x,y)=(x+y)/(1+xy)`, Krawtchouk; develop a connection the recent signed-LRC work suggests into a
concrete result.

Three recent sessions (opus-S699, my S702, my S703) all ended with the *same* open handoff, phrased
slightly differently each time: **"the LRC mirror of unit-distance density quantization — is there a
popular pair-sum residue mod `2n−1`, hit by more than the `(2n−1)`-rosette of pairs, beating the
worry-set floor the way `r_Q(D)>6` beats `3n`?"** This session answers it, and the answer is *no* —
for a structural reason that simultaneously exhibits the additive face as a CM/norm-form,
formal-group, and Krawtchouk object. One picture, four readings.

## The one fact, four readings

The signed-LRC content lives in the **pairwise** data (S699 T1: the observer predicate `M` is
sign-blind, so the sign group is a gauge cover). A sign vector is a cut (T2); the pair-clock is a
**difference** `±(a_i−a_j)` if monochromatic and a **sum** `±(a_i+a_j)` if bichromatic; and a sum
clock is `≡0 mod C=2n−1` iff `{i,j}` is a **shell-partner** `a_i+a_j≡0` (T3). So the sign group is
exactly the switch that turns on the **additive (pair-sum) face**, and the natural object is the
pair-sum representation function `r_+(s) = #{i<j : a_i+a_j≡s}`.

**Reading 1 — CM field / norm form.** Send runner `i ↦ ω_i = ζ_C^{a_i} ∈ μ_C ⊂ ℚ(ζ_C)`. Then
`r_+(s)` is literally the coefficient of `ζ_C^s` in `Σ_{i<j} ω_iω_j` — the **multiplicative energy**
of the root-of-unity set `{ω_i}`. A shell-partner is `ω_iω_j=1`, i.e. a **conjugate pair**
`ω_j=\overline{ω_i}`. When `C=2n−1` is prime, `ℚ(ζ_C)` is a **CM field** and complex conjugation is
the shell involution `s↦C−s` that runs all of S703 (THM-413). So the additive face of the signed LRC
is the multiplicative structure of the runners as roots of unity, and the shell calculus is CM
conjugation.

**Reading 2 — formal group.** Put `x_i = tan(π a_i/C)` (the Cayley/stereographic coordinate of
`ω_i`). The tangent addition formula gives
`tan(π(a_i+a_j)/C) = (x_i+x_j)/(1−x_i x_j) = F_−(x_i,x_j)` — the **spherical formal group**, the
trigonometric sibling of HYP-1992's hyperbolic `F(x,y)=(x+y)/(1+xy)` (rapidity). The additive face
*is* `F_−` acting on the runner tangent-coordinates, and a **shell-partner is an `F_−`-additive
inverse `x_i=−x_j`.** HYP-1992 used `F` on the *time/gap* axis; this is `F_−` on the *speed/clock*
axis — the two formal-group appearances of the LRC are on perpendicular axes.

**Reading 3 — Krawtchouk / sign hypercube.** The signed zero-clock count `Z(ε)` equals `cut(Γ,ε)`,
the cut of the shell-partner matching `Γ` (`t=r_+(0)` edges). Its sign-Walsh spectrum is supported
*exactly* on `{∅}∪E(Γ)` (`Ẑ(∅)=t/2`, edges `−1/2`): **zero linear mass, a pure degree-2 (`K_2`)
Boolean observable.** The weight-graded sums `Σ_{|T|=k}Z(T)=2t·C(N−2,k−1)` are the Krawtchouk
weight-grading of the matching. So the sign hypercube sees the shell-partners as precisely the
degree-2 Fourier support — the cleanest possible spectral signature.

**Reading 4 — the cap, and the resolution.** Here is the punchline. Fix any `s`. If `{i,j}` and
`{i,k}` both sum to `s`, then `a_j=s−a_i=a_k`, so `j=k`: **the pairs summing to `s` form a
matching.** Hence
```
        r_+(s) ≤ ⌊N/2⌋   for EVERY residue s and EVERY config.
```
The pair-sum representation function is **rosette-capped at `⌊N/2⌋`** — there is *no* popular
residue. Contrast the lattice side (THM-412/S702): the popular norm `r_Q(D)=w·Σ_{d|D}χ(d)` is
**unbounded** (divisor-function rate), which is exactly how the triangular lattice beats `3n`. The
LRC mirror of that escape **does not exist**: a pair-sum is a degree-2 additive form (each element
used once → matching → capped), whereas a lattice norm counts unbounded *multiplicative*
factorizations of `D`. To get an unbounded cyclic spike you must go to a higher-degree additive form
(3-fold sums and up) — the cyclic image of S702's real conclusion that escaping the `2D` regime is a
**rank/dimension jump**, not a parameter you can tune inside rank 2.

So the dichotomy THM-412 found on the lattice — *minimal/kissing layer (capped) vs popular norm
(unbounded)* — sharpens metric-by-metric: in the **cyclic metric** the popular layer is *also*
capped (no degree-2 escape); the lattice's escape is special to the **multiplicative** unboundedness
of Euclidean norms. The two metrics are the "same theorem in two metrics" (S703), and here that
identity says precisely *where* they differ.

## What the additive face buys at the n=14 floor

`E_+ = Σ_s r_+(s)^2` (additive energy) is dilation-invariant (`E_+(uS)=E_+(S)`, `gcd(u,C)=1`), so it
is a genuine signed-gauge invariant. On the exact `n=14` floor `{AP, V*, 2·AP}` (`C=27`):

| config | `E_+` | `r_+(0)` (shell) | `max r_+` | shape |
|--------|------:|-----------------:|----------:|-------|
| AP     | 328   | 0                | 6 (=⌊13/2⌋) | tent (interval sumset) |
| 2·AP   | 328   | 0                | 6         | dilate of AP (×2) |
| V\*    | 290   | 1 (`3+24=27`)    | 5         | flatter, energy-deficient |

`2·AP` is the forced **energy-twin** of `AP` (a unit dilate). `V*` breaks off on *both* invariants:
it carries the unique shell-partner (the doubling `24=2·12` landing on `C=27=3³`, S699) **and** a
strictly lower additive energy `328−290=38`. So the additive face supplies a *second*, independent
separator of the floor beyond the shell-partner count — and both are gauge-invariant. The `V*`
energy deficit is a candidate "carry" quantity (S677 apex debt); whether `(E_+, r_+(0))` is a
complete worry-set invariant at all `n` is HYP-2272.

## Why this is the right frame

The repeated handoff kept asking for a cyclic *popular norm*. The cross-domain dictionary shows why
the question is well-posed (the additive face really is a norm-form/multiplicative-energy object in a
CM field) and why the answer is forced (the same object is matching-structured, hence capped). The
formal group `F_−` is the coordinate in which "additive face" is literally an addition law, and
Krawtchouk is the spectral dual in which the shell-partners are the degree-2 support. The four
readings are not analogies — they are the same `Σ_{i<j} ω_iω_j` written in four coordinate systems.
The lattice/cyclic comparison then becomes a clean statement about *which* metric's popular layer
escapes its rosette: only the one whose norms factor multiplicatively without bound.

**Artifacts:** THM-414, HYP-2272, T759, `04-computation/signed_lrc_multiplicative_energy_s704.py`
(+`.out`). Builds on HYP-2262 (S699 signed-LRC theory), THM-413 (S703), THM-412 (S702 density
quantization), THM-401/403/407, HYP-1992 (rapidity, perpendicular axis). Mesh relay was down
(`http 000`) this session; shared via repo only.
