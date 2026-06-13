# THM-414 — The signed-LRC additive face: matching cap, degree-2 Krawtchouk structure, and cyclotomic/formal-group realization

**Status:** PROVED (elementary; every part verified by brute force, n≤7 exhaustive over all
2^{n-1} sign patterns, and the matching cap + dilation invariance over 4000 random configs).
**Source:** monad-explorer-2026-06-06-S704. Cross-domain (CM field / norm form / formal group /
Krawtchouk) development of the signed-LRC theory (S699 T1–T4, HYP-2262) and the repeated
S699/S702/S703 handoff ("the LRC mirror of unit-distance density quantization"). Builds on
THM-401 (shell modulus `C=2n−1`), THM-403 (`C`-th roots of unity), THM-413 (order-3 silent flip),
THM-412 (density quantization, the lattice side), distinct axis from HYP-1992 (rapidity).

---

## Setup

`N = n−1` runners with **distinct nonzero speed residues** `a_1,…,a_N ∈ (ℤ/C)^×`, `C = 2n−1`
(the generic case; the AP and its dilates and the `n=14` floor rows all qualify). A sign vector
`ε∈{±1}^N` is a 2-coloring (a **cut**) of the runners (S699 T2). The ordered **pair-clock** is
`c_{ij}(ε) = ε_i a_i − ε_j a_j (mod C)`. A pair `{i,j}` is a **shell-partner** iff
`a_i + a_j ≡ 0 (mod C)` (S699 T3 = THM-401).

Define the **pair-sum representation function** and **additive energy**
```
   r_+(s) = #{ i<j : a_i + a_j ≡ s (mod C) },        E_+ = Σ_{s∈ℤ/C} r_+(s)^2 .
```
`r_+(0)` is the number of shell-partner pairs.

---

## Part 1 (the matching cap) — PROVED

> **For every residue `s` and every config, the pairs `{i,j}` with `a_i+a_j≡s` are
> vertex-disjoint (form a matching); hence `r_+(s) ≤ ⌊N/2⌋`.**

**Proof.** Fix `s`. If `{i,j}` and `{i,k}` both sum to `s`, then `a_j ≡ s−a_i ≡ a_k`, so `j=k` by
distinctness of residues. Thus each runner lies in at most one pair summing to `s`; the pairs are
disjoint and number at most `⌊N/2⌋`. ∎

**Corollary (resolution of the "popular pair-sum" mirror — in the NEGATIVE).** The cyclic LRC
additive face admits **no popular-sum escape**: no residue is hit by more than `⌊N/2⌋` pairs, for
any config. This is the sharp contrast with the lattice side (THM-412 / S702), where the popular
norm `r_Q(D)=w·Σ_{d|D}χ(d)` is **unbounded** (divisor-function rate, Wigert). The repeatedly-asked
LRC mirror of "a popular norm `r_Q(D)>6` beating the kissing rosette `κ≤6`" *does not exist*: pair
sums are a degree-2 (each element used once) additive form, structurally matching-capped, whereas a
lattice norm counts unbounded multiplicative factorizations. Escaping the cap requires a
higher-degree additive form (3-fold sums, …) — the cyclic image of S702's rank/dimension jump.
Verified: `max_s r_+(s) ≤ ⌊N/2⌋` over 4000 random configs (0 failures); `=⌊N/2⌋` at the center for
the AP.

---

## Part 2 (degree-2 Krawtchouk structure of the signed zero-clock count) — PROVED

Let `Γ` be the **shell-partner graph** (vertices = runners, edges = shell-partner pairs). By Part 1
with `s=0`, `Γ` is a **partial matching**; set `t = |E(Γ)| = r_+(0)`. For a sign pattern `ε`, let
`Z(ε)` = number of pairs giving a zero clock (`c_{ij}(ε)≡0`).

> **(a)** A pair `{i,j}` is a zero-clock under `ε` iff it is a shell-partner **and** bichromatic
> (`ε_i≠ε_j`). Hence `Z(ε) = cut(Γ, ε)`.
>
> **(b)** The sign-Walsh (Fourier) spectrum of `Z` on the cube `{±1}^N` is supported **exactly** on
> `{∅} ∪ E(Γ)`:  `Ẑ(∅)=t/2`,  `Ẑ({i,j})=−1/2` for each shell edge, `0` otherwise. In particular
> `Z` has **zero linear (`K_1`) Fourier mass** — it is a pure degree-2 function.
>
> **(c)** Generating function and weight-grading (Krawtchouk weight sums):
> `Σ_ε x^{Z(ε)} = 2^{N−t}(1+x)^t`, so `Z` over the sign group is `2^{N−t}·Binomial(t,½)`; and
> `Σ_{|T|=k} Z(T) = 2t·C(N−2, k−1)` for the Hamming sphere of weight `k`.

**Proof.** (a) `ε_i a_i ≡ ε_j a_j`: if `ε_i=ε_j` then `a_i≡a_j` (false); if `ε_i=−ε_j` then
`a_i+a_j≡0` (shell-partner). (b) Write `ε_i = χ_i∈{±1}`; then
`Z(ε)=Σ_{\{i,j\}∈E(Γ)} (1−χ_iχ_j)/2`, whose Walsh expansion is read off term by term. (c) Each
matched edge independently is mono (2 of 4 sign choices, `Z+\!=0`) or bichromatic (`Z+\!=1`), giving
factor `(2+2x)`; the `N−2t` unmatched runners give `2^{N−2t}`. The sphere sum counts, per edge, the
`k`-subsets containing exactly one endpoint: `2·C(N−2,k−1)`. ∎

The Krawtchouk reading: `K_k(j;N)=Σ_i(−1)^iC(j,i)C(N−j,k−i)` is the Walsh transform of the
weight-`k` sphere indicator; by (b) the only nonzero zonal modes of `Z` are `K_0` and the degree-2
edge characters, so the signed zero-clock count is a **`K_2`-supported Boolean observable** — the
sign-hypercube dual of the shell-partner matching. `t` (the zero-clock capacity) is the unique
sign-gauge invariant it carries.

---

## Part 3 (cyclotomic / CM and formal-group realization) — PROVED

Map runner `i ↦ ω_i = ζ_C^{a_i} ∈ μ_C ⊂ ℚ(ζ_C)`.

> **(a)** Shell-partner `⟺ ω_iω_j = 1 ⟺ ω_j = \overline{ω_i}` (a **conjugate pair** of roots of
> unity). When `C` is prime, `ℚ(ζ_C)` is a CM field and complex conjugation is the shell involution
> `s↦C−s` (THM-413's fold).
>
> **(b)** `r_+(s)` equals the **multiplicative convolution** of the set `{ω_i}` (the coefficient of
> `ζ_C^s` in `Σ_{i<j} ω_iω_j`). So *the signed-LRC additive face is the multiplicative energy of the
> runners' roots of unity on `μ_C`.*
>
> **(c)** With the tangent (Cayley/stereographic) coordinate `x_i = tan(π a_i/C)`,
> `tan(π(a_i+a_j)/C) = (x_i+x_j)/(1 − x_i x_j) = F_−(x_i,x_j)`, the **spherical formal group** (the
> trigonometric sibling of the hyperbolic `F(x,y)=(x+y)/(1+xy)` of HYP-1992). Shell-partner
> `⟺ x_i = −x_j`, the `F_−`-additive inverse. So the additive face IS the formal group law acting on
> the runner tangent-coordinates, and the shell involution is `F_−`-negation.

**Proof.** (a) `ζ^{a_i}ζ^{a_j}=ζ^{a_i+a_j}=1 ⟺ a_i+a_j≡0`; conjugate since `|ω|=1`. (b) immediate
from `ω_iω_j=ζ^{a_i+a_j}`. (c) tangent addition formula; oddness of `tan`. ∎

---

## Secondary invariant (the dilation–energy split of the n=14 floor) — VERIFIED

`E_+` is invariant under dilation `S ↦ uS` (`gcd(u,C)=1`), since `r_+` is permuted by `s↦us`
(verified, 4000 trials, 0 failures). On the exact `n=14` floor `{AP, V*, 2·AP}` (C=27):
```
   AP   : E_+ = 328,  r_+(0)=0 (no shell-partner),  max r_+ = 6 at s=13,14,15  (tent profile)
   2·AP : E_+ = 328,  r_+(0)=0,                       max r_+ = 6              (= a dilate of AP, ×2)
   V*   : E_+ = 290,  r_+(0)=1 (shell-partner 3+24),  max r_+ = 5             (flatter)
```
So `2·AP` is the energy-twin of `AP` (forced, being a dilate), while `V*` carries **both** the
unique shell-partner (`r_+(0)=1`, S699/HYP-2262) **and** strictly lower additive energy. The signed
additive face thus carries two gauge-invariant separators of the floor: the shell-partner count
`r_+(0)` (degree-2 Krawtchouk support, Part 2) and the additive energy `E_+` (Part 1's `r_+`,
dilation-graded).

---

## Honest status

- **Proved & verified:** Part 1 (matching cap `r_+(s)≤⌊N/2⌋`, 4000-config check), its corollary
  (no popular-sum escape), Part 2 (a)(b)(c) (exhaustive over all sign patterns, n≤7), Part 3
  (a)(b)(c) (n=7 and the full n=14 floor). Dilation invariance of `E_+` (4000 trials).
- **Standard math, applied (the contribution = the mapping):** tangent addition / formal group,
  cyclotomic conjugation, Walsh/Krawtchouk duality. The *new* content is the dictionary tying the
  signed-LRC additive face to all four, and the **negative resolution** of the popular-pair-sum
  mirror via the matching cap.
- **Not claimed:** a proof of LRC(14). The deliverable is the additive-face structure theory and
  the resolution that the cyclic side has no popular-norm escape (the rank-2 ceiling), sharpening
  THM-412's lattice dichotomy into a metric-by-metric statement.
- **Open (→ HYP-2272):** is `E_+` (with `r_+(0)`) a complete invariant of the worry-set floor at
  every `n`? Does a higher-degree additive form (3-fold sums) realize the cyclic image of the rank
  jump? Is the `V*` energy deficit `328−290=38` a meaningful carry quantity (S677 apex debt)?

**Artifacts:** `04-computation/signed_lrc_multiplicative_energy_s704.py`,
`05-knowledge/results/signed_lrc_multiplicative_energy_s704.out`. Cross-refs: HYP-2262 (signed-LRC
theory), THM-413, THM-412, THM-401/403/407, HYP-1992 (rapidity, distinct axis), HYP-2272.
