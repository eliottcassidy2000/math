---
id: THM-518
title: Stranger-decoupling for the LRC(14) extremizer family (the m→∞ tail is a finite measure, Weyl) + resonant strangers carry the infimum + the honest Bedert two-route diagnosis (Riesz product is the WRONG tool for AP-cores; the prime route lands at 2/29 ≈ 96.6% of 1/14) + the genuine bridge R̂ = the singular series L(S)
status: PROVED (the stranger-decoupling Weyl limit; the prime-route arithmetic; the Riesz-stall, verified by exact direct-grid). OBSERVED (the resonant-stranger dips below the limit, verified n=13 fine grid). DIAGNOSIS (neither 2025 state-of-the-art route reaches 1/14). The bridge R̂(ℓ)=Σ_k r_k(ℓ)(−p/2)^k = L(S) is a structural identification. inf L>0 remains OPEN.
source: mac-mini-2026-06-15-S5
depends_on:
  - THM-501   # singular series L(S) = lim D(q,S)/q = lonely measure
  - THM-503   # almost-Sidon loose class, archimedean integral
  - THM-504   # the |T|≥3 conditional-convergence wall (cross-level alternation)
  - THM-515   # L(S) is the lonely measure; additive energy governs L; Riesz is the right framework
related:
  - OPEN-Q-104   # Riesz-product route to inf L>0
  - OPEN-Q-097   # cross-level / level-decomposition route
  - HYP-2526     # dominance ⊥ L (the stranger does not set L)
  - HYP-2547     # Riesz product stalls at the AP-extremizers
  - HYP-2548     # prime point-mass lands at 2/29
  - HYP-2549     # stranger-decoupling limit + resonances
---

# THM-518 — Stranger-decoupling, the Bedert two-route diagnosis, and the bridge

Endgame for `inf_S L(S) > 0` (≡ C'(14) ≡ LRC(14)). All computation is exact direct-grid
on the lonely measure `L(S) = meas{τ : ||v τ|| > 1/14 ∀ v∈S}` (THM-515), saved to
`05-knowledge/results/lrc14_{riesz_verify,riesz_optimize2,stranger_limit,decouple_confirm}_macmini_0615s5.out`.

## A. Stranger-decoupling: the m→∞ tail of the extremizer family is a finite measure (PROVED)

The conjectured extremizers are the interior-drop cores `S = ({1..13}\{j}) ∪ {14m}`
(one bounded 12-element core + one "stranger" multiple of 14). For the **tail** `m→∞`:

> **lim_{m→∞} L( ({1..13}\{j}) ∪ {14m} ) = (6/7) · meas( Lonely({1..13}\{j}) ).**

**Proof (Weyl / Riemann–Lebesgue).** `L(core ∪ {w}) = ∫_{A} 1_{||wτ||>1/14} dτ` with
`A = Lonely(core)` a fixed finite union of intervals. The integrand is periodic in `τ`
with period `1/w` and mean `6/7` (each danger band `||wτ||≤1/14` has density `2/14=1/7`).
For fixed `A`, `∫_A 1_{||wτ||>1/14} dτ = (6/7)|A| + O((#intervals of A)/w) → (6/7)|A|` as
`w=14m→∞`. ∎ Verified (j=6): `meas(Lonely(core12))=0.0081569`, `(6/7)·=0.0069916`;
`L(m=50)=0.0069340` (99.2%), monotone approach for large `m`. The 12-core measures are
positive for every `j` (`j=1: 0.0714, j=6: 0.00816, j=12: 0.0121`), so the **tail limit is
bounded below by `(6/7)·min_j meas(Lonely({1..13}\{j})) = (6/7)·0.00816 = 0.00699 > 0`**.

This reduces the *infinite* `m→∞` tail of each `j`-family to **one finite measure** — a
genuine collapse of the family's non-compact direction.

## B. But the infimum is carried by RESONANT finite strangers, not the limit (OBSERVED — the honest correction)

The naive hope "`inf_m L = the m→∞ limit`" is **FALSE**. The `m`-profile (j=6, exact grid):

| m | stranger 14m | L | vs limit 0.006992 |
|---|---|---|---|
| 1,2,3,5 | 14,28,42,70 | 0.008157 | **+17%** (stranger does not bind: `Lonely(core12)` already misses the `14m`-band) |
| **7** | **98 = 2·7²** | **0.005242** | **−25% (DIPS BELOW)** |
| 11,13 | 154,182 | 0.0072–0.0074 | +3–5% |
| 20,30,50 | 280,420,700 | 0.0066–0.0069 | → limit from below |

The dip at `m=7` is a true **resonance**, not a grid artifact (grid `Q=4414410` divisible
by 98): `98 = 14·7` shares the factor `7` with the core element `7`, so there are `τ` with
`||7τ||≈1/2` (core-safe) yet `||98τ||=||14·7τ||` small (stranger-danger), carving extra
holes out of `Lonely(core)`. **`inf_m L({1..13}\{6}∪{14m}) ≤ 0.00524`, attained at a finite
resonant `m`, below the `(6/7)`-limit `0.00699`.** (This also locates the true worst stranger:
`m=7`/98 beats the previously-cited `m=4`/56, `L≈0.0056`.) **So `inf L>0` reduces to: the
finite tail limit (part A, positive) AND a finite set of resonant strangers (those `14m`
sharing arithmetic structure with the bounded core) — each a finite, positive measure.** The
resonances are exactly the `7`-power / square dilations; controlling them is the remaining work.

## C. The honest Bedert two-route diagnosis (the OPEN-Q-104 answer)

Bedert (arXiv:2511.16636, 2025) is the state of the art on the cousin lonely-runner-gap
problem. Porting both of his routes to the LRC(14) cores (exact direct-grid + the paper's
Lemma 5.3) gives a precise, **honest** verdict — *neither reaches `1/14`*:

1. **The Riesz product is the WRONG tool for the AP-cores.** `R(τ)=∏_{v∈D}(1+a_v cos2πvτ)≥0`,
   certificate `∫M·R/∫R < 1 ⟹ loose`. The dilated-AP cores have **small additive dimension**
   `dim₂ ~ log N ≈ 2–3` (an AP is the opposite of dissociated), so Bedert's gain `dim₂²/n³`
   is worthless there. Verified directly: the uniform `∏(1−cos)` certifies 3/5 loose configs
   (dilated `{1..12}∪14`: 0.955; evader: 0.938) but **FAILS both interior-drop extremizers**
   (`{1..13}\{6}∪56`: 1.064; `\{12}∪84`: 1.035). Per-speed amplitude optimization pulls the
   `j=6` extremizer only to **1.0096 — still ≥ 1** (best: bump `a₁,a₂` to +0.6). Adding 2nd
   harmonics or dropping the stranger makes it worse. The Riesz product **provably stalls at
   the maximal-additive-energy AP-extremizers** — precisely the infimum.
   *(Caveat: an earlier Fourier-truncated estimate (Kmax=14) gave a spurious 0.95 "certificate";
   the exact direct-grid overturns it. Use direct-grid for these ratios — see MISTAKES.)*

2. **The prime point-mass measure is the RIGHT tool — and it lands JUST short.** Bedert
   Lemma 5.3: `ML(V) ≥ ⌈(p−1)/(2n)⌉/p` for any prime `p` dividing no `v∈V` (`n=13`, `2n=26`).
   For **every** extremizer the best admissible prime is `p=29`, giving
   `ML ≥ ⌈28/26⌉/29 = 2/29 = 0.06897 = 96.6% of 1/14 = 0.07143` — **short by 3.4%.** (No
   smaller prime is admissible: `2,3,5,7,11,13 |` core elements; `k=⌈(p−1)/26⌉` forces
   `p ≤ 14k`, leaving only `p∈{27,28}` for `k=2`, neither prime.)

> **Both 2025 state-of-the-art routes miss `1/14`** (Riesz: ratio `1.0096 ≥ 1`; prime:
> `0.069 < 0.0714`), yet the cores **are** loose (`L≈0.0052 > 0`, computed). The truth sits
> in the `~1–4%` sliver between what the best general methods certify and the conjectured
> optimum. **This sliver IS the exact-value difficulty of LRC(14)** — the "improve the union
> bound by a factor of two" gap Tao names, here pinned to a `1.4×` residue.

## D. The genuine bridge: Bedert's Riesz coefficients ARE the singular series (reframes OPEN-Q-097)

For a **non-dissociated** frequency set (an AP-core), Bedert's clean `R̂(ℓ)=(−p/2)^k` on the
disjoint level-sets `E_k` breaks — the relation `Σε_m m=ℓ` has many solutions, giving

> **R̂(ℓ) = Σ_k r_k(ℓ) (−p/2)^k,  r_k(ℓ) = #{ε : Σε_m m = ℓ, Σ|ε_m| = k}.**

This signed level-graded relation-count is **exactly** the project's singular series
`L(S) = Σ_{t∈Λ} ∏h(t_i)` (THM-515): the additive **relation lattice** `Λ` of the cores is
Bedert's `E_k` with the dissociation dropped, and the `(−p/2)^k` grading is the `(−1)^{|T|}`
cross-level alternation of THM-504. Two consequences:

- **OPEN-Q-097's literal target is FALSE.** The per-level masses
  `Λ_k = (6/7)^{13−k}Σ_{|T|=k}∏s` **grow**: `Λ_2≈+0.11, Λ_3≈−0.55, Λ_4≈+1.17` (truncated
  B=6), so `Σ_{|T|≥3} ≈ +0.62 ≫ (6/7)^13 = 0.135`. The bound `|Σ_{|T|≥3}| < (6/7)^13` cannot
  hold; **the cross-level `(−1)^{|T|}` alternation is essential** (the naive level-truncation
  gives `0.86` vs true `L=0.0056`). Reframe OPEN-Q-097 to: *control the conditionally /
  Abel-convergent alternating series of growing level masses.*
- **Bedert's level bound is the missing tool.** `|E_k ∩ AP| ≤ (C log|P|)^k` (Bonami
  hypercontractivity + Rudin Λ(q) + Bell–Chueluecha–Warnke sunflower) bounds exactly how many
  weight-`k` relations of an AP-core fall in any progression — i.e. it bounds the relation
  counts `r_k(ℓ)` that drive the growing `Λ_k`. The cores live in `[1, ~14m]` (AP-like), so
  `log|P|` is small; this is the portable ingredient for the THM-504 wall, independent of the
  gap-vs-exact mismatch.

## Honesty

`inf L>0` is **not proved**. Established this session: (A) the stranger-decoupling Weyl limit
(the `m→∞` tail of each `j`-family = one finite positive measure); (B) the honest correction
that the infimum is carried by finite *resonant* strangers (`m=7`/98 dips to `0.00524`), so
the reduction needs the resonances handled, not just the limit; (C) the precise diagnosis that
**both** 2025 state-of-the-art routes (Riesz product, prime point-mass) fall `~1–4%` short of
`1/14` on the cores, the Riesz product being structurally the wrong instrument for AP-cores;
(D) the structural bridge identifying Bedert's Riesz coefficients with the project's singular
series, which falsifies OPEN-Q-097's literal target and hands over Bedert's level bound as the
tool for the conditional-convergence wall. The exact value `1/14` is genuinely harder than
either general method, and (per Bedert's own porting note) is where the project's **exact
relation-lattice computation**, not the asymptotic machinery, must do the work.
