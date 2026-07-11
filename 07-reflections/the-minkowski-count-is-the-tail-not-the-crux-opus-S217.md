---
source: opus-2026-07-11-S217
status: SHARPENING — the "support-6 Minkowski successive-minima count" bounds only the DISSOCIATED TAIL of
  corr(E); it is provably insufficient (and covolume-blind) for the LOW-HEIGHT CORE that actually carries the
  correction. The core is finite-census + the PROVED far-element peel (THM-546); the genuine open frontier is
  the ungapped accounting (Plat↔Δ entanglement), a SIGNED object, not an absolute lattice count.
tags:
  - lrc14
  - minkowski
  - successive-minima
  - THM-546
  - THM-538
  - covolume
  - signed-cancellation
  - far-element-peel
---

# The Minkowski count is the tail, not the crux

**opus-2026-07-11-S217.** Tasked to attempt the support-6 Minkowski successive-minima count and to play the
angles against each other. Three concept scouts + an exact experiment
(`04-computation/lrc14_covolume_signed_vs_absolute_opus_S217.py`) settle what that count can and cannot do.
The headline is a correction to the target itself.

## The object

`meas(S7(E)) = M7(k) + corr(E)`, `corr(E) = Σ_{0≠n∈Λ°(E)} K(n)`, `Λ°(E) = {n∈ℤ^{k−1}: Σ nⱼeⱼ=0}` (rank
k−2), `K(n) = D7(n mod 7)/∏ⱼnⱼ`, `|Re D7| ≤ 0.1431`. The LRC(14)-S3 crux is `meas(S7(E)) ≤ cap_k`
(consec the maximizer); the binding row is k=8, where `corr(consec_8) = 0.303` against `margin_8 =
cap_8 − M7(8) = 0.357` (slack 0.054). The proposed closer: `|corr| ≤ Σ|K(n)| ≤ c₁⁶/(λ₁⋯λ₆)` on Λ°(E), the
successive-minima product `∏λᵢ ≍ covol(Λ°) = |e/gcd|₂`.

## Two proved walls, and a third my experiment adds

**Wall 1 — sign-blindness (F3, mac-mini Angle F, PROVED-lossy).** `Σ|K(n)|(AP,k=8) = 1.773 ≫ margin_8 =
0.357` (dangerous shapes 4–10×). The smallness of `corr` (0.303, not 1.773) is signed cancellation of
`Re D7` across cosets — a covolume/successive-minima count cannot see it. Any absolute bound loses ≥5×.

**Wall 2 — the absolute tail diverges harmonically (MISTAKE-078).** The free envelope `Σ 0.6973/|n|` over a
box grows without bound (7.42 at radius 10⁵); the box-truncated lattice sum is only *conditionally*
convergent. So even restricting to the tail, an absolute bound does not exist.

**Wall 3 — covolume misses the short relations (this session, EXACT).**

| family | meas(S7) | corr | covol | corr·covol |
|---|---|---|---|---|
| consec {0..7} | 0.3272 | **0.3027** | 11.83 | 3.58 |
| AP dilated 5· | 0.3272 | 0.3027 | 11.83 | 3.58 |
| stranger {0..6,**40**} | 0.1905 | 0.1660 | 41.1 | 6.83 |
| stranger {0..6,**400**} | 0.1969 | 0.1724 | 400.1 | **68.99** |
| wide dissociated | 0.0261 | 0.0016 | 243.2 | 0.39 |

`corr·covol` is **not bounded**: the stranger's `corr` stays ≈ 0.17 whether the far element is 40 or 400,
while covol scales with it. The reason is structural — the correction is carried by the **short** relation
`1+2−3 = 0` inside the `{1,2,3}` core, which is present at `λ₁ = √3` no matter how far the stranger sits.
Covolume (a global `∏λᵢ` quantity) is dominated by the far element and is blind to that short relation. So
**even setting signs aside, covol/successive-minima is the wrong ruler for the part of `corr` that matters.**

Corollary worth keeping: dilated APs have *identical* covol and corr (covol `= |e/gcd|₂` divides the
dilation out) — **scale-invariance (THM-531) is a covolume identity**, not a separate mechanism.

## The correct decomposition (the sharpening)

`corr(E)` splits by relation height, and the two pieces need opposite tools:

1. **Low-height core** (short relations, bounded `|n|` — the sub-AP content). This *carries* `corr`; consec
   maximizes it because consec is the densest low-height relation lattice (F2: AP minimizes covol ⇔
   maximizes short-relation count). It is **finite**: for a bounded core it is the exhaustive census (consec
   the verified argmax, 0 exceedances); far elements are stripped off by the **PROVED far-element peel**
   (THM-546: `|Δ_w| ≤ (6/49)·V(E')/w`, the multi-D sum collapsing to a *1-D Abel sum against the mod-7
   character*). This is the ONE place the signed structure was successfully captured — and it is rank-1.
2. **High-height tail** (dissociated, all `|n|` large). Here `corr → 0` (wide dissociated: 0.0016). This is
   where covolume/successive-minima genuinely applies and where LEM-022 lifts (its rank-1 separation count
   `N(K,M) ≤ 1+4KM/P(w)` is the two-coordinate case of the lattice-point count). But by Walls 1–2 the tail is
   still only *conditionally* convergent, so even here the count must be **signed**, not absolute.

**Net: the Minkowski successive-minima count is the tail, not the crux.** The crux (the low-height core) is
finite-census + THM-546; the tail is where the count lives, and there it must carry the `D7` sign.

## What this says the last gap actually is

Not "execute the absolute `c₁⁶/(λ₁⋯λ₆)` count" — that is proved insufficient on three independent grounds.
The remaining frontier is the **ungapped accounting**: wide families with no single peelable far element and
not a dilated AP (higher-Freiman-dimension GAPs). THM-546 closes the *gapped* regime (a peelable ratio ≳ 4·
core-span); dilated APs are THM-531; bounded spread is the census. The ungapped-wide-non-AP residual is held
by the computational **Plat↔Δ entanglement** (a wide base shrinks `Φ = p0 + p1/7`, compensating the larger
`Δ_w`; margin ≥ 0.22), which is **not yet a clean theorem**. That entanglement — a SIGNED, measure-side
statement — is the honest open object, and it is the multi-dimensional generalization of THM-546's 1-D Abel
peel, not of the absolute covolume count.

**Actionable:** the productive lift of LEM-022 is THM-546's *signed 1-D peel* generalized to a multi-D Abel
summation-by-parts against `D7` over Λ°(E) (the THM-504 prescription: Abel across the support filtration +
Pólya–Vinogradov on the mod-7 character). The successive-minima product enters only as the convergence/
sparsity input to that signed sum. The absolute Minkowski count should be retired as a *closer* (kept only as
a density input), and the ungapped Plat↔Δ entanglement promoted to the named open lemma.

→ THM-538 (identity/kernel), THM-546 (the proved signed peel — the template), THM-531 (dilation = covol
identity), THM-504 (Abel-across-support prescription), HYP-2607 (consec maximizes L_y — the dual of the core
extremality), LEM-022 (rank-1 tail case), the-signed-cancellation-is-support-6-not-t3-opus-S216 (prior).
