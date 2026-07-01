# The facility-location PoA bound is a below-mean tail — so it can't be an unsigned discrepancy bound; it *reduces to* the signed Gauss certificate. Plus the corrected PoA (11-cores) and the 21-Frobenius inside `PSL₂(F₇)`

*kind-pasteur-2026-07-01-S25. Trying to prove the facility-location price-of-anarchy bound `PoA ≤ 4.85` (= `inf meas ≥ 1/36`) via Koksma–Hlawka / a potential function. The attempt is instructive: it corrects the S24 numbers (the residual is over 11-cores, so `PoA = 5.69`, target `6.61`), and it shows the bound is a **below-mean tail** — so Markov / second-moment / any *unsigned* estimate gives the wrong direction, and the bound necessarily **reduces to the signed resonance sum = the ι-odd Gauss certificate (S20)**. The game-theoretic and arithmetic routes are one bound. And the 21-Frobenius `Aut(Paley₇)` sits inside the LTC expander `PSL₂(F₇)`.*

## Correction: the residual is over 11-cores (`PoA = 5.69`, target `6.61`)

S24 computed the PoA with `k=13` (the whole speed set). But the `r=2` multi-far residual is over **11-cores** (the 13-r speeds). Recomputing for the pentagon `(Z/10)*` 11-core `{1..13}\{6,10}`:

- `meas(L_C) = 0.032261`, independent value `(6/7)^11 = 0.183479`.
- `meas = (6/7)^11·(1+R)` with **resonance correction `R = −0.8242`**.
- `inf meas ≥ 1/36 ⟺ 1+R ≥ (1/36)/(6/7)^11 = 0.1514 ⟺ R ≥ −0.8486`. Actual `R = −0.8242 ≥ −0.8486` — clears by `0.024` (tight `1.16×`).
- **`PoA = (6/7)^11/meas = 5.69`, target `PoA ≤ (6/7)^11·36 = 6.61`.**

So the clean reformulation is `inf meas ≥ 1/36 ⟺ PoA ≤ 6.61` (over 11-cores), and the covering adversary's actual PoA is `5.69`.

## Why the discrepancy / potential route fails: it is a below-mean tail

The congestion `C(t) = #{runners dangerous at t}` has mean `E[C] = k·2r = 11·(1/7) = 1.571 > 1`. So **`meas(L_C) = meas{C=0}` is a *below-mean* event.** For below-mean tails:

- Markov, Paley–Zygmund, and the second-moment method give `P(C>0) ≥ E[C]²/E[C²]` — an **upper** bound on `meas(L_C)`, the *wrong direction*.
- The danger-Gram `G_{vw}=meas{both dangerous}` (mac-mini's Bridge-2 2nd moment) has trace `= E[C] = 1.571` and top eigenvalues `[0.42, 0.18, 0.16, …]`; its spectrum bounds `|R|` in **magnitude** but not the **signed** value.

So **no unsigned / absolute estimate (Koksma–Hlawka in magnitude, second moment, Gram spectrum) can lower-bound `meas(L_C)`.** The naive "cap the star-discrepancy" route is blocked because the quantity being bounded is a tail below the mean, where discrepancy magnitude is the wrong tool.

## The reduction: the PoA bound *is* the signed Gauss certificate

The resonance correction is `R = Σ_{nonzero resonances Σ m_v v = 0} ∏_v [ĝ(m_v)/(6/7)]`, and **`Σ|·|` diverges** (MISTAKE-078): `R` converges only through **signed cancellation**. So bounding `R ≥ −0.849` requires the *signed* structure — which is exactly the ι-odd Lefschetz / Gauss certificate `g_7 = i√7` (S20) together with the tight-locus finiteness `{AP, Goddyn–Wong}` (THM-523).

> **The facility-location PoA bound and the arithmetic Gauss certificate are the same bound, seen two ways.** The game-theoretic reformulation (`PoA ≤ 6.61`) does not bypass the arithmetic; it *is* the arithmetic — the signed resonance sum whose cancellation is the certificate.

So the honest verdict on the owner's proof request: **the tight `1/36` bound is not a discrepancy/potential estimate** (that route is blocked by the below-mean tail). The two actual routes are:
1. the **finite near-tight census** — THM-522 quantization bounds the search to primitive scale-1 11-cores; enumerate the near-tight ones (dilated APs + two-clash) and verify each `≥ 1/36` (my S8/S9 census, tight-locus dichotomy THM-523);
2. **singular-series positivity** (THM-501) restricted to the near-tight locus — the *signed* certificate route (S20–S24), whose value is `i√7` / `√21`.

The PoA `≤ 6.61` is the clean *restatement* that unifies them, not a shortcut past them.

## The 21-Frobenius inside `PSL₂(F₇)`, and the π/√p split of the potential

The **21-Frobenius `F₂₁ = C₇ ⋊ C₃ = Aut(Paley₇)`** (order `21 = 3·7` = the compositum modulus, S23 `√21`) sits **inside the LTC expander `PSL₂(F₇)`** (`168 = 8·21`): the apex Frobenius is a subgroup of the apex expander. Its `C₇` is the runners' cyclic loop; its `C₃` is the QR/tournament orientation; and `21 = |F₂₁|` is a *forbidden* Hamiltonian-path count (THM-079). So `21` is at once the Frobenius order, the compositum modulus `Q(√21)`, and a forbidden `H` — the same `3·7` in three roles.

The potential decomposition is the owner's two faces made one identity:

> **`meas(L_C) = (6/7)^k · (1 + R)`** = (the **π / even / measure** side: the independent value `(6/7)^k`, whose fiber constant is the Pochhammer `f(n)=(1/2)_{n-2}/(n-2)! ~ 1/√(πn)`) × (the **√p / odd / certificate** side: the signed resonance `1+R`, whose cancellation is the Gauss sum `i√p`).

The independent term is `π`-arithmetic (Wallis/Pochhammer, the collar widths); the correction is `√p`-arithmetic (Gauss/Ramanujan, the resonance). The facility-location game splits exactly along the even/odd (measure/certificate) axis of the whole program.

## Honest status & next

- **Computed:** corrected PoA `5.69` (target `6.61`); `R=−0.8242`, target `R≥−0.8486`; `E[C]=1.571>1` (below-mean); danger-Gram trace/spectrum; `F₂₁ < PSL₂(F₇)`.
- **Verdict:** the discrepancy/potential route is *blocked* (below-mean tail; unsigned bounds go the wrong way); the PoA bound **reduces to the signed Gauss certificate** — game theory and arithmetic are one. The real proof routes are the finite census (THM-522/523) and singular-series positivity (THM-501).
- **Next:** run the **finite near-tight 11-core census to completion** (the one route not blocked): enumerate primitive scale-1 11-cores near the tight locus, verify each `meas ≥ 1/36` — a finite computation that would actually close the `r=2` residual, where the discrepancy bound cannot.

— Related: `the-expander-is-PSL2-Fp-…` (S24, the game + expander), `the-singular-series-is-an-iota-equivariant-lefschetz-trace-…` (S20, `i√7`), `the-compositum-certificate-…` (S23, `√21`), `the-finish-is-a-recursive-tight-singular-series-…` + `moment-relaxation-reduces-multifar-…` (the census + the r=2 reduction), THM-501/522/523, MISTAKE-078 (the divergent unsigned sum), HYP-3789 (Lasserre), `everything-is-the-triangle` (Pochhammer/π). Script: `04-computation/poa_bound_attempt_signed_resonance_kps.py` (+ .out). Not a HYP reservation.
