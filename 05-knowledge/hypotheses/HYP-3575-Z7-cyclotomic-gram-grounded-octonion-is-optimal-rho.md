---
id: HYP-3575
title: Grounding klein-S7b/the owner's transitivity reframe -- rho_j>=c IS the spectral gap of the Z_7 autocorrelation Gram, which is PSD by Bochner (= an SOS) and SET-INDEPENDENT because the Z_7^* multiplier permutes the Fourier modes (the transitive group-average); the eigenvalues live in Q(cos 2pi/7) (the S75e Fejer-Bochner field). NEW: the apex-7 octonion/Fano set QR={1,2,4} is a PERFECT (planar) DIFFERENCE SET mod 7, so its autocorrelation is FLAT (lambda_k=2 for all k!=0) -- the MAXIMAL spectral gap -- meaning the octonion/Paley structure is the OPTIMAL rho_j (confirming klein-S6's relocation of the octonion from the dead b_1^- route to the descent's Z_7 core, now with the concrete role: octonion = flat autocorrelation = best-case floor). And CV(H)^2~2/n (THM-589, the S_n group-average, even-overlap doubling-2) is the PROVEN finite rehearsal of the same transitive-group-2nd-moment-Fourier-SOS object.
status: SYNTHESIS + verified grounding (Z_7 Gram PSD/Bochner, set-independence by Z_7^* orbit, QR{1,2,4} flat spectrum = perfect difference set, CV(H)~2/n). Extends klein-S7b HYP-3566; does NOT close rho_j>=c (the actual descended-core spectral gap over all configs is the remaining computation).
source: mac-mini-2026-06-29-S27
related:
  - HYP-3566  # klein-S7b: the floor is a transitivity problem; rho_j>=c = Z_7 cyclotomic reference-collapse (this grounds it)
  - HYP-3570  # mac-mini: ESSENTIAL x BOUNDED reframe (this is the BOUNDED half made concrete)
  - THM-589   # CV(H)^2~2/n -- the PROVEN S_n rehearsal
  - HYP-3535  # S75e Fejer-Bochner SOS (= the Z_7 cyclotomic Gram positivity)
  - THM-580   # 2-adic descent (peels to the Z_7 apex)
  - HYP-3547  # the apex-7 octonion/Fano (QR{1,2,4}) -- here it is the OPTIMAL rho_j
  - HYP-3563  # b_1^- (the DEAD octonion route klein-S6 refuted; relocated here)
results:
  - 04-computation/Z7_cyclotomic_gram_SOS_rehearsal_macmini_20260629.py
  - 05-knowledge/results/Z7_cyclotomic_gram_SOS_rehearsal_macmini_20260629.out
---

# HYP-3575 -- the Z_7 cyclotomic Gram, grounded; the octonion is the optimal rho_j

## The transitivity reframe, grounded
The owner/klein-S7b (HYP-3566): a TRANSITIVE group turns a 2nd moment into a single GROUP-AVERAGE whose
positive-definite certificate in the group-Fourier basis is an SOS; on `Z_7` that is the cyclotomic basis
= the S75e Fejer-Bochner SOS; so `rho_j>=c` = positivity of the `Z_7` cyclotomic Gram, set-independent by
transitivity. VERIFIED here on the `Z_7` autocorrelation Gram (`a(d)=#{x in S: x+d in S}`, circulant):
- **PSD = SOS:** every eigenvalue `lambda_k = |sum_{x in S} omega^{kx}|^2 >= 0` (Bochner). The Gram of a
  safe set is automatically a sum of squares; `rho_j>=c` is its SPECTRAL GAP (the smallest mode bounded
  below) -- the Fejer-Bochner minorant.
- **SET-INDEPENDENT:** the `Z_7^*` multiplier `u` sends `lambda_k -> lambda_{ku}`, permuting the nonzero
  modes, so the SPECTRUM is invariant across the orbit (VERIFIED: `{1,2,4},{3,5,6},...` all give the same
  multiset). The gap depends on the `Z_7` STRUCTURE, not the config -- the transitivity payoff.
- **The S75e field:** the eigenvalues are real cyclotomic values in `Q(cos 2pi/7)` (roots of
  `8x^3+4x^2-4x-1`), the totally-real cubic of the apex prime -- exactly the field of S75e's `F_7`.

## NEW: the octonion/Fano set is the OPTIMAL rho_j (a flat spectrum)
The QR set `{1,2,4}` mod 7 -- the Paley arc rule, the Fano line, the octonion triple (HYP-3547) -- is a
**perfect (planar) difference set** mod 7: its six nonzero differences are `{1,2,3,4,5,6}` each once. Hence
its autocorrelation is FLAT (`a(0)=3, a(d)=1` for `d!=0`), giving the FLAT spectrum
> `lambda_0 = 9`, `lambda_k = |S| - 1 = 2` for all `k != 0`  (VERIFIED),
the **MAXIMAL spectral gap** for a 3-element set (a non-difference-set like `{1,2,3}` drops to gap `0.31`).
So **the apex-7 octonion/Fano structure is the OPTIMAL `rho_j`** -- the flattest autocorrelation, the
best-conditioned cyclotomic Gram. This RELOCATES the octonion (klein-S6 refuted it in the `b_1^-` cycle
space; HYP-3566 relocated it to the descent's `Z_7` core) and gives it a CONCRETE role: the octonion is
the perfect-difference-set / flat-autocorrelation structure that MAXIMIZES the floor's spectral gap. The
two 7's (b_1^- coincidence vs apex arithmetic, klein-S6) are now reconciled: the apex octonion is real and
lives HERE, in the `Z_7` Gram, as the optimal case.

## The proven rehearsal: CV(H)~2/n is the same object on S_n
`CV(H)^2 = W(n)/n! - 1 ~ 2/n` (THM-589) is the `S_n` group-averaged 2nd moment, PROVED bounded/clean. Its
positivity certificate is the even-overlap parity (`W(n)`'s per-run `1+(-1)^m` -- the DOUBLING-2 / 2-adic
mechanism, klein-S5), which is the SAME doubling the 2-adic descent (THM-580) uses to peel `14=2.7` down to
the `Z_7` apex. So:
> `rho_j>=c` (LRC floor, `Z_7`-average, cyclotomic SOS, set-independent) and `CV(H)^2~2/n` (PROVEN,
> `S_n`-average, even-overlap SOS) are ONE object -- a transitive-group-averaged 2nd moment certified
> positive in the group-Fourier basis -- on `Z_7` vs `S_n`, linked by the shared doubling-2.
The rehearsal is the theorem; the floor is the same theorem on the descended `Z_7` apex (Heegner `h=1` =
the gentlest cyclotomic = the flattest Gram).

## The remaining computation (the actual close)
What is grounded: the MECHANISM (Z_7 Gram = PSD/SOS, set-independent, octonion-optimal) and the rehearsal
(`S_n`, proved). What remains: show the ACTUAL descended core (the odd part `O_j` at every binding level,
over ALL covering configs) is `Z_7`-cyclic-invariant (klein-S7b's test) so its Gram IS this `Z_7` Gram
with gap `>= c`; where it is NOT `Z_7`-invariant, the `Gamma_0(14)` congruence (HYP-3553) supplies the
missing transitive symmetry. Concrete next: run klein's `Z_7`-invariance test on the binding descended
cores; the octonion flat-spectrum is the target the gap must clear.

## What it buys
Turns the transitivity reframe from a proposal into a grounded mechanism: `rho_j>=c` is a Bochner-PSD
spectral-gap of a set-independent `Z_7` cyclotomic Gram (= S75e), with the octonion/Fano as the optimal
(flat-spectrum) case and `CV(H)~2/n` as the proven `S_n` rehearsal. The floor's open piece is now "the
descended core's Gram clears the octonion-optimal gap," a single spectral question on the apex.
