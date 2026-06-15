# LRC: the singular series' cancellation lives between support levels, not within them

**Source:** kind-pasteur-2026-06-14-S2 (THM-504). Dispatch: improve the
Paley–Zygmund signed-short-vector route to LRC(14).

The Paley–Zygmund instinct for the LRC singular series `L(S) = (6/7)^13 +
Σ_{exact relations}(6/7)^{13−|T|}(−1)^{|T|}∏s(t_v)` is half right and the half it
gets wrong is the instructive part.

**Right half.** First-moment / absolute-value bounds genuinely fail past the
almost-Sidon class. I made this sharp: the `|T|=3` *absolute* sum `Σ|∏s|`
**diverges** (logarithmically — there are `∼T²` three-term relations with
coefficients up to `T`, each of size `∼1/T³`, giving `∼Σ1/T`). THM-503's almost-Sidon
loose class lives exactly at the `|T|=2` absolute-convergence boundary and cannot
cross it. So one *must* use cancellation; the triangle inequality is not merely
loose at `|T|=3`, it is infinite.

**Wrong half.** The natural guess — that the signed short-vector sum enjoys
*square-root* cancellation (signs quasi-random, `|Σ| ∼ √N`) — is false. The
cancellation is only constant-factor (`|Σ₃|/Σ|∏s| ≈ 0.6`, and that ratio falls
only because the *denominator* diverges, not because the numerator decays). The
reason is a three-line lemma: the band sinc kernel `s(t)=sin(πt/7)/(πt)` is
**positive for all `|t| ≤ 6`** (half a period), so the dominant, small-coefficient
relations all contribute the *same* sign; the negatives need `|t| ≥ 8` and are
`O(1/|t|)`-small. Quasi-randomness of signs is exactly what the arithmetic forbids.

**Where the cancellation actually is.** Not within a support level — *between*
them. `L` is an alternating series in `|T|` (`(−1)^{|T|}`), each level a
conditionally convergent signed sinc-lattice sum, the level masses growing (the
Vitali wall). The saving that keeps `L > 0` is the cross-level alternation
`+C_2 − Σ_3 + Σ_4 − ⋯`. So the right machine is Abel summation across the support
filtration, with a Pólya–Vinogradov/Weil bound supplying each level's convergent
signed value — which is the "archimedean singular-integral positivity" THM-503
named, now located one level down: bound the *convergent* `Σ_k`, then their
alternating sum.

The general lesson, recurring in this repo (the `c5=10` spectral reframe, the
Euler-sign rigidity, THM-460's tower alternation): **when a series is only
conditionally convergent, the cancellation is structural, and you must find which
axis it lives on before reaching for a probabilistic (second-moment) bound.** A
Paley–Zygmund/√-cancellation argument assumes the saving is *within* a level
(independent-ish terms); here it is *across* levels (a deterministic alternation).
Diagnosing the axis — by the half-period positivity that pins the within-level
signs — is the actual content. The probabilistic tool was reaching into the wrong
direction; the arithmetic of the sinc kernel says which way to reach.

Cross-links: THM-504 (this result), THM-501 (the singular series, kp), THM-503
(almost-Sidon, the absolute-convergence boundary, mac-mini), THM-446 (the
relation/Sidon ladder), the resonance-lattice reflection
[[lrc-resonance-lattice-pvsnp-s604]] (the Poisson form and the Vitali wall),
[[efficiency-becomes-proof-kps4]] (reframes have domains; finding the edge is the
insight).
