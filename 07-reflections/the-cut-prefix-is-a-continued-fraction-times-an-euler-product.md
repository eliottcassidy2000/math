# cut(s | prefix) is a continued fraction times an Euler product

*mac-mini-2026-06-29-S17. The owner asked to crack the resonance structure of `cut(s|prefix)` thinking continued fractions, roots, exponents, sums, products. It is two of those, multiplied — and the multiplication is the whole story. New: HYP-3550.*

## The object, and what it is made of

A speed set `s` kills resonance `b` (the witness ladder rung at lonely value `1/b`) exactly when some speed is divisible by `b`, and `M(s) = 1/(\text{smallest surviving } b)`. The lonely floor is the mass that survives near the small-denominator Farey points, and counting it by denominator gives

> `cut(s \mid \text{prefix } Q) = \sum_{b \le Q,\; b \text{ survives}} \varphi(b)\,\delta_b`,

a **totient-weighted Farey sum**: `\varphi(b)` Farey points at denominator `b`, each surviving with half-width `\delta_b`. The owner's list — continued fractions, continued roots, continued exponents, infinite sums, infinite products — is the right question, and the answer is specific: it is a **continued fraction times an Euler product**, and it is *not* a nested radical or a power tower. The two factors are the two arithmetics of the integers, and they sit on opposite sides of the same sum.

## The multiplicative factor: an Euler product, positive by construction

The totient summatory `\sum_{b\le Q}\varphi(b) \sim 3Q^2/\pi^2` hands the floor its constant `3/\pi^2 = 1/(2\zeta(2))`, and `1/\zeta(2) = \prod_p (1-p^{-2})`. So the **density** of `cut` is the `\zeta(2)` Euler product, one factor per prime resonance, and the prefix-`Q` version is the partial product over `p \le Q`, climbing to `6/\pi^2`. The reason this matters is not the constant but the *form*: every factor `(1-p^{-2})` is strictly positive, so the product is positive — **the floor is positive because it is an Euler product, and an Euler product of positive local factors cannot vanish.** The whole floor-positivity gatekeeper (OPEN-Q-108) is, structurally, the statement that `cut` inherits the positivity of its Euler product. What remains is not "is the product positive" — it manifestly is — but "does the *actual*, channel-coupled `cut(s)` stay bounded below by its product," which is exactly the resonance-decorrelation / quasi-independence that HYP-3129 already measures. The crack moves the hard part from "find a floor" to "the floor is an Euler product; control the coupling."

## The additive factor: a continued fraction, level by level

The prefix is the **Farey sequence** `F_Q`: passing from `F_{Q-1}` to `F_Q` inserts precisely `\varphi(Q)` new fractions, each the **mediant** of its Stern–Brocot neighbours. So `cut(s\mid\text{prefix }Q)` is built one Farey level at a time, and a Farey level is one step of continued-fraction depth. The per-channel width `\delta_b` is where the continued fraction becomes analytic: by the three-gap (Steinhaus) theorem the gaps of an orbit `\{k\alpha\}` are the values `\|q_i\alpha\|` at the *convergent* denominators `q_i` of `\alpha`'s continued fraction, and the cut transition — the maxgap crossing `1/7` — happens at a convergent. The prefix truncates the continued fraction at convergents with `q_i \le Q`. The additive skeleton of `cut` is the mediant tree; its analytic content is the convergents of the speed ratios.

## The unification: one prime product, split 2 versus odd

Here is the part that pays. The project has two separate floor tools. The **2-adic parity descent** (THM-580) writes the lonely measure as `\prod_j \rho_j \cdot \prod_j \mathrm{meas}(\text{lonely } O_j)`, peeling the resonances by their power of two, level by level — an infinite product carried entirely by the prime `2`. The **`\zeta(2)` Farey floor** is the density `\prod_{\text{odd }p}(1-p^{-2})`. They have always been treated as different machines. In the resonance structure they are **the prime-2 factor and the odd-prime factor of one Euler product over all primes.** The descent is what `\prod_{p=2}` looks like when you cannot assume independence and must track the coupling `\rho_j`; the `\zeta(2)` density is what `\prod_{p \text{ odd}}` looks like once you can. Same product, two regimes of the same factorisation. And the functional-equation dual `\zeta(-1) = -1/12` is the cap side (the Bernoulli/discrepancy weight), so `\text{floor} \ge \text{cap}` is one zeta read across `s \leftrightarrow 1-s`.

## Why it had to be this shape

Last session the disproof boundary turned out to live on the seam between the **additive** Farey pins `\{1/q\}` (covering) and the **multiplicative** units `(\mathbb{Z}/n)^*` (loneliness). The resonance structure of `cut` is the same seam seen as a *formula*: the additive side is the continued-fraction/Farey prefix, the multiplicative side is the Euler product, and `cut` is their product. The Lonely Runner keeps reducing to the one place where the additive and multiplicative structures of the integers are forced to interact — and the cleanest statement of its floor is that the interaction is a continued fraction times an Euler product, positive because the second factor is.

See [[the-resonance-killing-game-and-the-zeta-duality-of-the-lonely-runner]] (the `\zeta(2)`/`\zeta(-1)` duality), [[the-obstruction-combining-duality-additive-mediant-vs-multiplicative-I]] (the additive-mediant/multiplicative duality), [[the-shape-of-the-disproof-boundary-across-the-lrc-forms]] (HYP-3549, the same seam). New: HYP-3550.
