# The duality web at apex 7: even/odd, ±, and add/mult are three axes of one structure — and max-H picks the additive pole while LRC's obstruction is the multiplicative one

*opus-2026-06-29. Owner: do the homology, then explore even/odd duality and how it connects to ±,
add/mult, etc. — extend. Merging mac-mini's three-pillars (HYP-3547), the add/mult gauge stack
(HYP-1984), klein's R-even/R-odd metagraph (THM-587), and my Fourier/free-monoid threads.*

## The homology piece (HYP-3544, first concrete result)
The metagraph 1-skeleton's **R-odd Betti numbers are nonzero**: `b₁⁻ = 1, 7, 119` (n=4,5,6), via Hopf
trace `tr(R|H₁)=tr(R|C₁)−SC+1` (`SC=#R-fixed classes`, matching `P_n(−1)`). So the R-odd (Borsuk–Ulam)
obstruction is *topologically real*, not just a spectral block. (`b₁⁻ ≠ NS`, so it's a genuine
higher-degree refinement; the full level-graded equivariant homology is klein's open deliverable.)

## Three dualities, one structure, hinged on apex 7 (verified)
Apex `7` is one of only **two primes (3, 7) that are simultaneously Mersenne, Heegner, and 3-mod-4**
(verified: 31,127,8191 are Mersenne but `h(Q(√−M))>1`). Its three arithmetic properties are three
*orthogonal* dualities, and they are mac-mini's three proof pillars:

| axis | duality | apex-7 property | pillar | tournament/LRC face |
|---|---|---|---|---|
| **EVEN/ODD** | `N = 2^h·odd_core` (HYP-1984) | **Mersenne** `7=2³−1` | 2-adic descent (THM-580) | `R`-even/`R`-odd metagraph; `14=2·7` peels to the apex-7 face |
| **POSITIVE/NEGATIVE** | real ⊕ imaginary; SOS vs obstruction | **3-mod-4** `χ(−1)=−1` | Borsuk–Ulam (THM-581) | Gauss sum `i√7` (imag/neg, `R`-odd); cap's negative R-odd residual |
| **ADD/MULT** | additive char ⊕ multiplicative char | **Heegner** `Q(√−7)`, h=1 | Fejér–Bochner SOS (HYP-3535) | character (mult) vs interval (add) Fourier |

And **Mersenne ∩ Heegner ∩ 3mod4 = {3,7} ⟹ LRC(6), LRC(14)** are the only doubly-even apex cases —
`28=T(7)=C(8,2)` is the 2nd even perfect number = the arc-count of the 8-tournament (Euclid–Euler:
perfect = `T(Mersenne)` = arc-hypercube dim, HYP-3546). The three axes *intersect* exactly at the LRC
frontier.

## My extension 1: H lives at the ADD/MULT interface (and so does the maximizer)
The two canonical structures on tournaments are dual operations:
- **MULTIPLICATIVE** (condensation / free monoid): `H(X⊕Y)=H(X)·H(Y)`, so `H = ∏ H(SC primes)`.
- **ADDITIVE** (OCF / odd cycles): `H = 1 + 2Σ_k α_k 2^k = I(Ω,2)`, a *sum* over disjoint cycle families.

> **`H` is simultaneously a product over the multiplicative primes and a sum over the additive cycle
> families.** This is *why* `H` resists a clean formula — the two structures don't align (the #P-hard
> residual). The **H-maximizer is multiplicatively IRREDUCIBLE (a single SC prime) but additively
> MAXIMAL (most cycle families)** — it wants additive richness packaged in a multiplicative atom.

## My extension 2: max-H picks the ADDITIVE pole; LRC's obstruction is the MULTIPLICATIVE pole
My Fourier duality gives the add/mult axis concretely:
- **ADDITIVE pole** = the **interval/Dirichlet kernel** = the **half-turn** tournament (max-H optimum,
  `R`-even, real, provable). The connection set `{1,…,(n−1)/2}` is an *additive* interval.
- **MULTIPLICATIVE pole** = the **character/Gauss sum** = **Paley** (`QR = {1,2,4}` mod 7 = the Fano /
  octonion *multiplication* rule, HYP-3547). It is `R`-odd, imaginary `i√7`, the **LRC obstruction**.

> So **max-H optimizes toward the *additive* (interval/Dirichlet) pole, while LRC's hardness lives on
> the *multiplicative* (character/Gauss/octonion) pole.** The maximizer is the additive one *precisely
> because* the multiplicative (Paley/octonion) one is too "pseudorandom" (flat spectrum) — and that same
> multiplicativity is the LRC sign-obstruction. The two famous problems sit on opposite poles of the
> add/mult axis, both hinged by `ε = χ(−1)` (the even/odd reversal sign = the Gauss-sum reality).

## The unified picture
`ε = χ(−1) = (−1)^{(p−1)/2} = (−1)^{C(p,2)}` is the **single hinge** of all three axes: it is the
2-adic/Mersenne parity face (even/odd), the Gauss-sum reality (±, real vs `i√p`), and the
character-at-`−1` (add/mult, whether the multiplicative character sees the additive `−1`). At apex 7
all three say "odd / imaginary / 3-mod-4," which is why LRC(14) is the hard-but-poised case: the
additive pillar (Dirichlet, SOS-provable) handles the bulk, but the multiplicative pillar (Gauss,
`i√7`, the octonion sectors) carries the irreducible `R`-odd obstruction — the same `b₁⁻>0` the metagraph
homology now exhibits.

## Status
- **Verified:** `b₁⁻ = 1,7,119` (metagraph R-odd homology); Mersenne∩Heegner∩3mod4 = {3,7}; perfect =
  `T(Mersenne)` = arc-count.
- **New (opus):** `H` at the add/mult interface (product-of-primes ⊕ sum-of-cycles); max-H = additive
  pole, LRC obstruction = multiplicative pole; `ε` as the common hinge of even/odd, ±, add/mult.
- **Open (shared):** the full R-odd equivariant homology = the LRC odd index (klein HYP-3544).

Related: mac-mini HYP-3547 (apex-7 three pillars), HYP-3546 (perfect=arc-count), HYP-1984 (add/mult gauge),
klein THM-587 (R-even/R-odd, signed cycle index), opus character-vs-interval + free-monoid + max-H
reflections, THM-580/581 (descent/BU pillars), OPEN-Q-108.
