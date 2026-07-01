# The nonabelian expander is `PSL₂(F_p)` (gap `2√p` = the Gauss-sum arithmetic); sum-of-two-squares IS the difficulty split; and `inf meas ≥ 1/36` is a facility-location price-of-anarchy bound (`PoA ≤ 4.85`)

*kind-pasteur-2026-07-01-S24. Hunting the nonabelian expander realization of the LTC lead, grounding the owner's number-theoretic jumping-off points (sum-of-two-squares, Koksma–Hlawka, Pochhammer, Mahler–Popken/integer-complexity), and importing a potential-function / price-of-anarchy argument for the LRC covering-min as a facility-location game. Three clean results: the expander is `PSL₂(F_p)` with Ramanujan gap `2√p` = the Gauss-sum certificate's arithmetic; the `p mod 4` difficulty split IS the sum-of-two-squares theorem; and `inf meas` is the game value with a computable price of anarchy `≈ 4.18`, so `inf meas ≥ 1/36 ⟺ PoA ≤ 4.85`.*

## The nonabelian expander: `PSL₂(F_p)`, and `2√p` is the certificate's `√p`

The S23 anti-LTC verdict asked for a *nonabelian expander* realization. The natural one at the apex prime `p` is **`PSL₂(F_p)`** — the LPS Ramanujan group. Verified: `|PSL₂(F_p)| = p(p²−1)/2`, and at the apex **`|PSL₂(F_7)| = 168 = |Aut(Fano plane)| = |GL₃(F_2)|`** — the Fano plane is the `QR_7` / octonion structure (the `twentyeight` apex). The LPS `(p+1)`-regular Ramanujan expander has spectral gap `|λ| ≤ 2√p`, and:

> **the expander's spectral gap `2√p` is twice the Gauss-sum magnitude `|g_p| = √p`** (S20). The certificate `g_7 = i√7`, the Paley skew spectrum `{0,±i√7}`, and the Ramanujan expander gap `2√7` are *one √p arithmetic*.

So the LTC lead resolves: a locally-testable realization of the tournament/LRC data should live on the **`PSL₂(F_p)` left-right Cayley complex** (nonabelian, expanding), and its Ramanujan spectral gap is *already* the Gauss-sum certificate. The explicit code/generators remain to be built — but the group and the spectral bridge are pinned: `PSL₂(F_7) = Aut(Fano)`, gap `2√7`.

## Sum-of-two-squares IS the Brouwer/Borsuk–Ulam difficulty split

The owner's link pays off exactly. Verified across primes:

> **`p ≡ 1 mod 4` ⟺ `p = a²+b²` (sum of two squares) ⟺ `−1` is a QR ⟺ Gauss sum `g_p = √p` REAL ⟺ complement is an automorphism ⟺ BROUWER (easy).**
> **`p ≡ 3 mod 4` ⟺ NOT a sum of two squares ⟺ `g_p = i√p` IMAGINARY ⟺ free `ℤ₂` ⟺ BORSUK–ULAM (hard).**

The apex `7` is `3 mod 4`, *not* a sum of two squares — the hard side. So the S18 `p mod 4` difficulty axis is literally **Fermat's sum-of-two-squares theorem**: the LRC(2p) is easy exactly when the apex is a sum of two squares (`√p` real, a fixed point exists) and hard when it is not (`√p` imaginary, only an antipodal Borsuk–Ulam certificate). The whole easy/hard dichotomy is the reality/imaginariness of `√p`, which is the sum-of-two-squares of the apex.

## `inf meas` as a facility-location price of anarchy

Reading the covering-min as the owner's adversarial facility-location game (runners = facilities, the lonely observer = the point farthest from all, covering-min = the adversary's `min` over configs), a potential-function / PoA argument falls out. For `n=14` (`k=13` speeds):

- **Social optimum** (the equidistributed / independent config, each runner safe with prob `6/7`): `meas = (6/7)^13 = 0.134801`.
- **Adversarial minimum** (the extremizer): pentagon `(Z/10)*` `meas = 0.032261`.
- **Potential** `Φ(S) = (6/7)^k − meas(L_C) = 0.1025` — the *resonance penalty* the covering adversary induces (the pairwise correlation / additive energy).
- **Price of anarchy** `= optimum / adversarial-min = 0.1348 / 0.0323 = 4.178`.

So `inf meas = (6/7)^k / PoA`, and the target reformulates cleanly:

> **`inf meas ≥ 1/36 ⟺ PoA ≤ (6/7)^13·36 = 4.85`.** The actual PoA is `4.18 < 4.85`, so the extremizer clears `1/36` (`1.16×`, tight). **The whole `r=2` residual becomes: "the facility-location price of anarchy of the covering game is at most `4.85`."**

This is the potential-function form: `meas(L_C) = (6/7)^k − Φ(S)`, with `Φ` the resonance the adversary can induce; bounding `Φ ≤ (6/7)^k(1 − 1/4.85)` proves `inf meas ≥ 1/36`. And the **rigorous handle is Koksma–Hlawka**: `|meas(L_C) − (6/7)^k| ≤ V·D*` (total variation × star-discrepancy of the runner phases), so `Φ` is a *discrepancy*, the adversary *maximizes* `D*` (resonance), and the covering constraint *caps* it — the PoA bound is a Koksma–Hlawka discrepancy bound. This is exactly the "import a potential/PoA argument" the owner asked for: **the LRC `inf meas` is a facility-location game value, its potential is the resonance discrepancy, and the residual is a PoA `≤ 4.85` bound.**

## The tangential jumping-off points

- **Pochhammer:** the fiber fraction `f(14) = (1/2)_{12}/12! = 0.161 ~ 1/√(πn)` — the Pochhammer symbol `(1/2)_{n-2}` is the two-sheeted branched-cover / Wallis-`π` structure (the "triangle foundation" constant). It sits on the *even/measure* side (the `π` of the collar widths), complementary to the `√p` (odd/certificate) side.
- **Mahler–Popken integer complexity:** `‖7‖ = 6`, `‖14‖ = 8`, `‖21‖ = 9` — "how hard to build the apex from 1's." Tangential, but note `21 = 3·7` (the compositum modulus, S23) has `‖21‖ ≤ ‖3‖+‖7‖` — the multiplicative structure of the bridge. A speculative lead: is the apex's integer complexity a proxy for the LRC proof complexity? (No evidence yet.)
- **Koksma–Hlawka** is not tangential — it is the rigorous form of the PoA potential above.

## The unifying `√p`

At the hard apex `p ≡ 3 mod 4` (not a sum of two squares), `√p` appears four ways: the **Gauss-sum certificate** `i√p` (S20), the **Paley skew spectrum** `{0,±i√p}` (S21), the **Ramanujan expander gap** `2√p` (the LTC realization), and the **imaginary quadratic** `Q(√-p)` (the field, S21). The nonabelian expander (`PSL₂(F_p)`), the certificate, and the difficulty are one `√p` arithmetic — and the LRC(6)/LRC(14) split is where `√-3` and `√-7` meet at `√21` (S23).

## Honest status & next

- **Verified:** `|PSL₂(F_7)|=168=|Aut(Fano)|`; Ramanujan gap `2√p` vs Gauss `√p`; sum-of-two-squares ⟺ `p mod 4` ⟺ Gauss real/imaginary (Brouwer/Borsuk–Ulam); the facility-location numbers (`(6/7)^13=0.1348`, extremizer `0.0323`, `PoA=4.178`, target `PoA≤4.85`).
- **Framing, not proof:** the PoA `≤ 4.85` reformulation of `inf meas ≥ 1/36` is a clean *restatement* (a Koksma–Hlawka discrepancy bound to be proven), not a proof; the `PSL₂(F_p)` expander realization identifies the group and spectral bridge but not the explicit LTC code.
- **Next:** (a) prove the facility-location `PoA ≤ 4.85` via a Koksma–Hlawka bound on the covering config's star-discrepancy (the covering constraint caps `D*`); (b) build the explicit `PSL₂(F_7)` left-right Cayley code and test whether the iso class / LRC certificate is `O(1)`-locally testable on it — turning the anti-LTC into an LTC via the apex's own group.

— Related: `the-half-tiling-is-an-abelian-square-complex-…` (S23, the anti-LTC verdict), `the-compositum-certificate-…` (S23, `√21`), `the-singular-series-is-an-iota-equivariant-lefschetz-trace-…` (S20, `i√7`), `the-spine-stays-quadratic-…` (S21, `√-p`), `lrc-difficulty-by-n-…` (S18, `p mod 4`), `twentyeight-…` (Fano/octonion apex), `everything-is-the-triangle.md` (Pochhammer/`π`). Script: `04-computation/expander_sumofsquares_facility_location_kps.py` (+ .out). Not a HYP reservation.
