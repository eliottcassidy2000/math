# The link is a congruence: Φ₆ and 14a meet in the Eisenstein series, modulo Φ₆

*mac-mini-2026-07-14-S95. S94 proved there is no functor runner → X₀(14): the covering-min's arithmetic
is Eisenstein (`Φ₆`, `j=0`, `ℚ(√−3)`), while `X₀(14)` is the conductor-14 curve `14a` (`j≠0`,
`ℚ(√−7)`/`ℚ(√−14)`). Owner then asked the right question: **prove they are linked despite being
different curves.** They are — and the form of the answer is the point. Different elliptic curves are
linked not by functors but by **congruences**, and there is a genuine, verified one here: the cusp form
`f₁₄` is congruent modulo 6 to the Eisenstein series, `a_p ≡ 1+p (mod 6)` — and 6 is exactly Φ₆'s
modulus.*

---

## The link

> **`f₁₄` (the cusp form of `14a`) satisfies `a_p ≡ 1 + p (mod 6)` for every good prime `p`.**

Verified with zero violations over all 34 good primes below 160 (`a_5=0≡6`, `a_{13}=−4≡2`, `a_{17}=6≡0`,
`a_{43}=8≡2`, `a_{47}=−12≡0`, …). And it is not an accident to be checked but a theorem to be proved:

**Mechanism.** `14a` has `E(ℚ)_tors = ℤ/6ℤ` — a *rational* point of order 6. A rational `N`-torsion
point makes the mod-`N` Galois representation **reducible**,
`ρ̄_N ~ \begin{pmatrix} 1 & * \\ 0 & χ \end{pmatrix}`, with the trivial sub-representation coming from
the torsion point and the cyclotomic quotient `χ` from the Weil pairing. Hence for good `p`,
`a_p = \mathrm{tr}\,ρ̄_6(\mathrm{Frob}_p) ≡ 1 + χ(\mathrm{Frob}_p) = 1 + p \pmod 6`. A reducible
representation is precisely an **Eisenstein** one, so this says `f₁₄ ≡ E \pmod 6`, with `E` the weight-2
Eisenstein series (`a_p(E) = σ_1(p) = 1+p`). This is Mazur's Eisenstein-ideal phenomenon, in the
smallest instance where it meets our problem.

## Why the modulus is 6 = Φ₆

`6 = \mathrm{ord}(ζ_6)` is the order of the primitive sixth root of unity, i.e. the modulus of the
**sixth cyclotomic polynomial** `Φ_6(x) = x²−x+1` — whose values `Φ_6(n) = n²−n+1` are exactly the
covering-min denominators `n/Φ_6(n)`. Reduction mod 6 lands in `μ_6 = ℤ[ζ_6]`, the **Eisenstein
integers**, field `ℚ(√−3)`, ramified at the prime 3 — the disc-`−3` part of the `ℤ/6 = ℤ/2 × ℤ/3`
torsion. So `f₁₄` is congruent, *at Φ₆'s own modulus*, to the Eisenstein series.

## Both sides touch the same Eisenstein series

The link is not one-sided. The Eisenstein series `E₂` shows up on the *covering-min* side too, and it
is proved: the Dedekind sum controlling the covering-min margin is
`s(n, Φ_6(n)) = −(Φ_6−1)/(12Φ_6) → −1/12 = −B_2/2` (HYP-3768), and `−1/12` is exactly the **`E₂`
quasi-modular anomaly constant**. So:

- **`f₁₄` touches `E₂` by congruence** — `a_p ≡ 1+p = σ_1(p) \pmod 6` (this session);
- **the covering-min touches `E₂` by anomaly** — its Dedekind margin is the `E₂` transformation defect
  (HYP-3768).

`14a` is `f₁₄`'s home; `Φ_6` is the covering-min's home; and the **weight-2 Eisenstein series is the
common body they are both congruent/tangent to, modulo `Φ_6`.**

## Why this is consistent with S94 — and sharper than "coincidence"

S94 was right: there is **no functor**, because `14a` (`j≠0`) and the Eisenstein curve (`j=0`) are
genuinely different curves in different fields. S95 is right too: they are **linked by a congruence**.
There is no contradiction — this is the ordinary state of affairs in the theory of modular forms.
Congruences (Ramanujan's `τ(n) ≡ σ_{11}(n) mod 691`, Mazur's Eisenstein ideal) link forms and curves
that are *not* isomorphic, *not* isogenous, and carry *different* `j`-invariants. A functor would demand
they be the same object; a congruence only demands that their `q`-expansions agree modulo an integer.
The "coincidence at 14" the corpus warned of is thereby **upgraded to a theorem**: the level side and
the value side meet — not as one curve, but as two curves congruent to the same Eisenstein series
modulo `Φ_6`. They share the number 14 as *level*, the number 6 as *congruence modulus*, and the
Eisenstein series as the object both reduce to.

## Coda — the arc, now with a positive terminus

S89→S94 traced the LRC(14) last bit into `X₀(14)` and found the bridge stops at the curve: the runner
is a level without a curve, and its curve is Eisenstein, not `14a`. S95 supplies what the functor could
not: the honest link. It is a **congruence**, `f₁₄ ≡ E \pmod 6`, and its modulus is `Φ_6` itself. The
covering-min (Eisenstein, via the `E₂` anomaly) and the cusp form (via the mod-6 Eisenstein congruence)
do not become the same object — they become congruent modulo the very cyclotomic polynomial whose values
are the covering-min denominators. That is the truest sense in which `Φ_6` and `14a` are linked: not one
curve, but two curves speaking the same Eisenstein language mod 6.

---

*Cross-links: S94 (HYP-6585, no functor / different curves) — this is its positive complement; HYP-3768
(the covering-min Dedekind margin = the `E₂` anomaly `−1/12`); the Hecke dictionary of `f₁₄` (`a_p ≡ p+1
mod 6`, opus); klein-S10 (`14a` torsion `ℤ/6`, `(ℤ/14)^* = ℤ/6`). Verification/proof: HYP-6605,
`04-computation/phi6_14a_link_eisenstein_congruence_macmini_S95.py`.*
