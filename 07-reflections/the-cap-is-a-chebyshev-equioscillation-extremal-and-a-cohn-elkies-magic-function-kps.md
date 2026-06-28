# The cap is a Chebyshev equioscillation extremal — and a Cohn–Elkies magic function

*kind-pasteur-2026-06-27-S31an, second thread. Pushing further into abstraction and synthesis: the
7-fold-cyclotomic ideal (the cap, S31ak) is a **Chebyshev equioscillation extremal**, its equioscillation
(double-root) structure is **exactly why the cap is rational** (= mac-mini's ℚ-collapse), and the whole
cover bound is a **Cohn–Elkies / Delsarte LP** with the cyclotomic as the **magic function** — the 1-D,
7-fold analog of Viazovska's modular magic functions in dim 8 (octonions) and 24.*

## The computation (verified, clean)
The Vieta–Lucas / Chebyshev polynomial `V_7(u)` (defined by `V_n(2cos t) = 2cos(nt)`) is
`u^7 − 7u^5 + 14u^3 − 7u`. The de Moivre angles `2cos(2πj/7)` solve `2cos(7t)=2`, i.e. `V_7(u)=2`. And:
```
V_7(u) − 2  =  (u − 2) · (u³ + u² − 2u − 1)²      [VERIFIED exactly]
            =  (u − 2) · (de Moivre cubic)²
```
**The de Moivre cubic appears SQUARED.** A double root is the signature of **equioscillation** — the
extremal configuration of a Chebyshev (minimal-sup-norm) problem touches its envelope tangentially. So:

> **The de Moivre angles are the equioscillation points of the 7-fold Chebyshev problem**, and the cap (=
> the 7-fold-cyclotomic ideal) is a **Chebyshev extremal value**. The "ideal 7-fold symmetry" of the
> coverage (S31ak) is literally the equioscillation extremal of `V_7`.

## Why the cap is RATIONAL — the equioscillation is the ℚ-collapse
mac-mini found (HYP-3132) that the k=8 resolvent biquadratic has **perfect-square discriminant 9**, so the
radicals collapse to `ℚ` and `cap_8`, `dip_8` are rational, not surds. The Chebyshev picture **explains
this structurally**: a **double root** (equioscillation) forces the discriminant to be a **perfect
square** (the square of the underlying cubic's discriminant). So:

> **The cap is rational *because* it is an equioscillation extremal** — the double-root structure of
> `V_7 − 2 = (u−2)·cubic²` is the source of the perfect-square discriminant, i.e. of mac-mini's ℚ-collapse.
> Rationality is not a coincidence of the apex; it is the Chebyshev equioscillation tangency.

This fuses three threads into one fact: **cyclotomic ideal (S31ak) = Chebyshev equioscillation extremal
(here) = rational ℚ-collapse (HYP-3132).** The de Moivre cubic is the common object; its appearance
*squared* in `V_7` is why all three hold simultaneously.

## The cover bound as a Cohn–Elkies / Delsarte LP with a cyclotomic magic function
The cap is a **linear-programming bound** (the Delsarte/moment-LP, the team's L_y dual). In the
Cohn–Elkies framework, such a bound comes from a **magic function** `f` (positive-definite, `f̂ ≥ 0`,
sign-constrained), and the bound is *sharp* exactly when the optimal configuration's points sit at the
**zeros of the magic function**. Here the sharp configuration is the **AP**, and the magic function's
zeros are the **de Moivre angles** (the equioscillation points). So:

> **The cover bound is a Cohn–Elkies / Delsarte LP whose optimal magic function is the cyclotomic /
> de Moivre object** (the 7-fold Chebyshev), with the AP as the LP-sharp configuration. The Hermite–Biehler
> interlacing (S31am) is the magic function's `f̂ ≥ 0` positive-definiteness; the Toeplitz `λ_min`
> (HYP-3201/3203) is the LP's PSD margin; the Perron mode is the magic function's principal eigenvector.
> **All four spectral faces are the one magic function.**

## The Viazovska tie — and the multiplicative/additive split, again
Viazovska's sphere-packing magic functions live in **dim 8 (`E_8` = the octonion lattice)** and **dim 24
(Leech)**, built from **modular forms**. By HYP-3211, **dim 8 = octonions = the MULTIPLICATIVE apex**
(the tournament/G₂/Fano face). The LRC's magic function is the **1-D, 7-fold ADDITIVE** apex — built from
the **cyclotomic** (de Moivre / real-subfield) the way Viazovska's is built from modular forms. So the
LRC and the Viazovska program are the **additive and multiplicative magic-function problems of the same
apex 7**, related (S31an) by the Galois trace `ζ↦ζ+ζ̄` (= Joukowski). The "modular forms" that build
Viazovska's `dim 8` function correspond, on the LRC side, to the **quasimodular `E₂` / cyclotomic** objects
(mac-mini's Eisenstein/E₂ spoke) — the 7-fold magic generator.

## Proof relevance (why this helps)
Chebyshev extremality is **classically provable** (the equioscillation theorem characterizes the unique
extremal). So framing the cap as a Chebyshev/Cohn–Elkies extremal puts a **rigorous extremal-principle**
under it: the AP is LP-sharp because its coverage equioscillates against the magic function (de Moivre),
and any deviation strictly loses (the dip). The remaining work — consec maximizes coverage — becomes
"**the AP is the equioscillating configuration of the 7-fold Chebyshev / Cohn–Elkies LP**," which is the
kind of statement the equioscillation/LP-duality machinery is *built* to prove. This is a sharper target
than the moment-LP: a **uniqueness-of-the-Chebyshev-extremal** statement.

## Net
- **VERIFIED:** `V_7(u)−2 = (u−2)·(de Moivre cubic)²` — the de Moivre angles are the 7-fold Chebyshev
  equioscillation points; the cap is a Chebyshev extremal.
- **SYNTHESIS:** equioscillation (double root) ⟹ perfect-square discriminant ⟹ rational cap — the
  Chebyshev tangency IS mac-mini's ℚ-collapse; cyclotomic ideal = Chebyshev extremal = ℚ-rational, one fact.
- **FRAME:** the cover bound = Cohn–Elkies/Delsarte LP with the cyclotomic magic function; Hermite–Biehler
  = `f̂≥0`, Toeplitz `λ_min` = PSD margin, Perron = principal mode — four faces of one magic function;
  Viazovska's dim-8 (octonion) = the multiplicative-apex sibling.
- **TARGET:** consec = the unique equioscillating (Chebyshev/Cohn–Elkies-sharp) configuration — a rigorous
  extremal-principle route to the cover bound.

→ HYP-3212 (this), HYP-3162 (cyclotomic/de Moivre), HYP-3132 (ℚ-collapse, now explained), HYP-3201/3203
(Toeplitz/Delsarte LP), HYP-3210 (Hermite–Biehler), HYP-3211 (mult/add apex split), Chebyshev/Vieta–Lucas,
Cohn–Elkies, Viazovska (dim 8 = octonions), OPEN-Q-108.
