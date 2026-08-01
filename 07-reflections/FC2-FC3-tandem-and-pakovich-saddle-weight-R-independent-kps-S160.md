---
source: kind-pasteur-2026-07-24-S160 (Opus 4.8)
status: RESULT. (1) FC(2) modular transversality scan D=2..6: ALL transversal (capacity=P) -- the diagnostic
  behaves IDENTICALLY for FC(2) (known-true) and FC(3), so it cannot distinguish them (both lean true; an
  isolated FC(3) counterexample would be invisible). (2) FC(3) pushed to D=6 (first 14 leaks independent, no
  overshoot). (3) Pakovich/saddle test: for a composition psi=R o W, the saddle weight C=exp(psi_{D-1}/(D psi_D))
  is INDEPENDENT of R -- it equals exp(depression-shift of W). Consequence (verified concretely): a composition
  inherits W's roots-of-unity cancellation and can only COARSEN it, never refine it, so no Pakovich composition
  closes the FC leak. Analytic obstruction complementing the (algebraic) transversality; both say FC(3) true
  in the natural family.
tags: [factorial-conjecture, transversality, pakovich, moment-problem, saddle-point, roots-of-unity, composition]
related: [kps-S154, kps-S155, kps-S157, kps-S159]
---

# FC(2)/FC(3) in tandem, and the Pakovich composition saddle weight

## 1. FC(2) transversality scan (control) -- D=2..6, all transversal
Antisymmetric family (odd moments auto-zero), factorial functional `L(x^a y^b)=a!b!`, leak-Jacobian rank mod
`p=2^31-1`:
| D | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|
| P | 1 | 3 | 5 | 8 | 11 |
| rank | 1 | 3 | 5 | 8 | 11 | (= P each) |
**All transversal (capacity=P, tower terminates).** FC(2) is believed *true*, and it behaves **identically** to
FC(3) (kps-S159) under this diagnostic. **So generic transversality cannot separate FC(2)-true from FC(3):** it
supports "FC(3) behaves like FC(2)" (leans true), but a hidden *isolated* FC(3) counterexample stays invisible.

## 2. FC(3) to D=6
`P=26`; tested `J=14` leaks (`f^41`): rank `= 14 = min(J,P)` -- **first 14 leaks independent, no overshoot**
(full `J>=P` would need `f^78`). Extends the fully-confirmed D=2,3,4,5 (kps-S159). No degree shows a rank drop.

## 3. Pakovich composition vs the saddle weight -- the R-independence
Framework (interpretation): vanishing-moment solutions of a polynomial moment problem factor as **`psi=R o W`**
(Pakovich); the moments `L(psi^m)=int psi^m e^{-|x|}` have a steepest-descent asymptotic whose dominant subleading
factor is the **saddle weight** `C=exp(psi_{D-1}/(D psi_D))` (the depression shift of the exponent). Moments can
*vanish* only if saddle contributions cancel -- i.e. only if `C` is a **root of unity** (the Conway-Jones
condition the seed exploits).

> **VERIFIED FACT.** For `psi=R o W`, `C=exp(psi_{D-1}/(D psi_D)) = exp( w_{deg W-1}/(deg W * w_{deg W}) )` --
> the depression shift of the **inner** map `W` alone. **It is independent of the outer map `R`.**
Checked exactly for `R = y, y^2+7y+4, 2y^3-y+9, y^5` over `W=5z^2-2z+3`: all give `C=e^{-0.2}`. A nontrivial
root-of-unity `C` requires a **complex** `W` (imaginary shift); real `W` gives `C=1` only.

**Consequence.** `R` is invisible to the saddle weight. So a composition `R o W` inherits its root-of-unity
cancellation from `W`, and `R` can only **coarsen** it, never refine it. Concretely, for the seed
`W=x+omega y+omega^2 z` and `R(y)=y^d` (`psi=W^d`): `L(psi^k)=(dk)!*[3|dk]`, so
- `d` coprime to 3 (`d=1,2,4`): leak pattern **preserved** (`3|k`, just reindexed);
- `3 | d` (`d=3,6`): phases collapse, **all** moments nonzero -- cancellation **destroyed**.
Never does `R` close the `3|k` leak that `W` carries. **No Pakovich composition `psi=R o W` closes the FC leak
via the saddle mechanism.**

## 4. Three independent lines, one verdict
- **Algebraic (transversality):** capacity=P through D=5 (full) and D=6 (partial); no overshoot. No counterexample
  as a family.
- **Analytic (saddle/Pakovich):** `C` is `R`-independent; composition inherits/coarsens `W`'s cancellation, never
  refines -- no composition closes the leak.
- **Structural (kps-S155 radial/angular):** `R` is the *radial* reparametrization (invisible to `C`), `W` the
  *angular* part (sets the roots of unity). FC is the radial restriction; the angular cancellation is exactly the
  roots-of-unity one, which leaks at multiples. Same split, same conclusion.

All three say **FC(3) holds in the natural (roots-of-unity / composition) family.** The single surviving opening
is an **isolated, non-composition, non-transversal** `psi` -- a genuine Kontsevich-Zagier period-rigidity
coincidence (kps-S154), which none of these three methods can reach.

## 5. Honest caveats
The saddle-weight/Pakovich reading in Sec 3 is an **interpretation** (I do not have the originating framework's
exact definitions); the *fact* `C` is `R`-independent, and the concrete `W^d` leak patterns, are verified and
stand on their own. FC(2)/FC(3) transversality is generic (isolated special points not excluded). The
convergence of three methods is strong evidence, not a proof; the proof lives in the period-rigidity gap.

Files: `/tmp/{fc2scan,fc3_d6,pakovich}.py`. Builds on kps-S154/S155/S157/S159.
