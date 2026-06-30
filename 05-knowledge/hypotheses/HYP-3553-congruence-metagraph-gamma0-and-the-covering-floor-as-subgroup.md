---
id: HYP-3553
title: The Gamma_0(N) <-> metagraph <-> LRC dictionary -- developing the four bullets. (1) Han-Lee's congruence Siegel 2nd moment is the LRC floor WITH the covering built in as the subgroup Gamma_0(N) (index psi(N)=N prod(1+1/p)), giving a SET-INDEPENDENT floor depending only on N via phi(N),psi(N),J2(N) -- cleaner than the totient sum. (2) The metagraph Burnside MASS FORMULA is the template for a closed-form LRC FIRST moment under congruence (what the union bound lacks). (3) The metagraph mu_2 (ordering pairs) is the finite rehearsal for the S_2 variance Var(N_R) (THM-579). (4) G_n(N) = the LEVEL-N congruence metagraph (Z/N-circulant tournaments, marked Z/N) is the tournament image of Gamma_0(N) / X_0(N); its mass = the dihedral Burnside (1,1,2,4,4,6 for N=3,5,7,9,11,13) with the Paley class as the CM/cusp point. Tournaments ARE PSL(2,Z)-modular (the-modular-tournament); the covering congruence is the level structure.
status: SYNTHESIS + verified arithmetic (psi/phi/J2; circulant-metagraph masses). The dictionary and the four bullets are a research PROGRAM; the exact congruence-floor constant is the Han-Lee computation (structure shown, constant not pinned). Not a proof.
source: mac-mini-2026-06-29-S20
related:
  - HYP-3561   # the metagraph is a finite Siegel transform (this adds the congruence subgroup)
  - HYP-3550   # the totient-sum / Euler-product floor (this replaces it by Gamma_0(N), set-independent)
  - THM-579    # the sheet-count CV / Var(N_R) 2nd moment (bullet 3)
  - THM-586    # Paley dihedral Burnside = the G_n(N) mass / CM point (bullet 4)
  - THM-585    # rotational/circulant tournaments (the level-N structure)
external: arXiv:2507.05905 (Han-Lee, congruence Siegel moments); Gamma_0(N), X_0(N), Dedekind psi
reflections:
  - the-modular-tournament.md          # tournaments = PSL(2,Z); X_0(p) genus controls OCR (kps S18e)
  - zeta2-governs-the-lonely-runner-floor.md
results:
  - 04-computation/congruence_metagraph_gamma0_floor_macmini_20260629.py
  - 05-knowledge/results/congruence_metagraph_gamma0_floor_macmini_20260629.out
---

# HYP-3553 -- the Gamma_0(N) <-> metagraph <-> LRC dictionary

## The dictionary
Building on [the-modular-tournament] (tournaments ARE PSL(2,Z): S=complement, ST=3-cycle, T=vertex-add;
X_0(p) genus controls the OCR residual) and HYP-3561 (the metagraph is a finite Siegel transform):
| Siegel / Gamma_0(N) (Han-Lee) | tournament / metagraph | LRC floor |
|---|---|---|
| `SL(2,Z)` | the tournament modular group | speed relabeling |
| **`Gamma_0(N)`** (index `psi(N)=N prod(1+1/p)`) | **`G_n(N)` level-N congruence metagraph** | the **covering mod N** |
| `X_0(N)` modular curve | the metagraph at level N | the covering nerve |
| primitive `(p,q)=(p0,q0) mod N` | tournament/resonance w/ covering residue | Farey point mod N |
| density `1/phi(N)`, `zeta(2)` norm | congruence orbit density, Burnside norm | surviving-resonance density |
| 1st moment (Siegel mass formula) | Burnside **mass formula** | `E[#lonely]` |
| 2nd moment (Rogers/Schmidt) | metagraph `mu_2` (ordering pairs) | `Var(N_R)` (THM-579) |

## Bullet 1 -- the floor with covering BUILT IN (set-independent)
The LRC floor `c_q >= 1/(2 zeta(2)) = 3/pi^2` is the Farey/`zeta(2)` density; currently the covering enters
as a totient SUM `sum phi(b) delta_b` (HYP-3550) or set-dependent inclusion-exclusion. Han-Lee's congruence
Siegel 2nd moment puts the covering INSIDE the average as the subgroup `Gamma_0(N)`: the floor becomes a
congruence 2nd moment over `Gamma_0(N)\H`, with the covering modulus `N` entering through `psi(N), phi(N),
J2(N)` ONLY -- **set-independent** (depends on `N`, not the speed set). Verified arithmetic at `N=14`:
`phi=6, psi=24=[SL(2,Z):Gamma_0(14)], J2=144`. The congruence per-class density `1/(zeta(2) J2(N))` summed
over the `phi(N)` surviving unit classes gives the structure; the exact constant is the Han-Lee evaluation.
This is "cleaner than the totient sum" and removes the set-dependence that blocks the uniform floor.

## Bullet 2 -- the mass formula = the missing first moment
The union bound (`floor >= sum - overlaps`) is an INEQUALITY with no clean first moment. The metagraph's
Burnside **mass formula** (`#classes = (1/n!) sum_sigma Fix(sigma)`, exact) is the template: an exact average
= the FIRST moment `E[#surviving lonely points]`. Han-Lee's Siegel 1st-moment formula (with congruence) is
its continuous twin. Replacing the union bound by a mass-formula first moment + a second-moment
concentration is the first/second-moment method -- and Han-Lee supplies both moments WITH congruence.

## Bullet 3 -- mu_2 (ordering pairs) = the shared second-moment engine
The metagraph's `mu_2` -- the count/correlation of ordered PAIRS of arcs (the 2-point function) -- is the
finite, exactly-computable rehearsal of the LRC `S_2` variance `Var(N_R)` (THM-579, the sheet-count CV) and
of Han-Lee's `int hat f^2` (the Rogers/Schmidt pair correlation). All three are PAIR second moments; the
metagraph one is bounded and clean (CV(H) ~ 0.5-0.6, HYP-3561), the testbed for the analytic `S_2` bound.

## Bullet 4 -- G_n(N), the tournament image of the covering congruence
**Definition.** `G_n(N)` = the LEVEL-N congruence metagraph = the iso classes of `Z/N`-circulant
tournaments (a marked `Z/N` structure, the tournament analog of `X_0(N)`'s marked cyclic `N`-subgroup),
under the multiplier group `(Z/N)*`. A circulant tournament on `N` (odd) vertices is a sign-choice over the
`(N-1)/2` pairs `{r, N-r}`; the multiplier group acts; **mass = the dihedral Burnside (THM-585/586)**.
Computed: `mass(G_n(N)) = 1,1,2,4,4,6` for `N=3,5,7,9,11,13` -- a tiny, finite-index structured subset of
the full metagraph (just as `Gamma_0(N)` is finite-index in `SL(2,Z)`). The **Paley class** (QR connection
set, `N` prime `=3 mod 4`) is the distinguished **CM / cusp point** -- the tournament image of the
`Gamma_0(N)` cusp, where the OCR-genus law (kps) lives.

## Creative extensions (metagraph structures for proofs)
- **Hecke operators.** Vertex-addition `T` (the modular `T`) is a Hecke-like raising operator `G_n -> G_{n+1}`;
  the metagraph recursion is its action. Diagonalizing it (eigen-tournaments) would give H-recursions, the
  analog of Hecke eigenforms; the OCR/`X_0(p)` genus predicts which eigenvalues are "rational."
- **Genus law as a closed-form predictor.** kps's OCR-denominator law (closed form iff `X_0(p)` genus 0/1)
  says WHICH metagraph invariants admit closed forms -- a free predictor for the LRC first-moment program
  (the floor constant should be "rational/modular" exactly when the covering modulus `N` has `X_0(N)` small genus).
- **Metagraph as expander.** The spectral-moment gap (HYP-3561) measures mixing; a spectral-gap bound would
  give typical-case concentration (most tournaments/covering sets behave like the mean) -- the qualitative
  content of the second-moment method.

## What it buys
A precise dictionary turning the LRC covering floor into a `Gamma_0(N)` congruence second moment (covering
as a subgroup, set-independent), a concrete new object `G_n(N)` (the level-N congruence metagraph, the
tournament image of `X_0(N)`), and a research program (first+second moment with congruence, Hecke
recursions, genus-law closed-form predictor). The exact congruence-floor constant is the next computation.
