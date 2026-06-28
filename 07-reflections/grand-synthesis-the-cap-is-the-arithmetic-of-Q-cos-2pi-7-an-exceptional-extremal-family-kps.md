# Grand synthesis: the LRC(14) cap is the arithmetic of ℚ(cos 2π/7) — an exceptional extremal family

*kind-pasteur-2026-06-27-S31an capstone. The owner asked for ever more creative, abstract synthesis.
Three threads this session — octonion/G₂, Chebyshev equioscillation, the cyclotomic discriminant — collapse
into ONE statement: **the LRC(14) cap/dip is the arithmetic of the real cyclotomic field `ℚ(cos 2π/7)`**
(the de Moivre cubic `x³+x²−2x−1`). Everything the team has found about the cap is a facet of this single
number field, and LRC(14) is the **additive, 1-D, 7-fold member of an exceptional extremal family** that
includes Viazovska's sphere packings. Below: the collapse, the meta-picture, and a batch of forward leads.*

## The collapse — one field, six facets
`F = ℚ(cos 2π/7)`, the totally-real cubic subfield of `ℚ(ζ₇)`, defined by the de Moivre cubic
`m(x)=x³+x²−2x−1`, roots `2cos(2πj/7) = {−1.8019, −0.4450, 1.2470}`, discriminant `disc(m)=49=7²`,
class number 1, Galois group `C₃`. Every cap-fact is a property of `F`:

| Facet | Statement | = property of `F` |
|---|---|---|
| **Cyclotomic ideal** (S31ak) | cap = the 7-fold-symmetric ideal; zeros at de Moivre angles | the roots of `m` |
| **Chebyshev** (HYP-3212) | `V₇(u)−2 = (u−2)·m(u)²`; cap = equioscillation extremal | `m` is the 7-fold Chebyshev factor |
| **Rationality** (HYP-3132) | cap, dip ∈ ℚ (disc 9 ℚ-collapse) | the **double root** `m²` ⟹ perfect square ⟹ `F` rationalizes |
| **Discriminant** (Thread 3) | deep denominators carry `7²` (`dip₈`, `m_P`: `v₇=2`) | `disc(m)=7²` |
| **Joukowski / trace** (S31am/3211) | `z+1/z` bridges the two maps | the **Galois trace** `ζ↦ζ+ζ̄=2cos`, `ℚ(ζ₇)→F` |
| **Magic function** (HYP-3212) | cover bound = Cohn–Elkies/Delsarte LP | the LP dual is the `F`-cyclotomic magic function |

**One field, `F=ℚ(cos 2π/7)`, is the arithmetic home of the cap.** The rationality (the only reason the
proof is *finite/solvable*) is the equioscillation double-root of `F`'s defining polynomial; the `7²` in
the denominators is `disc(F)`; the Joukowski bridge is `F`'s trace map.

## The meta-picture — an exceptional extremal family at the apex
LRC(14) is not an isolated hard case; it is the **additive, 1-D, 7-fold** member of a family of extremal
problems governed by **magic functions at exceptional structures**:

| problem | dimension / structure | magic function from | apex object |
|---|---|---|---|
| Viazovska sphere packing | dim 8 (`E₈` = octonion lattice) | modular forms | **octonions** (mult. apex) |
| Viazovska / CKMRV | dim 24 (Leech) | modular forms | Golay/`M₂₄` |
| **LRC(14)** | **1-D, 7 sectors** | **cyclotomic `ℚ(cos 2π/7)`** | **the de Moivre cubic** (add. apex) |

The **apex prime/dimension** carries an **exceptional structure** (Fano/octonions/`E₈` at the
multiplicative face; the cyclotomic cubic at the additive face — HYP-3211), and the extremal bound is a
**magic-function / equioscillation** statement. **LRC(14) is the cyclotomic (additive) sphere-packing
problem of the heptagon.** Viazovska's dim-8 (octonion) and LRC's 1-D-7-fold are the multiplicative and
additive magic-function problems of the *same* apex 7, related by the Galois trace.

This is the deepest answer yet to "why do the exceptional structures crowd into LRC(14)" (mac-mini's
four-faces): because the cap is a magic-function extremal at the apex 7, and 7 is the smallest prime with
nontrivial exceptional structure (Fano/octonion multiplicatively, cubic-cyclotomic additively). The
`14=2·7` stacks the 2-adic (the dyadic/`E₈`-doubling) on it.

## Forward creative leads (the next abstractions, flagged)
1. **⊕ Stark / L-values.** `F`'s units are the de Moivre angles; the **Stark conjecture** ties the
   `L'(0,χ)` of `F` to a regulator of these units. **Is the cap (or `m_P`) an `L`-value / regulator of
   `ℚ(cos 2π/7)`?** The `7²=disc` denominator is the conductor — a strong hint. This would give the cap a
   *closed arithmetic form*.
2. **⊕ Modular magic function.** Make the Cohn–Elkies magic function **explicit** as a (quasi)modular form
   for `Γ₀(7)` — the level-7 analog of Viazovska's level-1 construction. mac-mini's `E₂`/Eisenstein spoke
   is the generator. This would *construct* the LP-sharp certificate.
3. **○ Beraha / Tutte.** `B₇=2+2cos(2π/7)=3.247` is a Beraha number (chromatic-zero accumulation point of
   the Tutte polynomial). **Is the cap a Tutte/chromatic value of a 7-structured graph** (Heawood/Fano)?
   The cyclotomic appears in both chromatic roots and the cap.
4. **○ Mahler measure / Lehmer.** `m(x)` is a non-reciprocal cubic; its Mahler measure `M(m)=|−1.8019|·…`
   the dip ↔ a height/Mahler deviation from the Kronecker (cyclotomic, `M=1`) ideal.
5. **○ Subshift / transfer operator.** The empty-sector sequence as `x` sweeps is a Sturmian/rotation
   subshift (THM-536/485); the de Moivre angles = its spectral measure; the cap = a cylinder probability
   via the transfer operator's Perron eigenvalue (ties Perron HYP-3210 to dynamics).

## Net
- **SYNTHESIS:** the cap/dip of LRC(14) is the arithmetic of `F=ℚ(cos 2π/7)`; its six facets (cyclotomic
  ideal, Chebyshev equioscillation, rationality, `disc=7²`, Joukowski-trace, magic function) are one field.
- **META:** LRC(14) = the additive/cyclotomic member of the magic-function exceptional-extremal family
  (Viazovska dim-8/octonion = the multiplicative sibling), the heptagon's sphere-packing problem.
- **LEADS:** Stark/L-value closed form for the cap; explicit `Γ₀(7)` modular magic function; Beraha/Tutte;
  Mahler/Lehmer; subshift transfer operator. The first two are proof-relevant (a closed cap + a constructed
  certificate).

→ HYP-3213 (this), HYP-3211/3212 (octonion/Chebyshev), HYP-3162 (cyclotomic), HYP-3132 (ℚ-collapse),
HYP-3099/3210 (two maps/Joukowski), THM-577 (rational cap), Stark, Viazovska, Cohn–Elkies, Beraha, the
four-faces-of-14, OPEN-Q-108.
