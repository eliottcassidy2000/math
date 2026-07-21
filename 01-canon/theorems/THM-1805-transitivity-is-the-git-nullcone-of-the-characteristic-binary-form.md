---
id: THM-1805
title: "TRANSITIVITY IS THE GIT NULLCONE OF THE CHARACTERISTIC BINARY FORM — the representation-theoretic root of 'tournaments are in/transitivity itself'. The map T ↦ char_A(x), homogenized to the degree-n binary form char_A(x,y) = yⁿ·char_A(x/y), sends a tournament to a point of Sym^n(ℂ²), the SL₂-representation of binary n-ic forms. (1) TRANSITIVE ⟺ char_A = xⁿ = the maximally UNSTABLE binary form (a single root of multiplicity n) = the GIT NULLCONE — because a transitive tournament is strictly triangular, hence nilpotent (THM-895's λ=0 ⟺ transitive is exactly 'the characteristic form is nilpotent = nullcone'). (2) THE UNSTABLE LOCUS IS EXACTLY THE TRANSITIVE TOURNAMENTS, verified n = 3..6: by Hilbert–Mumford a binary n-ic is unstable iff it has a root of multiplicity > n/2, and the count of tournaments with max-root-multiplicity > n/2 is 6, 24, 120, 720 = n! = the number of labeled transitive tournaments EXACTLY, with NO non-transitive unstable form (every unstable one is acyclic, c₃ = 0). So the only tournament whose adjacency characteristic polynomial has a root of multiplicity exceeding n/2 is the transitive one, and there it is xⁿ. (3) THE TRACE MOMENTS ARE THE SL₂-INVARIANTS: tr(Aᵏ) = power sums of the roots = the coefficients of the form, so the moment-nullcone ladder (THM-1775) is literally the coefficient/invariant data of the characteristic form, and 'all moments vanish ⟺ transitive' is 'the form is in the nullcone'. (4) THE FIBERS FORGET THE PERMANENT: T ↦ char_A is not injective — its fibers are the co-spectral classes, exactly where H splits at n = 6 (THM-1780, H = 13 vs 17). The characteristic binary form is the SL₂-invariant shadow of the tournament; H is the #P datum it cannot see. (5) THE SKEW FORM char_S is the EVEN companion (spectrum ±iλ, so char_S(−x) = (−1)ⁿchar_S(x)), and the half-dictionary x ↦ 2x+1 (S = 2A−J+I) carries char_A to char_S — but the NULLCONE lives on the A-side, since tr(S²) = −n(n−1) ≠ 0 forces S never nilpotent: transitive is xⁿ for A yet char_S = x(x⁴+10x²+5) at n=5, the {0,½,1} vs {−1,0,1} asymmetry made spectral"
status: >
  (1) PROVED (transitive ⟺ nilpotent ⟺ char_A = xⁿ is classical; the GIT-nullcone identity is
  the definition of the nullcone for binary forms).
  (2) VERIFIED by exhaustive labeled census n = 3,4,5,6: unstable count = 6,24,120,720 = n!,
  and every unstable tournament has c₃ = 0 (transitive).  The general statement "only the
  transitive tournament has a char_A-root of multiplicity > n/2" is a CONJECTURE for n ≥ 7,
  strongly evidenced; not proved.  The max-root-multiplicity distributions are recorded
  (e.g. n=6: mult 1 on 28608, mult 2 on 2480, mult 3 on 960, mult 6 on 720).
  (3) The trace-moments = power-sums = invariants identification is Newton's identities, exact.
  (4) The fiber = co-spectral class statement is definitional; the split at n=6 is THM-1780.
  (5) char_S even is the skewness of S (exact); tr(S²) = −n(n−1) is a one-line count; the
  half-dictionary shift is THM-1555.
  This is a UNIFYING FRAME with one new verified fact (the unstable locus = transitive
  tournaments) and several exact identifications; it proves no open problem but places
  transitivity, the moment ladder, and H's #P-jump inside binary-form GIT.
source: kind-pasteur-2026-07-20-S128c131 (owner: representation theory of binary forms and how it relates to tournaments = in/transitivity itself)
depends_on:
  - THM-895     # λ = 0 ⟺ transitive (nilpotent)
  - THM-1775    # the moment-nullcone template
  - THM-1780    # H leaves the ladder at n=6 (the fiber forgets the permanent)
  - THM-1555    # the half-dictionary S = 2A − J + I
related: [THM-1725, THM-1620]
script: 04-computation/tournament_binary_form_git_kps_S128c131.py (+ .out)
---

# THM-1805 — transitivity is the GIT nullcone of the characteristic binary form

The owner asked for the representation theory of binary forms and its relation to tournaments,
"which are in/transitivity itself." Here is the exact statement of that relation.

## The map

`char_A(x) = det(xI − A)` homogenizes to the **binary form** `char_A(x,y) = yⁿ char_A(x/y)` of
degree `n`, a point of `Sym^n(ℂ²)` — the irreducible `SL₂`-representation of binary `n`-ic
forms. So

> **`T ↦ char_A`  is a map  {tournaments on `n` vertices} → {degree-`n` binary forms}.**

`SL₂` acts on the target; the tournament symmetries (`S_n` relabeling, complement `T ↦ T^op`)
act on the source and both fix `char_A` (relabeling conjugates `A`; complement transposes it).

## (1)–(2) The nullcone is transitivity

A binary form's `SL₂` **GIT nullcone** — the unstable locus, forms with an `SL₂`-orbit closure
containing `0` — is by Hilbert–Mumford exactly `{`a root of multiplicity `> n/2}`, and its
extreme point is `xⁿ` (one root, multiplicity `n`).

> **`T` transitive `⟺` `char_A = xⁿ` `⟺` `char_A` is the maximally unstable form.**

(`⟸`/`⟹`: transitive `A` is strictly triangular in the linear order, so nilpotent, `char_A =
xⁿ`; conversely `char_A = xⁿ` means `A` nilpotent means acyclic means transitive — THM-895.)

And the **whole** unstable locus is nothing but the transitive tournaments:

| `n` | unstable tournaments (max root mult `> n/2`) | `n!` | non-transitive unstable |
|---|---|---|---|
| 3 | 6 | 6 | none |
| 4 | 24 | 24 | none |
| 5 | 120 | 120 | none |
| 6 | 720 | 720 | none |

**Every unstable tournament is acyclic.** So the only tournament whose characteristic polynomial
has a root of multiplicity exceeding `n/2` is the transitive one, and there the form is `xⁿ`.
(Verified `n ≤ 6`; conjectured beyond.)

## (3) The moments are the invariants

`tr(Aᵏ) = Σ λᵢᵏ` are the power sums of the roots of `char_A`, hence — by Newton — the
**coefficients** of the binary form, i.e. the basic `SL₂`-covariant data. So the
moment-nullcone template (THM-1775) applied to the trace functional *is* the statement that a
binary form lies in the nullcone iff its coefficient sequence is that of `xⁿ`, and the detection
depth `n` is the number of coefficients of a degree-`n` form. **The rational floor of the
moment ladder is the coefficient map of `Sym^n`.**

## (4) The fibers forget the permanent

`T ↦ char_A` is far from injective: its fibers are the **co-spectral classes**. That is exactly
where `H` splits — at `n = 6`, the fiber over `(0,0,12,12,10,48)` carries `H = 13` and `H = 17`
(THM-1780). So:

> the characteristic binary form is the **`SL₂`-invariant shadow** of the tournament; `H` (the
> permanent, `#P`) is the datum living in the fiber, invisible to the form.

This is the binary-form reading of the whole ladder: `Sym^n`-invariants sit at the bottom
(rational floor), the `#P` permanent at the top (THM-1780), and the fiber is the gap between
them.

## (5) The even companion and the ½

The skew matrix `S = A − Aᵀ` has spectrum `±iλ` (skew-symmetric), so `char_S` is an **even** form
(`char_S(−x) = (−1)ⁿ char_S(x)`) — a binary form with the extra `ℤ/2` (Weyl) symmetry
`x ↦ −x`, i.e. really a form in `x²`. The half-dictionary `S = 2A − J + I` (THM-1555) is the
change of variable `x ↦ 2x+1` between `char_A` and `char_S`. But the **nullcone lives on the
`A`-side only**: `tr(S²) = −n(n−1) ≠ 0` forces `S` never nilpotent, so the transitive tournament
is `xⁿ` for `A` yet `char_S = x(x⁴+10x²+5)` at `n = 5`. This is the `{0,½,1}` vs `{−1,0,1}`
asymmetry (THM-1555) made spectral: negation-symmetry (odd/even, the sign world) forbids the
nullcone, which lives in the affine world where transitivity can degenerate to a point.

## The unification

Two of the project's binary forms now sit in one frame:

- **tournament** — `char_A`, degree `n`, nullcone = transitive = `xⁿ`;
- **TNC/GMC kernel** — `R(u)`, degree `D`, nullcone = one-sided = degenerate; its monodromy
  exponents (THM-1725) are the ramification of the associated map, i.e. `SL₂`-covariant data.

The moment-nullcone ladder is a **ladder of binary forms with their GIT nullcones**, and
"in/transitivity" is the `n`-vertex tournament's position relative to the nullcone of its own
characteristic form. The odd/even axis the owner keeps returning to is the `SL₂` **Weyl
involution** `x ↦ −x`; the roots of unity are the **semistable/polystable spectra** (the Paley
tournament's eigenvalues are Gauss sums — roots of unity averaged).

## Named next

- **Prove (2) for all `n`:** show no non-transitive tournament has a `char_A`-root of
  multiplicity `> n/2`. Likely via the Perron–Frobenius simplicity of the strongly-connected
  part plus a rank bound.
- **The strictly-semistable stratum** (max root multiplicity `= n/2`, `n` even): `960` at
  `n = 6`. These are the "boundary" tournaments — a distinguished class worth identifying (an
  eigenvalue of multiplicity exactly `n/2`).
- **Covariants as tournament invariants:** the discriminant of `char_A` (vanishing = repeated
  eigenvalue), the Hessian, and the low transvectants are `SL₂`-covariants; compute what
  tournament structure each detects, extending THM-1780's "the form forgets `H`" to "which
  covariants see which strata".
