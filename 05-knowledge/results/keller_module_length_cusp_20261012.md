# The cusp curve carries no constituent: the A_3-length of THM-1300's Keller module is three on its fibre stratification (HYP-9257's value four is refuted)

2026-10-12, opus-2026-10-12-S24 (third addendum). **PROVED AS A SKETCH on
cited standard facts:** the characteristic-cycle bookkeeping below.
**FINITE-EXACT (sympy):** the Jelonek surface's singular locus is exactly the
cusp curve, and the transversal slice is an ordinary cusp (multiplicity two,
Milnor number two). **CITED:** the Kashiwara–Dubson index formula for
characteristic cycles, the local Euler obstruction of a plane curve germ
(equal to its multiplicity), Riemann–Hilbert for direct images, and
THM-4621(9) (the three named constituents). **OPEN caveat:** constituents
supported on strata finer than the fibre-count stratification (isolated
points of the cusp curve, the `a = 0` locus) are not addressed; nothing in
THM-4621 suggests any. Self-audited; no independent peer audit.

## Setting (THM-1300, THM-4620, THM-4621)

`F : A^3 -> A^3` is THM-1300's Keller map, étale everywhere (`det JF = -2`),
generically three-to-one, with Jelonek surface `J = {c_3 = 0}`,
`c_3 = 27a^2c^2 - 18abc + 16a + b^3c - b^2`, over whose points one sheet stays
and two escape, and cusp curve `C = {(4/(27t^2), 4/(3t), t)} = Sing J` over
which all three sheets escape. The Keller module `V_F = K[x,y,z]` with the
`A_3`-action through the dual fields is the direct image `F_+ O`, and by
Riemann–Hilbert its de Rham complex is the perverse sheaf
`P = Rf_* C_(A^3)[3]`. HYP-9257 guessed `length(V_F) = 4`: the three
constituents of THM-4621(9) — `O = K[F]`, the minimal extension `IC(std)` of
the rank-two standard `S_3` local system on `A^3 \ J`, the delta module
`IC_J` along `J` — plus one more supported on `C`.

## 1. The local structure of J along C (FINITE-EXACT)

`c_3` and its gradient vanish identically on `C`; solving `c_3 = ∇c_3 = 0`
gives exactly `a = 4/(27c^2)`, `b = 4/(3c)`: `Sing J = C`. The slice of `J`
by a plane transversal to `C` at the point `t = 1` is a plane curve of
multiplicity two with tangent cone the double line `3(3u - v)^2` and Milnor
number two (the resultant of the two partials vanishes to order exactly two):
an ordinary cusp, `A_2`. So `J` is a cuspidal edge along `C`, and its local
Euler obstruction at points of `C` is `Eu_J = 2` (an irreducible plane curve
germ has `Eu` equal to its multiplicity; products with a smooth factor do not
change it).

## 2. The stalk Euler characteristics of P (PROVED)

Let `B` be a small ball around a point `x` of the target and `φ(x) = χ(f^{-1}(B))`,
the Euler characteristic of the stalk of `Rf_* C` (the étale map `f` is
locally trivial along the strata, so the stalk is the cohomology of
`f^{-1}(B)`). Since `f` is étale, `f^{-1}(B)` is an open 3-manifold.

* `x` in `U = A^3 \ J`: three disjoint balls, `φ = 3`.
* `x` in `J \ C`: `f^{-1}(B)` is the staying sheet's ball (`χ = 1`) together
  with the escaping pair, a connected double cover of `B \ J ≃ S^1`
  (transposition monodromy; `χ = 0`): `φ = 1`.
* `x` in `C`: `B \ J ≅ (ball in C^2 minus a cusp) × interval ≃ S^3 minus the
  trefoil`, `χ = 0`; its connected three-sheeted cover `X` (monodromy `S_3`,
  the braid-group representation of the trefoil) has `χ = 0`; the staying
  sheet over `B ∩ (J \ C)` is a surface `S ≃ S^1` (`χ = 0`) attached to `X`
  along a collar; Mayer–Vietoris gives `φ = 0`.

In the perverse normalization `φ_P = -φ`: `-3, -1, 0`.

## 3. The index formula and the multiplicities (PROVED on cited facts)

Write `CC(V_F) = m_3 [T*A^3] + m_J [T*J] + m_C [T*C]` (no point strata
occur in the fibre stratification). The Kashiwara–Dubson index formula
`φ_P(x) = Σ_S (-1)^{dim S} m_S Eu_{S̄}(x)` gives, successively,

    x in U:      -m_3 = -3                      =>  m_3 = 3,
    x in J \ C:  -3 + m_J · 1 = -1              =>  m_J = 2,
    x in C:      -3 + 2 · Eu_J(x) - m_C = 0     =>  m_C = 2·2 - 3 = 1.

The value `m_J = 2` agrees with THM-4621's étale-local form
`O ⊕ O[1/s] ⊕ j_+ L_(-1)` near `J \ C` (`T*J` multiplicities `0 + 1 + 1`),
which calibrates the sign convention.

Now distribute `m_C = 1` among the constituents. `CC` is additive in exact
sequences and every multiplicity is a non-negative integer.

* `O` contributes `0` to `T*C`.
* `IC(std)` contributes some `a_C ≥ 0`.
* `IC_J` (the minimal extension of a rank-one local system `L` on `J \ C`;
  generically the delta module `O[1/s]/O`) contributes `b_C ≥ 1`: in the
  transversal slice `IC_J` is the IC complex of the cusp `Γ` with coefficients
  `L|_Γ`, and the index formula on the slice gives `b_C = Eu_Γ(0) + χ(IC_x)`
  with `χ(IC_x) = -1` (one branch, trivial monodromy) or `0` (nontrivial
  monodromy), i.e. `b_C = 1` or `2`.
* Any constituent supported on `C` is the minimal extension of a local system
  of rank `r ≥ 1` on an open subset of the smooth curve `C` and contributes
  `r ≥ 1` to `T*C`.

Hence `1 = m_C = a_C + b_C + Σ r ≥ 0 + 1 + Σ r`, so `a_C = 0`, `b_C = 1`
(the monodromy of `L` around `C` is trivial), and **no constituent is
supported on `C`**. With THM-4621(9) (`length ≥ 3`, the three constituents
accounting for `m_3 = 1 + 2`, `m_J = 0 + 1 + 1`):

> **length_(A_3)(V_F) = 3 on the stratification `U, J \ C, C`.**

The guess of HYP-9257 — that the drop of the fibre count from one to zero
along `C` must be paid by a constituent supported on `C` — is refuted: the
drop is paid by the cuspidal-edge geometry of `J` itself (`Eu_J = 2`), which
raises the `T*C`-multiplicity of the delta module along `J` to exactly one.

## 4. What is not covered, and the decisive remaining check

* The stratification used is the fibre-count stratification. A constituent
  supported on a finer stratum (an isolated point of `C`, or the curve
  `J ∩ {a = 0}` if the local structure of `J` or of the escaping sheets
  changes there) would not be detected by the three equations above. The
  `a = 0` stratum is not on `C` (`a = 4/(27t^2) ≠ 0`); THM-4621 (audit A)
  reports the same fibre counts `3/1/0` there. A full answer computes `φ`
  at a point of `J ∩ {a = 0}` and the local type of `J` there.
* The three standard inputs — the index formula, `Eu` of a curve germ, and
  the stalk-as-cohomology-of-`f^{-1}(B)` for an étale algebraic map — are
  cited, not re-proved; the topological computation of `φ(C) = 0` is the
  only new input, and it is robust: every piece has Euler characteristic
  zero.
* A `b`-function or characteristic-cycle computation in Singular/Macaulay2
  would confirm `m_C = 1` directly; the repo's toolkit lacks them.

## 5. Reproduction

    python3 -B 04-computation/experiments/keller_cusp_length_20261012.py

The [script](../../04-computation/experiments/keller_cusp_length_20261012.py)
checks the singular locus, the slice multiplicity, tangent cone and Milnor
number, and performs the integer bookkeeping of section 3 as assertions.
