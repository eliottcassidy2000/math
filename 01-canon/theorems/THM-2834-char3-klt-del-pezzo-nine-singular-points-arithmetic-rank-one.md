---
id: THM-2834
title: "A char-3 klt del Pezzo with nine tame singular points and arithmetic Picard rank one"
status: >
  PROVED (quasi-smoothness, tame klt singularity census, Fano index, Method-A
  cap, exact-family decoupling, corrected intersection lattice, and
  Frobenius-orbit rank argument) + VERIFIED-EXACT (irreducibility/separability,
  strata counts, corrected Gram/Frobenius matrices, and point counts over
  F_{3^k}, k = 1..4, matching the predicted Frobenius eigenvalues).
  X = X_14 in P(2,2,7,7) over F_3 is a klt del Pezzo surface with NINE
  singular geometric points (7 x A_1 + 2 x 1/7(1,1), all tame), Fano index 4,
  arithmetic Picard rank rho(X/F_3) = 1, geometric rank 7.
source: mac-mini-2026-07-28-S171 (external open-problem raid; Epoch
  FrontierMath "Surface with a High Number of Singularities", char 3, > 7 sings)
depends_on: []
related: []
script: 04-computation/delpezzo_char3_P2277_ninesing_verify_macmini_S171.py
output: 05-knowledge/results/delpezzo_char3_P2277_macmini_S171.out
script_sha256: 96c8d345d8768de1d9ea63a216aba3a6650e3f9b44c31a388813508ca2232a0c
output_sha256: e93e82425bf18293f724c62917ebde46eb47bb33f9827a83dc13cfd5ae6be90f
hash_basis: LF-normalized bytes
---

# THM-2834 — X_14 ⊂ P(2,2,7,7)/F_3: nine tame singular points, ρ_arith = 1

> **Corrected lattice audit (MISTAKE-318).**  The first version incorrectly
> wrote `C_i ≡ D_i` and consequently gave a positive diagonal intersection
> matrix.  In fact `div(ell_i)=C_i+D_i=2H`, so
> `D_i=2H-C_i`.  The repaired intersection matrix and negative-cycle
> Frobenius action below still prove, rather than merely suggest,
> `rho(X_bar)=7` and `rho(X/F_3)=1`.

## The surface

    X = { F = 0 } ⊂ P(2,2,7,7),   F = A(x0,x1) + x2^2 + x3^2,
    A(x0,x1) = x0^7 + x0^2 x1^5 + 2 x1^7   over F_3,
    a(t) = t^7 + t^2 + 2 irreducible (hence separable) over F_3.

* Well-formed (all weight triples coprime), weights prime to 3 (tame).
* Quasi-smooth: `∇F = (A_0, A_1, 2x_2, 2x_3)`; vanishing forces `x2=x3=0` and,
  by Euler (`7·A = x0 A_0 + x1 A_1`, `7 ≡ 1 ≠ 0` in `F_3`), a common root of
  the separable pair `(a, a')` — none; at `(1:0)`, `A_0 = 7x0^6 ≠ 0`.
* Fano: `-K_X = O_X(18 - 14) = O_X(4)` ample, index 4.  klt: all
  singularities are tame cyclic quotients.
* Singular census = ambient strata only: 7 geometric points on the (2,2)-line
  (roots of `a`; each `1/2(1,1) = A_1`), 2 geometric points on the (7,7)-line
  (`x2^2 + x3^2 = 0`, rational over `F_9`; each `1/7(2,2) ≅ 1/7(1,1)`).
  **Nine** singular points, exceeding the requested `> 7`.

## Picard rank

Work over `k=Fbar_3`, put `H=[O_X(1)]`, choose `iota^2=-1`, and factor

    A = product_{i=1}^7 ell_i,       x2^2+x3^2 = s t,
    C_i = V(ell_i,s),                D_i = V(ell_i,t).

The three divisor identities that control the calculation are

    div(ell_i)=C_i+D_i=2H,           div(s)=sum_i C_i=7H,
    div(t)=sum_i D_i=7H.                                      (1)

Thus `D_i=2H-C_i`; in particular, `C_i` and `D_i` are **not** numerically
equivalent.  Weighted intersection theory gives

    H^2 = 14/(2·2·7·7)=1/14,         H.C_i=2·7/(2·2·7·7)=1/14. (2)

For completeness, all local intersections can be read on the tame quotient
covers.  Distinct `C_i,C_j` meet only at the `1/7(1,1)` point `s=0`, so their
intersection is `1/7`.  The pair `C_i,D_i` meets only at the corresponding
`A_1=1/2(1,1)` point, with intersection `1/2`; `C_i,D_j` are disjoint for
`i!=j`.  Intersecting `C_i+D_i=2H` with `C_i` now gives

    C_i^2=-5/14,        C_i.C_j=1/7 (i!=j),
    C_i.D_i=1/2,        C_i.D_j=0 (i!=j).                     (3)

The Gram matrix of `C_1,...,C_7` has eigenvalue `1/2` on the all-ones line
and eigenvalue `-1/2` on its six-dimensional sum-zero complement.  Its
determinant is `1/128`, so these seven curve classes are independent, while
`H=(sum C_i)/7`.

There cannot be further geometric classes.  The function field is rational:
on `x1!=0` its degree-zero subfield is

    k(x0/x1, s^2/x1^7),

because `st/x1^7=-a(x0/x1)`.  On the minimal resolution `pi:Y->X`, the seven
`A_1` points contribute seven disjoint `(-2)`-curves.  Each `1/7(1,1)` point
contributes one `(-7)`-curve with discrepancy `-5/7`.  Hence

    K_X^2=16H^2=8/7,
    K_Y^2=8/7+2·(-7)(5/7)^2=-6,
    rho(Y)=10-K_Y^2=16.

The nine exceptional curves are independent, so
`rho(X_k)=rho(Y)-9=7`.  Combined with the nonsingular Gram matrix, this proves

    NS(X_k) tensor Q = <C_1,...,C_7>_Q,              rho(X_k)=7.             (4)

Finally, arithmetic Frobenius cycles the seven roots and interchanges `s,t`.
After cyclic indexing,

    Frob(C_i)=D_{i+1}=2H-C_{i+1}.

Writing `v_i=C_i-H` gives `sum v_i=0` and the clean action

    Frob(v_i)=-v_{i+1}.                                             (5)

Thus the six sum-zero eigenvalues are `-zeta_7^j`, `j=1,...,6`, while `H`
is fixed.  None of the six is `1`, so

    rho(X/F_3)=rank NS(X_k)^Frob=1.                                 (6)

The exact matrix audit checks preservation of (3), Frobenius order `14`, and
fixed-space dimension `1`.  Independent point counts give
`#X(F_{3^k})=q^2+c_kq+1` with `c_k=2,0,2,0` for `k=1,...,4`, exactly the
trace of `{1} union {-zeta_7^j:j=1,...,6}`.  More generally the trace is
`1+(-1)^(k+1)` when `7` does not divide `k`.

## Two structural obstructions (proved, recorded for reuse)

1. **Method-A cap.**  For a smooth `Y = CI(2,2) ⊂ P^4` with a diagonal `±1`
   involution `g` free in codimension 1, signature `(3,2)` is forced and
   `Y ∩ P^1_-` is empty on a smooth `Y`: at a common root `p` of the two
   restricted binary quadratics `B_i = ℓ_p m_i`, both gradients are
   proportional to `dℓ_p` (the `x_+`-partials vanish at `x_+ = 0`), so `Y` is
   singular at `p`.  Hence `|Y ∩ Fix(g)| <= 4`: quotients of smooth CI(2,2)
   by diagonal involutions can never reach 8 singular points.
2. **Exact-family decoupling.**  For degree `14` in the particular weights
   `(2,2,7,7)`, the equation
   `2(a+b)+7(c+d)=14` has no solution with both `a+b>0` and `c+d>0`.
   Consequently every such hypersurface has `F=A(x0,x1)+C(x2,x3)`.  Over the
   algebraic closure, a quasi-smooth member has seven distinct linear factors
   in `A` and two distinct linear factors in `C`, producing the curve lattice
   above.  Therefore this **specific hypersurface family** cannot have
   geometric Picard rank one; the construction achieves arithmetic rank one.
   No broader four-weight decoupling assertion is made here.

## Submission data (Epoch, Method B)

    Weights: [2, 2, 7, 7]
    F = x0^7 + x0^2*x1^5 + 2*x1^7 + x2^2 + x3^2   (over ZZ/3)
