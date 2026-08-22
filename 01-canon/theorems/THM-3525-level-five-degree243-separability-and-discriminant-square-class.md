---
id: THM-3525
title: "Level-five degree-243 separability gate and discriminant square class"
status: >
  PROVED + VERIFIED-EXACT + SPLIT-REPRESENTATION AUDITED.  For the fixed
  sporadic Keller map, the generic fifth x-coordinate eliminant has full
  degree 243 and is separable.  Equivalently, over a splitting field the 81
  fifth-stage cubic blocks are cubic, separable, and pairwise coprime.  An
  exact good fibre at (1,1,1) over F_251 is independently reconstructed by
  nested coefficient algebras with Newton interpolation and by FLINT regular
  matrices with multiplicative Fourier inversion; both give the same
  244-coefficient hash.  The odd-block recursion is therefore lawful and
  gives [Delta_5]=[-2 R_5].  No irreducibility or image-equation claim for
  R_5, fifth nonproperness component, all-level law, arbitrary-map theorem,
  or general Jacobian-conjecture conclusion follows.
source: codex/fixed-level-five-degree243/2026-08-16
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2582-odd-block-discriminant-tower-and-composite-jelonek-square-class
  - THM-3498-level-four-old-boundary-cancellation-and-degree81-discriminant-gate
  - THM-3506-fixed-keller-five-face-norm-transform-and-271-99-boundary
related:
  - THM-3504-level-four-sporadic-keller-image-prime-and-four-component-nonproperness
  - THM-3521-fixed-R5-finite-sheet-unit-and-next-old-L-clearing
  - THM-3522-fixed-keller-five-face-renewal-propagation
  - THM-3523-fixed-R6-finite-sheet-unit-and-next-old-L-clearing
scripts:
  - 04-computation/keller_level_five_degree243_finite_field_probe_independent_20260816.py
  - 04-computation/keller_level_five_degree243_fourier_flint_independent_audit_20260816.py
outputs:
  - 05-knowledge/results/keller_level_five_degree243_finite_field_probe_independent_20260816.out
  - 05-knowledge/results/keller_level_five_degree243_fourier_flint_independent_audit_20260816.out
script_sha256:
  - 25e5bdcb299bd3291125a2bc5a45208f245ea53af39c37df36f8f0946c16115d
  - 044f535a19a7142da23b541c8dddfd9bccfeae20d0e895bc0e0384905f12423a
output_sha256:
  - 91cf818499aeb3e0ccbb9b71ed7ead06983fbbac9c1fb4e6dd5904e5312d703e
  - a5cabf0f78538ff7719ae203f6be3916a6b918ed13542ce9c38d125b361366a5
semantic_sha256:
  - 1397fb3a3173e8dfbe867f4fe4c2d527ef33da3c94fb60fcb227f8c60b9d15b7
  - 0f44c226329aedf3a3c232dc44fcd228cfce9587483d92f47090baac2f2be7ea
coefficient_sha256: 912f32ec0b9b375d9db2ba71d7fdf224456c86862871ae0f8f92bac5038f00ab
hash_basis: LF-normalized bytes for files; ascending F_251 coefficient ledger
---

# THM-3525 -- the fifth eliminant has all 243 distinct generic roots

**PROVED + VERIFIED-EXACT + SPLIT-REPRESENTATION AUDITED.**

Retain the fixed sporadic Keller map `F:C^3->C^3`, its inverse cubic chart,
and the cleared norms

```text
G=L^43 N(J),                 R_5=L^271 N(G).             (1)
```

THM-3498 proves that the fourth x-eliminant is generically separable and

```text
[Delta_4]=[2G].                                      (2)
```

The remaining hypothesis for the next THM-2582 odd-block step was generic
separability and coprimality of the 81 fifth-stage cubic blocks.  This theorem
closes that gate without expanding `R_5` or the characteristic-zero degree-243
eliminant.

## 1. The rank-81 index algebra and fifth cubic core

Let `K=Q(a,b,c)`.  Starting at the generic target `q_0=(a,b,c)`, iterate the
proved inverse chart four times.  If `A_i` is the resulting finite etale
algebra and `q_i=(x_i,y_i,z_i)` its universal inverse point, then

```text
dim_K A_i=3^i,                    i=1,2,3,4.             (3)
```

Write, at `q_4`,

```text
L_4=L(q_4),
T_4=4-3y_4z_4,
E_5(X)=L_4 X^3+T_4 X-2z_4.                              (4)
```

The monic inverse cubic is obtained from (4) by dividing by `L_4`.  The
unnormalized fifth x-eliminant is

```text
P_5(X)=Norm_(A_4/K)(E_5(X)).                            (5)
```

It has degree at most

```text
3 dim_K A_4=3*81=243.                                  (6)
```

Multiplying each block by its nonzero leading coefficient changes the full
eliminant only by a nonzero scalar in `K`; it does not change its roots,
degree, separability, or discriminant square class.

## 2. A lawful good fibre over F_251

Specialize the target to

```text
(a,b,c)=(1,1,1) over F_251.                            (7)
```

The first exact route represents every algebra element by nested coefficient
tuples.  It constructs the complete tower

```text
1 -> 3 -> 9 -> 27 -> 81,                               (8)
```

checks `F(q_i)=q_(i-1)` at all four stages, and certifies every inverse-chart
denominator, leading `L_i`, and cubic derivative as a unit.  It evaluates
the absolute norm in (5) by four successive relative cubic determinants at
the `244` points `X=0,...,243`, then uses exact consecutive-point Newton
interpolation.

The result satisfies

```text
deg P_5=243,
gcd(P_5,P_5')=1 in F_251[X].                            (9)
```

Literal `81 by 81` multiplication determinants independently agree with the
recursive norm at the held-out points

```text
(X,P_5(X))=(0,179),(244,143),(249,220).                (10)
```

Two hostile controls are exact.  Interpolation from only the first `243`
nodes has degree at most `242` and fails at `X=243`.  Replacing the constant
term `-2z_4` by `+2z_4` changes the value at zero from `179` to

```text
-179=72 mod 251,                                       (11)
```

as required by the odd dimension `81`.

## 3. Independent regular-matrix/Fourier reconstruction

The second route shares no tuple representation, Gaussian elimination, or
additive interpolation with the first.  Every algebra element is its FLINT
regular multiplication matrix.  It evaluates the literal `81 by 81`
determinant at all `250` nonzero elements of `F_251`.

Since `deg P_5<=243<250`, multiplicative Fourier inversion is alias-free:

```text
c_k=-sum_(x in F_251^*) P_5(x)x^(-k),       0<=k<250.  (12)
```

The recovered coefficients in degrees `244,...,249` are all zero, the
degree-243 coefficient is nonzero, and the unused `X=0` determinant equals
the recovered constant coefficient.  Direct Horner evaluation reproduces
all `250` determinant values.  Both routes give the same ascending
244-coefficient digest

```text
912f32ec0b9b375d9db2ba71d7fdf224456c86862871ae0f8f92bac5038f00ab. (13)
```

FLINT independently gives the squarefree factor-degree multiset

```text
1^2, 2^6, 3^4, 4^4, 6^4, 9, 12^3, 24, 36^3,          (14)
```

where every factor exponent is one and the degrees sum to `243`.  The fibre
is not irreducible; irreducibility is neither needed nor inferred.

## 4. Good reduction proves the generic gate

All denominators and derivatives used in (7)--(8) are units modulo `251`.
Thus (9) is a lawful specialization of the generic construction (5), not a
polynomial formed after a singular chart collapse.  Its nonzero leading
coefficient and nonzero discriminant prove that the corresponding
characteristic-zero rational functions are not identically zero.

Over an algebraic closure of `K`, the finite etale algebra `A_4` splits into
`81` points.  Equation (5) becomes the product of their `81` conjugate cubic
blocks.  Full degree and squarefreeness imply simultaneously that

```text
every block has degree 3;
every block is separable;
distinct blocks are pairwise coprime.                  (15)
```

These are open conditions, so one lawful good fibre proves the generic
characteristic-zero statement.  This is exactly the missing applicability
gate in THM-2582; it does not assert that (14) is the generic factorization.

## 5. The fifth discriminant square class

THM-2582 gives, for the odd block tower,

```text
[Delta_(r+1)]=[N(Delta_r)][Delta_1],
[Delta_1]=[-L].                                        (16)
```

Apply (16) using (2), the cubic degree of `N`, and (1):

```text
[Delta_5]
 =[N(2G)][-L]
 =[2^3 N(G)][-L]
 =[-8 L N(G)]
 =[-8 R_5/L^270]
 =[-2R_5].                                             (17)
```

Here `L^270=(L^135)^2`, and `8/2=4` is a square.  Both the minus sign and the
constant class `[2]` are therefore retained.  Equation (17) is an equality
in `K^*/K^{*2}`.  It does not imply that `R_5` is irreducible or squarefree.

## 6. Exact boundary and reproduction

This theorem proves the fixed-map generic degree-243/separability gate,
pairwise coprimality of the 81 cubic blocks, and the square class (17).  It
does **not** prove:

- that `R_5` is irreducible, squarefree, or the equation of `F(V(G))`;
- a fifth nonproperness component or `S_(F^5)=V(LHJGR_5)`;
- a bounded-degree squarefree specialization of `R_5` itself;
- a degree-729 gate, all-level induction, arbitrary-map theorem, or
  counterexample to the Jacobian conjecture.

Run from the repository root:

```text
python -B 04-computation/keller_level_five_degree243_finite_field_probe_independent_20260816.py
python -B -O 04-computation/keller_level_five_degree243_finite_field_probe_independent_20260816.py
python -B 04-computation/keller_level_five_degree243_fourier_flint_independent_audit_20260816.py
python -B -O 04-computation/keller_level_five_degree243_fourier_flint_independent_audit_20260816.py
```

Each normal/optimized pair byte-matches its stored transcript. **QED.**
