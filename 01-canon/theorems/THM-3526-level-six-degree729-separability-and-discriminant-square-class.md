---
id: THM-3526
title: "Level-six degree-729 separability gate and discriminant square class"
status: >
  PROVED + VERIFIED-EXACT + SPLIT-REPRESENTATION AUDITED.  For the fixed
  sporadic Keller map, the generic sixth x-coordinate eliminant has full
  degree 729 and is separable.  Equivalently, over a splitting field the 243
  sixth-stage cubic blocks are cubic, separable, and pairwise coprime.  At the
  lawful target (1,1,1) over F_733, recursive coefficient tuples with symbolic
  polynomial norms and literal FLINT regular matrices with 732-point
  multiplicative Fourier inversion independently recover the same
  730-coefficient polynomial.  The odd-block recursion is therefore lawful
  and gives [Delta_6]=[2 R_6].  No irreducibility or image-equation claim for
  R_6, sixth nonproperness component, all-level law, arbitrary-map theorem,
  or general Jacobian-conjecture conclusion follows.
source: codex/fixed-level-six-degree729/2026-08-16
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2582-odd-block-discriminant-tower-and-composite-jelonek-square-class
  - THM-3521-fixed-R5-finite-sheet-unit-and-next-old-L-clearing
  - THM-3525-level-five-degree243-separability-and-discriminant-square-class
related:
  - THM-3498-level-four-old-boundary-cancellation-and-degree81-discriminant-gate
  - THM-3522-fixed-keller-five-face-renewal-propagation
  - THM-3523-fixed-R6-finite-sheet-unit-and-next-old-L-clearing
scripts:
  - 04-computation/keller_level_six_degree729_recursive_tuple_probe_20260816.py
  - 04-computation/keller_level_six_degree729_fourier_flint_independent_audit_20260816.py
outputs:
  - 05-knowledge/results/keller_level_six_degree729_recursive_tuple_probe_20260816.out
  - 05-knowledge/results/keller_level_six_degree729_fourier_flint_independent_audit_20260816.out
script_sha256:
  - ad5d57b124f37b23fe9541c61bc4d919106b2db6c87a38632fabf3e7c076b8de
  - 5ae0a538af7a65f4f5371c889e4119c9b2674636b28d59ce74f882f153d606f3
output_sha256:
  - 7d152cdbf720012f1dd162bfbe603ae43691717b8fb7550a2770e3d67016eba1
  - 65df39e35ed031d0f237c79be8c504caa23d9187599df2505b94e4a78d71235a
semantic_sha256:
  - 8009c86f1c8f290829df2ba8332dc2b09929b08cbd55376f48a13acd8c2c427c
  - d3eedf7368e7b98681ca41b76529c49403a351f79f2eb858e69df95362dc3518
coefficient_sha256: 7aba23e306b00b14b8c60c34f9762ba8b35aecac111065058dfe9d4b3f1ecd51
hash_basis: LF-normalized bytes for files; ascending F_733 coefficient ledger
---

# THM-3526 -- the sixth eliminant has all 729 distinct generic roots

**PROVED + VERIFIED-EXACT + SPLIT-REPRESENTATION AUDITED.**

Retain the fixed sporadic Keller map `F:C^3->C^3`, its inverse cubic chart,
and the cleared norms

```text
R_5=L^271 N(G),              R_6=L^1699 N(R_5).       (1)
```

THM-3525 proves that the fifth x-eliminant is generically separable and

```text
[Delta_5]=[-2R_5].                                    (2)
```

The remaining hypothesis for the next THM-2582 odd-block step is generic
separability and coprimality of the 243 sixth-stage cubic blocks.  This
theorem closes that gate without expanding `R_6` or the characteristic-zero
degree-729 eliminant.

## 1. The rank-243 index algebra and sixth cubic core

Let `K=Q(a,b,c)`.  Starting at the generic target `q_0=(a,b,c)`, iterate the
proved inverse chart five times.  If `A_i` is the resulting finite etale
algebra and `q_i=(x_i,y_i,z_i)` its universal inverse point, then

```text
dim_K A_i=3^i,                    i=1,...,5.           (3)
```

At `q_5`, write

```text
L_5=L(q_5),
T_5=4-3y_5z_5,
E_6(X)=L_5X^3+T_5X-2z_5.                              (4)
```

The unnormalized sixth x-eliminant is

```text
P_6(X)=Norm_(A_5/K)(E_6(X)).                          (5)
```

It has degree at most

```text
3 dim_K A_5=3*243=729.                                (6)
```

Multiplying the conjugate monic blocks by their nonzero leading
coefficients changes the full eliminant only by a nonzero element of `K`.
It does not change its roots, degree, separability, or discriminant square
class.

## 2. A lawful good fibre over F_733

Specialize the target to

```text
(a,b,c)=(1,1,1) over F_733.                           (7)
```

The first exact route represents elements by nested cubic coefficient
tuples.  Recursive `3 by 3` adjugates construct the complete tower

```text
1 -> 3 -> 9 -> 27 -> 81 -> 243.                       (8)
```

All five graph identities `F(q_i)=q_(i-1)` hold.  At each level the leading
`L_i`, cubic derivative, y-chart denominator, and `x_i^3` have nonzero
absolute norm.  Their exact ledger is

```text
(level,dim,Norm(L),Norm(derivative),Norm(y-den),Norm(x^3))
(1,  3, 25,593,330,455)
(2,  9,663,178, 54,347)
(3, 27,172,118,287,557)
(4, 81,192, 42,634, 44)
(5,243,511, 85,465,723).                              (9)
```

The terminal leading coefficient has norm `364`.  Keeping `X` symbolic and
taking five successive polynomial-valued relative determinants gives the
exact degree ledger

```text
3 -> 9 -> 27 -> 81 -> 243 -> 729.                    (10)
```

For the resulting polynomial in `F_733[X]`, FLINT gives

```text
deg P_6=729,
gcd(P_6,P_6')=1.                                     (11)
```

Its 42 irreducible factors all occur with exponent one.  Seven separately
evaluated recursive norms agree with Horner evaluation; at `X=730`, the
literal flattened `243 by 243` determinant is `74`.  Flipping the constant
term from `-2z_5` to `+2z_5` negates the value at zero, as the rank 243 is
odd, and deleting the leading term fails at `X=730`.

## 3. Independent regular-matrix/Fourier reconstruction

The second route shares neither nested coefficient tuples, recursive
adjugate inversion, polynomial-valued relative norms, nor additive
interpolation with the first.  Every algebra element is represented by its
literal FLINT regular multiplication matrix.  It independently rebuilds the
five inverse graphs and obtains the leading/derivative norm pairs

```text
(25,593),(663,178),(172,118),(192,42),(511,85),       (12)
```

matching the corresponding columns of (9).

It then evaluates the literal `243 by 243` determinant in (5) at all `732`
nonzero elements of `F_733`.  Since

```text
deg P_6<=729<732=|F_733^*|,                           (13)
```

multiplicative Fourier inversion is alias-free.  Moreover
`1/732=-1 mod 733`, so for `0<=k<732`,

```text
c_k=-sum_(x in F_733^*) P_6(x)x^(-k).                (14)
```

The recovered coefficients in degrees `730` and `731` are both zero, the
degree-729 coefficient is `364`, and the unused `X=0` determinant equals the
recovered constant coefficient.  Direct Horner evaluation reproduces all
732 determinant values.  The sample ledger

```text
(X,P_6(X))=(1,696),(2,491),(729,331),
             (730,74),(731,17),(732,607)              (15)
```

agrees with the recursive route.  Most decisively, both representations give
the same ascending 730-coefficient digest

```text
7aba23e306b00b14b8c60c34f9762ba8b35aecac111065058dfe9d4b3f1ecd51. (16)
```

The regular-matrix route independently gives derivative gcd one and the
same squarefree factor-degree ledger.

## 4. Good reduction proves the generic gate

Every denominator and cubic derivative used in (7)--(9) is a unit modulo
733.  Thus (11) is a lawful specialization of the generic construction (5),
not a polynomial created by a singular chart collapse.  Its nonzero leading
coefficient and nonzero discriminant prove that the corresponding
characteristic-zero rational functions are not identically zero.

Over an algebraic closure of `K`, the finite etale algebra `A_5` splits into
243 points.  Equation (5) becomes the product of their 243 conjugate cubic
blocks.  Full degree and squarefreeness imply simultaneously that

```text
every block has degree 3;
every block is separable;
distinct blocks are pairwise coprime.                 (17)
```

These are open conditions.  One lawful good fibre therefore proves the
generic characteristic-zero statement.  The factor degrees seen over
`F_733` are not promoted to a generic factorization.

## 5. The sixth discriminant square class

THM-2582 gives, for an odd block tower,

```text
[Delta_(r+1)]=[N(Delta_r)][Delta_1],
[Delta_1]=[-L].                                       (18)
```

Apply (18) using (1)--(2) and the cubic degree of `N`:

```text
[Delta_6]
 =[N(-2R_5)][-L]
 =[(-2)^3N(R_5)][-L]
 =[8LN(R_5)]
 =[8R_6/L^1698]
 =[2R_6].                                             (19)
```

Here `L^1698=(L^849)^2`, and `8/2=4` is a square.  Both signs in the first
line are load-bearing: their product makes the level-six class positive.
Equation (19) is an equality in `K^*/K^{*2}`.  It implies neither
irreducibility nor squarefreeness of `R_6` itself.

## 6. Exact boundary and reproduction

This theorem proves the fixed-map generic degree-729/separability gate,
pairwise coprimality of the 243 cubic blocks, and the square class (19).  It
does **not** prove:

- that `R_6` is irreducible, squarefree, or an image equation;
- a sixth nonproperness component or a formula for `S_(F^6)`;
- finite-sheet nonvanishing at every later norm rung;
- an unconditional all-level discriminant induction (THM-3528 subsequently
  proves raw polynomial complete packets, not separability);
- an arbitrary-map theorem, classification of Keller maps, `JC(2)`,
  `DC(2)`, or a general Jacobian-conjecture counterexample.

Run from the repository root:

```text
python -B 04-computation/keller_level_six_degree729_recursive_tuple_probe_20260816.py
python -B -O 04-computation/keller_level_six_degree729_recursive_tuple_probe_20260816.py
python -B 04-computation/keller_level_six_degree729_fourier_flint_independent_audit_20260816.py
python -B -O 04-computation/keller_level_six_degree729_fourier_flint_independent_audit_20260816.py
```

Each normal/optimized pair matches its stored transcript line-for-line.
**QED.**
