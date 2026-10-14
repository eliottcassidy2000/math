# A critical-free toy clutch with a logarithmic Jacobian obstruction

**Date:** 2026-10-14. **Status:** PROVED exact counterexample to the toy all-degree critical-point extrapolation; complete affine toy classification; exact symbolic checks. **Novelty unconfirmed.** This does not refute the fixed accessory-pair hypothesis H-J and does not produce a Keller pair.

## 1. Scope and source

The starting family is

\[
P_B(x,z)=(xz^2+B(x)z)^2+2xz+x,
\]

the specialization `V=x`, `A=2x` of the clutch family in THM-3225 and Section 6 of `euclid_rigid_hinges_and_clutch_growth_20261014.md`. The original fixed accessory pairs of THM-3225 have `deg V=16` and `deg A=8`; the present toy pair is not among them. Therefore the example below disproves an all-degree extrapolation to this toy pair, not the original H-J assertion about those accessory pairs.

The high-exponent theorem for `B_N=1+x^N` remains intact. The new point is that one exceptional affine clutch can absorb the entire critical resultant while leaving a different, exact obstruction to a polynomial Jacobian mate.

## 2. Complete half-slope calculation

For any scalar `b`, put

\[
B_b=b+x/2,\quad t=2z+1,\quad u=xt+2b.
\]

The affine coordinate change from `(x,z)` to `(x,t)` gives

\[
P_{B_b}=-2b+u+\frac{u^2(t-1)^2}{16}.
\]

Define

\[
Q=1+\frac{u(t-1)^2}{8}.
\]

Differentiating with `x,t` as independent coordinates gives

\[
(P_{B_b})_x=tQ,\qquad
(P_{B_b})_t=xQ+\frac{u^2(t-1)}8.
\]

If `t != 0` and the first derivative vanishes, then `Q=0`. This forces both `u != 0` and `t != 1`, so the second derivative cannot vanish. Consequently all critical points have `t=0`. There `u=2b` and

\[
(P_{B_b})_t=\frac{b+4}{4}x-\frac{b^2}{2}.
\]

### Theorem 2.1: exact escape and its exception

- If `b != -4`, the polynomial `P_(b+x/2)` has precisely one critical point,
  \[
  (x,z)=\left(\frac{2b^2}{b+4},-\frac12\right).
  \]
  It is Morse: in the original coordinates its Hessian determinant is `-(b+4)^2/4`.
- If `b=-4`, the polynomial
  \[
  P_*(x,z)=\bigl(xz^2+(x/2-4)z\bigr)^2+2xz+x
  \]
  has **no finite critical point**. At `t=0`, its remaining derivative is the nonzero constant `-8`.

As `b` tends to `-4`, the unique critical point goes to infinity with `z=-1/2`. This is a concrete loss of finite critical mass at an exceptional leading-coefficient resonance, not a hidden new finite critical point.

## 3. Exact resultant and all affine toy clutches

Write `d=B'` and

\[
\begin{aligned}
F_B={}&B^4(1-2d)
+x(4B^3d^2-2B^3d+4B^2d-4B^2)\\
&+x^2(-12Bd^2+12Bd+4)+8x^3(d-1)^3.
\end{aligned}
\]

Direct substitution into the inherited universal resultant gives

\[
K_B=-8x^3F_B,\qquad
\operatorname{Res}_y(R_1,R_2)=x^3K_B,
\]

where

\[
R_1=2y(y+B)(2y+B)+2x^2,
\quad R_2=x^3+x^2y+(2xB'-B)y^2(y+B).
\]

In the half-slope family this simplifies exactly to

\[
F_{b+x/2}=x\bigl[-2b^2+(b+4)x\bigr],
\quad
K_{b+x/2}=8x^4\bigl[2b^2-(b+4)x\bigr].
\]

Thus

\[
K_{x/2-4}=256x^4,
\quad
\operatorname{Res}_y(R_1,R_2)=256x^7.
\]

Every root of this resultant is boundary-supported. This independently agrees with the direct gradient proof.

### Corollary 3.1: unique critical-free affine toy clutch

Among all affine polynomials `B=b+e x`, the sole critical-free member is `B=x/2-4`.

**Proof.** If `b=0`, then `(0,-1/2)` is a critical point. Suppose `b != 0`. When `e=0`, `F_B` has nonzero constant coefficient `b^4` and degree three, with leading coefficient `-8`, so it has a root away from zero. When `e != 0,1/2`, its constant coefficient is `b^4(1-2e)` and its degree-four leading coefficient is `e^4(2e-1)`; again it has a nonzero root. The exact gradient-resultant equivalence on `x != 0` then supplies a critical point. The case `e=1/2` was completely solved above. ∎

For a genuine degree `n >= 2` clutch, the top term has degree `5n-1` in `F_B`, with coefficient `2n(2n-1) lc(B)^5`; hence `deg K_B=5n+2`. But the affine resonance invalidates a proposed argument that used this degree formula for every positive degree. At slope `1/2`, the degree drops from the generic affine value seven to five, and at intercept `-4` it drops again to four. The boundary order is then exactly four and absorbs everything.

## 4. The exceptional submersion still has no polynomial Jacobian mate

For `B=x/2-4`, set `u=xt-8` as above. The fibre at value eight factors as

\[
P_*-8=u\left(1+\frac{u(t-1)^2}{16}\right).
\]

It contains the smooth closed component

\[
C:\ xt=8\cong\mathbb G_m.
\]

On `C`, the derivatives in `(x,t)` are `P_x=t` and `P_t=x`.

### Theorem 4.1: logarithmic obstruction

There is no polynomial `H(x,z)` for which `Jac_(x,z)(P_*,H)` is a nonzero constant.

**Proof.** Suppose the Jacobian equals `c != 0`. Restrict `H` to the hyperbola by putting

\[
h(t)=H\left(8/t,(t-1)/2\right)\in k[t,t^{-1}].
\]

The chain rule and `t=2z+1` give, on `C`,

\[
c=2(tH_t-xH_x)=2t h'(t).
\]

For every Laurent polynomial, the constant coefficient of `t h'(t)` is zero. It cannot equal the nonzero constant `c/2`. ∎

Equivalently, a Jacobian mate would require the differential `(c/2) dt/t` to have a Laurent-polynomial primitive. Its nonzero logarithmic residue forbids that primitive. This is an algebraic residue obstruction, not a claim that absence of critical points implies existence of a mate.

Combining Corollary 3.1 with Theorem 4.1 shows that **every affine toy clutch has no polynomial constant-nonzero-Jacobian mate**. Exactly one of those exclusions requires an obstruction other than a finite critical point.

The same argument already forbids a mate for every half-slope member with `b != 0`: the fibre at value `-2b` contains `xt=-2b`, on which the same Euler derivative appears. When `b=0`, the explicit critical point suffices.

## 5. Mechanism for the fixed accessory problem

This example gives a genuine warning and a replacement test, not a proof or disproof of the original H-J.

1. Boundary-avoiding `B` and a large generic resultant degree do not by themselves force interior critical points. Exceptional leading terms and boundary jets must be tracked together.
2. Complete boundary support need not produce a Keller candidate: inspect the exceptional fibres and restrict the Hamiltonian derivation to each component. A nonzero logarithmic residue can obstruct a mate even when the gradient ideal is the unit ideal.
3. For the actual THM-3225 accessory pairs, `deg V=16`, `deg A=8`. The inherited resultant has degree exactly `96` for every `deg B <= 8`, because `64V^6` is the unique top term. For `deg B=n >= 10`, its unique top term is `-4AB^3J_B^2V`, of degree `5n+54`; its leading coefficient is nonzero because `2n-16 != 0`. Degree nine is a genuine resonance wall: the terms `-4AB^3J_B^2V` and `8B^3J_BV^3` tie at degree 99, with combined coefficient `16 beta^4 v^3(v-a beta)`, where `a=lc(A)`, `v=lc(V)`, `beta=lc(B)`. Thus degree 99 holds off the wall `beta=v/a`; the degree on that wall is not settled here. These statements follow from the displayed universal formula and do not assert a boundary multiplicity bound.
4. At each multiple `T` boundary root, the coefficient from THM-3225 remains independent of the higher clutch jets: for multiplicity `m >= 3` and `B(alpha)=b != 0`, the order remains `3m-1`, with leading coefficient `16m(m-1)b^5v^3/(m-2)`. The remaining all-degree difficulty is concentrated at the simple `S` root. The affine two-jet cap cannot simply be reused for arbitrary polynomial `B`.

A useful strengthened target is therefore: for each fixed accessory pair and every admissible `B`, either an interior critical point survives, **or** a directly computable period/residue on a special-fibre component excludes a polynomial Jacobian mate. The toy exception proves why these are logically different routes. The strengthened target remains OPEN for the fixed accessory pairs.

## 6. General fibrewise residue criterion

Let `P` be a polynomial over a characteristic-zero field, and let `C` be a smooth irreducible component of one of its fibres. The Hamiltonian derivation

`D_P=P_x partial_z - P_z partial_x`

preserves the ideal of `C` and restricts to a derivation on its coordinate ring. A polynomial mate with `Jac(P,H)=c != 0` must satisfy `D_P(H|C)=c`.

If `C` is isomorphic to `G_m` with coordinate `t` and the restricted derivation is `lambda t partial_t`, `lambda != 0`, no such mate exists: the constant coefficient of the Euler derivative of a Laurent polynomial is zero. More invariantly, a mate would give `d(H|C)=c omega`, where `omega(D_P)=1`. A nonzero residue of this differential, or a nonzero period over the complex numbers, forbids a regular primitive.

The logarithmic hypothesis matters. A general nonzero derivation `h(t) partial_t` on `G_m` does not automatically obstruct a mate. If it is nowhere zero, then `h(t)=lambda t^r` for an integer `r`. The component equation is solvable in Laurent polynomials exactly when `r != 1`, by

`q(t)=c t^(1-r)/(lambda(1-r)) + constant`.

For `r=1`, the required primitive is logarithmic and does not belong to the coordinate ring. Solvability of this restricted equation when `r != 1` is only the removal of this componentwise obstruction; it does not construct a global polynomial Jacobian mate. If a more general `h` has a zero on `G_m`, evaluation at that zero already obstructs a nonzero constant image.

## 7. Reproduction and limitations

The companion `toy_clutch_exception_20261014.py` verifies the universal resultant independently in formal symbols; the specialized factorization; the normal form and gradient identities; the Hessian; the affine leading coefficients; and the actual accessory-degree envelope. It also verifies the logarithmic restriction identity for generic polynomial monomials through total degree eight. The latter finite check illustrates the identity; the Laurent-polynomial proof above holds in all degrees.

No remote or live checkout was modified by this worker. No theorem number is reserved. No assertion about arbitrary polynomial maps, the truth of the planar Jacobian conjecture, or prior literature novelty is made.
