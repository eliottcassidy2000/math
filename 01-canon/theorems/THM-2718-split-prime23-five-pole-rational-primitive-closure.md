---
id: THM-2718
title: "Split prime-23 five-pole rational-primitive closure"
status: >
  PROVED + INDEPENDENTLY HOSTILE-AUDITED.  On the B_0!=0, y not
  identically zero, lambda!=0 prime-23 scale chart of the split polynomial
  exact-square-prefix degree-twenty-two even-Faber bank, no physical Keller
  trajectory exists.  THM-2713 gives five distinct poles of the chosen-sheet
  q-function, while the third Faber flux and THM-2214's rational-primitive
  lemma allow q=A_src/U to have at most two pole points on the source P1.
  The eleven odd Faber seeds, B_0=0, y=0, lambda=0, the broader split branch,
  and JC(2) remain open.
source: thm2704-hostile-audit-2026-07-28-five-pole-primitive
audit: lrc-narrow-debt-queue-2026-07-28-five-pole-primitive
depends_on:
  - THM-2129-quartic-faber-three-coefficient-boundary-classification
  - THM-2214-nonsplit-terminal-quartic-spectral-curve-closure-through-degree-ten
  - THM-2713-split-prime23-component-divisor-budget-and-perfect-power-normal-form
related:
  - THM-2411-degree-twenty-two-first-flux-pole-divisor-square-class-reduction
  - THM-2704-split-even-prime23-generic-genus-eighty-nine
---

# THM-2718 -- five old poles do not fit through one rational primitive

**PROVED + INDEPENDENTLY HOSTILE-AUDITED.**  The former namespace collision
was repaired by moving the later coordinate-first claimant to THM-2740.

## 1. Exact statement

Let `C` be the weighted prime-23 curve of THM-2713 arising from a polynomial
Keller pair on the following precise chart:

```text
split polynomial exact-square prefix;
reduced mate degree 22;
target-translated even-Faber bank E22+B_0 E14+C_0 E10+D_0 E6+E_0 E2;
all eleven odd Faber coefficients zero;
B_0!=0,                    y not identically zero;
Phi_Q=lambda!=0,           Psi_Q=W_0 in C.                 (1)
```

Choose `rho in C*` with `rho^2=B_0` and use the scaled coordinates

```text
t=rho/y,              v=u/y^2,              zeta=Z/y^3.   (2)
```

Then there is no physical polynomial Keller trajectory on `C`.

Equivalently, the response curve

```text
F2=0,                   zeta F1^4=eta t^23, eta!=0         (3)
```

may have abstract special members with rational normalization, but none can
be the trajectory of a Keller pair satisfying (1).  This statement closes
only the chart (1).  It does not include `B_0=0`, `y=0`, `lambda=0`, any odd
Faber seed, the full split degree-twenty-two branch, or `JC(2)`.

## 2. The split sheet makes the third flux rational

Use the polynomial exact-square-prefix notation

```text
P=H^2+L,
H=V_src z^2+B_src z+C_src,       L=A_src z+E_src,           (4)
```

with `V_src,A_src,B_src,C_src,E_src in C[x]` and `V_src!=0`.
On the split deck, `V_src` is a square in `C(x)`.  In fact one may choose

```text
V_src=U^2,                         U in C[x]\{0}.            (5)
```

Indeed, write a rational square root as `r/s` in lowest terms.  From
`r^2=s^2 V_src`, every irreducible factor of `s` would divide `r`, so `s` is
a unit.  Absorb its constant into `r`.

On one split sheet the depressed-quartic coordinates satisfy

```text
w=Uz+B_src/(2U),
q_phys=A_src/U,                    T=q_phys^2.               (6)
```

Here `A_src` is not identically zero: the chosen-sheet first-flux equation
has nonzero right side because `lambda!=0`, hence `q_phys` is not the zero
function.

Let `Phi_Q,Psi_Q,R_Q` be the three exact Faber observables of THM-2129.  Its
Hamiltonian identity, after the change from `z` to `w`, is

```text
J_(x,w)(P,Q)
 =(w^2+p/4) Phi_Q' + w Psi_Q' + R_Q'
 =kappa/U,                         kappa in C*.              (7)
```

The first two observables in (1) are constants.  Comparing coefficients in
`w` therefore gives

```text
R_Q'=kappa/U,                       U R_Q'=kappa.             (8)
```

No nonsplit parity is used here.  In the nonsplit theorem parity supplied
the stronger value `Phi_Q=0`; equation (8) needs only `Phi_Q'=0`, which is
still true for the split constant `Phi_Q=lambda`.

All quantities in `R_Q` are rational on the chosen split sheet: the centered
quartic coordinates `d_ctr,s_ctr,T,q_phys` lie in `C(x)` and every Faber
coefficient is constant.  Thus

```text
R_Q in C(x),
```

so (8) is an honest rational-primitive equation over the source line.

## 3. The rational-primitive pole-capacity lemma

We use the elementary lemma proved in THM-2214, Section 7.3.

> Let `U in C[x]\{0}`, `S in C(x)`, and `U S'=kappa in C*`.  Then either
>
> ```text
> U is constant and S is affine-linear,                         (9a)
> ```
>
> or, after translating `X=x-xi`,
>
> ```text
> U=u_0 X^m,             S=s_0+s_1 X^(1-m),
> m>=2,                  u_0 s_1!=0.                           (9b)
> ```

For completeness, if nonconstant `U` has degree `D` and `h` distinct roots,
then a rational primitive of `1/U` has map degree `D-h`.  Its fibre over its
value at infinity has multiplicity `D-1`; hence `D-1<=D-h`, so `h=1`.
Direct integration gives (9b).  A simple root would instead create the
forbidden logarithmic term.

Apply the lemma to `(U,S)=(U,R_Q)` in (8).  It follows that `U` is constant
or has only one distinct finite root.  Consequently the quotient

```text
q_phys=A_src/U,                    A_src in C[x],             (10)
```

has at most two distinct pole points on the projective source line: the one
possible finite root of `U` and `x=infinity`.  Cancellation in (10) can only
remove poles, never create another pole point.

## 4. A physical trajectory must dominate the prime-23 curve

The physical coefficient functions define on the generic point

```text
gamma_0:P1_x ---> C,                 x |-> [1:t:v:zeta].     (11)
```

Because `C` is projective, this rational map extends across the omitted
source points.  It cannot be constant.  To see this without a genericity
assumption, suppose its value were constant.  The generic image lies in
`h=1`.  Since `rho!=0`, equations (2) would successively make

```text
t constant => y constant;
v constant => u=v y^2 constant;
zeta constant => Z=zeta y^3 constant.                       (12)
```

The chosen-sheet reconstruction from the first flux is

```text
q_phys=chi*t^5/F1(t,v,zeta),             chi in C*,          (13)
```

so `q_phys` and `T=q_phys^2` would be constant as well.  Hence
`d_ctr=u/T`, `s_ctr=y/11`, and every input of the exact third observable
would be constant.  This would make `R_Q` constant, contradicting (8).

THM-2713 proves that `C` is geometrically integral for every parameter value
in (3).  Let `nu:X->C` be its normalization.  The nonconstant map (11) from
the normal source line therefore factors through a nonconstant, hence finite
and surjective, morphism

```text
gamma:P1_x -> X.                                           (14)
```

Riemann--Hurwitz already makes (14) impossible if `g(X)>0`.  Thus any
remaining physical trajectory would force

```text
X=P1.                                                       (15)
```

## 5. Five old pole fibres contradict the primitive lemma

THM-2713 gives on `X` five distinct old normalization points

```text
O=O_1+...+O_5
```

and, for the normalized chosen-sheet function `q_C`,

```text
div(q_C)=5N-3O,                   ord_(O_i)(q_C)=-3.        (16)
```

The physical `q_phys` in (6) differs from `q_C` only by the nonzero constant
in (13), so it has the same pole divisor.  Each fibre
`gamma^(-1)(O_i)` is nonempty, and these five fibres are pairwise disjoint.
The single source point `x=infinity` belongs to at most one of them.  Hence
four of the old points admit distinct finite preimages `p_i`.

If `e_i` is the ramification index of (14) at such a preimage, pullback of
(16) gives

```text
ord_(p_i)(q_phys)=-3e_i<0.                               (17)
```

Using (10) at the finite point `p_i`,

```text
ord_(p_i)(U)>ord_(p_i)(A_src)>=0.                        (18)
```

Thus `U` has at least four distinct finite roots.  This contradicts (9),
which permits at most one.  Therefore (14), and hence the physical Keller
trajectory, cannot exist.

## 6. What closes and what does not

The proof combines two statements on different spaces:

```text
normalization X:       q_C has five distinct old pole points;
source P1_x:           U R_Q'=kappa gives U at most one finite root. (19)
```

Surjectivity of the physical map is the bridge.  It is essential: the
five-pole divisor is not by itself an abstract rationality obstruction, and
the rational-primitive lemma is not an equation intrinsic to every abstract
member of (3).

The conclusion is exactly

```text
{physical Keller trajectories on chart (1)}=empty.          (20)
```

The following remain outside (20):

1. the `B_0=0` coefficient chart, where `rho^2=B_0` cannot define (2);
2. `y=0`, which is not represented by the signed scale `t=rho/y`;
3. `lambda=0`, where the prime-23 chosen-sheet relation degenerates;
4. every one of the eleven odd Faber seeds allowed on a split deck;
5. a proof that every abstract curve (3) has positive genus; and
6. the full split degree-twenty-two branch, `JC(2)`, or `DC(2)`.

No new computation is used.  The curve integrality and divisor input are the
proved THM-2713 package; the differential identity and rational-primitive
lemma are the proved THM-2129/THM-2214 package.  The new step is the
source-to-normalization pole-capacity comparison (14)--(18).
