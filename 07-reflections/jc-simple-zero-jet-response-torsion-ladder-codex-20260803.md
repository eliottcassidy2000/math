# Simple-zero deformations preserve the exact torsion ladder and replace the scalar bridge by a jet response

**Status:** PROVED + EXACT-SYMBOLIC CONTROLS; canonized in
[THM-3326](../01-canon/theorems/THM-3326-linear-in-z-unit-response-trichotomy-and-jet-torsion.md).
It extends the mechanism of THM-3318.

## Theorem-ready statement

Let `K` be a characteristic-zero field, `R=K[x,z]`, `r>=2`,
`lambda in K*`, and `f in K[x]`.  Put

```text
P=f(x)+lambda*x^r*z,             b=f(0),
p=P-b,                           phi=f-b=x*u(x),
D_P=P_x*partial_z-P_z*partial_x,
C_P=R/D_P(R).
```

Assume `c=f'(0)=u(0)` is nonzero.  Then `(P_x,P_z)=R`.  There is a unique
polynomial `A in K[z][x]` of `x`-degree less than `r` such that

```text
A*P_x = 1 mod x^r.
```

With

```text
B=(1-A*P_x)/(lambda*x^r),        m=A_x+B_z,
mu=[m] in C_P,                   theta=[1] in C_P,
```

the row `(A,B)` is a polynomial Bezout row and `mu` is the canonical
divergence class.  Its exact annihilators are

```text
Ann_K[P](theta)=((P-b)^(r-1)),
Ann_K[P](mu)=((P-b)^r).                                  (1)
```

There is a unique `F(T) in K[T]`, `deg F<r-1`, satisfying the truncated
interpolation law

```text
F(phi(x)) = (r-1)*u(x)/f'(x) mod x^(r-1).                (2)
```

It has `F(0)=r-1`, and the marked bridge is

```text
(P-b)*mu = -F(P-b)*theta.                                (3)
```

Consequently the abstract cyclic block does not change under any such
deformation:

```text
K[P]*mu ~= K[T]/(T^r),
K[P]*theta=(P-b)K[P]*mu ~= K[T]/(T^(r-1)).               (4)
```

In the normalized case requested in the research prompt, `f(0)=0`, simply
replace every `P-b` above by `P`.

The response polynomial has an intrinsic formal-coordinate description.  If
`psi` is the compositional inverse of `phi`, then

```text
F(T) = (r-1)*T*psi'(T)/psi(T) mod T^(r-1).               (5)
```

Thus the scalar relation from THM-3318 survives literally if and only if

```text
f(x)=b+c*x mod x^r.                                      (6)
```

If the first nonlinear term is `a*x^d`, with `2<=d<=r-1`, then

```text
F(T)=(r-1)*(1-(d-1)*a*c^(-d)*T^(d-1)+higher terms).      (7)
```

The first bridge deformation therefore records exactly the first nonlinear
jet of `f`, while `lambda` does not enter `F`.

Moreover, `F` is a complete truncated jet passport once the linear
coefficient `c` is retained.  If `alpha=psi'(0)=c^(-1)`, then

```text
T*d/dT log(psi(T)/(alpha*T)) = F(T)/(r-1)-1 mod T^(r-1). (JP)
```

Termwise integration reconstructs `psi mod T^r`, and hence `phi mod x^r`,
from `(c,F)`.  The abstract module in `(4)` forgets this jet, while its marked
canonical generators and bridge recover it.

## Proof audit

### 1. Unimodularity and the canonical row

Since

```text
P_x=f'(x)+lambda*r*x^(r-1)*z,       P_z=lambda*x^r,
```

`P_x` is a unit modulo `x^r`: its constant coefficient is `c`.  Truncating
its formal inverse gives the unique `A` above.  Hence `B` is polynomial and

```text
A*P_x+B*P_z=1.
```

Any two Bezout rows differ by `(k*P_z,-k*P_x)`.  Their divergences differ by
`-D_P(k)`, so `[A_x+B_z]` is independent of the row.

### 2. Localization computes the complete primitive fibres

Inverting `x` gives an exact coordinate presentation

```text
S=R[x^(-1)]=K[P,x,x^(-1)],
z=(P-f(x))/(lambda*x^r),
D_P=-lambda*x^r*partial_x                 (P held fixed),
ker(D_P:S->S)=K[P].                                      (8)
```

The kernel assertion follows coefficientwise for Laurent polynomials and is
where characteristic zero is used.  Define

```text
h=-A/(lambda*x^r),
Q0=x^(1-r)/(lambda*(r-1)).                               (9)
```

The row correction by `h` is

```text
(A+h*P_z, B-h*P_x)=(0, 1/(lambda*x^r)).                 (10)
```

The right side is `(Q0_z,-Q0_x)`, hence its divergence is zero.  With the
sign convention in `(8)`, this gives

```text
D_P(h)=m,                        D_P(Q0)=1.               (11)
```

Equation `(8)` then classifies all localized primitives as

```text
{q:D_P(q)=m}=h+K[P],             {q:D_P(q)=1}=Q0+K[P].  (12)
```

### 3. Pole orders prove both annihilator equalities

Write

```text
U=u(x)+lambda*x^(r-1)*z,          p=P-b=x*U.
```

Here `U(0,z)=c` and `A(0,z)=c^(-1)`.  The displayed killing primitives are

```text
p^(r-1)*Q0=U^(r-1)/(lambda*(r-1)) in R,
p^r*h=-A*U^r/lambda in R.                              (13)
```

For sharpness, write a nonzero polynomial `G(p)=p^e*H(p)` with `H(0)!=0`.
If `G(p)theta=0` or `G(p)mu=0`, equation `(12)` forces a polynomial primitive
to equal, respectively,

```text
G(p)*Q0+L(P),                    G(p)*h+L(P),
```

for some `L in K[P]`.  Their exact `x`-orders are `e+1-r` and `e-r`, with
nonzero leading coefficients because `c`, `H(0)`, and `A(0)` are nonzero.
No polynomial `L(P)` cancels a negative Laurent power.  The thresholds in
`(13)` are therefore sharp, proving `(1)`.

### 4. Truncated composition produces the deformed bridge

Substitution `T -> phi(x)` is an automorphism

```text
K[T]/(T^(r-1))  -->  K[x]/(x^(r-1))
```

because `phi'(0)=c!=0`; this proves existence and uniqueness in `(2)`.  Also,
modulo `x^(r-1)`, one has

```text
A*U = u/f'.                                             (14)
```

Therefore

```text
g=p*h+F(p)*Q0
 =(F(p)/(r-1)-A*U)/(lambda*x^(r-1))                    (15)
```

is polynomial.  Applying `D_P` and using `(11)` gives

```text
D_P(g)=p*m+F(p),
```

which is exactly `(3)`.  Evaluating `(2)` at zero gives `F(0)=r-1`.  Hence
`F` is a unit modulo `p^(r-1)`, so `(3)` and `(1)` imply `(4)`.

For `(5)`, put `x=psi(T)` in `u/f'=phi/(x*phi')`; the result is
`T*psi'/psi`.  Its logarithmic derivative is constant modulo `T^(r-1)`
exactly when `psi(T)=c^(-1)T mod T^r`, equivalently `(6)`.  Expanding the
inverse at the first nonlinear term gives `(7)`.

### 5. The jet response is a composition cocycle

The formal elasticity

```text
E_psi(T)=T*psi'(T)/psi(T)=F(T)/(r-1)
```

obeys an exact chain rule.  For composable simple-zero formal coordinates,

```text
E_(psi1 o psi2)(T)=E_psi1(psi2(T))*E_psi2(T).           (CC)
```

Thus successive nonlinear changes of the `x`-coordinate do not add bridge
errors; their response passports multiply after pullback.  Equation `(JP)`
also proves the reconstruction claim: divide the coefficient of `T^j` in
`F/(r-1)-1` by `j`, exponentiate, multiply by `c^(-1)T`, and compositionally
invert.  This recovers `f-b mod x^r`.  It is a precise source/target/sidecar
connection: the cyclic torsion block is the target, the source jet maps to it
through `(2)`, the abstract module destroys that jet, and the marked bridge
polynomial is an exact sidecar that restores it.

## Exact examples and hostile boundaries

For `b=0`, `c=1`, and `f=x+a*x^2`, the first bridges are

```text
r=2: F(T)=1,
r=3: F(T)=2*(1-a*T),
r=4: F(T)=3*(1-a*T+3*a^2*T^2).
```

For `r=5` and `f=x+a*x^3`,

```text
F(T)=4*(1-2*a*T^2).
```

At `r=3`, `f=x+a*x^2`, the exact canonical representative from the truncated
Bezout row is

```text
m=2*(-2*a^2*x+5*a+6*lambda*x*z),
P*mu=-2*(1-a*P)*theta.
```

This is the cheapest counterexample to the claim that the THM-3318 bridge
always remains scalar.  The stronger survivor is the unit-polynomial bridge
`(3)` and the unchanged annihilator ladder `(1)`.

The hypotheses have real boundaries.

- If `f'(0)=0`, both `P_x` and `P_z` vanish on `x=0`; the gradient is not
  unimodular, so the canonical Bezout divergence class is not defined.
- At `r=1`, the unit primitive is logarithmic rather than Laurent-polynomial,
  and this finite torsion ladder disappears.
- In positive characteristic the Laurent derivative has a larger kernel and
  `(r-1)^(-1)` may not exist; the characteristic-zero proof does not extend.
- These are explicit non-Keller first coordinates: `theta!=0` proves there is
  no polynomial `Q` with `Jac(P,Q)=1`.  The result is a cokernel/extension
  mechanism, not a case or counterexample of the planar Jacobian conjecture.

## Reproduction

Run

```text
python 04-computation/jc_linear_z_jet_torsion_response_trichotomy_20260803.py
```

The script checks the Bezout identity, both localized primitives, sharp lower
poles, polynomial killing powers, bridge polynomiality and sign, four
symbolic deformation families, six exact hostile jets through `r=7`, and the
multi-root controls used by THM-3326.
