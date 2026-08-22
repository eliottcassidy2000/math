---
id: THM-3326
title: "Linear-in-z unit-response trichotomy and jet-deformed torsion"
status: >
  PROVED + GENERIC-SYMBOLIC + FINITE-EXACT HOSTILE-AUDITED.  For every
  gradient-unimodular P=f(x)+g(x)z over a characteristic-zero field, the unit
  response class has a complete root-support trichotomy.  If g is constant it
  vanishes.  If g has one geometric root of multiplicity d>=2, its exact
  annihilator is ((P-f(a))^(d-1)); the canonical divergence class has
  annihilator ((P-f(a))^d), and their marked bridge has nonzero constant term.
  Together with f'(a), it exactly records the (d-1)-jet of f through the
  logarithmic derivative of its inverse formal coordinate.  If g has at least two geometric roots, the
  unit response generates a free K[P]-submodule.  The latter obstruction is
  the nonexistence of a rational primitive for dx/g, not gradient singularity.
audit: >
  The exact companion checks four symbolic deformation families, six hostile
  simple-zero jets through d=7, both sharp lower poles and killing powers,
  the marked bridge and signs, and three smooth two-root Bezout rows with
  nonzero residues.  Normal and optimized runs byte-match; the source has no
  assertion nodes or floating literals.  The rational-primitive and
  annihilator classifications are uniform proofs, not bounded inferences.
source: root/second-creative-jacobian-lrc/2026-08-03
depends_on:
  - THM-3318-hamiltonian-divergence-torsion-ladder-for-x-plus-xr-z
related:
  - THM-2063-one-fiber-linear-planar-keller-pairs
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
script: 04-computation/jc_linear_z_jet_torsion_response_trichotomy_20260803.py
output: 05-knowledge/results/jc_linear_z_jet_torsion_response_trichotomy_20260803.out
script_sha256: f1ef6ddf8c3e4fa4ab488e582eeb715335f8b1747c4f8bf918c677e6d72fee0b
output_sha256: 9c559b93894d4f1c00127ab36c6cb2e916104b58f9d73c73a275e2751cd4aa0f
hash_basis: LF-normalized bytes
---

# THM-3326 -- linear-in-`z` unit-response trichotomy and jet-deformed torsion

**PROVED + GENERIC-SYMBOLIC + FINITE-EXACT HOSTILE-AUDITED.**

## 1. Complete trichotomy

Let `K` be a characteristic-zero field, `R=K[x,z]`, and

```text
P=f(x)+g(x)z,                  0!=g in K[x],
D_P(q)=P_x q_z-P_z q_x,
C_P=R/D_P(R),                 theta=[1] in C_P.          (1)
```

Assume the gradient ideal `(P_x,P_z)` is the unit ideal.  Then exactly one of
the following occurs.

1. **Constant coefficient.**  If `g in K*`, then `theta=0`.

2. **One geometric root.**  If `g` has one root over an algebraic closure,
   then

   ```text
   g=lambda*(x-a)^d,          a in K, d>=2,              (2)
   ```

   and `f'(a)!=0`.  Put `b=f(a)` and `p=P-b`.  For the canonical divergence
   class `mu` defined below,

   ```text
   Ann_(K[P])(theta)=(p^(d-1)),
   Ann_(K[P])(mu)=(p^d).                                  (3)
   ```

   There is a unique `F(T) in K[T]`, `deg F<d-1`, with `F(0)=d-1` such that

   ```text
   p mu=-F(p)theta.                                       (4)
   ```

   Hence

   ```text
   K[P]mu ~= K[T]/(T^d),
   K[P]theta=p K[P]mu ~= K[T]/(T^(d-1)).                 (5)
   ```

3. **At least two geometric roots.**  Then

   ```text
   Ann_(K[P])(theta)=(0),                                (6)
   ```

   so `K[P]theta` is a free rank-one `K[P]`-module.

Thus the number of deleted vertical divisors is a load-bearing response
sidecar: one repeated divisor gives finite torsion, while two or more give a
torsion-free obstruction.  In every nonconstant case `theta!=0`, so `P` has
no polynomial Jacobian mate.

THM-2063 already proves that a Keller pair with one component affine in a
source variable is tame.  The new content here is the exact cokernel module,
its root-support trichotomy, and the marked jet response; this is not a new
case of `JC(2)`.

## 2. Gradient unimodularity determines the allowed roots

Over an algebraic closure, let `alpha` be a root of `g`.  On the vertical line
`x=alpha`,

```text
P_z=0,                    P_x=f'(alpha)+g'(alpha)z.       (7)
```

If `g'(alpha)!=0`, one choice of `z` makes `(7)` zero.  If both derivatives
vanish, every `z` is a common gradient zero.  Therefore

```text
(P_x,P_z)=R
 iff for every root alpha of g,
        g'(alpha)=0 and f'(alpha)!=0.                    (8)
```

In particular every root of a nonconstant `g` is repeated.  If there is only
one geometric root, `(2)` follows; the root descends to `K` because the
coefficient of `x^(d-1)` is `-lambda*d*a`.

For constant `g`, the polynomial `-x/g` satisfies

```text
D_P(-x/g)=1,                                             (9)
```

proving the first branch.

## 3. Rational-primitive obstruction for multiple roots

Suppose a nonzero `G(T) in K[T]` and `Q in R` satisfy

```text
D_P(Q)=G(P).                                             (10)
```

In the rational coordinate field `K(P)(x)`, solve

```text
z=(P-f(x))/g(x),                D_P=-g(x) partial_x      (11)
```

with `P` held fixed.  Equation `(10)` makes

```text
d/dx(-Q/G(P))=1/g(x).                                    (12)
```

The following elementary lemma is therefore decisive.

> For a nonconstant polynomial `g` over a characteristic-zero field,
> `1/g` has a rational primitive after any constant-field extension if and
> only if `g=lambda*(x-a)^d` with `d>=2`.

To prove the lemma, factor over an algebraic closure as

```text
g=lambda product_(i=1)^k (x-a_i)^(m_i),
d=sum_i m_i.                                             (13)
```

A rational derivative has zero residue, so no `m_i` can equal one.  If
`H'=1/g`, then `H` has pole order exactly `m_i-1` at `a_i`; its reduced
denominator has degree

```text
sum_i(m_i-1)=d-k.                                        (14)
```

After subtracting its value at infinity,

```text
H=c x^(1-d)+O(x^(-d)),                  c!=0,            (15)
```

so `H` has order `d-1` at infinity.  A proper rational function with
denominator degree `d-k` has infinity order at most `d-k`.  Hence
`d-1<=d-k`, forcing `k=1`.  Conversely,

```text
-1/(lambda*(d-1)*(x-a)^(d-1))                           (16)
```

is a primitive in the one-root case.  Applying the lemma to `(12)` proves
`(6)` whenever `g` has multiple geometric roots.

This obstruction can coexist with a smooth gradient.  For

```text
f=x,
g=lambda*x^r*(1+c x)^s,             r,s>=2,             (17)
```

one has `g|(g')^2` and the explicit polynomial Bezout row

```text
A=1-g'z,                  B=((g')^2/g)z^2,
A P_x+B P_z=1.                                             (18)
```

Yet `g` has two roots and `(6)` applies.  At `lambda=c=1`, `r=s=2`,

```text
g=x^2(1+x)^2,
A_x+B_z=(6+20x+20x^2)z,                                  (19)
1/g=1/x^2-2/x+1/(x+1)^2+2/(x+1).                         (20)
```

The residues `-2,+2` exhibit the logarithmic obstruction directly.

## 4. The one-root canonical row and localized fibres

Translate `x-a` to `x`, so `(2)` becomes `g=lambda*x^d`.  Write

```text
phi(x)=f(x)-b=x u(x),              u(0)=c=f'(0)!=0,
p=P-b=x U,
U=u(x)+lambda*x^(d-1)z.                                  (21)
```

Since

```text
P_x=f'(x)+lambda*d*x^(d-1)z,
P_z=lambda*x^d,                                          (22)
```

there is a unique `A in K[z][x]`, of `x`-degree below `d`, satisfying

```text
A P_x=1 mod x^d.                                         (23)
```

Define

```text
B=(1-A P_x)/(lambda*x^d),
m=A_x+B_z,
mu=[m] in C_P.                                           (24)
```

Then `(A,B)` is a polynomial Bezout row.  Any two such rows differ by
`(hP_z,-hP_x)`, whose divergence is `-D_P(h)`, so `mu` is canonical.

After inverting `x`,

```text
S=R[x^(-1)]=K[P,x,x^(-1)],
D_P=-lambda*x^d partial_x,
ker(D_P:S->S)=K[P].                                      (25)
```

Put

```text
h=-A/(lambda*x^d),
Q_0=x^(1-d)/(lambda*(d-1)).                              (26)
```

The row correction is

```text
(A+hP_z,B-hP_x)=(0,1/(lambda*x^d))=(Q_(0,z),-Q_(0,x)),  (27)
```

and therefore

```text
D_P(h)=m,                    D_P(Q_0)=1.                 (28)
```

Equation `(25)` classifies the complete localized primitive fibres as

```text
{q:D_P(q)=m}=h+K[P],          {q:D_P(q)=1}=Q_0+K[P].    (29)
```

## 5. Sharp annihilators

The displayed killing primitives are polynomial:

```text
p^(d-1)Q_0=U^(d-1)/(lambda*(d-1)),
p^d h=-A U^d/lambda.                                    (30)
```

For sharpness, write a nonzero polynomial response as
`G(p)=p^e H(p)`, `H(0)!=0`.  By `(29)`, any polynomial primitive must equal

```text
G(p)Q_0+L(P)             or             G(p)h+L(P),     (31)
```

with `L in K[P]`.  Along `x=0`, the exact orders of these two expressions are

```text
e+1-d,                               e-d,                (32)
```

and their leading coefficients are nonzero because `U(0)=c` and
`A(0)=c^(-1)`.  A polynomial `L(P)` has no negative `x`-order and cannot
cancel either leading pole.  Equations `(30)--(32)` prove both equalities in
`(3)`.

## 6. The bridge is an inverse-coordinate jet passport

Composition by `phi` is an automorphism of
`K[T]/(T^(d-1))`.  Thus there is a unique `F`, `deg F<d-1`, satisfying

```text
F(phi(x))=(d-1)u(x)/f'(x) mod x^(d-1).                  (33)
```

Modulo `x^(d-1)`, equations `(21)--(23)` give `A U=u/f'`.  Therefore

```text
k=p h+F(p)Q_0
 =(F(p)/(d-1)-A U)/(lambda*x^(d-1)) in R.                (34)
```

Applying `D_P` yields

```text
D_P(k)=p m+F(p),                                         (35)
```

which is exactly `(4)`.  Also `F(0)=d-1`, so `F` is a unit modulo
`p^(d-1)`; `(3)--(4)` imply `(5)`.

Let `psi` be the compositional inverse of `phi`.  Substituting
`x=psi(T)` in `(33)` gives the intrinsic formula

```text
F(T)=(d-1) T psi'(T)/psi(T) mod T^(d-1).                (36)
```

Thus the abstract cyclic module forgets the nonlinear jet, but the marked
bridge together with `c=f'(0)` recovers it.  Indeed

```text
T d/dT log(psi(T)/(c^(-1)T))=F(T)/(d-1)-1 mod T^(d-1), (37)
```

so `(c,F)` reconstructs `psi mod T^d` and hence `phi mod x^d`.  The bridge is
scalar exactly when `f(x)=b+c x mod x^d`.

For the first nonlinear hostile, `d=3`, `f=x+a x^2`,

```text
F(T)=2(1-aT),                    p mu=-2(1-a p)theta.    (38)
```

The response elasticity `E_psi=T psi'/psi=F/(d-1)` also obeys

```text
E_(psi_1 o psi_2)(T)=E_(psi_1)(psi_2(T)) E_(psi_2)(T),  (39)
```

so successive formal coordinate changes compose multiplicatively after
pullback.

## 7. Scope and audit

The companion checks `(23)--(35)`, both sharp lower poles, four symbolic
families, six hostile jets through `d=7`, and three multi-root rows of the
form `(17)`.  The two execution modes are identical.

The hypotheses are sharp for this proof.  If `f'(a)=0` in the one-root case,
the gradient is not unimodular and the canonical divergence class is not
defined.  At `d=1`, `dx/g` is logarithmic and the unit response is
torsion-free.  In positive characteristic the Laurent derivative has a
larger kernel and `(d-1)^(-1)` can fail.  None of these excluded first
coordinates is a Keller counterexample: the theorem diagnoses why a mate
does not exist.
