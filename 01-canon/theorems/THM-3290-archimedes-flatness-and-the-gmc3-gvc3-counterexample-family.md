---
id: THM-3290
title: "Archimedes flatness and the GMC(3)/GVC(3) counterexample family"
status: >
  PROVED + VERIFIED-EXACT; INDEPENDENT IMMUTABLE AUDIT PENDING.
  A supplied three-variable Gaussian-moment counterexample is proved, not
  merely checked, and generalized.  With `rho=t^2+xy` the operator
  `Delta=4d_x d_y+d_t^2` is exactly the Laplacian of `rho`, so `Delta^k` on
  degree-`2k` forms is the Gaussian moment functional.  On the sphere `rho=1`
  the configuration collapses to ONE complex variable via
  `xC_nu = rho^(2nu+1) - t^2 A^(2nu)`, and the spherical average becomes a
  single Taylor coefficient
  `<x^(2delta) R_nu^k> = C(k-nu, k-delta) * 2^(2k)(2k)!/(4k+1)!!`.
  The vanishing is a collision of two orders: the `t`-antiderivative is FLAT
  to order `2k+1` at `a=1` while the prefactor has degree only `k-nu`.  This
  proves `<R_nu^k>=0` exactly for `k>=nu` (sharp, with sub-threshold value
  `(-1)^k C(nu-1,k) F_(2k)(1)`), yields an infinite family `P_nu=R_nu^nu`,
  `Q_nu=x^(2nu)`, `D=Delta^((4nu+2)nu)`, and reproduces the supplied closed
  form `2^(8m+1)(6m+1)!(2m)!(12m+3)!!/(4m+1)!!` at `nu=1`.  The mechanism
  consumes the Archimedes coordinate `t`, which is why it has no two-variable
  analogue.
audit: >
  The exact companion compares THREE independent routes: direct three-variable
  polynomial algebra with the operator applied literally, the proved closed
  form, and an explicit `t`-expansion of the `z`-constant term.  It verifies
  the supplied object verbatim (degree 12, 23 terms, key identity,
  `Delta^(6m)(P^m)=0` and the closed form for `m=1..5`), the dictionary
  `L(rho^k)=(2k+1)!!`, twelve cross-route agreements, the sharp threshold and
  sub-threshold closed form for `nu<=5`, the family for `nu<=5` and `m<=6`,
  the underlying hypergeometric identity to `m=40`, that the vanishing is
  created by the `t`-average and not termwise, and a hostile control showing
  the flatness/degree collision is sharp.  Normal and `-O` replay are
  byte-identical.  A concurrent session (boxeph, 2026-08-03) independently
  verified the same finite instances by full sympy expansion and recorded the
  two all-`m` statements as open; this theorem proves both, and the two
  implementations agree on every shared instance.  Independent immutable audit
  of the proof itself is pending.
source: death-star-gvc3-counterexample-2026-08-03
depends_on: []
related:
  - 05-knowledge/results/gvc3-delta6-counterexample-verification-boxeph.md
  - THM-2022-gmc2-frobenius-lowest-balanced-face
  - THM-1435-zhao-vc-witness-transport-machinery-and-the-closed-shortcut
  - THM-1490-the-gaussian-moment-counterexample-verified-proved-shortened-and-obstructed
external:
  - "The nu=1 member was supplied to this session by the owner.  Its stated
    arXiv identifier (2606.17854) is WRONG: that identifier is Ajwani--Gajjala
    --Raman--Ray, *Counterexamples to Wegner's Conjecture for Rectangles*
    (cs.CG), which contains none of this material.  Provenance of the nu=1
    object is therefore UNRESOLVED; see section 6.  Independently, GMC is known
    false in every dimension >= 3 -- see the CORE-PAPERS entry for Long,
    *Small Counterexamples to the Gaussian Moments Conjecture*
    (arXiv:2607.18186), whose smallest three-variable example has degree 4 and
    five terms and is therefore SMALLER than the nu=1 object proved here."
script: 04-computation/archimedes_flatness_gmc3_family_thm3290.py
output: 05-knowledge/results/archimedes_flatness_gmc3_family_thm3290.out
script_sha256: 69bb5dcf03ca1acf15542fe214d4916d26fc8f8f04f104becc12a300798c1ec0
output_sha256: bdf22d1896c8d994177fe547050298ea9f8e6a8e9d5cca5f12e96e9c19e385b9
hash_basis: LF-normalized bytes
---

# THM-3290 -- Archimedes flatness and the GMC(3)/GVC(3) counterexample family

**PROVED + VERIFIED-EXACT; INDEPENDENT IMMUTABLE AUDIT PENDING.**

A three-variable counterexample to the Generalized Vanishing Conjecture was
supplied as a bare list of polynomials plus a closed form.  Checking it is
finite work; this theorem instead identifies *why* it works, proves it for all
`m`, and turns the mechanism into an infinite family with a closed form and a
sharp threshold.

## 1. The dictionary: this is a Gaussian moment

Work in `Q[x,y,t]` and put

```text
rho=t^2+xy,          Delta=4 d_x d_y + d_t^2.                     (1)
```

`Delta` is exactly the Laplacian of the quadratic form `rho`: substituting
`x=u+iv`, `y=u-iv` gives `Delta=d_u^2+d_v^2+d_t^2` and `rho=u^2+v^2+t^2`.
Define the moment functional

```text
L(f)=(exp(Delta/2) f)(0)= sum_(k>=0) Delta^k f(0) / (2^k k!).      (2)
```

For `f` homogeneous of degree `2k` this collapses to `L(f)=Delta^k f/(2^k k!)`,
and `L(f)=E[f(G)]` for a standard Gaussian `G` in `(u,v,t)`.  The companion
checks `L(rho^k)=(2k+1)!!`, the `n=3` radial moments.  Hence

```text
Delta^(6m)(P^m)=0   <=>   L(P^m)=0   <=>   E[P(G)^m]=0.           (3)
```

**A GVC statement about `Delta^6` in three variables is a GMC(3) statement.**
That is the first thing to record: the supplied object belongs to the
repository's Gaussian-moments lane, not to a separate operator lane.

## 2. Radial-spherical split and the Archimedes coordinate

For `f` homogeneous of degree `2k`,

```text
L(f)=(2k+1)!! * <f>,                                              (4)
```

where `<f>` is the mean of `f` over the unit sphere `u^2+v^2+t^2=1`.  Put
`z=u+iv`, so on the sphere `z zbar = 1-t^2`.  By **Archimedes' theorem** the
uniform measure on `S^2` is the product of the uniform measure in `t` on
`[-1,1]` with the uniform measure in the argument of `z`.  Substituting
`zbar=(1-t^2)/z` turns any polynomial into a Laurent polynomial in `z` over
`Q[t]`, and the argument average extracts its `z`-constant term:

```text
<f> = (1/2) int_(-1)^(1) CT_z(f) dt.                              (5)
```

## 3. The collapse to one variable

Let

```text
A=rho+x^2,      x C_nu = rho^(2nu+1) - t^2 A^(2nu),               (6)
```

for `nu>=1`.  The right side of `(6)` is divisible by `x`: modulo `x` both
`rho` and `A` reduce to `t^2`, so it reduces to `t^(4nu+2)-t^2 t^(4nu)=0`.
Hence `C_nu` is a polynomial, homogeneous of degree `4nu+1`.  Set

```text
R_nu = A C_nu^2,      deg R_nu = 8nu+4.                           (7)
```

At `nu=1`, `(6)` is the supplied identity `xC=rho^3-t^2A^2` and `R_1` is the
supplied `P=AC^2` of degree 12 with 23 terms; the companion checks that the
supplied `C=y rho^2-2x t^2 rho-x^3 t^2` is literally `(rho^3-t^2A^2)/x`.

On the sphere `rho=1` we have `x=z`, `A=1+z^2`, and `(6)` reads
`zC_nu=1-t^2A^(2nu)`, so

```text
R_nu = A (1-t^2 A^(2nu))^2 / z^2.                                 (8)
```

Everything is now a function of `w:=z^2` alone.  Writing `a:=1+w=A`,

```text
CT_z( x^(2delta) R_nu^k ) = [w^(k-delta)] ( a^k (1-t^2 a^(2nu))^(2k) ).  (9)
```

## 4. The flatness collision

Substituting `u=a^nu t` in `(5)` and `(9)`,

```text
int_0^1 (1-a^(2nu) t^2)^(2k) dt = a^(-nu) F_(2k)(a^nu),          (10)

F_n(s) := int_0^s (1-u^2)^n du,                                  (11)
```

so that

```text
<x^(2delta) R_nu^k> = [w^(k-delta)] ( a^(k-nu) F_(2k)(a^nu) ).   (12)
```

Write `H(a):=F_(2k)(a^nu)`.  Then

```text
H'(a)=nu a^(nu-1) (1-a^(2nu))^(2k),                              (13)
```

which vanishes to order `2k` at `a=1`.  Hence

```text
H(a)-H(1)  vanishes to order 2k+1 at a=1.                        (14)
```

Since `[w^j]g(1+w)=g^(j)(1)/j!`, split `(12)` as
`H(1) a^(k-nu) + a^(k-nu)(H(a)-H(1))`.  By Leibniz every derivative of order
`<= k-delta <= k <= 2k` of the second summand vanishes at `a=1`, using `(14)`.
Therefore only the first summand survives and

```text
<x^(2delta) R_nu^k> = H(1) * C(k-nu, k-delta)
                    = C(k-nu, k-delta) * 2^(2k)(2k)!/(4k+1)!!,   (15)
```

with `C(.,.)` the generalized binomial coefficient and
`F_(2k)(1)=int_0^1(1-u^2)^(2k)du=2^(2k)(2k)!/(4k+1)!!` the Beta integral.
**Equation `(15)` is the whole theorem.**

**Two orders in collision.**  `C(k-nu,k-delta)` vanishes precisely when the top
`k-nu` is a nonnegative integer strictly below the bottom `k-delta`.  The
vanishing is therefore *not* a cancellation among terms: it is the statement
that a polynomial prefactor of degree `k-nu` is annihilated by a derivative of
strictly higher order `k-delta`, because the `t`-antiderivative `F` is flat
there.  The companion confirms directly that the `z`-constant term is a nonzero
polynomial in `t` whose `t`-average vanishes.

## 5. Consequences

**(a) Sharp threshold.**  Taking `delta=0` in `(15)`:

```text
<R_nu^k> = 0   <=>   k >= nu,                                    (16)
```

and below threshold the exact value is

```text
<R_nu^k> = (-1)^k C(nu-1,k) * 2^(2k)(2k)!/(4k+1)!!   (1<=k<nu).  (17)
```

Verified against direct polynomial algebra: `<R_2^1>=-8/15`,
`<R_3^1>=-16/15`, `<R_3^2>=128/315`.

**(b) The supplied object, proved for all `m`.**  At `nu=1` every `k>=1`
satisfies `(16)`, so `L(P^m)=0` for all `m>=1`, i.e. `Delta^(6m)(P^m)=0`.
With `delta=1`, `C(m-1,m-1)=1`, so

```text
L(x^2 P^m) = (12m+3)!! * 2^(2m)(2m)!/(4m+1)!!,                   (18)

Delta^(6m+1)(x^2 P^m) = 2^(6m+1)(6m+1)! * L(x^2 P^m)
   = 2^(8m+1)(6m+1)!(2m)!(12m+3)!!/(4m+1)!!.                     (19)
```

`(19)` is exactly the supplied closed form, now proved rather than sampled.
It is an `m>=1` statement: at `m=0` the formula would read `6` while
`Delta(x^2)=0`, because `x` is an isotropic linear form and `x^2` is harmonic.

**(c) An infinite family.**  `R_nu` itself is not a counterexample for `nu>=2`,
by `(17)`.  Its `nu`-th power is: put

```text
P_nu = R_nu^nu,   Q_nu = x^(2nu),   D_nu = Delta^((4nu+2)nu).     (20)
```

Every exponent occurring is `nu m >= nu`, so `(16)` gives `D_nu^m(P_nu^m)=0`
for all `m>=1`, while `(15)` with `delta=nu` gives
`<Q_nu P_nu^m> = 2^(2 nu m)(2 nu m)!/(4 nu m+1)!! != 0`, hence
`Delta^nu D_nu^m(Q_nu P_nu^m) != 0` and a fortiori `D_nu^m(Q_nu P_nu^m) != 0`.
The first members are

```text
nu=1: deg P=12,  Q=x^2,  D=Delta^6;
nu=2: deg P=40,  Q=x^4,  D=Delta^20;
nu=3: deg P=84,  Q=x^6,  D=Delta^42.                             (21)
```

So GVC fails in three variables for `Delta^6`, `Delta^20`, `Delta^42`, ... .

## 6. Scope, and what is emphatically NOT proved

- **No Jacobian consequence.**  Zhao's equivalence with the Jacobian Conjecture
  is for the plain Laplacian `Delta`, not its powers.  Nothing here touches
  `JC`, `JC(2)`, `DC(2)`, or THM-1435's VC-witness dimension bracket
  `5 <= vcwd <= ~20`, which concerns `Delta` itself with homogeneous
  Hessian-nilpotent `P`.  The `P` here is homogeneous, but the hypothesis in
  play is `Delta^(6m)(P^m)=0`, not `Delta^m(P^m)=0`.
- **No GMC(2) tension.**  `GMC(2)` is PROVED (THM-2022 via NC2).  The mechanism
  above cannot descend: it consumes the Archimedes coordinate `t`, whose
  antiderivative supplies the order-`2k+1` flatness in `(14)`.  In two
  variables the sphere is a circle, `(5)` has no `t`-integral, and `<f>` is the
  bare `z`-constant term -- there is no integration left to be flat.  **The
  dimension boundary of GMC is the appearance of an Archimedes coordinate.**
  This is a structural reading of the known `2`/`3` boundary, not a new proof
  of either side.
- **The threshold is about the `x`-direction, not about witness degree.**
  `(15)` covers only witnesses `Q=x^(2delta)`, and there the gate is
  `delta>=nu`.  It must NOT be read as "a witness of degree `<2nu` cannot
  exist".  A hostile monomial sweep at `nu=2` (where `P_2=R_2^2` has degree 40)
  finds that `x^2` indeed fails, exactly as `(15)` predicts, but `y^2`, `xy`
  and `t^2` all give `L(Q P_2) != 0`.  So the minimal witness degree is `2` for
  both `nu=1` and `nu=2` and does not grow with `nu`; the earlier guess that it
  would was refuted by this probe before it was recorded anywhere.  Witnesses
  odd in `t` always vanish, since `P_nu` is even in `t`.
- **No minimality claim.**  Long's `arXiv:2607.18186` reports a three-variable
  GMC counterexample of degree `4` with five terms; the `nu=1` object here has
  degree `12` and 23 terms and is therefore not minimal.  The contribution is
  the proof, the mechanism, the closed form `(15)`, and the family, not size.
- **Provenance is unresolved.**  The identifier supplied with the `nu=1`
  object, `arXiv:2606.17854`, is Ajwani--Gajjala--Raman--Ray, *Counterexamples
  to Wegner's Conjecture for Rectangles*, and contains none of this material.
  No priority claim is made or available for the `nu=1` object.  What this file
  claims as its own is sections 2--5.

## 7. Exact evidence

Run

```text
python 04-computation/archimedes_flatness_gmc3_family_thm3290.py
python -O 04-computation/archimedes_flatness_gmc3_family_thm3290.py
```

and compare LF-normalized bytes with the declared output.  Three independent
routes are cross-checked: literal polynomial algebra with the operator applied
term by term, the closed form `(15)`, and an explicit `t`-expansion of `(9)`.
Exact integer and rational arithmetic only; no floating point, random sampling,
imported executable, or assertion-sensitive test.

**QED.**
