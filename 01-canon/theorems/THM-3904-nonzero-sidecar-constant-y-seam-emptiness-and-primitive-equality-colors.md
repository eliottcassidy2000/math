---
id: THM-3904
title: "Nonzero-sidecar constant-y seam emptiness and primitive equality colors"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In the THM-3881
  residual, the f!=0 lane with T,f in k[x] is empty, without using the
  source address.  On every positive common-degree seam, removing the gcd of
  the two leading colors leaves, up to units and interchange, one square and
  a times one square.  The latter is a necessary passport only: THM-3902 has
  positive two-jet lifts, so positive equality and JC(2) remain OPEN.
source: root / post-THM-3899 constant-seam and primitive-color audit, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS.  The 3,103-gate primary expands the full
  residual, checks every strict-degree row, primewise gcd allocation, both
  first-response boundaries, the exact even-quartic discriminant, and positive
  shells.  A separate 734-gate implementation reconstructs the residual from
  THM-3881, checks 513 valuation allocations and 2,352 nonzero linear pairs
  over F_7, and retains hostile positive-degree controls.  Normal, optimized,
  and frozen streams match for both paths.
depends_on:
  - THM-3881-cusp-ideal-residual-transport-rank-two-matrix-factorization
  - THM-3899-nonzero-sidecar-y-degree-tariff-and-equianharmonic-equality-colors
related:
  - THM-3897-f-zero-residual-all-degree-global-emptiness
  - THM-3901-nonzero-sidecar-osculating-response-fan
  - THM-3902-nonzero-sidecar-equality-color-two-jet-response
script: 04-computation/jc2_nonzero_sidecar_constant_y_seam_and_primitive_colors_thm3904.py
output: 05-knowledge/results/jc2_nonzero_sidecar_constant_y_seam_and_primitive_colors_thm3904.out
script_sha256: 0bea345474e7c3149cb0aab7c3dece7b4570131dd3be33514a2470a20c7bead9
output_sha256: 8e2b5c38a67e96dd8cc2331616bcbc7717cf72f68d848d99f8abaa0b972600b8
semantic_sha256: 5739aa3e6870c6ece7e10d5a15109e9c25c6a41ebecdc4567c950e8f1a47caec
independent_audit_script: 04-computation/jc2_nonzero_sidecar_constant_y_seam_independent_audit_thm3904.py
independent_audit_output: 05-knowledge/results/jc2_nonzero_sidecar_constant_y_seam_independent_audit_thm3904.out
independent_audit_script_sha256: b8f04f650f51e537262c27ee8b8fe773951ce1773e6bdc5db9a0eaddd02b079d
independent_audit_output_sha256: 622b9ce979a6168bc92056cea52cbc542198c97c455125a90e57df458d93871b
independent_audit_semantic_sha256: b79a8f3123e1b2675d30da2620973e2482a341c1a4cf32948c0ca55448ed5e3d
hash_basis: raw LF bytes
---

# THM-3904 -- the constant-y nonzero-sidecar seam is empty

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `k` of characteristic zero.  In `D=k[x,y]` put

```text
a=x+1,                       L=9x+4,
F=15x^2+15x+4,              K=y^2-F,
P=aL^2,                     r=aT+Kf,
A=KT+aPf,                   B=Pf^2-T^2.                  (1)
```

The THM-3881 residual is

```text
S(T,f)=L^4+2(3A+3P+r^2)L^2f+(8A+6P+3r^2)B.              (2)
```

The new closure and the retained positive-degree passport are as follows.

1. If `T=u(x)`, `f=v(x)!=0`, and `S(T,f)=G^2` in `D`, then a contradiction
   follows.  This assertion does not require `T(0,0)=4f(0,0)`.
2. If `deg_y(T)=deg_y(f)=n>=1`, the primitive leading colors have the exact
   square allocation in Section 3.  It is necessary, not sufficient.

## 1. The complete constant-y seam

Let `T=u(x)` and `f=v(x)!=0`.  As a polynomial in `K`, `(2)` has the form

```text
S=A_2K^2+A_1K+A_0,
A_2=v^2{L^2v(2+3av)-3u^2}.                              (3)
```

The top coefficient cannot vanish.  Indeed, `A_2=0` would give

```text
3u^2=L^2v(2+3av).                                       (4)
```

The factor `2+3av` is nonzero because its value at `a=0` is `2`.  Equation
`(4)` forces `L|u`; write `u=Lw`.  Then

```text
3w^2=v(2+3av).                                           (5)
```

If `u=0`, its right side is already nonzero.  Otherwise the left side has
even degree, whereas for `d=deg(v)>=0` the right side has degree `2d+1`.
This contradiction includes constant nonzero `v`.

Consequently `A_2!=0`, so `S` has y-degree four.  Write a hypothetical square
root as

```text
G=alpha*y^2+beta*y+gamma,              alpha!=0.         (6)
```

The residual is even in `y`; its y-cubic coefficient is zero.  Equation
`2alpha*beta=0` forces `beta=0`, so `G` is linear in `y^2`, equivalently in
`K`.  Hence the discriminant of `(3)` must vanish.  Direct expansion gives

```text
Disc_K(S)
 =-4(2+3av){L^2v(1+2av)-2u^2}^3.                       (7)
```

Again `2+3av` is nonzero, so `(7)` forces

```text
2u^2=L^2v(1+2av).                                       (8)
```

Again `L|u`; writing `u=Lw` gives `2w^2=v(1+2av)`, whose two sides have even
and odd degrees.  This final contradiction proves the x-only lane empty.  No
source address, square-root sign, or specialization was used.  The affine
translation `K=y^2-F` preserves the quadratic discriminant used in `(7)`.

## 2. Leading equation on positive equality

For context, let `n>=1` and write

```text
T=u*y^n+u_1*y^(n-1)+...,
f=v*y^n+v_1*y^(n-1)+...,
G=gamma*y^(2n+2)+gamma_1*y^(2n+1)+....                   (9)
```

The leading equation is

```text
gamma^2=3v^2(aL^2v^2-u^2).                              (10)
```

Its parenthesis is nonzero by the odd `a`-valuation.  Primewise UFD
valuations in `(10)` give `v|gamma`.  Choose `rho,i in k` with
`rho^2=3`, `i^2=-1`, and write

```text
gamma=rho*v*z,                 z^2+u^2=a(Lv)^2.          (11)
```

## 3. The gcd sidecar and primitive colors

Let

```text
c=gcd(z,u),        z=c*z_0,        u=c*u_0,
gcd(z_0,u_0)=1.                                             (12)
```

For every irreducible prime `pi`, `(11)` gives

```text
2 ord_pi(c) <= ord_pi(a)+2 ord_pi(Lv).                   (13)
```

For `pi!=a` this immediately gives `ord_pi(c)<=ord_pi(Lv)`; for `pi=a`,
integrality gives the same inequality after taking the floor.  Thus

```text
c | Lv.                                                   (14)
```

Put `C=Lv/c`.  The two primitive colors are coprime, since a common divisor
of `z_0+i*u_0` and `z_0-i*u_0` divides both `2z_0` and `2i*u_0`.  Their
product is

```text
(z_0+i*u_0)(z_0-i*u_0)=aC^2.                            (15)
```

Unique factorization therefore gives, up to interchanging the colors, a unit
`mu in k^*`, an `epsilon in {0,1}`, and coprime `R_-,R_+ in k[x]` such that

```text
z_0+i*u_0 = mu*a^epsilon*R_-^2,
z_0-i*u_0 = mu^(-1)*a^(1-epsilon)*R_+^2,
R_-*R_+   is associated to C.                            (16)
```

After absorbing a unit into a square root, the last association may be made
an equality.  Formula `(16)` is the complete primitive allocation: recording
only the product loses the common gcd sidecar `c`.

THM-3902 supplies the exact first two response jets below `(11)` and explicit
positive lifts.  Thus `(16)` cannot be promoted to equality-lane emptiness.

## 4. Exact degree partition and boundary

Combining the present closure with current canon gives

```text
f=0                              T=0 only                 THM-3897
f!=0, T,f in k[x]                empty                    THM-3904
f!=0, deg_y(f)>=1, deg_y(T)<deg_y(f) empty                THM-3899
f!=0, deg_y(T)=deg_y(f)>=1       two-jet passport         THM-3902
f!=0, deg_y(T)>deg_y(f)          response fan             THM-3901
```

This is a complete degree-regime filtration, not a square classification.
Positive common degree, strict `T`-leading degree after its response tariffs,
a polynomial Keller atlas, and `JC(2)` remain **OPEN**.

## 5. Exact replay

Run from the repository root:

```bash
python3 -B 04-computation/jc2_nonzero_sidecar_constant_y_seam_and_primitive_colors_thm3904.py
python3 -B -O 04-computation/jc2_nonzero_sidecar_constant_y_seam_and_primitive_colors_thm3904.py
python3 -B 04-computation/jc2_nonzero_sidecar_constant_y_seam_independent_audit_thm3904.py
python3 -B -O 04-computation/jc2_nonzero_sidecar_constant_y_seam_independent_audit_thm3904.py
```

The primary and independent streams must byte-match their frozen raw-LF
outputs.  The finite-field atlas is a hostile control, not a proof substitute;
the proof is the address-free UFD/discriminant argument above.  **QED.**
