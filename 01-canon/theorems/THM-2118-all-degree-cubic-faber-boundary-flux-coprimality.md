---
id: THM-2118
title: "All-degree cubic Faber boundary-flux coprimality and cubic source-fiber closure"
status: >
  PROVED. For every n>=1 not divisible by three, THM-2084's centered
  boundary polynomial e_n and flux polynomial phi_n have no common complex
  root. A translated coefficient formula identifies e_n with a terminating
  Gauss hypergeometric polynomial; an exact first-order boundary-flux identity
  and its regular-singular ODE prove coprimality in every degree. Uniform
  valuation arguments then remove every finite pole of the cubic depressed
  coordinates. Combined with THM-2102's power-free-face theorem, this proves
  that every planar Keller pair having an output-pencil member cubic along a
  linear source fiber is a polynomial automorphism. This does not prove JC(2);
  the next source-fiber frontier is quartic.
source: codex-2026-07-22-JC2-cubic-Faber-all-degree
depends_on:
  - THM-2084
  - THM-2102
related:
  - THM-2110
  - THM-2113
script: 04-computation/jc2_cubic_faber_boundary_flux_coprimality_codex_20260722.py
output: 05-knowledge/results/jc2_cubic_faber_boundary_flux_coprimality_codex_20260722.out
script_sha256: TO_BE_FILLED
output_sha256: TO_BE_FILLED
hash_basis: repository blobs with LF line endings
---

# THM-2118 -- all-degree cubic boundary-flux coprimality and closure

## 1. Statement and Faber notation

For an integer `m>=1` not divisible by three, retain THM-2084's cubic
Faber polynomials

```text
P_0=z^3+pz+q,
E_m(z;p,q)=Pol_z P_0^(m/3),
Phi_m(p,q)=3[z^-1]P_0^(m/3),
R_m(p,q)=3[z^-2]P_0^(m/3).                         (1)
```

Here the fractional power is the unique formal Laurent branch with leading
term `z^m`. On the centered boundary `P_0(1)=0`, put

```text
e_n(a)=E_n(1;a,-1-a),
phi_n(a)=Phi_n(a,-1-a).                              (2)
```

We prove two assertions.

> **Boundary-flux theorem.** If `n>=1` and `3` does not divide `n`, then
> ```text
> gcd_C[a](e_n,phi_n)=1.                             (3)
> ```

> **Cubic source-fiber theorem.** Let `(P,Q)` be a planar complex polynomial
> Keller pair. If some nonconstant member of its output pencil has degree
> three along a linear source-fiber direction, then `(P,Q)` is a polynomial
> automorphism.

The second assertion is a source-fiber theorem, not a bound on generic cover
degree. It closes the cubic stratum left open by THM-2084 and THM-2110, but it
does not prove the planar Jacobian conjecture.

## 2. Translation gives one quadratic coefficient sequence

Fix `n` and write

```text
alpha=n/3,                    b=a+3,
T(u)=1+3u+b u^2.                                    (4)
```

On the boundary in (2), translation by `t=z-1` gives

```text
z^3+az-(1+a)=t(t^2+3t+b)
             =t^3 T(t^-1).                          (5)
```

Translation preserves the polynomial part at infinity. It also preserves the
coefficient of `z^-1`: terms `z^-k` with `k>=2` remain of order at most
`t^-2`, while `z^-1=t^-1+O(t^-2)`. Hence (5) gives the exact pair

```text
e_n(a)   =[u^n]T(u)^alpha,
phi_n(a)=3[u^(n+1)]T(u)^alpha.                       (6)
```

Expanding the first coefficient yields

```text
e_n(a)=sum_(0<=j<=floor(n/2))
  binom(alpha,n-j) binom(n-j,j) 3^(n-2j) b^j.        (7)
```

Its constant coefficient in `b` is

```text
3^n binom(n/3,n) != 0,                               (8)
```

because `n/3` is not an integer. In particular `e_n` is not the zero
polynomial. The ratio of consecutive coefficients in (7) gives the terminating
hypergeometric form

```text
e_n(a)=3^n binom(n/3,n)
  * 2F1(-n/2,(1-n)/2;1-2n/3;4(a+3)/9).              (9)
```

One upper parameter is a nonpositive integer, and the lower parameter is
never an integer when `3` does not divide `n`, so (9) is an ordinary
polynomial identity with no analytic continuation issue.

## 3. The boundary-flux identity and hypergeometric ODE

Let `C_j` denote the coefficient of `b^j` in (7), extended by zero outside
its displayed range, and let `D_j` be the coefficient of `b^j` in
`[u^(n+1)]T(u)^alpha`. Directly from the two binomial formulas,

```text
3(n+1)D_j
 =2(n-2j+2)C_(j-1)+3(3j-2n)C_j.                    (10)
```

The right side is exactly the coefficient of `b^j` in

```text
2n(b-3)e_n-b(4b-9) d e_n/db.
```

Since `b=a+3`, equation (10) proves the first-order identity

```text
(n+1)phi_n
 =2na e_n-(4a^2+15a+9)e_n'.                        (11)
```

The Gauss equation for (9), after the affine change
`x=4(a+3)/9`, is

```text
(a+3)(4a+3)e_n''
 -(2n-3)(2a+3)e_n'+n(n-1)e_n=0.                    (12)
```

Both equations are polynomial identities over `Q[a]`.

## 4. Coprimality, including both singular points

Suppose `e_n(a_0)=phi_n(a_0)=0`, and set

```text
Delta(a)=(a+3)(4a+3).                                (13)
```

If `Delta(a_0)!=0`, equation (11) gives `e_n'(a_0)=0`. A root of
multiplicity `r>=2` makes the `Delta e_n''` term in (12) have exact order
`r-2`, while the other two terms have larger order. This is impossible.

It remains to exclude the regular-singular points. If `e_n` vanished to the
positive integer order `r` at `a=-3`, the coefficient of order `r-1` in
(12) would give

```text
-9(r-1)+3(2n-3)=0,                  r=2n/3.           (14)
```

This is not an integer because `3` does not divide `n`. At `a=-3/4`, the
same lowest-order comparison gives

```text
9(r-1)-(3/2)(2n-3)=0,               r=n/3+1/2.       (15)
```

The right side is never an integer: `(2n+3)/6` has odd numerator. Thus
`e_n` has no root at either zero of `Delta`, and a common root is impossible.
This proves (3).

The excluded degrees are a necessary boundary control. If `3|n`, then
`P_0^(n/3)` is already a polynomial, so both `e_n` and `phi_n` vanish
identically on `P_0(1)=0`.

## 5. Uniform removal of the centering poles

We now apply (3) to THM-2084's all-degree cubic normal form. After a linear
target change, write the cubic pencil member and a reduced mate as

```text
P=A(x)y^3+B(x)y^2+C(x)y+D(x),          A!=0,
J(P,Q)=kappa in C*,
n=min_H deg_y(Q-H(P)).                                (16)
```

THM-2084 proves `3` does not divide `n`, writes `A=U^3`, and over `C(x)`
introduces

```text
h=B/(3U^2),       z=Uy+h,
p=C/U-3h^2,       q=D-h^3-ph,
P=z^3+pz+q.                                            (17)
```

Modulo a polynomial target shear, its Faber normal form is

```text
Q=sum_(m<=n, 3 not dividing m) c_m E_m,       c_n!=0, (18)
```

with constant coefficients. For the corresponding linear combinations
`Phi=sum c_m Phi_m` and `R=sum c_m R_m`, the Keller equation is

```text
Phi'=0,                         R'=kappa/U.            (19)
```

Let `v` be a finite-place valuation and suppose `v(h)=-H<0`. Since the
original constant coefficient

```text
D=h^3+ph+q                                             (20)
```

is polynomial, put `rho=-v(p)` and split into three exhaustive regimes.

### 5.1. The regime `rho>2H`

Here `q~-ph`, so a monomial `p^i q^j` in `Phi_m`, where
`2i+3j=m+1`, has pole order

```text
(i+j)rho+jH.                                           (21)
```

Increasing `j` by two lowers (21) by `rho-2H>0`. Thus the unique largest
term has `j=0` for odd `m` and `j=1` for even `m`; its binomial coefficient
is nonzero because `m/3` is not an integer. Its order is

```text
((m+1)/2)rho                    if m is odd,
(m/2)rho+H                      if m is even.           (22)
```

At `m=n` this strictly exceeds every lower Faber order. Hence `Phi` has a
pole, contradicting its constancy in (19).

### 5.2. The regime `rho<2H`

Now `q~-h^3`. Every term of `E_m(h;p,q)` has pole order at most `mH`, and
the coefficient at that order is

```text
E_m(1;0,-1)
 =sum_(j=0)^floor(m/3) (-1)^j binom(m/3,j)
 =(-1)^r binom(m/3-1,r) !=0,           r=floor(m/3).   (23)
```

The top summand `c_n E_n(h)` therefore has the unique pole order `nH` in
`Q(x,0)`; all lower representatives have order below `nH`. This contradicts
the polynomiality of the original reduced mate at `y=0`.

### 5.3. The balanced regime `rho=2H`

For suitable nonzero leading coefficients,

```text
p~a h^2,                         q~-(1+a)h^3.           (24)
```

The unique top orders in polynomiality of `Q(x,0)` and constancy of `Phi`
respectively force

```text
e_n(a)=0,                         phi_n(a)=0.            (25)
```

This contradicts (3). All three regimes are empty, so

```text
h in C[x].                                               (26)
```

## 6. Uniform removal of the remaining depressed-coefficient poles

Equation (17) and polynomiality of the original coefficients now give

```text
Up in C[x],                         q+hp in C[x].        (27)
```

If `p` had a finite pole of order `rho>0`, then `q=O(p)` there.

If `n` is odd, `Phi_n` has the unique term of maximal ordinary `(p,q)`
degree

```text
constant * p^((n+1)/2).                                (28)
```

It has strictly larger pole order than every other top term and every lower
`Phi_m`, contradicting `Phi'=0`. If `n` is even, `E_n(h)` instead has the
unique maximal ordinary-degree term

```text
binom(n/3,n/2) p^(n/2),                                (29)
```

whose coefficient is nonzero. Since `h` is regular and `q=O(p)`, it
strictly dominates every other term in `Q(x,0)`, again a contradiction.
Thus

```text
p,q in C[x].                                             (30)
```

Now `R` is polynomial, so (19) makes `kappa/U=R'` polynomial. A polynomial
whose reciprocal times a nonzero constant is polynomial is a unit. Hence

```text
U in C*,                                                (31)
```

and `(x,y)->(x,z=Uy+h)` is an honest triangular polynomial coordinate.

## 7. The power-free face closes the cubic stratum

It remains only to use the new connection to THM-2102. In polynomial
coordinates we have

```text
P=z^3+p(x)z+q(x).                                       (32)
```

If `p,q` are both constant, then

```text
J(P,Q)=-(3z^2+p) Q_x                                   (33)
```

cannot be a nonzero constant. Thus at least one of `p,q` is nonconstant.
Set the degree of a constant or zero entry to zero, put

```text
a=deg p,             b=deg q,
w_x=6,               w_z=max(3a,2b)>0.                (34)
```

The leading `w`-face of (32) is `z^3` together with every one of

```text
p_a x^a z                 when w_z=3a,
q_b x^b                   when w_z=2b,                 (35)
```

that occurs. At least one displayed extra term occurs.

This face is power-free in the exact sense of THM-2102. Indeed, if it were
`lambda H^m` with `m>=2`, its `z`-degree three would force `m=3` and
`deg_z H=1`. Write `H=c(x)z+r(x)`. Comparison of the leading `z^3`
coefficient makes `c(x)` a nonzero constant, while the absent `z^2`
coefficient gives `r(x)=0`. Then `lambda H^3` is a scalar multiple of
`z^3`, contradicting (35).

THM-2102 now makes `P` triangular and the Keller pair a polynomial
automorphism. This proves the cubic source-fiber theorem.

The all-degree upper-Newton primitive problem from THM-2084 remains a valid
question about arbitrary Faber combinations, but it is no longer required to
close the cubic Keller stratum: the polynomial depressed cubic exposes a
power-free positive-weight face first. The next unresolved source-fiber
degree is four. `JC(2)` and `DC(2)` remain open. QED.

## 8. Exact referee and validity boundary

The companion script works over exact rational polynomial rings. For every
`1<=n<=80` prime to three it independently checks:

- the translated coefficient formulas (6)--(7) against the original Faber
  sums;
- the first-order identity (11) and ODE (12);
- constant gcd of `e_n` and `phi_n`;
- the nonzero coefficients used in both unbalanced valuation regimes and in
  the odd/even depressed-pole split.

It also checks the degree-thirteen boundary resultant from THM-2110, the
necessary common-zero control at degrees divisible by three, and representative
power-free depressed-cubic faces. Checks use explicit exceptions rather than
Python `assert`, so normal and optimized runs execute the same referee. The
proof above, not the finite sweep, supplies every all-degree quantifier.

There is no intrinsic binary relation in this argument: the faithful carrier
is the ordered coefficient/valuation filtration `(n;H;rho;e_n;phi_n;R_n)`.
Forcing a tournament would discard the parity and pole-order sidecars that
make Sections 5--7 work.
