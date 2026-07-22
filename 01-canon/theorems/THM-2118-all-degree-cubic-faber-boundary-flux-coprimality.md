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
script_sha256: a4b4ae5243c8c693c4264929b0b6e77be91ac80a5660bf99ed61a12260c71e85
output_sha256: 0018a5b61aa5ab99ccd8cada63b46bcccf76303f39b2127516795fa7072598a2
hash_basis: repository blobs with LF line endings
---

# THM-2118 -- all-degree cubic boundary-flux coprimality and closure

## 1. The two centered polynomials

Retain THM-2084's depressed cubic

```text
P_a(z)=z^3+a z-(a+1),                                  (1)
```

so `P_a(1)=0`. For `n>=1`, `3 not|n`, let

```text
E_n(z;a)=Pol_z P_a(z)^(n/3),
e_n(a)=E_n(1;a),
phi_n(a)=3 [z^(-1)] P_a(z)^(n/3).                     (2)
```

These are exactly the boundary and flux polynomials in THM-2084. The theorem
is

```text
gcd_(C[a])(e_n,phi_n)=1              for every 3 not|n. (3)
```

Put

```text
s=a+3,                   Q_s(x)=1+3x+s x^2.           (4)
```

The translation from `a` to `s` preserves polynomial gcds.

## 2. Lagrange inversion exposes adjacent coefficients

Let `f(w)=w+O(w^(-1))` be the formal inverse at infinity defined by

```text
P_a(f(w))=w^3.                                          (5)
```

The standard Faber generating identity is

```text
f'(w)/(f(w)-z)=sum_(m>=0) E_m(z;a) w^(-m-1).           (6)
```

Integrating in `w` and taking `z=1`, `t=w^(-1)`, gives

```text
-log(t(f(1/t)-1))=sum_(m>=1) e_m(a)t^m/m.              (7)
```

Let `r_1,r_2` be the other roots of `P_a`. The numbers

```text
A=1-r_1,             B=1-r_2
```

satisfy `A+B=3`, `AB=s`. If

```text
g=t(f(1/t)-1),                  u=t/g,                  (8)
```

then the factored cubic equation becomes

```text
g(g+At)(g+Bt)=1,
g=Q_s(u)^(-1/3),                u=t Q_s(u)^(1/3).       (9)
```

Lagrange inversion applied to (7)--(9) yields the exact finite coefficient
formula

```text
e_n(s)=[x^n] Q_s(x)^(n/3).                             (10)
```

The flux is the adjacent coefficient:

```text
phi_n(s)=3[x^(n+1)] Q_s(x)^(n/3).                      (11)
```

One quick normalization check for (11) is to differentiate the polynomial
part directly:

```text
d/da P_a(z)^(m/3)=(m/3)(z-1)P_a(z)^((m-3)/3).
```

Evaluating its polynomial part at `z=1` leaves precisely the first omitted
coefficient, so

```text
e_m'(a)=(m/9)phi_(m-3)(a),              m>=4.          (12)
```

Apply (12) to `m=n+3` and differentiate (10) to obtain (11), including its
factor three.

## 3. The boundary-flux identity

Put `F=Q_s^(n/3)`. The coefficient recurrence obtained from
`Q_s F_x=(n/3)(3+2sx)F` and differentiation in `s` gives

```text
(n+1)phi_n
 =s(9-4s)e_n'+2n(s-3)e_n.                             (13)
```

Equivalently, the coefficient of `x^n` in

```text
s(9-4s)F_s+2n(s-3)F-3F_x
```

is zero. This is a formal coefficient identity, not an asymptotic relation.

Consequently, if `e_n(s_0)=phi_n(s_0)=0` and
`s_0(9-4s_0)!=0`, then

```text
e_n'(s_0)=0.                                           (14)
```

Thus every common zero away from `0,9/4` would be a double zero of `e_n`.

## 4. The hypergeometric ODE

Expanding (10) gives

```text
e_n(s)=d_0 * 2F1(-n/2,(1-n)/2;1-2n/3;4s/9),           (15)
d_0=3^n binom(n/3,n)!=0.
```

The hypergeometric series terminates at degree `floor(n/2)`. Gauss's equation
becomes

```text
s(9-4s)e_n''
 +[9-6n+(4n-6)s]e_n'
 -n(n-1)e_n=0.                                        (16)
```

At every ordinary point `s notin {0,9/4}`, a solution of (16) is determined
by its value and first derivative. Equations (14) and (16) would therefore
force the nonzero polynomial `e_n` to vanish identically, a contradiction.

It remains to exclude the two singular points. Formula (10) gives

```text
e_n(0)=3^n binom(n/3,n),
e_n(9/4)=(3/2)^n binom(2n/3,n).                        (17)
```

For an integer `n>=1`, `binom(x,n)` vanishes exactly when
`x in {0,1,...,n-1}`. Since `3 not|n`, both `n/3` and `2n/3` are nonintegers.
Hence both values in (17) are nonzero. There is no common zero anywhere in
`C`, proving (3). QED.

## 5. Frontier effect

THM-2084 isolated two all-degree laws behind the finite calculations:

1. the balanced centered-pole noncollision `gcd(e_n,phi_n)=1`;
2. the upper-Newton primitive noncollision for the pair `(Phi_n,R_n)` at
   polynomial infinity.

The coprimality argument in Sections 1--4 closes the first law for every
admissible degree. By itself it does not show that every finite pole is
balanced, control the depressed-coefficient poles, or separate lower Faber
representatives. Sections 8--10 supply those additional valuation arguments
and then invoke THM-2102, closing the cubic source-fiber stratum without the
upper-Newton law. That law remains an independent question about arbitrary
Faber combinations. General proper-power descent and planar JC(2) remain open.

## 6. Exact referee and assumption challenge

The companion performs two independent exact checks.

- A SymPy/Rational path reconstructs `E_n` and `phi_n` from the original
  THM-2084 finite binomial sums, independently reconstructs (10)--(11), and
  checks (13), (16), (17), and the gcd through `n=80`.
- A standard-library `Fraction` path generates the two adjacent coefficient
  polynomials independently and runs rational Euclidean gcds through `n=200`.

The challenged assumption was that THM-2084's low-degree resultants were
inherently finite. The correct object is instead the adjacent-coefficient
pair of one quadratic fractional power; the ODE supplies the all-degree
rigidity. Tournament Analysis is not natural for this referee: there is no
intrinsic pairwise binary relation or tie Hamiltonian path, and forcing one
would erase the derivative and singular-point data used by the proof.

## 7. Independent translated-coefficient proof and cubic closure

The first proof establishes the all-degree coprimality theorem. The following
independent coefficient derivation then supplies the uniform valuation and
power-free-face steps needed for cubic source-fiber closure. Equation labels
restart within this proof.

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

### 7.1. Translation gives one quadratic coefficient sequence

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

### 7.2. The boundary-flux identity and hypergeometric ODE

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

### 7.3. Coprimality, including both singular points

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

## 8. Uniform removal of the centering poles

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

### 8.1. The regime `rho>2H`

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

### 8.2. The regime `rho<2H`

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

### 8.3. The balanced regime `rho=2H`

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

## 9. Uniform removal of the remaining depressed-coefficient poles

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

## 10. The power-free face closes the cubic stratum

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

## 11. Exact referee and validity boundary

The companion is exactly the referee described in Section 6: over rational
polynomial rings it checks the translated formulas, differential identities,
singular values, and gcd through `n=80`, with an independent `Fraction` gcd
census through `n=200`. It does **not** computationally certify the valuation
splits or the power-free-face application in Sections 8--10; those are
all-degree paper arguments. Its validations use Python `assert`, so only the
ordinary run is evidentiary; an optimized transcript is not an independent
check.

There is no intrinsic binary relation in this argument: the faithful carrier
is the ordered coefficient/valuation filtration `(n;H;rho;e_n;phi_n;R_n)`.
Forcing a tournament would discard the parity and pole-order sidecars that
make Sections 8--10 work.
