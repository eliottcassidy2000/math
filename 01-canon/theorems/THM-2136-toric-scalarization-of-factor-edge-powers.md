---
id: THM-2136
title: "Toric scalarization of factor-edge powers and the Hermite compatibility budget"
status: >
  PROVED. Every positive-weight homogeneous irreducible in C[x,y] is an axis
  or a primitive toric binomial. On the normalization of any such divisor,
  a homogeneous factor quotient is a scalar monomial. Consequently every
  THM-2134 factor-edge polynomial has the exact form
  P(T)=t^A R(t^B T), with A,B integral and R in C[T]. Its residue-field
  m/gcd(m,n)-th-power condition is equivalent to R being that same power;
  the scalar exponent condition is automatic. Across distinct non-axis toric
  factors, the coefficients are leading Hermite jets. Their exact
  interpolation map has rank min(K+1,sum(v_i+1)), so low degree creates
  explicit compatibility relations while sufficiently high degree permits
  arbitrary local data.
  Many factors alone therefore give no contradiction and no planar-JC proof.
source: codex-2026-07-22-JC2-toric-factor-edge-scalarization
depends_on:
  - THM-2134
related:
  - THM-2102
  - THM-2132
---

# THM-2136 -- toric scalarization of factor-edge powers

Let `w=(w_1,w_2)` be a positive integral weight. Normalize it by

```text
g_w=gcd(w_1,w_2),        a=w_1/g_w,        b=w_2/g_w,
gcd(a,b)=1.                                           (1)
```

This theorem answers the first compatibility question left by THM-2134.
The answer has two parts. For one factor, all function-field coefficients
collapse to constant coefficients after a monomial change of variable. For
several factors, the missing compatibility is exactly a finite Hermite-jet
problem; it is not supplied merely by having many factors.

## 1. Classification of weighted-homogeneous irreducibles

Every nonconstant `w`-homogeneous irreducible in `C[x,y]`, up to a nonzero
scalar, is exactly one of

```text
x,
y,
x^b-lambda y^a,                  lambda in C^*.       (2)
```

### Proof

Let `pi` be such an irreducible. If `x|pi` or `y|pi`, irreducibility gives
the first or second case. Suppose neither variable divides `pi`. There is
then a monomial of `pi` with `x`-exponent zero and one with `y`-exponent
zero. If the common weighted degree is `q`, these two monomials show

```text
q=g_w a b ell                                      (3)
```

for some positive integer `ell`. Every solution of

```text
a i+b j=a b ell
```

has `b|i` and `a|j`. Hence, with

```text
X=x^b,                    Y=y^a,                     (4)
```

one has

```text
pi=H(X,Y)                                             (5)
```

for an ordinary homogeneous binary form `H` of degree `ell`. Over `C`, the
form `H` is a product of `ell` linear forms. The factors `X` and `Y` are
excluded by the assumption on the axes. Irreducibility therefore forces
`ell=1`, which gives the third form in (2).

For completeness, each binomial in (2) is irreducible. Choose `mu in C^*`
with `mu^b=lambda` and consider

```text
C[x,y] -> C[t],          x |-> mu t^a,       y |-> t^b. (6)
```

Its kernel contains `x^b-lambda y^a`. Divide any kernel element by this
binomial as a polynomial in `x`. The remainder has `x`-degree less than `b`.
After (6), two remainder monomials `x^i y^j` and `x^(i')y^(j')` with
`0<=i,i'<b` have the same `t`-exponent only if

```text
a(i-i')=b(j'-j).
```

Coprimality makes `b|(i-i')`, hence `i=i'` and then `j=j'`. Thus no
distinct remainder monomials cancel, the remainder is zero, and the kernel
is exactly the principal binomial ideal. Its quotient embeds in the domain
`C[t]`, so the binomial is prime and therefore irreducible. QED.

## 2. A homogeneous quotient restricts to one scalar monomial

For an irreducible `pi` in (2), use the following normalization data:

| `pi` | normalization of `pi=0` | `d_pi=deg_w(pi)` | `omega_pi` |
|---|---|---:|---:|
| `x` | `x=0, y=t` | `w_1` | `w_2` |
| `y` | `x=t, y=0` | `w_2` | `w_1` |
| `x^b-lambda y^a` | `x=mu t^a, y=t^b`, `mu^b=lambda` | `g_w a b` | `g_w` |

Let `A!=0` be `w`-homogeneous of degree `q`, and put

```text
v=v_pi(A).                                            (7)
```

Then there is a unique nonzero scalar `gamma` such that, on the normalization
in the table,

```text
(A/pi^v)|_(pi=0)
   =gamma t^E,
E=(q-v d_pi)/omega_pi in Z_(>=0).                    (8)
```

For an axis factor, divide out the exact axis power. Modulo that axis, a
homogeneous polynomial has one possible surviving pure monomial, giving
(8). For a binomial factor, the quotient is not divisible by `pi`, so its
image under (6) is nonzero. Every one of its monomials has `t`-exponent

```text
a i+b j=(q-v d_pi)/g_w,                              (9)
```

and their nonzero sum is therefore one scalar times that monomial. This also
proves the integrality and nonnegativity asserted in (8).

## 3. Scalarization of a THM-2134 factor edge

Use all notation of THM-2134. In particular, `h` has weighted degree `d`,

```text
e=v_pi(h),
E_pi={delta>=0 : v_pi(f_delta)+tau delta=me},
r_0=gcd(E_pi minus {0}),                              (10)
```

and the factor-edge polynomial is

```text
P_pi(T)=sum_(delta in E_pi) p_delta T^(delta/r_0)
          in kappa_pi[T].                            (11)
```

The normalization in Section 2 identifies `kappa_pi` with `C(t)`. Put

```text
K_pi={k>=0 : k r_0 in E_pi}.                         (12)
```

For `k in K_pi`, the edge identities give

```text
deg_w(f_(k r_0))=md-k r_0,
v_pi(f_(k r_0))=me-tau k r_0.                        (13)
```

Applying (8) yields

```text
p_(k r_0)=gamma_k t^(A_pi+k B_pi),       gamma_k in C^*, (14)
```

where

```text
A_pi=m(d-e d_pi)/omega_pi,
B_pi=r_0(tau d_pi-1)/omega_pi.                       (15)
```

Both numbers in (15) are integers. The first is the exponent in (8) for
`h^m/pi^(me)`. For the second, (14) first shows

```text
A_pi+k B_pi in Z                    for every k in K_pi. (16)
```

The positive elements of `K_pi` have gcd one by the definition of `r_0`.
Choose integers `c_k` with `sum c_k k=1`. Since `A_pi` is integral,

```text
B_pi=sum c_k(k B_pi) in Z.                           (17)
```

The integer `B_pi` may be negative; every exponent actually occurring in
(14) is nevertheless nonnegative by (8).

Define the constant-coefficient edge polynomial

```text
R_pi(S)=sum_(k in K_pi) gamma_k S^k in C[S].         (18)
```

Equations (14)--(18) give the exact scalarization

```text
P_pi(T)=t^A_pi R_pi(t^B_pi T)             in C(t)[T]. (19)
```

There is no coefficient function left after this monomial gauge.

## 4. The residue-field power is exactly a scalar polynomial power

Let

```text
g_0=gcd(m,n),                       M=m/g_0.          (20)
```

Then

```text
P_pi(T) is an Mth power in C(t)[T]
      iff
R_pi(S) is an Mth power in C[S].                     (21)
```

More generally, if `A,B in Z`, `R in C[S]` has `R(0)!=0`, and

```text
P(T)=t^A R(t^B T),                                   (22)
```

then `P` is an `M`th power in `C(t)[T]` if and only if

```text
M|A                 and              R=R_0^M
                                      for some R_0 in C[S]. (23)
```

Indeed, the change `T=t^(-B)S` is an automorphism over `C(t)` and sends
`P` to `t^A R(S)`. Unique factorization in `C(t)[S]` says that every
nonconstant irreducible multiplicity of `R` must be divisible by `M`.
Because `C` is algebraically closed, its constant factor has an `M`th root,
so `R=R_0^M`. The remaining unit is `t^A`; its `t`-adic valuation is divisible
by `M` exactly when `M|A`. The converse follows from the explicit root

```text
t^(A/M) R_0(t^B T).                                  (24)
```

In the factor-edge setting, the scalar exponent condition is automatic.
Apply (8) to `u=h/pi^e`. There are `beta in C^*` and

```text
H_pi=(d-e d_pi)/omega_pi in Z_(>=0)                  (25)
```

such that `bar(u)=beta t^H_pi`. Therefore

```text
A_pi=m H_pi,                                         (26)
```

and `M|A_pi`. This proves (21).

If THM-2134 gives `P_pi=Q_pi^M` and

```text
R_pi=S_pi^M,                                         (27)
```

then the root can be chosen in the scalarized form

```text
Q_pi(T)=t^(A_pi/M) S_pi(t^B_pi T),                   (28)
```

with

```text
S_pi(0)=beta^g_0.                                    (29)
```

Thus the local root is unique after its constant coefficient is normalized,
but (28) still belongs to the normalization of one divisor.

## 5. Distinct toric factors are a Hermite-jet problem

The remaining cross-factor issue can be stated exactly. Keep the normalized
weights (1), set `X=x^b`, `Y=y^a`, and let `A` be any nonzero `w`-homogeneous
polynomial of weighted degree `q`. There are unique integers

```text
0<=r<b,              0<=s<a,              K>=0       (30)
```

and a binary form `H_A(X,Y)` of ordinary degree `K` such that

```text
A=x^r y^s H_A(X,Y),
q/g_w=a r+b s+a b K.                                 (31)
```

Writing

```text
H_A(X,Y)=Y^K Phi_A(X/Y),             deg Phi_A<=K,   (32)
```

identifies the coefficients of `A` with those of `Phi_A`.

Fix distinct nonzero numbers `lambda_1,...,lambda_N`, put

```text
pi_j=X-lambda_j Y,
v_j=v_(pi_j)(A)=ord_(lambda_j)(Phi_A),                (33)
```

and choose `mu_j^b=lambda_j`. If

```text
c_j=Phi_A^(v_j)(lambda_j)/v_j!,                       (34)
```

then the scalar in (8) is exactly

```text
gamma_j=mu_j^r c_j.                                  (35)
```

Indeed, divide (31)--(32) by

```text
pi_j^v_j=Y^v_j(X/Y-lambda_j)^v_j
```

and substitute `x=mu_j t^a`, `y=t^b`. Formula (35) is the surviving leading
Taylor coefficient; all remaining powers of `t` are exactly the exponent in
(8).

For fixed nonnegative integers `v_1,...,v_N`, define the confluent evaluation
map

```text
J_K : C[Z]_(<=K)
   -> direct_sum_(j=1)^N C[epsilon]/(epsilon^(v_j+1)),
Phi |-> (Phi(lambda_j+epsilon) mod epsilon^(v_j+1))_j. (36)
```

Put

```text
B=sum_(j=1)^N (v_j+1).                               (37)
```

Then

```text
rank(J_K)=min(K+1,B),
codim(im J_K)=max(0,B-K-1).                          (38)
```

If `K<B`, a polynomial in the kernel is divisible by

```text
product_j (Z-lambda_j)^(v_j+1),                      (39)
```

whose degree is `B`; hence the kernel is zero. If `K>=B-1`, the Chinese
remainder theorem modulo the pairwise coprime powers in (39) makes `J_K`
surjective. These two observations prove (38), including the boundary
`K=B-1`.

In particular, arbitrary prescribed nonzero leading jets

```text
Phi(Z)=c_j(Z-lambda_j)^v_j
          modulo (Z-lambda_j)^(v_j+1)                (40)
```

can be realized simultaneously whenever

```text
K>=sum_j(v_j+1)-1.                                   (41)
```

For the leading scalars alone, put

```text
S=sum_j v_j,             D(Z)=product_j(Z-lambda_j)^v_j. (41a)
```

Because the `v_j` are the actual vanishing orders of `Phi_A`, one has
`S<=K` and `Phi_A=D Psi` with `deg Psi<=K-S`. If

```text
D_j(Z)=D(Z)/(Z-lambda_j)^v_j,
```

then

```text
c_j=D_j(lambda_j) Psi(lambda_j),                     (41b)
```

and every factor `D_j(lambda_j)` is nonzero. Ordinary evaluation at the
distinct `lambda_j` therefore gives the leading-scalar map the exact rank

```text
min(K-S+1,N)                                         (41c)
```

inside `C^N`, and the exact codimension

```text
max(0,N-K+S-1)=max(0,B-K-1).                         (41d)
```

Thus below (41) the leading scalars satisfy exactly `B-K-1` independent
linear compatibility relations. The harmless factors `mu_j^r` in (35) are
nonzero diagonal rescalings, so (41c)--(41d), not the full-jet rank in (38),
apply to the normalized edge scalars `gamma_j`.

## 6. What high factor count does and does not prove

For several non-axis toric factors sharing fixed face, edge index, and
monomial gauge data, equations (36)--(41) give a genuine interpolation
invoice. If one additionally claims that a chosen coefficient from the local
roots (28) is the restriction of one common global `w`-homogeneous root face,
then those coefficients must land in the Hermite image for that face's
weighted-degree budget `K`. THM-2134 by itself makes no such cross-factor
identification. Under this additional common-face claim, if

```text
K<sum_j(v_j+1)-1,                                    (42)
```

there are explicit linear relations to check. This is the exact
high-factor-count sidecar missing from a bare residue-field power statement.

It is not an automatic obstruction:

1. when (41) holds for the one common coefficient face under consideration,
   Hermite interpolation realizes arbitrary local leading jets on that face;
   realizing whole independently chosen local `M`th-power root polynomials
   requires the threshold and a common edge-index/gauge identification
   separately for every coefficient face;
2. when (42) holds, rank deficiency alone does not show that the nonlinear
   `M`th-power locus misses the Hermite image; and
3. a genuine global coarsened power satisfies every local power and every
   interpolation relation, regardless of the number of factors.

Moreover different factors can have different slopes `tau`, lattices `r_0`,
edge supports, and monomial gauges `B_pi`. Scalarization supplies no canonical
identification among those data. A contradiction therefore requires a fixed
cohort of factors, its actual weighted-degree budget, and a proved nonzero
Hermite or subresultant obstruction. THM-2136 exposes that finite target but
does not manufacture it and does not prove planar JC.
