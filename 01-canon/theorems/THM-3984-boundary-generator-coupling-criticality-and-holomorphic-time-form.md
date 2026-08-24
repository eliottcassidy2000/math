---
id: THM-3984
title: "Boundary-generator couplings force criticality and a holomorphic time form"
status: >
  PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT. In every
  height n>=2, the first-coordinate family x+h(x^n t)+alpha*p+beta*y is a
  submersion exactly when alpha=beta=0; every nonzero linear coupling to
  either boundary generator creates an affine critical point. More
  generally x+h(x^n t)+c*p^a*y^b is critical for every a,b>=0 with
  a+b>0 and c!=0. At height two the generic-fibre time form for the
  leading-x linear-y family is a
  nonzero holomorphic differential for every degree of h. Its exact
  infinity ledger gives genus 2 for constant h and
  (4d+1-gcd(3,d))/2 for degree d>=1, so no finite algebraic cover creates
  a rational mate. These closures concern the displayed first-coordinate
  cells only; nonlinear multiweight p/y sums and JC(2) remain open.
source: jc-extra-debt-local + jc-mixed-generator-submersion / post-THM-3983 coupling lane, 2026-08-24
depends_on:
  - THM-3973-exact-volume-simple-cubic-determinantal-affine-plane-completion
related:
  - THM-3982-polynomial-shear-submersion-rational-exactness-and-two-color-image
  - THM-3983-coordinate-boundary-constancy-and-rational-place-budget
script: 04-computation/jc2_boundary_generator_coupling_criticality_thm3984.py
output: 05-knowledge/results/jc2_boundary_generator_coupling_criticality_thm3984.out
script_sha256: f7c8d8acd0e2c7821716eee48671cb3a6f6f22122c7ab968f3be7ffec88bf6b1
output_sha256: 62ebe106accb79f579b85d81dc53bbaa4830c3c368bb8d5beb5614884cb0dc7c
semantic_sha256: d9598c0dedda6f5637e168c97050a9ca0a5eebeec0245ede653d7283972849d9
hash_basis: raw LF bytes
---

# THM-3984 -- the first boundary-generator couplings are all critical

**PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT.** Work over
an algebraically closed field `k` of characteristic zero. For every `n>=2`
put

```text
u=x^n t,                 r=u(u+1),                 s=u^2(u+1),
p=r/x^n=(1+x^n t)t,      y=s/x^(n+1)=x^(n-1)(1+x^n t)t^2.   (1)
```

Thus `p,y` are the two boundary generators of the height-`n` completion
`B_n` from THM-3973. The following three statements hold.

1. For every polynomial `h in k[u]`,

   ```text
   A=x+h(u)+alpha*p+beta*y                              (2)
   ```

   is a source-plane submersion if and only if
   `alpha=beta=0`.

2. More generally, for every `a,b>=0` with `a+b>0`, every polynomial `h`,
   and every `c!=0`,

   ```text
   A_(a,b)=x+h(u)+c p^a y^b                             (3)
   ```

   has an affine critical point with `x!=0`. Consequently `(3)` is never
   the first entry of a polynomial Keller pair (or of a Darboux pair whose
   second entry is regular at that point). No general rational-mate claim is
   made in this panel.

3. At height `n=2`, for `c!=0` and arbitrary `h`, the rational time form on
   the geometric generic fibre of `A=x+h(u)+cy` is a nonzero holomorphic
   differential. Put `d=0` when `h` is constant (including zero), and
   otherwise put `d=deg h`. Its genus and number of source-plane places at
   infinity are

   ```text
   d=0:       genus=2,                    r_infinity=3;
   d>=1:      genus=(4d+1-gcd(3,d))/2,    r_infinity=3+gcd(3,d). (4)
   ```

   The time form remains nonexact on every finite algebraic cover. Thus the
   rational-mate obstruction is intrinsic to the generic fibre, independent
   of the special critical fibres supplied by Statements 1--2.

## 1. Critical equations for the mixed linear cell

On `x!=0`, `(x,u)` are source coordinates. Write `g=h'(u)` and
`N=n+2`. Differentiating `(2)` at fixed `u` and at fixed `x` shows that a
critical point is exactly a common zero of

```text
F=x^N-n alpha r x-(n+1)beta s,
G=g x^(N-1)+alpha(2u+1)x+beta*u(3u+2).                  (5)
```

### 1.1 Both boundary coefficients are nonzero

Assume `alpha*beta!=0` and let

```text
R(u)=Res_x(F,G).                                         (6)
```

The endpoint and degree ledger is

```text
ord_(u=0) R=2,
[u^2]R=(-1)^n(n-1)alpha^N beta,
R(-1)=beta^N,                                            (7)
```

and, if `g!=0` has degree `e` and leading coefficient `ell`,

```text
deg R=N e+3(N-1),
lc R=(-1)^(N+1)(N-1)^(N-1) beta^(N-1) ell^N.            (8)
```

If `g=0`, the replacement row is

```text
deg R=2N,                         lc R=(3beta)^N.         (9)
```

Here is a direct derivation of the two delicate seams. At `u=0`, the only
common root of the specialized rows can come from `x=0`. The unique nearby
root of `G` is

```text
x=-(2beta/alpha)u+O(u^2),
```

and substitution into `F` gives `(n-1)beta*u^2+O(u^3)`.
The remaining local resultant factor is a unit; the Sylvester sign and its
value give the coefficient in `(7)`. At `u=-1`, one has `F=x^N` and
`G(0)=beta`, proving the last row of `(7)`.

For `(8)`, the `N` roots of `F` at `u=infinity` satisfy

```text
x_i ~ xi_i u^(3/N),              xi_i^N=(n+1)beta.       (10)
```

If `g!=0`, the row `g x_i^(N-1)` strictly dominates the other two terms in
`G`; multiplying over the roots and using their product gives `(8)`. If
`g=0`, the term `beta s'~3beta u^2` dominates and gives `(9)`.

Every degree in `(8)--(9)` is greater than two. Hence `R/u^2` is
nonconstant. Its constant term is nonzero by the exact order in `(7)`, and
`R(-1)!=0`, so `R` has a root `u_0` outside `{0,-1}`. At such a root
`F(0,u_0)=-(n+1)beta s(u_0)!=0`; every common `x_0` is nonzero and gives a
genuine affine critical point with `t_0=u_0/x_0^n`.

### 1.2 The two coefficient axes

The specializations `beta=0` and `alpha=0` lose powers of `x`, so they
must not be inferred by specializing `(6)`.

If `beta=0,alpha!=0`, divide the two rows by their common `x`. The remaining
equations are

```text
x^(n+1)=n alpha r,                 g x^n=-alpha(2u+1).    (11)
```

Their consecutive-power compatibility is

```text
E_p=(-n)^n r^n g^(n+1)+alpha(2u+1)^(n+1)=0.             (12)
```

For nonzero `g`, its first term has strictly larger degree than the second;
moreover `E_p(0)` and `E_p(-1)` are nonzero. Choose a root outside the two
colors. If `g` is nonzero there, the consecutive powers reconstruct `x`.
If `g` vanishes, `(12)` forces `u=-1/2`, and any solution of
`x^(n+1)=-n alpha/4` works. The same explicit choice handles `g=0`
identically.

If `alpha=0,beta!=0`, the critical equations are

```text
x^(n+2)=(n+1)beta s,             g x^(n+1)=-beta s'.     (13)
```

After stripping their forced `u^(n+2)` compatibility factor, the residual
is

```text
E_y=(n+1)^(n+1)u^n(u+1)^(n+1)g^(n+2)
      -(-1)^(n+2)beta(3u+2)^(n+2).                       (14)
```

It is nonconstant and nonzero at `u=0,-1`. At a root with `g!=0`, the two
consecutive powers reconstruct `x`. If `g=0` at the root, `(14)` forces
`u=-2/3`, where `s'=0`; choose any nonzero root of
`x^(n+2)=(n+1)beta s(-2/3)`. Again this also handles `g=0` identically.

Finally, when `alpha=beta=0`, one has `A=x+h(x^n t)` and the exact Euler
identity

```text
x A_x-n t A_t=x.                                         (15)
```

A critical point would first force `x=0`, but there `A_x=1`. Thus the
uncoupled polynomial shear is a submersion. Sections 1.1--1.2 prove the
if-and-only-if assertion in Statement 1.

## 2. Every single boundary monomial is critical

For `(3)`, keep arbitrary `n>=2` and set

```text
w=na+(n+1)b,                 M=w+1,
A_0=a+2b,                    S=a+b,
psi=u^A_0(u+1)^S.                                           (16)
```

Then `p^a y^b=x^-w psi`. On `x!=0`, criticality is exactly

```text
x^M=wc psi,
g x^(M-1)=-c psi'.                                       (17)
```

The consecutive-power compatibility has forced factors from both colors.
Indeed, raising the second row of `(17)` to the `M`-th power and the first
to the `w`-th power gives

```text
(-c psi')^M-(wc psi)^w g^M=0.
```

Since

```text
psi'=u^(A_0-1)(u+1)^(S-1)((2a+3b)u+A_0),
```

removing the common nonzero scalar and the forced color powers gives the
following residual address polynomial:

```text
E_(a,b)=w^w u^(M-A_0)(u+1)^(M-S)g^M
         -(-1)^M c((2a+3b)u+A_0)^M.                      (18)
```

Both stripped exponents are positive:

```text
M-A_0=(n-1)(a+b)+1,
M-S=(n-1)a+nb+1.                                         (19)
```

If `g!=0` has degree `e`, the first term of `(18)` exceeds the second in
degree by

```text
(n-2)(a+b)+1+M e>0.                                      (20)
```

Thus `E_(a,b)` is nonconstant. Moreover

```text
E_(a,b)(0)!=0,                         E_(a,b)(-1)!=0.    (21)
```

Choose a root `u_0` outside the two colors. If `g(u_0)!=0`, compatibility
of the two nonzero consecutive powers in `(17)` reconstructs a unique
nonzero `x_0`. If `g(u_0)=0`, equation `(18)` forces

```text
u_0=-A_0/(2a+3b),
```

which is distinct from both colors `-1,0` in characteristic zero; there
`psi'(u_0)=0`, and one may choose any nonzero `M`-th root from the first
row of `(17)`. When `g=0`
identically, start directly with this same address. In every case
`t_0=u_0/x_0^n` is finite and `(x_0,t_0)` is critical.

For `c=0`, identity `(15)` proves submersivity. Thus the sharper endpoint
is

```text
x+h(x^n t)+c p^a y^b is a submersion iff c=0             (22)
```

for every fixed `(a,b)!=(0,0)`. This includes every pure power of `p` or
`y` and every multiplicative mixed monomial.

The leading coefficient of `x` can harmlessly be any nonzero scalar: divide
the whole first coordinate by that scalar and apply `(22)`.

## 3. The height-two generic time form

Now take `n=2`, the monomial `(a,b)=(0,1)`, and `c!=0`; let the target
parameter `q` be transcendental over `k`, and set
`K=k(q)`. In `K(x,u)` the generic fibre is

```text
F=x^4+(h(u)-q)x^3+c u^2(u+1)=0.                          (23)
```

It is geometrically integral. Indeed, in `k(x,t)` use the `t`-adic
valuation. Its restriction to `k(A)` is trivial, and reduction sends

```text
A=x+h(x^2t)+cy                         to              x+h(0). (24)
```

Thus the residue identifies `k(A)` with the full field `k(x)`. If `E` were
a finite algebraic intermediate field between `k(A)` and `k(x,t)`, the
restricted valuation on `E` would still be trivial: an element algebraic
over a trivially valued field cannot have positive or negative value.
Reduction would inject `E` into `k(x)` over the already-surjective copy of
`k(A)`, forcing `E=k(A)`. Hence `k(A)` is relatively algebraically closed
in `k(x,t)`. Characteristic zero makes the extension regular, proving the
claim.

Since

```text
dx wedge dt=x^-2 dx wedge du,              F_x=x^3 A_x  on F=0,
```

a rational mate `J(A,Q)=1` would restrict to a primitive of the time form

```text
omega=x du/F_x=du/[x(4x+3(h-q))].                        (25)
```

We now compute its complete divisor on the smooth projective generic fibre.

### 3.1 Finite normalization places

There are exactly two points with `x=0`. At `P_0=(0,0)`, the local leading
equation is

```text
(h(0)-q)x^3+c u^2=0,
```

so the unique normalization branch has

```text
ord(x)=2,             ord(u)=3,             ord(omega)=0. (26)
```

At `P_-=(0,-1)`, the factor `u+1` is simple and `x` is a parameter:

```text
ord(x)=1,             ord(u+1)=3,           ord(omega)=1. (27)
```

The geometric generic source fibre is smooth by generic smoothness. At
every other finite point `x!=0`, either `F_x` or `F_u` is a unit. Using

```text
omega=x du/F_x=-x dx/F_u                                 (28)
```

in the appropriate chart shows that `omega` is a unit. Thus `(26)--(27)`
exhaust all finite zeros and poles.

### 3.2 Complete infinity ledger

Use the convention for `d` in Statement 3 and, for `d>=1`, write `h_d`
for the leading coefficient of `h`.
If `d=0`, the unique infinity face is

```text
x^4+c u^3=0.
```

Coprimality of `4,3` gives one normalization place with

```text
ord(x)=-3,             ord(u)=-4,            ord(omega)=1. (29)
```

Suppose `d>=1` and put `gamma=gcd(3,d)`. The Newton boundary has exactly
two unbounded faces. The first is

```text
x^4+h_d u^d x^3=0,
```

whose primitive edge gives one place `P_1` with

```text
ord_(P_1)(x)=-d,        ord_(P_1)(u)=-1,
ord_(P_1)(omega)=2d-2.                                  (30)
```

The second face is

```text
h_d u^d x^3+c u^3=0.                                    (31)
```

Its edge has lattice length `gamma`; the primitive binomial has exactly
`gamma` distinct roots in characteristic zero. Hence it gives precisely
`gamma` places `P_(2,j)`, with

```text
ord(u)=-3/gamma,        ord(x)=(d-3)/gamma,
ord(omega)=2d/gamma-1.                                  (32)
```

These faces are complete. Indeed, at a place where `u` has a pole, the
three possible dominant values come from `x^4`, `h_d u^d x^3`, and
`c u^3`. Pairwise balance of the first and third is defeated by the middle
term. The other two balances are exactly `(30)` and `(32)`. If `x` has a
pole while `u` is finite, `x^4` is uniquely dominant, which is impossible.
Thus no infinity branch is omitted. The face binomials are primitive and
squarefree, so no address collision or hidden multiplicity occurs.

The source coordinates satisfy `t=u/x^2`. Consequently the source generic
fibre omits `P_0,P_-` and all the infinity-face places. This proves the
place counts in `(4)`.

### 3.3 Divisor, genus, and nonexactness

Equations `(26)--(32)` give

```text
d=0:
  div(omega)=P_-+P_infinity;

d>=1:
  div(omega)=P_-+(2d-2)P_1
              +sum_(j=1)^gamma (2d/gamma-1)P_(2,j).      (33)
```

Every coefficient is nonnegative, so `omega` is a nonzero holomorphic
differential. The displayed divisor degrees are respectively

```text
2                                      if d=0,
4d-1-gamma                             if d>=1.           (34)
```

Since a canonical divisor has degree `2g-2`, equation `(34)` gives exactly
the genus ladder `(4)`. Equations `(33)--(34)` also prove there are no
unlisted zeros.

A nonzero holomorphic differential on a smooth projective curve cannot be
exact: if `omega=dQ`, every pole of `Q` would create a pole of `dQ`, so `Q`
would be globally regular and hence constant. More strongly, on a finite
cover a point of ramification index `e` changes an order `v` to

```text
e v+(e-1)>=0.                                             (35)
```

The pullback remains nonzero holomorphic and the same argument applies.
Therefore `omega` in `(25)` has no rational primitive, even after finite
algebraic extension.

## 4. Structural consequence and scope

The three panels identify a sharp leading-`x` boundary-coupling wall.
Polynomial shears `x+h(u)` are submersions, but adding either generator
linearly, adding both linearly, or adding any one monomial `p^a y^b` always
produces an affine critical point. At height two, the linear-y failure is
also visible on the generic fibre as a holomorphic time-form obstruction of
explicitly growing genus. THM-3983 supplies the complementary place-budget
view: the restriction to `D` is `h(-1)+c y`, and the actual fibre pays both
positive genus and the listed infinity places.

What remains includes sums of two or more nonlinear boundary weights,
arbitrary polynomials in `p,y`, and cells without a nonzero leading-`x`
term. The latter distinction is essential: at height two the adjacent family
`a p+c y^m` with `a*c!=0,m>=2` has no leading `x` and lies outside this
theorem. No general first coordinate in `B_n`, and no consequence for
unrestricted `JC(2)`, is claimed.

**QED candidate.**

## Reproduction

```bash
python3 04-computation/jc2_boundary_generator_coupling_criticality_thm3984.py
python3 -O 04-computation/jc2_boundary_generator_coupling_criticality_thm3984.py
```
