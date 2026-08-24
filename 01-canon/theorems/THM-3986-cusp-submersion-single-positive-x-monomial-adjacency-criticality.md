---
id: THM-3986
title: "Every single positive-x monomial adjacent to the cusp submersion is critical"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. At height
  two, alpha*p+gamma*y^m is a source-plane submersion for every m>=2 and
  nonzero alpha,gamma, but adjoining any one nonzero monomial
  lambda*x*p^r*y^s makes it critical, uniformly for all r,s>=0. The
  logarithmic critical matrix has a universal one-dimensional kernel. Its
  common-power compatibility map is surjective; a complete color and
  vanishing-kernel ledger isolates three regular u=0 families and one
  L=N family, all of whose tuned invalid roots leave a valid off-color
  root. The three low-L rows x,xp,xy close separately. Criticality excludes
  regular and polynomial mates only; no rational-mate obstruction or
  Jacobian counterexample is claimed.
source: jc-mixed-generator-submersion / post-THM-3984 cusp-adjacency lane, 2026-08-24
audit: >
  PASS (audit-thm3986-all-m + root / jc-cohn3709, 2026-08-24). Two
  independent audits rederived the logarithmic kernel, all compatibility
  exponents and merged-address orders, the three tuned u=0 repairs, the
  all-m L=N kernel-root parametrization and its simple-root escape, and the
  separate x, xp, and xy branches. They also checked the positive base wall
  and the critical-versus-rational-mate scope. Normal, optimized, and frozen
  outputs byte-match at CHECKS=42281; hashes agree.
depends_on:
  - THM-3973-exact-volume-simple-cubic-determinantal-affine-plane-completion
related:
  - THM-3984-boundary-generator-coupling-criticality-and-holomorphic-time-form
  - THM-3985-cusp-plane-rational-time-residue-and-height-two-mixed-submersion
script: 04-computation/jc2_cusp_submersion_positive_x_monomial_adjacency_thm3986.py
output: 05-knowledge/results/jc2_cusp_submersion_positive_x_monomial_adjacency_thm3986.out
script_sha256: 237cf03e0a000d2cb8160cac8f4167260ad220ef025615b1e48b663972562cec
output_sha256: ff65f73eeaebb1ad0193794d55ed1222630cb4afa2724831bc07fcdf4e67040f
semantic_sha256: 869bca81569a6fa2467e7d4d80d6a724be0c5d37ecc3982c99e55583019f493c
hash_basis: raw LF bytes
---

# THM-3986 -- every one-monomial positive-x adjacency is critical

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero. In the height-two
member of THM-3973 put

```text
u=x^2t,                 rho=u(u+1),              sigma=u^2(u+1),
p=x^-2 rho=t(1+u),      y=x^-3 sigma=xt^2(1+u).              (1)
```

For every `m>=2`, `r,s>=0`, and `alpha*gamma*lambda!=0`, the base row

```text
A_m=alpha*p+gamma*y^m                                      (2)
```

is an affine-plane submersion, whereas its single positive-`x` monomial
adjacency

```text
A_(m;r,s)=alpha*p+gamma*y^m+lambda*x*p^r*y^s               (3)
```

has an affine critical point with `x!=0`. Thus the discontinuity at
`lambda=0` is exact: every one-monomial adjacency in `(3)` destroys the
entire all-degree submersion family `(2)`.

Consequently no polynomial `Q in k[x,t]`, and hence no regular completion
coordinate, can satisfy `J(A_(m;r,s),Q) in k*`. This implication is only a
**regular-mate** obstruction. A rational `Q` may have a pole at a critical
point, so no rational-mate conclusion and no Jacobian counterexample is
claimed here.

## 1. The unperturbed family is genuinely positive

The exact height-two bracket is

```text
J_(x,t)(p,y)=-t p.                                        (4)
```

Since `J(y^m,y)=0`, a critical point of `(2)` would satisfy `tp=0`.
On `t=0`, however, `(A_m)_t=alpha`. On the other component of `p=0` one
has `u=-1`, with `x,t!=0` and `y=0`; because `m>=2`,

```text
dA_m=alpha dp,                    p_x=2xt^2!=0.            (5)
```

Thus `(2)` has no affine critical point. This short argument is included so
that the positive endpoint does not depend on the still-under-audit status of
the related THM-3985.

## 2. A universal logarithmic kernel

Work first off the two colors `u=0,-1`, where `(x,u)` are source coordinates
and every monomial below is nonzero. Write

```text
P=alpha*p,       V=gamma*y^m,       W=lambda*x*p^r*y^s,
w=2r+3s-1,       L=w-2=2r+3s-3,    N=3m-2.               (6)
```

After multiplying the `u`-derivative by `u(u+1)`, the two critical equations
are

```text
2P+3mV+wW=0,
(2u+1)P+m(3u+2)V+((w+1)u+r+2s)W=0.                       (7)
```

The cross product of the two coefficient rows is the unexpectedly simple
universal kernel

```text
P:V:W=m(3u+2-r):-(2u+s+1):m.                             (8)
```

Put

```text
C_r=3u+2-r,                  D_s=2u+s+1.                 (9)
```

Comparing the first two and the first and third coordinates in `(8)` gives

```text
x^N=-(m gamma/alpha) u^(2m-1)(u+1)^(m-1) C_r/D_s,
x^L=(lambda/alpha) C_r u^(r+2s-1)(u+1)^(r+s-1).          (10)
```

Conversely, at an address outside

```text
u=0,       u=-1,       C_r=0,       D_s=0,               (11)
```

any common nonzero solution `x` of `(10)` realizes `(8)`, hence both rows
in `(7)` vanish. It gives a genuine source point `t=u/x^2`.

## 3. The positive-L compatibility map

Assume `L>0` and let `delta=gcd(N,L)`. Two prescribed nonzero powers
`x^N=X`, `x^L=Y` have a common solution over `k` if and only if

```text
X^(L/delta)=Y^(N/delta).                                  (12)
```

For `(10)`, this is the fibre equation `H_(m;r,s)(u)=kappa`, where `kappa`
is a nonzero coefficient scalar and

```text
H_(m;r,s)=
 u^((mr+s-3m+1)/delta)
 (u+1)^((1-mr-s)/delta)
 C_r^((2r+3s-3m-1)/delta)
 D_s^(-L/delta).                                         (13)
```

All four exponents are integers because their numerators arise directly as
integer combinations of `N` and `L` in `(12)`. Their sum is

```text
-2N/delta.                                                (14)
```

Thus `H(u)` has infinity exponent `-2N/delta`, equivalently a zero of order
`2N/delta` at infinity. It is nonconstant. When `s=1`, the factor `D_s`
merges with `u+1`, but the combined order is still nonzero. Hence the
nonconstant rational map `H:P1->P1` is surjective and every nonzero fibre is
nonempty.

What remains is not existence but **address debt**: could every preimage of
the required `kappa` lie in `(11)`? The exact merged-order ledger says no.

* At `u=-1` the order never vanishes. At the root of `D_s` it is `-L/delta`,
  except for the already merged `s=1` row. The roots of `C_r` and `D_s`
  never coincide, since that would require `2r+3s=1`.
* The only regular `u=0` rows are

  ```text
  (r,s)=(2,m-1),       (1,2m-1),       (0,3m-1).          (15)
  ```

* A root of `C_r` is regular precisely when `L=N`, or

  ```text
  2r+3s=3m+1.                                             (16)
  ```

Every other forbidden address is a strict zero or pole of `H`, and is
therefore absent from the nonzero fibre.

### 3.1 Repairing the three regular color rows

After absorbing a harmless nonzero scalar into the fibre value, the three
rows in `(15)` reduce respectively to

```text
(u+1)(2u+m)=kappa,
(u+1)[2(u+m)]^2=kappa(3u+1),
(u+1)(2u+3m)^3=kappa(3u+2)^2.                             (17)
```

If `u=0` lies in the required fibre, the three tuned values of `kappa` are

```text
m,                  4m^2,                  27m^3/4.       (18)
```

The cleared equations in `(17)` have degrees `2,3,4`, and their derivatives
at the tuned root are

```text
m+2,                8m(1-m),              54m^2(1-m).    (19)
```

They are nonzero for `m>=2`. Thus `u=0` is only a simple root and at least
one further root remains. Directly in `(17)`, the other addresses in `(11)`
make exactly one side zero, so every such further root is valid. In the first
row this is especially visible from

```text
(u+1)(2u+m)-m=u(2u+m+2).                                  (20)
```

### 3.2 Repairing every regular kernel-zero row

The solutions of `(16)` have a unique parametrization

```text
r=3q+2,       s=m-1-2q,       0<=q<=floor((m-1)/2).      (21)
```

The case `q=0` is the first row of `(15)`. For `q>=1`, the forbidden root of
`C_r` is `u=q`, and `(13)` simplifies up to a nonzero scalar to

```text
H(u)=u^q/[(u+1)^(q+1)(2u+s+1)].                          (22)
```

If the required fibre is tuned to contain `u=q`, its logarithmic derivative
there is

```text
(H'/H)(q)=q/q-(q+1)/(q+1)-2/(2q+s+1)=-2/m!=0.            (23)
```

So the invalid root is simple. The cleared fibre equation

```text
u^q=kappa (u+1)^(q+1)(2u+s+1)                            (24)
```

has degree `q+2`; it therefore has another root. Its other forbidden
addresses are strict zeros or poles in `(22)`, so that root is valid. This
closes the extra kernel-zero seam that first appears for `m>2`.

Sections 3.1--3.2 prove criticality for every `L>0`.

## 4. The three low-L rows

For nonnegative `r,s`, the inequality `L<=0` leaves exactly

```text
(r,s)=(0,0),                 (1,0),                 (0,1),
correction=x,                correction=xp,         correction=xy. (25)
```

### 4.1 The x row

For `(r,s)=(0,0)`, the kernel `(8)` gives

```text
x^N=-(m gamma/alpha)u^(2m-1)(u+1)^(m-1)(3u+2)/(2u+1),
x^3=(alpha/lambda)u(u+1)/(3u+2).                         (26)
```

Because `N=3m-2` is coprime to `3`, their common-power condition is again a
single rational fibre. After discarding its nonzero coefficient, its factor
exponents on

```text
u,          u+1,          3u+2,          2u+1
```

are

```text
3m-1,       -1,           3m+1,          -3.             (27)
```

All are nonzero, and the infinity exponent is `2N`. Surjectivity therefore
gives a fibre point away from all four invalid addresses, and `(26)`
reconstructs a common nonzero `x`.

### 4.2 The xp row

Here

```text
A=p(alpha+lambda x)+gamma y^m.                            (28)
```

At either color `u=0,-1`, the choice

```text
x=-alpha/lambda                                           (29)
```

annihilates both derivatives. The `y^m` row vanishes to first order because
`m>=2`. This gives an explicit affine critical point.

### 4.3 The xy row

Put `f(u)=rho(alpha+lambda u)`. Then

```text
A=x^-2 f(u)+gamma*x^(-3m)*sigma(u)^m.                     (30)
```

Eliminating the two summands from the critical equations cancels `m` and
gives

```text
3f' sigma-2f sigma'=0
iff alpha-2lambda u-3lambda u^2=0.                        (31)
```

The quadratic never vanishes at `u=0`. It vanishes at `u=-1` only when
`lambda=alpha`; in that tuned case its other root is `u=1/3`. Thus it always
has a root away from both colors. At such a root `f!=0`: a simultaneous
zero of `alpha+lambda u` and `(31)` would force `lambda=alpha,u=-1`. Hence

```text
x^N=-3m gamma sigma^m/(2f)                                (32)
```

reconstructs a nonzero `x`. This closes the last low row and proves the
theorem.

## 5. Exact companion and scope boundary

The companion script verifies without Python `assert`:

* the source bracket and base submersion controls for `2<=m<=10`;
* the symbolic kernel, the four compatibility exponents, and the infinity
  sign;
* `42,000+` exact address-order checks over `2<=m<=12` and
  `0<=r,s<=3m+4`, recovering exactly `(15)--(16)`;
* all three color repairs, the uniform simple-root identity `(23)`, the
  `x,xp,xy` rows, and independent binomial-resultant controls.

Normal and optimized runs are required to byte-match the frozen output.

The proved cell contains exactly one appended monomial `x*p^r*y^s` and no
additional `h(u)` or coordinated sum. Several positive-`x` monomials can
cancel one another's logarithmic principal symbols, so arbitrary sums,
the full completion `B_2`, rational mates, and JC(2) all remain open.
