# Sharp finite recognition of the planar source depth modules

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
This report proves a bounded-support recognition statement, not a degree
bound for a hypothetical Keller pair. Planar JC remains open.

## 1. Inheritance and the new consumer

Use the source normalization of
[THM-4308 / source-normal-bracket-hasse-truncation-through-row-eight](../../01-canon/theorems/THM-4308-source-normal-bracket-hasse-truncation-through-row-eight.md):

```text
u=x^2*t, p=t*(1+x^2*t), y=x*t*p,
P_d=span_K{x^a*u^b*p^c*y^e : a,b,c,e>=0, a+b<=d}.
```

Here `K` has characteristic zero and `d>=0`. By that theorem and
[THM-3989 / cusp-log-laurent-conductor-and-nondividing-depth-reduction](../../01-canon/theorems/THM-3989-cusp-log-laurent-conductor-and-nondividing-depth-reduction.md),
this is the actual depth module in the stated cusp chart.

The complete principal polynomial image is inherited from
[THM-4369 / source-packet-pascal-circuit-kernel-and-boundary-basis](../../01-canon/theorems/THM-4369-source-packet-pascal-circuit-kernel-and-boundary-basis.md),
with its bound `R=ell+d-ceil(ell/2)`. It is not claimed as a new result.
THM-4364 and THM-4368 provide the simplex annihilators and their independent
bank ranks. The new consumers below identify that bank with the *full*
truncated depth annihilator and prove an optimal finite recognition bound
when the polynomial's support is fixed.

The hostile is `t^T`: a finite prefix can lie in a projected depth module
for substantially longer than its own support, although the polynomial
itself does not belong to the depth module. The corrected near miss is to
confuse arbitrarily extendable prefixes with a fixed polynomial. The needed
sidecar is the declared support cutoff, including its zero later rows.

## 2. Complete diagonal description

For an integer `ell`, set

```text
s=max(0,ceil(ell/2)), rho=max(0,ceil(ell/3)),
D=ell+d-s, L=D-rho,
F_ell(z)=sum_(n>=s) (-1)^(n-s) [x^(2n-ell)t^n]F * z^(n-s).
```

Assume the necessary row caps `deg_x F_n<=n+d`. Then only `ell>=-d`
occurs, and the last allowed row on the diagonal is `ell+d`.
The exact diagonal image, including its integral version, is

```text
V_(ell,d)=(1-z)^rho K[z]_(<=L),                         (1)
```

where a negative polynomial degree bound means zero. Diagonal images are
independent direct summands: every generator lies on one entire diagonal.
Thus a finite polynomial `F` belongs to `P_d` exactly when each `F_ell`
belongs to (1).

For completeness, a short proof also resolves the small intercepts absent
from THM-4369. A source monomial has diagonal `ell=2c+3e-a`, and its signed
trace is

```text
(-1)^(n0-s) z^(n0-s)(1-z)^(c+e), n0=b+c+2e.
```

Its degree is `ell+a+b-s<=D`, and `c+e>=rho`. For `ell>=2,d=0`, put
`f=floor(ell/2)`. The solutions of `2c+3e=ell` give precisely the Bernstein
polynomials

```text
(-1)^j (1-z)^rho z^j(1-z)^(f-rho-j), 0<=j<=f-rho.
```

After removing the common factor, these are an integral basis of all
polynomials of degree at most `f-rho`: their coefficient matrix is triangular
with diagonal entries `+/-1`. Multiplication by `u` preserves `ell`,
increases depth by one, and multiplies the trace by `-z`. Since `P_(d-1)`
and `u*P_(d-1)` lie in `P_d`, induction proves (1). At `ell=1`, the first
nonzero module is generated at depth one by `xp`, with trace `1-z`; depth
zero is zero. At `ell<=0`, the initial generator is `x^(-ell)` at depth
`-ell`, with constant trace; multiplication by `u` proves the assertion.

## 3. The simplex bank is the full truncated annihilator

For a projection through row `m`, the ambient diagonal has rows

```text
s<=n<=min(m,ell+d).
```

Its dimension is `a=max(0,min(m,ell+d)-s+1)`. Formula (1), whose monomial
multiplication matrix has a leading unit triangular block, gives

```text
rank pi_m(V_(ell,d)) = min(a,max(0,L+1)).                (2)
```

For `ell>=2` and `s<=m<=ell+d`, the full annihilator has the consecutive
simplex basis

```text
L_(m,ell,q)=sum_(n=s)^m (-1)^(n-s) binom(m+q-n,q)F_(n,2n-ell),
max(0,ell+d-m)<=q<=rho-1.                              (3)
```

Each row annihilates by THM-4364. Their number is exactly the codimension
in (2), and THM-4368 supplies a unit minor on consecutive orders. Therefore
they form the entire annihilator, also over the integers. The small
intercepts are governed directly by (1): `ell<=0` has no extra constraints;
`ell=1` has the same formula with `rho=1`, including its zero depth-zero
module. When `m>ell+d`, use (3) at `m=ell+d` and the row caps; omitted
coordinates above that endpoint must be zero.

In particular the complete matrix ranks in the canon follow from (2):

| Row m | Depth d | Ambient dimension | Image rank | Nullity |
|---:|---:|---:|---:|---:|
| 8 | 2 | 63 | 51 | 12 |
| 8 | 3 | 72 | 63 | 9 |
| 14 | 2 | 150 | 108 | 42 |
| 14 | 3 | 165 | 129 | 36 |
| 15 | 2 | 168 | 119 | 49 |
| 15 | 3 | 184 | 142 | 42 |

These are ranks of the unrestricted depth projection. They must not be
confused with its restriction to a bracket-compatible parameter family.

## 4. Exact finite recognition theorem

**Theorem.** Let `T>=1`, `d>=0`, and let

```text
F(x,t)=sum_(n=0)^T F_n(x)t^n, deg F_n<=n+d.
N(T,d)=floor(4T/3)+d+1.
```

Then

```text
F in P_d iff pi_N(F) in pi_N(P_d).                     (4)
```

The universal cutoff is optimal: `F=t^T` passes the projected test through
row `N(T,d)-1` but fails at row `N(T,d)` and does not belong to `P_d`.
For `T=0`, every polynomial satisfying the cap is already in `P_d`.

The zero coefficients of `F` after row `T` are part of the input to (4).
Allowing them to vary would change the theorem.

### Proof of the upper bound

Fix a nonzero diagonal with `rho>=1` and write its polynomial as `f(z)`,
of degree at most `M=min(T,ell+d)-s`. Put `L=D-rho` as above. Membership
means that

```text
g(z)=f(z)/(1-z)^rho
```

is a polynomial of degree at most `L`. A projected test through relative
degree `R` says exactly that `g_(L+1),...,g_R` vanish. If `M>=rho-1`,
testing through `R=D=L+rho` is sufficient by (1). If `M<rho-1`, the first
`M+1` such equations already force `f=0`. Indeed, their square coefficient
matrix, with `k=M+1<=rho` and `B=L+1>=0`, is

```text
H_(i,j)=binom(B+i-j+rho-1,rho-1), 0<=i,j<k,             (5)
```

where coefficients with `B+i-j<0` mean zero. Its determinant is

```text
det H = product_(i=0)^(k-1)
          (rho+B-1-i)! i! / ((rho-1-i)! (B+i)!),         (6)
```

which is a positive integer, hence nonzero in `K`. One proof of (6) is
induction by the Desnanot-Jacobi identity: write the determinant as `D_k(B)`;
the contiguous minors give

```text
D_k(B) D_(k-2)(B)
 =D_(k-1)(B)^2-D_(k-1)(B-1)D_(k-1)(B+1).
```

The displayed product satisfies this identity by cancellation of factorials;
`D_0=1`, `D_1(B)=binom(B+rho-1,rho-1)`, and `D_k(0)=1` supply the boundary
values. All divisors used in the induction are nonzero for `k<=rho`.

Thus the diagonal is decided by absolute row

```text
n_ell=ell+d-rho+min(M+1,rho).                           (7)
```

For `ell>=1`, its value is at most

```text
T+d+1+floor(ell/2)-ceil(ell/3)
 <=T+d+1+floor(ell/6)
 <=T+d+1+floor(T/3)=N(T,d),                             (8)
```

because `ell<=2T`. Intercepts `ell<=0` impose only the already assumed
caps. This proves (4). The `ell=1,d=0` zero module has `L=-1,B=0`; the
same triangular matrix argument includes it.

### Sharpness

For `F=t^T`, the only diagonal is `ell=2T`, with `s=T`,
`rho=ceil(2T/3)`, `D=T+d`, `L=T+d-rho`. Its signed polynomial is `1`.
Since `rho>=1`, it does not belong to (1). But the source image contains

```text
(1-z)^rho * sum_(j=0)^L binom(rho+j-1,j)z^j
 =1-binom(rho+L,L+1)z^(L+1)+O(z^(L+2)).                (9)
```

The first mismatch is at absolute row

```text
T+L+1=2T+d-ceil(2T/3)+1=N(T,d),                        (10)
```

with a nonzero integer coefficient. Formula (1) supplies an actual source
polynomial with the trace (9), so every earlier projected test passes.
This is a non-vacuous sharpness construction, not just absence of a known
annihilator.

## 5. The same sharp theorem at every completion height

The height parameter below indexes *surfaces*, not the dimension in the
Jacobian conjecture. For `h>=2`, use the family of
[THM-3973 / exact-volume-simple-cubic-determinantal-affine-plane-completion](../../01-canon/theorems/THM-3973-exact-volume-simple-cubic-determinantal-affine-plane-completion.md):

```text
u=x^h*t, p=t*(1+u), y=x^(h-1)*t^2*(1+u),
B_h=K[x,u,p,y],
P_(h,d)=B_h intersect {F: deg_x F_j <= (h-1)j+d for every j}.
```

At `h=2`, this is exactly the inherited depth module `P_d`. At other
heights, the definition is the stated weighted degree filtration; no
unproved identification with a different geometric pole depth is needed.
THM-3973, Section 5, proves the entire negative graded piece, with grading
`wt(x)=1, wt(t)=-h`, to be

```text
(B_h)_(-ell)=x^(-ell)u^ceil(ell/h)(1+u)^ceil(ell/(h+1))K[u]
                                                        for ell>0. (11)
```

For nonpositive `ell`, the piece is `x^(-ell)K[u]`. The signed diagonal
description (1) therefore holds with

```text
ell=h*j-r, s=max(0,ceil(ell/h)),
rho=max(0,ceil(ell/(h+1))), D=ell+d-s.
```

Consequently, for a polynomial of t-degree at most `T>=1` satisfying these
caps, the optimal universal recognition row is

```text
N_h(T,d)=floor(h^2*T/(h+1))+d+1.                       (12)
```

The argument through (7) is unchanged. For the upper bound, `ell<=h*T`
and

```text
n_ell <= T+d+1+ell-ceil(ell/h)-ceil(ell/(h+1))
       <= T+d+1+ell*(1-1/h-1/(h+1))
       <= h^2*T/(h+1)+d+1.
```

Taking the floor proves (12). The coefficient `1-1/h-1/(h+1)` is positive
for `h>=2`, so the continuous upper bound has its maximum at `ell=h*T`;
no monotonicity of the stepped ceiling expression is assumed.
For sharpness, `t^T` has
that intercept, `s=T`, `rho=ceil(h*T/(h+1))`, and the first mismatch in
(9) is at `h*T+d-rho+1=N_h(T,d)`. This proves the all-height formula.

This extension uses the actual full graded ring (11). It does not assume
that every element in the intersection filtration has an expression using
only individually capped generator monomials at other heights.

Characteristic zero is essential. At height two, `T=3,d=0`, the actual
source polynomial

```text
p^3-3y^2=t^3-3x^4*t^5-2x^6*t^6
```

equals `t^3` through the asserted cutoff `N=5` in characteristic three,
although the constant signed diagonal of `t^3` is not divisible by
`(1-z)^2`. It fails at row six. This is an actual-source hostile to an
unqualified characteristic-free version of the recognition theorem.

## 6. Connection to the requested S-matrix paper

[Li, arXiv:2608.29750v1, Sections 2 and 4](https://arxiv.org/html/2608.29750v1)
combines a defect budget with sufficiently many measurements and a Gram
projection. The concrete transfer here is to retain the whole diagonal
measurement vector and prove its exact kernel, instead of testing a selected
scalar debt. The replacement for a real energy contradiction is the
nonzero *integer* determinant (6). Positivity is used to prove that this
fixed integer is nonzero, not to order arbitrary complex source coefficients.

Source: raw source-normal depth coefficients. Target: truncated polynomial
division by `(1-z)^rho`. Map: signed diagonal extraction and inverse
binomial transform. Preserved predicate: exact projected depth membership.
Lost coordinate: later coefficients unless the support cutoff is supplied.
Sidecar: `T` and its prescribed zero tail. Decisive hostile: `t^T`, whose
first detection row is exactly (10).

No part of the S-matrix theorem is a mathematical dependency of (4).
The bridge is an explicitly implemented proof mechanism.

## 7. Scope and next use

For candidate source coordinates `A,C` with fixed t-degrees `T_A,T_C`,
membership in `P_2,P_3` is completely decidable at rows

```text
floor(4T_A/3)+3, floor(4T_C/3)+4,
```

provided their later rows are prescribed to be zero. This removes the
all-row depth ambiguity for a *given finite candidate*. It supplies neither
a bound on `T_A,T_C`, nor the bracket equations, nor a reason that a partial
source solution terminates. Arbitrary chart entry and planar JC remain open.

The [independent audit](planar_jc_long_20260906_memory_depth_audit.md)
accepts the proof, all-height extension, exact full bank, sharpness, and
characteristic boundary. Its separate source-generator/Hasse-jet method
does not import the producer or its reciprocal-series recognition engine.
It performs **828,684** live exact checks; the
[producer](../../04-computation/planar_jc_long_20260906_depth.py) performs
**45,048**. Both normal and optimized runs reproduce their respective
frozen outputs exactly. The finite universes and controls are declared
in the sources and audit; the universal quantifiers rest on the proofs.

```text
python3 -B 04-computation/planar_jc_long_20260906_depth.py
python3 -B -O 04-computation/planar_jc_long_20260906_depth.py
python3 -B 04-computation/planar_jc_long_20260906_memory_depth_audit.py
python3 -B -O 04-computation/planar_jc_long_20260906_memory_depth_audit.py
```

The session manifest pins the source/output bytes. The inherited principal
image remains credited to THM-4369 and the all-height grading to THM-3973.
