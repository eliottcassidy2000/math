---
id: THM-3873
title: "The first non-graph triangular parabola forces a projective companion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  A first genuinely
  non-graph polynomial A1
  branch in the marked-root universe has normalization
  (A,C)=(r^2,6(r+r^4)).  For every polynomial transverse quotient, a unique
  formal-profile first mismatch produces a distinct reduced companion whose
  affine normalization misses at least two projective points.  The N=0 row
  is treated separately; for N>=1 the only degree resonance is d=N-1 and
  has a nonzero perfect-square leading coefficient.
source: jc_sparse_direct_search / post-THM-3866--3870 non-graph A1 frontier, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root, 2026-08-23).  The audit reconstructed
  the marked profile in the normalization UFD; checked that formal monic
  division is valid by even/odd T-adic remainder separation; rederived the
  N=0 signs and all N>=1 degree regimes; and verified that the selected
  reduced factor has positive T-degree and is neither vertical nor the
  marked parabola.  The assertion-free companion checks
  the triangular normalization and inverse, fixed-coordinate non-graph
  seams, formal quotient coefficients q_0,...,q_6, the exceptional N=0
  residual, and every degree regime around d=N-1 for N=1,...,6.  Normal and
  optimized runs byte-match the frozen 244-gate transcript.
related:
  - THM-3866-all-polynomial-graph-branches-force-projective-companion
  - THM-3870-vertical-axis-polynomial-profile-first-mismatch-companion
script: 04-computation/jc2_first_nongraph_parabola_thm3873.py
output: 05-knowledge/results/jc2_first_nongraph_parabola_thm3873.out
script_sha256: ad764c9b37f4a3b0f5a02073ed601814410f276dba49e212793afaee0f890b22
output_sha256: ba110dce7516cc33765506c298b9c5d1c63dbb653787f3d9aaa33ed5b5fbfc25
semantic_sha256: 9edbf07d9ffd0ac57f86bf28012d5a83e56959f8740453e17b34d02b5a5ce275
hash_basis: raw LF bytes
---

# THM-3873 -- the first non-graph A1 branch has no polynomial escape

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an algebraically closed field `k` of
characteristic zero.  Put

```text
Delta_b=-27A^2b^2+8AC^3-54ACb+9C^2-54b.                         (1)
```

Define the triangular coordinate, branch, and marked profile

```text
T=C/6-A^2,                     C=6(T+A^2),
F=A-T^2,
b0=6A+8A^2T.                                                   (2)
```

Suppose `b in k[A,C]=k[A,T]` and `V(Delta_b)` contains `Gamma=V(F)`.
Then there is a unique `Q in k[A,T]` with

```text
b=b0+FQ.                                                        (3)
```

For every polynomial `Q`, the reduced branch packet contains an irreducible
component `Xi` different from `Gamma` whose affine normalization omits at
least two distinct points of its smooth projective normalization.  One lies
over finite `A=0,T=infinity` -- equivalently `A=0,C=infinity` -- and another
lies over `A=infinity`.  No degree bound is imposed on `Q`.

The branch `Gamma` is a polynomial `A1`, but in the fixed coordinates
`(A,C)` it is neither a graph `C=g(A)` nor a graph `A=h(C)`.  Thus this is the
first all-degree closure beyond the forward, reverse, and vertical graph
families of THM-3866 and THM-3870.

## 1. Exact normalization and fixed-coordinate non-graph seams

The triangular change in `(2)` is a polynomial automorphism.  The quotient
is

```text
k[A,C]/(F)=k[A,T]/(A-T^2)=k[T].                                (4)
```

Writing its normalization parameter as `r` gives the mutually inverse data

```text
A=r^2,
C=6(r+r^4),
T=C/6-A^2=r.                                                    (5)
```

Equivalently, the fixed-coordinate equation is

```text
36A-(C-6A^2)^2=0.                                               (6)
```

This is genuinely outside both graph theorems.  The two parameters `r` and
`-r` have the same `A`, but their `C`-values differ by `12r`, so `C` is not
a polynomial in `A`.  Conversely, `deg_r A=2` and `deg_r C=4`; an identity
`A=h(C)` with nonconstant `h` would give `2=4 deg h`, which is impossible.

## 2. The marked profile and the full transverse universe

The cusp identity is

```text
A^2 Delta_b=27(P^3-u^2),
P=1+(2/3)AC,                    u=1+AC+A^2b.                     (7)
```

On `(5)`,

```text
P=(1+2r^3)^2.                                                   (8)
```

If `Delta_b` vanishes on `Gamma`, the UFD law in `k[r]` and `u(0)=1`
force `u=(1+2r^3)^3`.  Solving `(7)` gives

```text
b|_Gamma=2r^2(3+4r^3)
          =(6A+8A^2T)|_Gamma=b0|_Gamma.                        (9)
```

The kernel of `k[A,T] -> k[r]` is the principal ideal `(F)`, so `(9)` is
equivalent to the unique full transverse form `(3)`.  This is not a selected
ansatz: it is every polynomial profile containing this branch.

## 3. Formal division and its coefficient law

In `k[T][[A]]`, substitute `C=6(T+A^2)` into the unique binomial profile

```text
b_*={(1+(2/3)AC)^(3/2)-1-AC}/A^2
   =sum_(j>=0) beta_j A^j C^(j+2),                              (10)

beta_j=binom(3/2,j+2)(2/3)^(j+2) !=0.                          (11)
```

On `F=0`, one has `A=T^2`, and the marked-root calculation `(8)--(9)` gives
`b_*|_F=b0|_F`.  Division by the monic quadratic `T^2-A` therefore has zero
remainder.  Explicitly, a possible remainder has the form
`rho_0(A)+T rho_1(A)`; after `A=T^2` its even and odd parts vanish separately
in `k[[T]]`.  Hence there is a unique formal quotient

```text
Q_*=(b_*-b0)/F=sum_(j>=0)q_j(T)A^j,             q_j in k[T].   (12)
```

Its degrees and leading coefficients are uniform:

```text
deg_T q_j=j,
lc_T(q_j)=gamma_j=-beta_j 6^(j+2) !=0.                         (13)
```

To prove `(13)`, let `n_j(T)` be the coefficient of `A^j` in `b_*-b0`.
In `C=6(T+A^2)`, the unique top contribution comes from the `j`th summand
of `(10)`, while every use of an `A^2` term lowers the `T`-degree by three.
Thus

```text
deg_T n_j=j+2,              lc_T(n_j)=beta_j6^(j+2).           (14)
```

The coefficient of `A^j` in `FQ_*=(A-T^2)Q_*` is

```text
n_j=q_(j-1)-T^2q_j,                   q_(-1)=0.                 (15)
```

Polynomiality was proved by the remainder argument; `(14)--(15)` now give
`(13)` inductively.  The first values are

```text
q_0=-6,
q_1=4T,
q_2=-6T^2,
q_3=12T^3+6,
q_4=-28T^4-12T.                                                (16)
```

In particular `Q_*` is not a polynomial in `A,T`.

## 4. First mismatch and the truncation residual

Every polynomial `Q` in `(3)` has a unique first mismatch `N>=0`:

```text
Q=Q_<(N)+A^N H,
Q_<(N)=sum_(j=0)^(N-1)q_j(T)A^j,
h(T)=H(0,T) != q_N(T).                                         (17)
```

The empty sum is used at `N=0`, and `H` may be zero for an exact finite
truncation.  Put

```text
b_N=b0+FQ_<(N),
u_N=1+AC+A^2b_N.                                                (18)
```

The formal tail in `(12)` and the exact quadratic response in `b` give

```text
Delta_(b_N)=A^N F R_N,                                         (19)

R_N(0,T)=54q_N(T).                                              (20)
```

Here `R_N` is polynomial.  Indeed, `F` divides `Delta_(b_N)` because
`b_N|_F=b0|_F`, and the completed response shows that its polynomial
quotient is divisible by `A^N` coefficientwise.

The row `N=0` is exceptional but explicit:

```text
R_0=108(16A^6-24A^3-16AT-3),
deg_T R_0=1,                    lc_T(R_0)=-1728A.               (21)
```

For every `N>=1`, put `gamma=gamma_(N-1)`.  The top term of `b_N` comes
from `-T^2 A^(N-1)q_(N-1)`, so

```text
deg_T b_N=N+1,                 lc_T(b_N)=-gamma A^(N-1),
deg_T u_N=N+1,                 lc_T(u_N)=-gamma A^(N+1).        (22)
```

The square term in `(1)` is then the unique top term.  Dividing by the
degree-two factor `F`, whose `T`-leading coefficient is `-1`, gives

```text
deg_T R_N=2N,
lc_T(R_N)=27gamma^2 A^N.                                       (23)
```

## 5. Exact nonlinear response

For every perturbation `delta b`,

```text
Delta_(b_N+delta b)-Delta_(b_N)
=-54u_N delta b-27A^2(delta b)^2.                               (24)
```

Use `delta b=A^N F H` from `(17)`.  Then

```text
Delta_b=A^N F S,

S=R_N-54u_NH-27A^(N+2)F H^2,                                 (25)

S(0,T)=54[q_N(T)-h(T)] !=0.                                   (26)
```

Thus `A` does not divide `S`, so the first mismatch cannot disappear into
vertical multiplicity.  A possible additional factor `F` in `S` is harmless:
the projective contact below does not lie on `F`.

## 6. Every degree regime

### 6.1. The exceptional row N=0

If `H=0`, `(20)--(21)` give

```text
deg_T S=1>0=deg_T S(0,T).                                      (27)
```

If `H!=0`, put `d=deg_T H`.  In `(25)`, the last term has unique degree
`2d+2`; the first two terms have degrees at most `1` and `d+1`.  Therefore

```text
deg_T S=2d+2>d>=deg_T S(0,T).                                 (28)
```

### 6.2. Rows N>=1 away from resonance

Assume `N>=1`.  If `H=0`, `(20)` and `(23)` give the strict drop
`2N>N`.  Now let `H!=0` and `d=deg_T H`.  By `(20)` and `(26)`,

```text
m:=deg_T S(0,T) <= max(N,d).                                   (29)
```

If `d<N-1`, the `R_N` term has unique degree `2N`; the other degrees are at
most `N+1+d<=2N-1` and `2d+2<=2N-2`.  Hence

```text
deg_T S=2N>max(N,d)>=m.                                        (30)
```

If `d>N-1`, then `d>=N`; the last term has unique degree `2d+2`, while the
other degrees are `2N` and at most `N+1+d<2d+2`.  Thus

```text
deg_T S=2d+2>max(N,d)>=m.                                      (31)
```

### 6.3. Resonance d=N-1

Let `tau(A)=lc_T(H)`.  All three terms in `(25)` may have degree `2N`.
Using `(22)--(23)` and `lc_T(F)=-1`, their combined leading coefficient is

```text
27A^N[gamma^2+2A gamma tau+A^2tau^2]
=27A^N(gamma+A tau)^2.                                        (32)
```

The bracket has nonzero constant term `gamma`, so it cannot vanish in
`k[A]`.  Consequently

```text
deg_T S=2N>N>=m.                                                (33)
```

Equations `(27)--(33)` prove in every case that

```text
D:=deg_T S > deg_T S(0,T)=m.                                  (34)
```

## 7. Multiplicity-safe projective companion

Homogenize `S` to degree `D` in `[T:Z]`.  Its special fibre is

```text
mathcal S(0;T,Z)=S_0^h(T,Z)Z^(D-m),                            (35)
```

where `S_0^h(1,0)` is the nonzero leading coefficient of `S(0,T)`.  Thus
the total residual reaches

```text
P_0=(A=0,[T:Z]=[1:0]).                                         (36)
```

Factor `S=c product_i G_i^(e_i)` with distinct irreducible `G_i` and all
multiplicities retained.  Actual-degree homogenization selects a reduced
factor `Xi=V(G_i)` whose closure passes through `P_0`.  It is not vertical,
by `(26)`.  It is not the marked branch, because the degree-two
homogenization

```text
F^h=AZ^2-T^2                                                     (37)
```

has value `-1` at `P_0`.  Nor is it a boundary component.  Hence `Xi`
dominates the `A`-line.  Explicitly, a selected factor of `T`-degree zero
would lie in `k[A]`; over the algebraically closed field, its closure could
pass through `A=0` only if it were the already excluded factor `A`.

Normalize the closure of `Xi` in `P1_A x P1_T`.  A point over `P_0` is
absent from the affine normalization because `T/Z` has a pole.  Surjectivity
to `P1_A` supplies another absent point over `A=infinity`, where `A` has a
pole.  The two points are distinct because their `A`-images are zero and
infinity.  Factor repetition changes only contact multiplicity and does not
affect this reduced-component argument.

## 8. Exact scope and next frontier

THM-3866 and THM-3870 close every polynomial graph in the fixed carrier.
This theorem closes the first coordinate-rectifiable but fixed-coordinate
non-graph branch and every one of its polynomial transverse deformations.
The coordinate `T=C/6-A^2` is one elementary triangular generator in a
Jung--van der Kulk rectification.  The argument does **not** yet transport
formally under an arbitrary source shear: such a shear changes both the
distinguished `A`-adic base and the fibre degree used in `(34)`.  A general
orbit theorem therefore needs a weighted degree/contact sidecar, rather
than bare invariance of `Delta_b`.

It does not close all polynomial `A1` branches: higher triangular depth in
an Abhyankar--Moh reduction, and marked-root pairs with more complicated
degree allocation, remain open.  No planar Jacobian counterexample and no
all-`A1` theorem is claimed.

## 9. Exact replay

Run

```bash
python3 04-computation/jc2_first_nongraph_parabola_thm3873.py
python3 -O 04-computation/jc2_first_nongraph_parabola_thm3873.py
```

Both commands must byte-match
`05-knowledge/results/jc2_first_nongraph_parabola_thm3873.out`.  The
assertion-free companion performs 244 exact gates.  It checks `(4)--(9)`,
constructs `q_0,...,q_6`, and replays the exceptional row `N=0` plus all
three degree regimes around `d=N-1` through `N=6`.  Equations `(12)--(34)`
are the all-degree proof; the finite grid is an adversarial control, not its
logical substitute.  Raw-LF SHA-256 hashes are recorded in the metadata.
