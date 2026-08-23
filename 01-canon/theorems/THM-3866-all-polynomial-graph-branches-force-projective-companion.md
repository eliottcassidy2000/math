---
id: THM-3866
title: "All polynomial graph branches force a projective companion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Every polynomial
  graph component of the
  depressed-cubic branch, with an arbitrary polynomial transverse quotient
  Q(A,C), has a distinct reduced companion whose affine normalization misses
  at least two projective points.  The proof compares Q with the unique
  formal binomial quotient at its first A-adic mismatch.  All three possible
  C-degree regimes have a strict finite-base degree drop; the only apparent
  resonance has a nonzero perfect-square leading coefficient.
source: jc_sparse_direct_search / first-mismatch completion of THM-3859 and THM-3863, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root, 2026-08-23).  The audit rederived the
  marked-root classification and checked completed monic division against
  the coefficient recursion for arbitrary polynomial s(A).  It verified the
  formal-to-polynomial A^N F divisibility, exact nonlinear response, every
  degree regime including N=0, T=0, and the d=N square resonance.  It also
  followed the factorization with multiplicities: the selected reduced
  factor through finite-base C-infinity is neither vertical, the marked
  graph, nor boundary, and its normalized closure has a distinct point over
  A-infinity.  The companion checks generic
  polynomial s(A)-jets through order four, the universal nonlinear profile
  perturbation, and 141 exact first-mismatch gates for N=0,...,4 in the
  d<N, d=N, and d>N regimes, including the resonant leading square.  Normal
  and optimized runs byte-match the frozen 141-gate transcript and both
  recorded hashes.
related:
  - THM-3859-marked-root-polynomial-graph-companion-puncture-obstruction
  - THM-3863-finite-binomial-hensel-peels-force-projective-companion-contact
script: 04-computation/jc2_all_polynomial_graph_first_mismatch_thm3866.py
output: 05-knowledge/results/jc2_all_polynomial_graph_first_mismatch_thm3866.out
script_sha256: c099664f9108c2a35b74213a0bf3c5ed53d5ddfe161bad2d321a90b42caeef15
output_sha256: 1377d14c155c9c2131366b044b0fe3c1d27f90350443f4887b4c0f7b3986f946
semantic_sha256: 759fa0f273f8340f83a648ce992b5018b8a61c5d19559532611001567d36fd7c
hash_basis: raw LF bytes
---

# THM-3866 -- every polynomial graph has a projective companion

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `k` of characteristic zero.  For
`b=b(A,C) in k[A,C]`, put

```text
Delta_b=-27A^2b^2+8AC^3-54ACb+9C^2-54b.                         (1)
```

Suppose `V(Delta_b)` contains a polynomial graph

```text
Gamma: F=C-g(A)=0,                  g in k[A].                   (2)
```

Then the reduced branch packet contains an irreducible component `Xi`
different from `Gamma` such that the affine normalization of `Xi` omits at
least two distinct points of its smooth projective normalization.  One
omitted point lies over

```text
A=0,                         C=infinity,                         (3)
```

and another lies over `A=infinity`.

This closes **all polynomial transverse quotients `Q(A,C)`** for graph
branches.  No degree bound is imposed on `g`, `b`, or `Q`.  It does not close
non-graph polynomial-curve components, and it does not close the vertical
axis `A=0`.

## 1. Coordinate-safe marked-root form

The cusp identity is

```text
A^2 Delta_b=27(P^3-u^2),
P=1+(2/3)AC,                    u=1+AC+A^2b.                     (4)
```

Restrict to the graph `(2)`, whose coordinate ring is the UFD `k[A]`.
Exactly as in THM-3859, `P^3=u^2` has a unique marked root congruent to one
at `A=0`.  It gives a unique `s=s(A) in k[A]` with

```text
g=6s(1+As),
b0=2s^2(3+4As),
F=C-6s(1+As),
b=b0+QF,                         Q in k[A,C].                    (5)
```

For completeness, the UFD step is short.  Write `P=z^2,u=z^3`; then

```text
Ag=(3/2)(z-1)(z+1),
A^2b(A,g)=(1/2)(z-1)^2(2z+1).                                  (6)
```

The special value is `z(0)=1`, so `z=1+2As`, and `(5)` follows.  Conversely
`(5)` gives this square/cube identity on `F=0`.  Since `F` is monic in `C`,
the quotient `Q` in `(5)` is unique.  This makes the later first mismatch
intrinsic to the fixed graph coordinate `A`; no choice of a representative
modulo `F` remains.

## 2. The unique formal comparator for arbitrary polynomial s(A)

In the complete ring `k[C][[A]]`, define

```text
b_*={(1+(2/3)AC)^(3/2)-1-AC}/A^2
   =sum_(j>=0) beta_j A^j C^(j+2),                              (7)

beta_j=binom(3/2,j+2)(2/3)^(j+2) !=0.                          (8)
```

It is the unique formal profile for which `Delta=0` and `u=1 mod A`.
On the polynomial graph in `(5)`,

```text
1+(2/3)AC=(1+2As(A))^2,
b_*|_(F=0)=b0.                                                   (9)
```

Monic division by `F` therefore gives a unique formal quotient

```text
Q_*=(b_*-b0)/F=sum_(j>=0)q_j(C)A^j       in k[C][[A]].          (10)
```

Although `s` may have arbitrary degree, every coefficient obeys the uniform
law

```text
q_j in k[C],                 deg_C q_j=j+1,
lc_C(q_j)=beta_j.                                                (11)
```

To see this without choosing coordinates on the graph again, write

```text
g(A)=sum_(i>=0)g_i A^i,        b0(A)=sum_(i>=0)c_i A^i.          (12)
```

The coefficient of `A^j` in `FQ_*=b_*-b0` is

```text
(C-g_0)q_j
=beta_j C^(j+2)-c_j+sum_(i=1)^j g_i q_(j-i).                   (13)
```

The right side vanishes at `C=g_0` by `(9)`, so division is polynomial.  Its
unique top term is `beta_j C^(j+2)`, proving `(11)` inductively.  In
particular `Q_*` has infinitely many nonzero `A`-coefficients and is not a
polynomial.

## 3. First mismatch and exact nonlinear response

Expand the polynomial quotient from `(5)` in powers of `A`.  Since `Q_*` is
not polynomial, there is a unique least `N>=0` at which it differs from
`(10)`.  Thus

```text
Q=Q_<(N)+A^N T,
Q_<(N)=sum_(j=0)^(N-1)q_j(C)A^j,
t(C)=T(0,C) != q_N(C).                                         (14)
```

The empty sum is used for `N=0`; `T` itself may be zero when `Q` is exactly
the finite truncation `Q_<(N)`.

Put

```text
b_N=b0+FQ_<(N),
u_N=1+AC+A^2b_N.                                                 (15)
```

The all-depth calculation of THM-3863, repeated here, gives

```text
Delta_(b_N)=A^N F R_N,                                         (16)

R_N(0,C)=54q_N(C).                                              (17)
```

For `N>=1`,

```text
deg_C R_N=2N+1,
lc_C(R_N)=-27 beta_(N-1)^2 A^N,                                (18)

deg_C u_N=N+1,
lc_C(u_N)=beta_(N-1)A^(N+1).                                   (19)
```

For `N=0`, the base companion has degree two, leading coefficient `8A`,
and its special fibre has degree one.

Here is the promised derivation, including the arbitrary-`s(A)` scope.  In
`k[C][[A]]` write

```text
H_N=sum_(j>=N)q_j(C)A^(j-N),
b_*=b_N+A^N F H_N.                                             (19a)
```

Since `Delta_(b_*)=0`, the exact quadratic response formula below shows
that `Delta_(b_N)` is divisible by `A^N F` in the complete ring.  Independently,
`F` divides `Delta_(b_N)` already in `k[A,C]`, because `b_N|_(F=0)=b0`.
Write `Delta_(b_N)=F L` there.  Monicity of `F` makes this polynomial
quotient agree with the completed quotient, so every coefficient of
`L in k[A,C]` is divisible by `A^N`; hence `R_N=L/A^N` is polynomial.
The completed response identity now gives `R_N(0,C)=54q_N(C)`, proving
`(17)`.  For `N>=1`, the unique top `C`-term of `b_N` is

```text
beta_(N-1) A^(N-1) C^(N+1).                                   (19b)
```

This gives `(19)`.  In `(1)`, the square `-27A^2b_N^2` then has the unique
top degree `2N+2`, so division by `A^N F` gives `(18)`.  For `N=0`, direct
use of `b_N=b0(A)` makes `8AC^3` the top term before division by `F`, and
specialization at `A=0` gives the asserted degree-one residual.  Thus
`(16)--(19)` require no bound or constancy assumption on `s`.

The exact response of `(1)` to any perturbation `delta b` is

```text
Delta_(b_N+delta b)-Delta_(b_N)
=-54u_N delta b-27A^2(delta b)^2.                               (20)
```

Use `delta b=A^N F T`.  Equations `(16)` and `(20)` give

```text
Delta_b=A^N F S,

S=R_N-54u_NT-27A^(N+2)FT^2,                                   (21)

S(0,C)=54(q_N-t) !=0.                                          (22)
```

Thus neither `A` nor the entire special fibre is a factor of `S`.

## 4. The three degree regimes and the resonant square

First suppose `N>=1` and `T!=0`; put

```text
d=deg_C T,                         tau(A)=lc_C(T).               (23)
```

The special degree satisfies

```text
m:=deg_C S(0,C) <= max(N+1,d).                                  (24)
```

There are three exhaustive generic regimes.

### 4.1. Lower mismatch degree d<N

The `R_N` term in `(21)` has degree `2N+1`, while

```text
deg_C(u_NT)<=N+1+d<=2N,
deg_C(FT^2)=2d+1<=2N-1.                                        (25)
```

Therefore

```text
deg_C S=2N+1>max(N+1,d)>=m.                                    (26)
```

### 4.2. Higher mismatch degree d>N

The last term of `(21)` has unique degree `2d+1`; the other two degrees are
`2N+1` and at most `N+1+d`.  Since `d>N`,

```text
deg_C S=2d+1>d=max(N+1,d)>=m.                                  (27)
```

### 4.3. Resonance d=N

All three terms may reach degree `2N+1`.  Their combined leading
coefficient is not an uncontrolled cancellation but the exact square

```text
-27A^N[beta_(N-1)^2
       +2A beta_(N-1)tau+A^2tau^2]
=-27A^N(beta_(N-1)+A tau)^2.                                   (28)
```

Its bracket has nonzero constant term `beta_(N-1)`, so `(28)` is not the
zero polynomial.  Hence

```text
deg_C S=2N+1>N+1>=m.                                            (29)
```

If `T=0`, then `S=R_N` and `(17)--(18)` give the same strict drop directly.

Finally let `N=0`.  If `T=0` or `deg_C T=0`, the base quadratic retains
degree two while `(22)` has degree one.  If `d>=1`, the last term in `(21)`
has unique degree `2d+1`, while `(22)` has degree at most `d`.  Thus the
strict inequality

```text
D:=deg_C S > deg_C S(0,C)=m                                    (30)
```

holds in every first-mismatch case, including all boundaries omitted by the
notation in `(23)`.

## 5. Degree drop selects a distinct bad component

Homogenize `S` in the `C` coordinate to its generic degree `D`, using
`[C:Z]` on `P1_C`.  If `S_0^h` is the degree-`m` homogenization of the
nonzero polynomial `S(0,C)`, then `(30)` gives

```text
mathcal S(0;C,Z)=S_0^h(C,Z) Z^(D-m).                            (31)
```

The value `S_0^h(1,0)` is its nonzero leading coefficient, while the extra
factor `Z^(D-m)` vanishes there.  Hence the total projective branch has
positive contact at

```text
P_0=(A=0,[1:0]).                                                 (32)
```

Factor with all multiplicities retained,

```text
S=c product_i G_i^(e_i),                                       (33)
```

where the `G_i` are distinct irreducibles.  Individual homogenization at
the actual `C`-degrees turns `(33)` into the projective product, so at least
one reduced component `Xi=V(G_i)` has closure through `P_0`.

This component is not vertical: `(22)` says `A` does not divide `S`, and an
irreducible curve supported over `A=0` would have defining factor `A`.
It is not the marked graph `F=0`, either, because

```text
F(0,C)=C-6s(0)                                                   (34)
```

has nonzero value at `[1:0]` after degree-one homogenization.  It is also
not a component contained in the boundary `Z=0`, since `G_i` is an affine
factor homogenized at its actual degree.  Therefore `Xi` is a distinct
companion and dominates the `A`-line.

Close `Xi` in `P1_A x P1_C` and normalize:

```text
nu: Xi_tilde -> Xi_bar.                                        (35)
```

At least one point `p_0` of `Xi_tilde` lies over `(32)`.  The pullback of
`C/Z` has a pole there, so `p_0` is absent from the affine normalization.
Since `Xi` dominates `A`, the projective morphism
`Xi_tilde -> P1_A` is nonconstant and hence surjective.  Choose a point
`p_infinity` over `A=infinity`; the function `A` has a pole there, so this
point is also absent.  Their `A`-images are respectively zero and infinity,
so

```text
p_0 != p_infinity.                                               (36)
```

This proves the theorem for reducible and nonreduced total packets as well
as irreducible ones.  Multiplicity in `(33)` changes only the contact count,
not the selected reduced component or its two deleted points.

## 6. Scope and structural meaning

THM-3859 closed the one-variable quotient `Q=q(A)` by a quadratic
discriminant.  THM-3863 showed that every canonical finite truncation of the
binomial quotient retains projective contact.  The first-mismatch argument
above closes the entire polynomial space between and beyond those two
families: any deviation either lies below the truncation degree, above it,
or on the perfect-square resonance `(28)`.

The common invariant with THM-3860 is a finite-divisor Newton edge.  There,
a vanishing multiplier in `M(S-f)=s` forces a rational normal coordinate to
have negative valuation.  Here, the strict degree drop `(30)` forces a
polynomial companion sheet to reach `C=infinity` over finite `A=0`.  In both
cases, local formal solvability does not make the normal coordinate globally
effective.

What remains is sharply different:

```text
- polynomial A1 components not expressible as graphs over A;
- the vertical branch A=0, equivalently b(0,C)=C^2/6;
- constructions outside this depressed-cubic carrier.           (37)
```

No planar Jacobian counterexample and no all-carrier no-go is claimed.

## 7. Exact replay

Run

```bash
python3 04-computation/jc2_all_polynomial_graph_first_mismatch_thm3866.py
python3 -O 04-computation/jc2_all_polynomial_graph_first_mismatch_thm3866.py
```

Both commands must byte-match
`05-knowledge/results/jc2_all_polynomial_graph_first_mismatch_thm3866.out`.
The assertion-free companion performs 141 gates.  Its finite universe is
generic `s(A)`-jets through order four and, in a hostile fixed graph, every
regime `d=N-1,N,N+1,N+2` for `N=0,...,4`.  The all-degree proof is
`(13)` and `(20)--(30)`; the finite grid is a replay and adversarial control,
not its logical substitute.  Raw-LF SHA-256 hashes are recorded in the
metadata.
