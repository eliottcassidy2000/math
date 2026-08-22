---
id: THM-3342
title: "Sublinear deadline excess is impossible for fair critical-run extractors"
status: >
  PROVED + FINITE-EXACT HOSTILE-AUDITED.  If a deterministic exactly fair
  unknown-bias coin extractor has a pathwise deadline T(n) on every stream
  whose initial constant run has length n, then T(n)-n is not o(n).  In
  particular no one extractor has T(n)<=n+D for a fixed D.  Equivalently,
  every feasible deadline has limsup (T(n)-n)/n>0.  This is a per-extractor
  nonattainment theorem: it supplies no uniform positive lower bound on that
  limsup and therefore does not prove that the infimum C* exceeds one.
source: death-star-2026-07-30-laneB, audited-and-promoted-codex-2026-08-12
depends_on:
  - THM-2966-spine-normal-form-for-critical-run-fair-extractors
related:
  - THM-2225-dyadic-critical-run-extractors-and-cyclic-checksum-shell-bisection
  - THM-3009-archimedean-floor-for-balanced-block-extractors
  - THM-3340-single-donor-cyclic-rotation-proves-all-pointwise-AMM-floors
script: 04-computation/amm12592_sublinear_excess_impossible_laneB_deathstar.py
output: 05-knowledge/results/amm12592_sublinear_excess_impossible_laneB_deathstar.out
script_sha256: 3000852d5403792bcfe77ece03b99de6c9519393cae7e419b018eb5aa90bfd08
output_sha256: 61e6bb3eb7d150539fe041691d9a974a941e7962f0ce5877672c06dd56eb34e7
hash_basis: working-tree bytes (LF)
external:
  - "F. Carlson, Uber Potenzreihen mit ganzzahligen Koeffizienten, Math. Z. 9 (1921), 1--13."
  - "P. Fatou, Series trigonometriques et series de Taylor, Acta Math. 30 (1906), 335--400."
  - "L. Kronecker, Zwei Satze uber Gleichungen mit ganzzahligen Coefficienten, J. reine angew. Math. 53 (1857), 173--175."
  - "G. Szego, Uber Potenzreihen mit endlich vielen verschiedenen Koeffizienten, Sitzungsber. Preuss. Akad. Wiss. (1922), 88--91."
---

# THM-3342 -- no sublinear deadline excess

Let independent bits satisfy

```text
P(X_i=0)=p,       P(X_i=1)=q=1-p,       0<p<1,
```

and let `n` be the length of the maximal constant initial run.  An extractor
is deterministic and exactly fair if its heads probability is `1/2` for every
`p`.  A deadline `T` is pathwise: every nonconstant stream of critical value
`n` must stop by flip `T(n)`.

## 1. Statement and quantifiers

**Theorem.** No exactly fair extractor admits a deadline

```text
T(n)=n+o(n).                                                       (1)
```

Consequently every feasible deadline satisfies

```text
limsup_(n->infinity) (T(n)-n)/n > 0.                              (2)
```

This includes the bounded-additive case `T(n)<=n+D`.  The quantifier in (2)
is over `n` after fixing one extractor.  The theorem does **not** give a
constant `epsilon>0` valid for all extractors.  Thus it proves nonattainment
of asymptotic slope one, but it does not prove

```text
C*:=inf{C: some D and some extractor have T(n)<=Cn+D} > 1.        (3)
```

The feasible slopes could still accumulate at one through different
extractors.

## 2. Integer window series

THM-2966 gives the exact spine identity.  Put

```text
d_m=T(m)-m-1.
```

Fairness forbids stopping on a constant prefix, so `d_m>=0`.  There are
integer Bernstein vectors

```text
0<=w_(m,k),v_(m,k)<=binom(d_m,k)
```

and decided-tree polynomials `W_m,V_m` such that

```text
sum_(m>=1) p^m q W_m(p)+sum_(m>=1) q^m p V_m(p)=1/2.              (4)
```

Regard the second summand as a power series in `u=q`.  Define

```text
P_m(p)=q W_m(p)=sum_j a_m(j)p^j,
Q_m(u)=(1-u) V_m(1-u)=sum_j b_m(j)u^j,

F(p)=sum_(m>=1)p^m P_m(p)=sum_(N>=0) f_N p^N,
G(u)=sum_(m>=1)u^m Q_m(u)=sum_(N>=0) g_N u^N.                    (5)
```

All `a_m(j),b_m(j),f_N,g_N` are integers.  Expanding the Bernstein basis
gives the elementary bound

```text
|a_m(j)|,|b_m(j)| <= 2*3^(d_m).                                  (6)
```

For instance

```text
f_N=sum_(m<=N) a_m(N-m),
|f_N|<=2(N+1)3^(max_(m<=N)d_m),                                  (7)
```

and similarly for `g_N`.  Under (1), `d_m=o(m)`, so (7) is
subexponential.  Hence `F` and `G` have radius of convergence at least one.
Equation (4) becomes

```text
F(p)+G(1-p)=1/2,                  0<p<1.                          (8)
```

## 3. Pólya--Carlson rationality

The right side of

```text
F(p)=1/2-G(1-p)                                                     (9)
```

is analytic on `|1-p|<1`; it glues to the power series for `F` on the lens
`|p|<1, |1-p|<1`.  Thus `F` analytically continues across an open arc of its
unit circle.  The same holds for `G`.

If `F` has radius greater than one, its integer coefficients tend to zero
exponentially and are eventually zero, so `F` is a polynomial.  Otherwise
its radius is exactly one.  The Pólya--Carlson theorem says that an integer
power series of radius one is rational or has the unit circle as a natural
boundary.  The continuation above excludes the latter.  Therefore `F` is
rational; identically, so is `G`.

Fatou's rational-series lemma writes a rational integer power series as
`A/B` with coprime `A,B in Z[p]` and `B(0)=1`.  Every zero `rho` of `B` has
`|rho|>=1`.  If `deg B=r`, the reciprocal polynomial

```text
B*(z)=z^r B(1/z)
```

is monic integral and all its roots have modulus at most one.  Their product
is the nonzero integer given by the leading coefficient of `B`, so every
root has modulus exactly one.  Kronecker's theorem now makes every pole of
`F` a root of unity.  The same conclusion holds for `G`, with arbitrary pole
multiplicity allowed.

## 4. Complementary charts force the sixth roots

Identity (8) is now an identity of rational functions.  If `z` is a pole of
`F`, then `1-z` is a pole of `G`.  Both are roots of unity, hence

```text
|z|=|1-z|=1.
```

The two unit circles meet only at

```text
zeta_6=(1+i*sqrt(3))/2,       conjugate(zeta_6),                   (10)
```

the roots of

```text
Phi_6(p)=p^2-p+1,             Phi_6(1-p)=Phi_6(p).                (11)
```

Thus every pole of either series belongs to this conjugate pair.

## 5. The integer-polynomial endgame

Separate the polynomial and proper parts:

```text
F(p)=A(p)+R(p)/Phi_6(p)^e,
G(u)=B(u)+S(u)/Phi_6(u)^h,                                    (12)
```

where the rational terms may be zero and their numerators have degree less
than their denominators.  Initially `A,B` are rational polynomials.  We show
they are integral.

The Taylor coefficients of a proper fraction with denominator `Phi_6^e`
have the form

```text
c_N=U(N)zeta_6^(-N)+conjugate(U(N))zeta_6^N,                    (13)
```

where `U` is a polynomial of degree at most `e-1`.  On each residue class
modulo six, (13) is therefore a polynomial in `(N-r)/6`.  For all sufficiently
large `N`, it equals the integer coefficient `f_N`, because the polynomial
part in (12) has ended.  A polynomial which is integer-valued on a tail of
the integers is integer-valued on the whole integer lattice: its finite
differences are integers on the tail, and the Newton expansion

```text
H(t)=sum_j Delta^j H(T) binom(t-T,j)                              (14)
```

propagates this backwards.  Hence every Taylor coefficient of the proper
part is integral.  Subtracting it from `F` shows `A in Z[p]`; likewise
`B in Z[u]`.

Insert (12) into (8) and use (11):

```text
1/2-A(p)-B(1-p)
 =R(p)/Phi_6(p)^e+S(1-p)/Phi_6(p)^h.                             (15)
```

The left side is a polynomial.  After a common denominator, the right side
is proper and tends to zero at infinity.  Both sides must vanish identically.
At `p=0` this gives

```text
A(0)+B(1)=1/2,                                                   (16)
```

but the left side is an integer.  This contradiction proves the theorem.

## 6. Exact hostile audit and boundaries

Run

```bash
python 04-computation/amm12592_sublinear_excess_impossible_laneB_deathstar.py
python -O 04-computation/amm12592_sublinear_excess_impossible_laneB_deathstar.py
```

Both runs pass.  The nine exact checks cover the decided-tree polynomial
box through depth three, Vandermonde refinement, coefficient bounds, the
sixth-root intersection, the period-six `1/Phi_6` coefficients and
quasi-polynomial tails, a von Neumann positive control for (4), hostile
bounded-depth periodic searches, and Newton backward propagation.  The
classical analytic theorems are cited dependencies, not machine-checked.

The mechanism ends sharply at sublinear depth.  If `d_m` grows linearly,
the coefficient bound in (7) is exponential and the radius can drop below
one, so Pólya--Carlson no longer applies on the unit circle.  A quantitative
uniform lower bound on (2), and hence any proof that `C*>1`, requires a new
capacity or transport argument.  QED.
