---
id: THM-3587
title: "Consecutive Keller fibre factor, toric-carrier, and coefficient-span gates"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For a complex
  planar Keller map, put T=P_y Q_x, so T+kappa=P_x Q_y.  A hypothetical
  nonautomorphism has multiplicity factor count at least four across these
  two consecutive fibres, distinct-factor count at least four, and a
  distinct-factor equality cell only of balanced type (2,2).  A scalar
  matching carrier is tame when all selected carrier levels are irreducible,
  and every positive toric-monomial matching carrier is tame in all degrees.
  Coefficient-matrix span at most two is tame.  Span three either has a
  rank-one constant observer or gauges to the binary constant-Hessian
  equation; in that trace gauge, a hypothetical nonautomorphism has the
  stronger distinct-factor floor five.  Balanced span-four maps and the
  general constant-nonzero binary Hessian problem remain OPEN.  No proof or
  counterexample to JC(2) is claimed.
source: boxeph consecutive-factorization/integrability session, 2026-08-20--21
audit: >
  PASS.  An independent hostile audit checked the UFD root-count argument,
  arbitrary factor multiplicities in the unbalanced equality cell, the
  toric exponent transport, the scalar-carrier degree lemma, Frobenius
  orientation, and all three branches of the rank-three factor-five proof.
  Ordinary, optimized, and stored exact-companion transcripts are identical.
depends_on:
  - THM-2063-one-fiber-linear-planar-keller-pairs
related:
  - THM-3367-berggren-spinor-pencil-hessian-gauge-and-affine-line-keller-closure
  - THM-3548-planar-keller-conductance-shadow-gates
  - THM-3551-one-ray-planar-jacobian-mate-no-go
  - HYP-8905-binary-symmetric-jc2-subcase-and-separate-descent-programs
script: 04-computation/keller_consecutive_factor_integrability_scout.py
output: 05-knowledge/results/keller_consecutive_factor_integrability_scout.out
reflection: 07-reflections/keller-consecutive-factorization-integrability-plaquette-boxeph-20260820.md
script_sha256: 4d155585f3537cf98b0e496ee0936c8196e3837ffaca9e49573a3308105646c9
output_sha256: 56565498e51f0df6dd16c7b6fc194f358e5087406b197f650ce2ebb10688adf6
hash_basis: raw LF bytes
---

# THM-3587 -- consecutive Keller fibre factor, toric-carrier, and coefficient-span gates

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over
`C`.  Let

```text
F=(P,Q),                 Jac(P,Q)=kappa in C*,
a=P_x,  b=P_y,           c=Q_x,  d=Q_y,
T=bc.                                                        (1)
```

Then

```text
ad=T+kappa,             a_y=b_x,          c_y=d_x.          (2)
```

The theorem gives necessary gates for a hypothetical nonautomorphism.  It
does not close the balanced full-span cell and does not prove `JC(2)`.

## 1. Statements

For a nonzero polynomial `f`, let `Omega(f)` count irreducible factors with
multiplicity and let `omega(f)` count distinct nonassociate irreducible
factors.

### A. Consecutive-factor floor and balanced equality

If `F` is not a polynomial automorphism, then

```text
Omega(T)>=2,                 Omega(T+kappa)>=2,
Omega(T)+Omega(T+kappa)>=4,                              (3)

omega(T)+omega(T+kappa)>=4.                              (4)
```

Moreover, equality in `(4)` forces

```text
omega(T)=omega(T+kappa)=2.                               (5)
```

The multiplicity floor `(3)` and distinct-factor floor `(4)` are different.
There are tame Keller maps with total `Omega=4` and total `omega=3`.

### B. Scalar and toric matching carriers

Suppose

```text
T=h(R),                 h in C[t],  R in C[x,y].         (6)
```

If every polynomial `R-r` is irreducible for every distinct root `r` of

```text
h(t)(h(t)+kappa),                                          (7)
```

then `F` is tame.  The corresponding statement for `T+kappa` follows by
swapping the two source coordinates, which exchanges the matching products
and changes the sign of `kappa`.

More strongly, let `p,q` be coprime positive integers.  If either matching
product is a polynomial in the positive toric monomial

```text
z=x^p y^q,                                                (8)
```

then `F` is tame, with no degree bound and no selected-level irreducibility
assumption.

### C. Coefficient-span observer and Hessian router

Let `S_F` be the span in `M_2(C)` of the nonconstant coefficient matrices of
`DF`.  Equivalently it is the span of `DF(z)-DF(z_0)` as `z` varies.  Then

```text
dim S_F<=2             implies             F is tame.    (9)
```

If `dim S_F=3`, let `W` generate the one-dimensional Frobenius annihilator
for

```text
<W,M>=tr(W^T M).                                          (10)
```

A singular `W` again makes `F` tame.  If `W` is invertible, the source change
`z -> W^T z` gives a constant-trace gauge

```text
F=(lambda x+H_y,       lambda y-H_x),
det DF=lambda^2+det Hess(H).                              (11)
```

This is a router to the binary constant-Hessian equation, not an
invertibility equivalence with `grad H`.

In the gauge `(11)`, redefine `T=P_yQ_x` and let the transformed nonzero
Jacobian again be denoted `kappa`.  A hypothetical nonautomorphism satisfies

```text
omega(T)+omega(T+kappa)>=5.                              (12)
```

The gauge qualification is load-bearing: the individual product `P_yQ_x`
is not invariant under arbitrary source changes.

## 2. Unit edges and the irreducible-carrier lemma

Equation `(2)` gives the four unit cross ideals

```text
(a,b)=(a,c)=(d,b)=(d,c)=C[x,y].                         (13)
```

If any edge is constant, `F` is tame.  For example, if `a=alpha!=0`, then
`P=alpha x+p(y)`.  In coordinates `(U,V)=(P,y)`, the determinant equation is
`Q_V=kappa/alpha`, so

```text
Q=(kappa/alpha)V+H(U).                                   (14)
```

If `a=0`, then `b,c` are units and the other triangular form follows.  This
is THM-2063 in its constant-derivative boundary.  Hence all four edges of a
hypothetical nonautomorphism are nonzero nonunits.

We use the following elementary lemma repeatedly.

**Irreducible-carrier lemma.**  Let `r in C[x,y]` be irreducible and let
`A,B in C[t]` be nonconstant.  If

```text
A(r)_y=B(r)_x,                                            (15)
```

then `r` is affine linear, and a nonzero linear combination of `A(r)` and
`B(r)` is constant.

Write `A'=gA_1,B'=gB_1` with `gcd(A_1,B_1)=1`.  Equation `(15)` becomes

```text
A_1(r)r_y=B_1(r)r_x.                                     (16)
```

Bezout remains Bezout after substituting `r`.  Thus `A_1(r)|r_x` and
`B_1(r)|r_y`.  Neither partial vanishes under the hypotheses.  If `A_1` or
`B_1` were nonconstant, its substitution would have degree at least `deg r`,
strictly larger than the corresponding partial.  Hence `A_1,B_1` are
constants, `r_x,r_y` are proportional, and `r` lies in a one-linear-form
ring.  Irreducibility over `C` makes it affine linear.  Also `A',B'` are
proportional, giving the claimed constant combination.  THM-2063 then closes
any Keller map to which the lemma applies.

## 3. Proof of the consecutive-factor theorem

The two factors on each side of `(2)` are nonunits, so `(3)` is immediate.
For `(4)--(5)`, suppose one matching side has only one distinct prime.  After
a source swap if necessary, write

```text
T+kappa=u r^n,                                            (17)
```

with `u in C*`, `r` irreducible, and `n>=2`.  If
`rho_1,...,rho_n` are the distinct nonzero roots of `u z^n=kappa`, then

```text
T=u product_(j=1)^n (r-rho_j).                           (18)
```

The levels are pairwise coprime, so `omega(T)>=n`.

If the total distinct count is at most three, then `n=2` and each level in
`(18)` has one distinct prime.  It is genuinely irreducible: if
`r-rho_j=v f^m` with `m>1`, then `r=vf^m+rho_j` splits over `C` into
nonconstant shifts of `f`, contradicting irreducibility of `r`.  Every edge
now lies in `C[r]`; the first curl and the irreducible-carrier lemma make the
map tame.  This proves `(4)`.

For equality four with distribution `(1,3)`, equation `(18)` gives `n<=3`.
If `n=3`, the same proper-power argument makes all three levels irreducible,
and the carrier lemma closes the map.  It remains to treat `n=2`, where one
opposite level has two distinct primes with arbitrary multiplicities and the
other is irreducible.

After a constant target gauge and rescaling, the relevant entries satisfy

```text
BC=R(R+2),             B_x=R_y,            C_y=R_x,     (19)
```

and the other two entries are scalar multiples of `R+1`.  Choose the sign so
that `R+2` is irreducible.  If `R+2` belongs to `C`, write

```text
B=S,                   C=E(R+2),            R=SE.       (20)
```

The two curls reduce exactly to

```text
(R+2)E_y=S E_x.                                           (21)
```

Because `gcd(R+2,S)=1`, the left factor divides `E_x`.  If `S` is
nonconstant, then

```text
deg(R+2)=deg S+deg E>deg E_x,
```

so `E` is constant.  If `S` is constant, `B` is a constant edge.  If `R+2`
belongs to `B`, write

```text
B=S(R+2),              C=E,                 R=SE.       (22)
```

The symmetric identity

```text
(R+2)S_x=E S_y                                             (23)
```

gives the same alternatives.  A nonconstant surviving allocation is
whole-level, so `R_x,R_y` are proportional and an edge `R` differs by a
constant from an edge `R+1`.  THM-2063 closes it.  This proof retains every
factor multiplicity and establishes `(5)`.

## 4. Proof of the scalar and toric gates

If `h` in `(6)` is constant, then `T` is constant.  A nonzero `T` makes
`b,c` units, while `T=0` makes `b=0` or `c=0`; either way there is a constant
edge.  Now assume `h` is nonconstant.  Under `(6)--(7)`, UFD factorization
puts all four edges in `C[R]`.  Choose a selected root `r_0`; the irreducible
polynomial `R-r_0` generates the same univariate ring.  The first curl and
the irreducible-carrier lemma give a constant directional derivative,
proving the scalar statement.

For the toric theorem, it suffices to suppose `T=h(z)` for `(8)`.  If `h` is
constant and nonzero, `b,c` are units; if `h=0`, then `b=0` or `c=0`.
Thus a constant edge appears in either case.  If `h(0)` is neither `0` nor
`-kappa`, all selected roots are nonzero and every binomial `z-rho` is
irreducible by primitivity of `(p,q)`.  Thus all entries lie in `C[z]`.  The
first curl is

```text
q x A'(z)=p y B'(z).                                     (24)
```

The exponent rays on the two sides are disjoint, so `A'=B'=0` and there is a
constant edge.

Suppose next that `h(0)=0`.  Extract every coordinate factor from the zero
level:

```text
a=A(z),                 d=D(z),
b=x^r y^s B(z),         c=x^u y^v C(z),
B(0)C(0)!=0.                                               (25)
```

The curls force positive integers `delta,epsilon` with

```text
r=p delta+1,           s=q delta-1,
u=p epsilon-1,         v=q epsilon+1,                    (26)

qA'=z^(delta-1)[(p delta+1)B+pzB'],
pD'=z^(epsilon-1)[(q epsilon+1)C+qzC'].                  (27)
```

Put `K=delta+deg B` and `L=epsilon+deg C`.  The leading ratios are

```text
lc(A)/lc(B)=(pK+1)/(qK),
lc(D)/lc(C)=(qL+1)/(pL).                                 (28)
```

Top cancellation in `ad-bc` would make their product one, but

```text
(pK+1)(qL+1)-pqKL=pK+qL+1>0.                            (29)
```

This is impossible.  If `h(0)=-kappa`, a source swap turns the old
`T+kappa` into the new zero-level matching product and changes `kappa` to
`-kappa`.  A source swap also proves the alternative where `T+kappa`, rather
than `T`, is initially toric.

## 5. Proof of the coefficient-span gates

Use the Frobenius pairing `(10)`.  If a nonzero rank-one

```text
W=u v^T                                                    (30)
```

lies in `S_F^perp`, then

```text
<W,DF>=u^T DF v                                           (31)
```

is constant.  Invertible source and target completions make `(31)` a
constant derivative entry, and THM-2063 applies.

Every complex matrix subspace of dimension at least two contains a nonzero
singular matrix: restrict the quadratic determinant to any projective line.
If `dim S_F<=2`, its annihilator has dimension at least two and therefore
contains such a rank-one matrix.  This proves `(9)`.

If `dim S_F=3`, its annihilator is generated by one matrix `W`.  A
nonautomorphism forces `W` invertible.  For

```text
G(z)=F(W^T z),                                             (32)
```

cyclicity of trace gives

```text
tr DG=tr(DF W^T)=tr(W^T DF)=constant.                    (33)
```

The transformed determinant is `kappa_tr=kappa_old det W in C*`; from this
point relabel `kappa_tr` as `kappa`.  Thus every `kappa` in the trace gauge,
including `(34)--(35)`, denotes the transformed determinant.

Subtracting `lambda id`, with twice `lambda` equal to this trace, leaves a
divergence-free polynomial vector field.  Polynomial exactness on `C^2`
gives `G-lambda id=(H_y,-H_x)`, proving `(11)`.

## 6. The rank-three distinct-factor floor five

In `(11)`, put

```text
R=H_xy,                  mu=det Hess(H)=kappa-lambda^2.
```

Then

```text
a=lambda+R,        b=H_yy,        c=-H_xx,       d=lambda-R,
H_xx H_yy=R^2+mu,                                          (34)

T=-(R^2+mu),             T+kappa=lambda^2-R^2.            (35)
```

The global theorem gives total distinct-factor count at least four and
balances equality four.  We exclude that equality in three cases.

If `lambda mu!=0`, the four selected `R`-values in `(35)` are distinct;
otherwise `kappa=0`.  Equality four gives one distinct prime per selected
level.  If one level were a proper power `u f^m`, any other level would be
`u f^m+constant` and would split into `m` distinct shifts of `f`.  Hence all
selected levels are irreducible.  The scalar-carrier theorem closes the map.

If `lambda=0`, then `mu=kappa!=0`.  Equality four says that `R` has two
distinct primes and each of `R-alpha,R+alpha`, where `alpha^2=-mu`, has one.
The preceding proper-power comparison makes both opposite levels
irreducible.  Since `b,c` are nonunits, each takes one whole level.  The curl
`b_x=a_y` makes `R_x,R_y` proportional, producing a constant derivative and
a tame map.

It remains to treat `mu=0`, where `lambda!=0`.  Equality four makes `R` have
exactly two distinct primes.  Put

```text
B=b=H_yy,                 C=-c=H_xx.                     (36)
```

Then

```text
BC=R^2,                  B_x=R_y,          C_y=R_x.      (37)
```

Write

```text
R=rho f^m g^n,
B=beta R L,              C=beta^(-1)R/L,
L=f^r g^s,               -m<=r<=m,  -n<=s<=n.           (38)
```

The two curls express `R_y` in two ways and give the rational Burgers
equation

```text
L_y=beta L L_x.                                             (39)
```

If `L` is constant, `R_x,R_y` are proportional and THM-2063 applies.  If a
prime `h` among `f,g` is mixed and has nonzero `L`-exponent `e`, then

```text
ord_h(L_y)=e-1,             ord_h(L L_x)=2e-1.           (40)
```

Equation `(39)` would force `e=0`, a contradiction.  Thus every mixed prime
has exponent zero in `L`.  One mixed and one pure prime leaves a nonzero
power of a one-variable prime in `L`, impossible in `(39)`.  Two pure primes
in the same variable make `(39)` force `L` constant.

The only remaining case has, after relabelling, `f=f(x)` and `g=g(y)`.
Irreducibility over `C` makes both affine linear.  Exponent comparison in
`(39)` forces `(r,s)=(1,-1)`.  The original curls then give

```text
beta(m+1)f_x=n g_y,          beta m f_x=(n+1)g_y,        (41)
```

which implies `mn=(m+1)(n+1)`, impossible for positive `m,n`.  This closes
the last equality-four branch and proves `(12)` without invoking a general
nonhomogeneous Hessian-zero theorem.

## 7. Sharp controls, exact companion, and open boundary

The tame map

```text
v=x+y,                  Q=y+v^2/2,             P=v+Q
```

has

```text
DF=[v+1,v+2; v,v+1],    T=v(v+2),              T+1=(v+1)^2.  (42)
```

It has total `Omega=4`, total `omega=3`, and the explicit inverse
`v=P-Q,y=Q-v^2/2,x=v-y`.

The determinant-one matrix

```text
[xy+1,xy+2; xy,xy+1]                                      (43)
```

has the same consecutive-factor pattern but both curl defects equal `x-y`.
It is the load-bearing hostile showing that factorization alone is not
integrability.

Two exact tame shear words occupy the balanced `(2,2)` cell.  The first has
coefficient span three and the Frobenius rank-one observer
`W=[-1,0;1,0]`, reading `c-a=-1`.  The second has coefficient span four,
with a displayed coefficient minor `-1/32`.  In both cases
`gcd(T_x,T_y)=1`, but both `T` have critical points.  Thus
gradient-coprime is neither the unit-gradient-ideal condition nor a
counterexample certificate.

The exact companion verifies 53 identities and controls, including ordinary
and optimized byte-identical replay, an independent SymPy/python-FLINT
gradient-gcd computation, the multiplicity-safe allocation identities, the
toric gap `(29)`, the rank-three matching formulas, and explicit inverses.
The companion and detailed reflection are linked in the metadata.

What remains OPEN is substantial.  A smallest derivative-table survivor may
still have balanced factor pattern `(2,2)` and coefficient span four; the
three-shear control proves that this passport is populated but not that it
contains a nonautomorphism.  In span three, `(12)` is only a stronger
factor-count gate inside the constant-trace gauge.  THM-3367 closes its
affine-line coefficient image and top homogeneous ray, while the general
constant-nonzero binary Hessian problem in HYP-8905 remains open.  No
nonproperness, collision, inverse, or `JC(2)` conclusion follows from the
factor counts alone.
