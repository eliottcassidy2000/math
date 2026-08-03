# The Hesse moment kernel closes support four in the degree-four cyclic eigenspace

**Status:** PROVED algebraic mixed-moment formula + FINITE-EXACT projective
support-`<=4` exclusion, pending independent immutable audit; canonical
statement
[THM-3321](../01-canon/theorems/THM-3321-hesse-moment-kernel-and-cyclic-quartic-support-four-exclusion.md).
Support five remains OPEN.  This note does not prove `FC(3)` and does not
treat non-eigenvectors.

The exact artifact is
[`degree_four_support4_macaulay_scout_20260803.py`](../04-computation/degree_four_support4_macaulay_scout_20260803.py),
with frozen
[`output`](../05-knowledge/results/degree_four_support4_macaulay_scout_20260803.out).

## Inheritance pass

The closest proved mechanism is
[THM-3310](../01-canon/theorems/THM-3310-degree-four-cyclic-eigenspace-on-the-triangle.md).
It puts

```text
z=s1+omega s2+omega^2 s3,       zbar=s1+omega^2 s2+omega s3
```

and identifies the degree-four `omega`-eigenspace with

```text
g=A zbar+B z^2+C z zbar^2+D z^3 zbar+E zbar^4.              (1)
```

It proves the mixed-moment selection rule and excludes every coefficient
support of size at most three.  Its corrected near miss is also the immediate
frontier: chained affine resultants become too large on support four, while an
affine normalization risks losing projective infinity.

The canonical hostile example is **MISTAKE-363** in
[`MISTAKES.md`](../01-canon/MISTAKES.md): reducing a simplex moment of total
degree `D` modulo `p` is invalid unless `p>D+2`.  Here the largest tested
moment is `m=21` and `deg(g)<=4`, so the exact guard is

```text
p>4m+2=86.                                                   (2)
```

The least-used relevant sidecar is
[HYP-9079](../05-knowledge/hypotheses/HYP-9079-cyclic-eigenspace-selection-rule-for-the-factorial-conjecture.md):
cyclic character forces only the moments with `3|m` to survive.  It is a
selection rule, not an exclusion, but it makes the seven generators
`M_3,M_6,...,M_21` the correct initial moment ideal.

The live concept board was:

| object | representation / operation | question or lost datum |
|---|---|---|
| mixed moments `mu(a,b)` | bivariate coefficient kernel | can the table become a closed sequence? |
| invariant triangle quotient | `r=z zbar,u=z^3,v=zbar^3`, `uv=r^3` | does its torus preserve the integral? |
| coefficient vector | projective `P^4` | how many coefficients may lawfully be normalized? |
| moment conditions | homogeneous ideal in coefficient space | can emptiness replace resultant elimination? |
| support-four boundary | five coordinate `P^3`s | are infinity and lower supports included? |
| support-five interior | complement of the coordinate boundary | what is the smallest exact next certificate? |

The inheritance change is therefore methodological but not cosmetic: retain
THM-3310's monomial coordinates and selection rule, replace its affine
resultant chain by one homogeneous graded-rank question, and build the moment
entries from a closed kernel rather than a barycentric table.

## 1. A closed form for every mixed moment

Let `lambda` range over `1,omega,omega^2`.  For `N=a+b`, polarization of the
Dirichlet simplex formula gives

```text
< (Xz+Yzbar)^N >
  = 2 N!/(N+2)! h_N(X lambda+Y lambda^(-1) : lambda^3=1).     (3)
```

The generating function of the complete homogeneous polynomials on the right
is

```text
product_(lambda^3=1) 1/(1-X lambda-Y lambda^(-1)).            (4)
```

The denominator is the Hesse cubic

```text
product_(lambda^3=1)(1-X lambda-Y lambda^(-1))
  =1-X^3-Y^3-3XY.                                             (5)
```

Put

```text
C(a,b)=[X^aY^b] 1/(1-X^3-Y^3-3XY).                           (6)
```

Comparing the coefficient of `X^aY^b` in (3) yields the exact formula

```text
mu(a,b)=<z^a zbar^b>
       =2 a! b! C(a,b)/(a+b+2)!,                              (7)

C(a,b)=sum multinomial(i+j+k;i,j,k) 3^k,                      (8)
```

where the sum is over `i,j,k>=0` satisfying

```text
3i+k=a,          3j+k=b.                                     (9)
```

Equations (7)--(9) are a proof, not an interpolation.  They immediately show
that `mu(a,b)=0` exactly when `a-b` is nonzero modulo three; in the surviving
congruence class the sum is positive.  They also turn the mixed-moment array
into an efficiently computable sequence:

```text
C(a,b)=C(a-3,b)+C(a,b-3)+3C(a-1,b-1)+[a=b=0],                (10)
```

with negative indices zero.  Thus either a one-dimensional finite sum (8) or
a constant-work lattice recurrence (10) replaces the six-index barycentric
expansion.  This is the requested closed-form/sequence bridge inside the
factorial-conjecture lane, not merely an implementation speedup.

The script multiplies the three cyclotomic factors in (5) exactly, compares
(7) against an independent exact six-index Dirichlet expansion for all 66
pairs `a+b<=10`, recovers THM-3310's named moments, and checks (10) at all
7,225 lattice points `0<=a,b<=84`.  These are controls for the proved identity,
not the source of its proof.

### Connection passport

- **Source:** the uniform Dirichlet moment functional on the triangle.
- **Target:** coefficients of the rational Hesse kernel `(1-X^3-Y^3-3XY)^-1`.
- **Map:** polarization (3), followed by coefficient extraction.
- **Preserved predicate:** the exact value and cyclic selection rule of every
  monomial moment `mu(a,b)`.
- **Destroyed information:** the pointwise location in the triangle and any
  positivity argument for a general complex `g`.
- **Needed sidecar:** the factorial denominator in (7), especially guard (2)
  under modular reduction.
- **Cheapest hostile test:** compare with direct barycentric expansion near all
  three residue classes, including the negative-index boundaries of (10).

## 2. The toric reframe is useful, but its torus is not a symmetry

Set

```text
r=z zbar,        u=z^3,        v=zbar^3,        uv=r^3.       (11)
```

Multiplying (1) by `z` gives the lower-dimensional invariant presentation

```text
z g=A r+B u+C r^2+D r u+E r v.                               (12)
```

This places the five coefficients on the affine toric surface (11), and it is
the best current representation for attacking full support.  It does **not**
license the proposed Hesse-torus coefficient normalization.

Indeed, the formal torus `z->t z`, `zbar->t^-1 zbar` gives the five basis
weights

```text
(-1,2,-1,2,-4).                                               (13)
```

The pure terms of `M_3` are all nonzero, with torus weights

```text
(-3,6,-3,6,-12)                                               (14)
```

and coefficients

```text
(1/10, 1/28, 37/1540, 137/10010, 1/91).                      (15)
```

Since (14) contains three different characters attached to nonzero terms,
`M_3` is not torus-covariant.  The torus preserves the equation `uv=r^3` but
destroys the simplex moment functional.  When `t^3=1`, all five coefficients
in (1) acquire the same scalar and the action is projectively trivial.

The exact audit of the named actions is therefore:

- global scaling of `(A,B,C,D,E)` is a lawful projective normalization that
  we use;
- the cyclic permutation of the triangle acts projectively trivially on this
  eigenspace;
- a reflection exchanges the `omega` and `omega^2` eigenspaces rather than
  permuting these five coordinates internally; and
- complex conjugation is anti-linear, not a second algebraic torus.

Consequently the five coordinate hyperplanes must all be checked.  Setting,
say, both `A=1` and `B=1` would be an unjustified quotient loss.

## 3. Homogeneous Macaulay closure of every support-four chart

Write `M_m(A,B,C,D,E)=<g^m>`.  For each coefficient `J` let

```text
I_J=(M_3,M_6,M_9,M_12,M_15,M_18,M_21)|_(J=0)                 (16)
```

in the four remaining coefficient variables.  All generators are homogeneous
of coefficient degree `m`.  At target degree 21, the Macaulay map is

```text
direct_sum_(m=3,6,...,21) R_(21-m)  --multiply by M_m--> R_21. (17)
```

There are 2,926 rows and

```text
dim R_21=binomial(24,3)=2,024                               (18)
```

columns.  The script builds (17) exactly modulo `p=101`.  Guard (2) holds, so
every denominator through `(4*21+2)!` is a unit.  For **each** of
`J=A,B,C,D,E`, the rank is 2,024.  The approach to closure is also stable:

```text
target degree:       18    19    20    21    22
rank deficiency:     79    46     7     0     0              (19)
```

The degree-22 matrix is freshly rebuilt, rather than inferred by multiplying
the degree-21 result.  It is an implementation-level degree-constancy guard.
The rank engine has a known-empty control `(x0,x1,x2,x3)`, which has full cubic
rank, and the hostile nonempty control `(x0,x1,x2)`, which misses exactly
`x3^3`.

The rank claim is backed by displayed maximal minors, not only a Boolean rank
return.  Order rows by `M_3,...,M_21` and then by the script's lexicographic
weak compositions; order all degree-21 columns by the same convention.  In
every deletion chart, the leading columns of the reduced transpose select the
same 2,024 zero-based rows

```text
0-1690, 2146-2310, 2315-2320, 2601-2716, 2821-2866.          (20a)
```

The determinants of those five square submatrices modulo 101 are

```text
deleted coefficient:        A    B    C    D    E
determinant mod 101:        27   67   42   50   82.          (20b)
```

The script reconstructs each selected square matrix and asks FLINT for its
determinant separately from the preceding rank call.

### Why one guarded modular rank proves the characteristic-zero statement

All entries of (17) lie in the localization `Z_(101)`.  Full column rank modulo
101 is witnessed by the explicit `2024 x 2024` minors (20a)--(20b).  Each same
minor is therefore a nonzero rational number, so (17) has full column rank
over `Q`.  Hence

```text
R_21 subset I_J.                                               (20)
```

In particular every pure power of degree 21 belongs to the ideal:

```text
x0^21,x1^21,x2^21,x3^21 in I_J.                              (21)
```

Any common zero over an algebraic closure of `Q` must therefore have all four
coordinates zero.  There is no projective point on the hyperplane `J=0`.
This includes points at infinity and every lower-support boundary; no
saturation, affine chart, or separate boundary resultant is being hidden.

Every vector of coefficient support at most four lies on at least one of the
five hyperplanes.  Thus the finite-exact conclusion is:

```text
No nonzero degree-four cyclic eigenvector of coefficient support <=4
satisfies all simplex moment conditions.                              (22)
```

In fact the seven conditions through `M_21` already exclude it.  This extends
THM-3310's support-`<=3` frontier by one complete support layer.

The projective construction also repairs an attractive but unsound shortcut.
An affine modular certificate `1 in I mod p` need not lift by itself: a
characteristic-zero solution can acquire denominators and escape into the
projective boundary on reduction.  The homogeneous full-degree certificate
(20), with pure powers (21), is exactly the missing infinity sidecar.

## 4. What support five now looks like

The only remaining coefficient chart has

```text
A B C D E !=0.                                                 (23)
```

We use global projective scaling; the named Hesse torus supplies no further
normalization.  No classification of all continuous moment symmetries is
claimed.
The toric representation (11)--(12), the closed kernel (7), and the
support-four boundary certificate are now the three reusable inputs.

A direct degree-21 homogeneous Macaulay map in five variables has

```text
13,972 rows, 12,650 columns,
4,755,220 raw sparse incidences before modular cancellations.          (24)
```

It is structurally incapable of full column rank.  For five forms of degrees
`3,6,9,12,15` in five variables, the degree-21 quotient is at least the
complete-intersection coefficient `1705`.  The `M_18` row block can remove at
most `dim R_3-1=34` dimensions: multiplication by `M_3` supplies one universal
Koszul overlap.  `M_21` removes at most one more.  Hence

```text
dim (R/I)_21 >=1670,             rank Mac_21 <=10980.        (25)
```

This repairs the first proposed support-five experiment: no degree-21
full-rank certificate can exist.  The formal product

```text
product_(d=3,6,...,21)(1-t^d)/(1-t)^5                      (26)
```

has coefficient `39` at degree `28` and first becomes nonpositive at degree
`29`.  Degree `29` is therefore the first full-rank candidate not ruled out by
this count, with a raw `66486 x 40920` map; this is not a generic-Hilbert-series
claim.  A degree-21 sparse run can test only whether the Hesse ideal reaches
the maximal possible rank `10980`.  Promotion needs a sufficient-degree
guarded certificate or a different exact saturation/affine-chart argument.

The fallback lower-dimensional slice is the common-radial-factor locus

```text
AD=BC,
g=(1+lambda r)(A zbar+B z^2)+E zbar^4,                        (25)
```

but it should be attacked only after (24): (25) is a genuine support-five
sublocus, yet proving it empty would not settle the generic chart.

## 5. Scope and stopping certificate

**PROVED:** the rational Hesse kernel (5), closed moment formulas (7)--(10),
and the torus non-covariance witnessed by (14)--(15).

**FINITE-EXACT:** at `p=101`, under the explicit `101>86` denominator guard,
all five degree-21 Macaulay maps have full column rank; the nonzero-minor lift
therefore proves the characteristic-zero projective exclusion (22).

**OPEN:** full coefficient support five; `FC(3)` outside this degree-four
cyclic eigenspace; all non-eigenvector cases; and any claim that (12) makes the
moment functional torus-invariant.

The stopping certificate is constructive: support four is closed by (20), and
support five has the exact entry oracle (7)--(10), the universal obstruction
(25), the first degree-29 candidate not excluded by that count, and a complete
boundary certificate.
