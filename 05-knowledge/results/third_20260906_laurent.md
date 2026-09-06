# Forced characteristic factors from the lower carry, and endpoint-39 detection

**Status: PROVED all-height divisibility and endpoint-39 detection;
[independent analytic and resultant audit accepted](third_20260906_laurent_audit.md).** No theorem ID is reserved. The structural result
below applies to every channel height. Its first consumer independently
certifies the actual family:

```text
support (-39,2g-39,3g-39), g>=20, gcd(g,39)=1:
first nonzero constant-term moment is exactly g or 2g.
```

All three coefficients may be arbitrary nonzero complex numbers. Both
alternatives occur for every allowed g. The general all-height positivity
of the remaining characteristic factors, and general trinomial two-rung
separation, remain **OPEN**.

## 1. Inheritance, novelty and the concept board

The closest proved mechanism is the
[endpoint-33 exact carried certificate](second_20260906_laurent.md),
with the [all-height quotient degree interface](overnight7_20260906_laurent_midpoint_transport.md),
Section 5. The real-root input for the first row is
[THM-4436, complete factorial-row roots](../../01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md).
The new factor theorem itself is algebraic and does not use that real-root
input.

The canonical hostile is the
[actual all-integer phase cone separator](continuing1_20260906_laurent_cone_separator.md).
The corrected near misses are the fixed-row amplitude recipe, repaired by
[mixing four original windows](continuing2_20260906_laurent_sparse_amplitude.md),
and scalar anchoring, whose
[effective two-anchor tail](continuing2_20260906_effective_anchor.md)
still leaves a continuous finite-phase question. Those are scoped inputs,
not uniform noncancellation theorems. The least-used sidecar here is the
**negative-integer parameter degeneration together with its inverse carry**.

The board contains six objects: complete return fibers; characteristic
coefficients; coefficient-field degenerations; the lower carry; the
symmetric Bezout matrix; and the minimal phase/derivative observer cone.
The first two connect to the fourth through the third. The cheap Bezout
probe below rejects a proposed coefficientwise matrix-positivity explanation;
the degeneration produces the positive structural result instead. The
amplitude repairs retain their original branch and field restrictions and
are not dependencies of the present proof.

At the inherited base, targeted searches for endpoint39, the literal
support, forced Laurent
characteristic factors, negative-integer degenerations and small-root carry
orders found no existing statement of either result below. Endpoint33 is
inherited, not counted again. No external-priority claim is made.

After this proof and its certificate were frozen, incoming commits
`28f41c8461` and `01c8b6887b` supplied an independently proved and audited
[endpoint-39 certificate](long_frontier_sep06_endpoint39.md) for exactly
the same actual family and normalized quotient. Its 258 positive shifted
characteristic coefficients establish the same noncancellation conclusion.
The connection is exact: the six characteristic polynomials there are the
products of our forced factors and the six residual polynomials here.
Our all-height degeneration theorem reduces this endpoint certificate
to 208 positive residual coefficients, and reduces the largest required
interpolation bank from 73 nodes to 58. The endpoint closure is therefore
independently reproduced; the additional result here is the uniform
factor structure and its certificate compression.

## 2. The all-height normalized rows

Fix an integer h>=1. Use x as a formal parameter and define

```text
p_(h,x)(t)=sum_(j=0)^h
  [(2h+1)! (x+h)_(h-j)/((3h-3j)!(1+2j)!)] t^j,

q_(h,x)(t)=sum_(e=-1)^(2h)
  [(2x+2h)_(2h-e)/((6h-3e)!(2+2e)!)] t^e.                (1)
```

The first polynomial is monic. Its constant coefficient is
`p_0=(2h+1)! (x+h)_h/(3h)!`. Initially invert t over
`Q(x)[t]/(p_(h,x))`. The inverse-carry denominator cancels exactly:

```text
q_(-1)/p_0
 = [2^(h+1)(3h)!/((2h+1)!(6h+3)!)]
     *x*product_(a=0)^(h-1)(2x+2a+1).                       (2)
```

Splitting `(2x+2h)_(2h+1)` into even and odd factors proves (2).
Consequently the remainder

```text
R_(h,x)(t)=q_(h,x)(t) mod p_(h,x)(t)
         =sum_(j=0)^(h-1) R_j(x)t^j
```

belongs to `Q[x,t]`, even though t is not invertible over the entire
unspecialized ring `Q[x,t]/(p)`. Giving x and t weight one yields
`deg_x R_j<=2h-j`. Let T_(h,x) be its multiplication matrix in the
basis `1,t,...,t^(h-1)`, and write

```text
C_(h,x)(z)=det(zI-T_(h,x))
         =z^h+sum_(k=1)^h c_(h,k)(x)z^(h-k).                 (3)
```

The entry in row i, column j has degree at most `2h+j-i`. Index sums
cancel in each principal minor, giving the inherited bound

```text
deg_x c_(h,k)<=2hk.                                         (4)
```

The polynomial matrix T, not an undefined value of t^-1 at t=0, is the
object specialized at the boundary parameters in the next section. The
negative-x points are algebraic boundary probes, not claimed instances
of the positive-charge source hypotheses.

## 3. All-height forced-factor theorem

For `1<=k<=h`, put

```text
D_(h,k)(x)=product_(r=2)^h (x+r)^max(0,min(k-h+r,r-1)).         (5)
```

**Theorem.** In `Q[x]`, `D_(h,k)` divides `c_(h,k)` for every h,k.
Equivalently,

```text
D_(h,k)=product_(j=1)^k (x+h-k+j)^j,              k<h,
D_(h,h)=product_(r=2)^h (x+r)^(r-1).                         (6)
```

The empty products are one. Thus the residual coefficient
`b_(h,k)=c_(h,k)/D_(h,k)` is a rational polynomial of degree at most

```text
2hk-k(k+1)/2,                   k<h,
2h^2-h(h-1)/2,                  k=h.                         (7)
```

This is a divisibility statement. Exact multiplicity of these factors
for all heights is not asserted; the exact orders through h=6 below are
**FINITE-EXACT** controls. Neither divisibility nor a smaller degree bound
alone proves positivity of the residual polynomials.

### Proof by the actual small-root response

Fix `2<=r<=h` and put `delta=x+r`. Directly from (1),

```text
p_(h,-r)(t)=t^r a(t),       a(0)!=0,
[p_(h,-r+delta)]_0 = A delta+O(delta^2),       A!=0.          (8)
```

For j<r each p_j has a simple factor delta, while p_r(-r) is nonzero.
The first assertion uses `(h-r)_(h-j)=0` exactly when j<r; the simple
parameter zero follows because that falling factorial has one zero factor.
Likewise,

```text
[q_(h,-r+delta)]_-1 = B delta+O(delta^2),       B!=0,
[q_(h,-r)]_e=0 for 0<=e<2r.                                (9)
```

The first nonzero regular coefficient at delta=0 has index e=2r.
All these facts follow from the same single-zero-factor calculation, now
with `(2h-2r)_(2h-e)`.

To see the small roots explicitly, set `delta=epsilon^r` and
`t=epsilon v`. Then

```text
epsilon^-r p_(h,-r+epsilon^r)(epsilon v)
 = a(0)v^r+A+O(epsilon).                                   (10)
```

The polynomial `a(0)v^r+A` has r distinct nonzero complex roots.
At each one, coefficientwise implicit lifting gives a formal series
`v_i(epsilon)` with that nonzero constant term. This is elementary formal
recursion: the coefficient of the next unknown term is multiplied by
`r*a(0)*v_i(0)^(r-1)!=0`. Thus the r small roots have form
`t_i=epsilon v_i(epsilon)`.

At these roots the **inverse carry** in (9) dominates every regular term:

```text
q_(h,-r+epsilon^r)(t_i)
 = epsilon^(r-1) [B/v_i(0)+O(epsilon)].                      (11)
```

Indeed the carry is `O(epsilon^r)/t_i`; every regular term of index
less than 2r is `O(epsilon^r)*t_i^e`, and every remaining term is
`O(t_i^(2r))`. The leading coefficient in (11) is nonzero. Since
q=R modulo p over the fraction field, (11) also describes the eigenvalues
of the **polynomial** multiplication operator R along these branches.
No value of q at the singular pair `(delta,t)=(0,0)` is used.

For completeness, the r small roots form a monic factor of p over
`C[[delta]]`: at delta=0 the factors `t^r` and a(t) are coprime,
so their coefficients lift uniquely one order at a time. The needed
linear correction equation is solvable by the Bezout identity for these
two coprime polynomials. The two factors remain coprime over `C[[delta]]`.
Chinese remainders split the multiplication operator R into a small block
of size r and a complementary block of size h-r, both with entries in
`C[[delta]]`.

Write the small block's characteristic coefficients as s_j(delta),
`j=1,...,r`. They belong to `C[[delta]]`, while (11) makes each
j-fold eigenvalue product divisible by `epsilon^(j(r-1))`. Therefore

```text
ord_delta s_j >= ceil(j(r-1)/r)
              = min(j,r-1).                                (12)
```

Possible cancellations can only increase this order. The complementary
block has integral formal coefficients, with no negative delta powers.
In the product of the two characteristic polynomials, a contribution to
c_(h,k) must use at least `max(0,k-(h-r))` small eigenvalues. Applying
(12) gives

```text
ord_(x=-r) c_(h,k) >= max(0,min(k-h+r,r-1)).                  (13)
```

The distinct rational linear factors x+r are coprime. Taking their product
proves (5), and summing the exponents gives (6)--(7).

### Carry-deletion boundary

If the e=-1 term is deliberately deleted from (1), its r small response
values instead have order at least delta: the regular terms in (9) give
`O(epsilon^r)`, and the first surviving regular term gives
`O(epsilon^(2r))`. This argument also applies to r=1. Consequently the
no-carry multiplication determinant is divisible by

```text
product_(r=1)^h (x+r)^r.                                    (14)
```

Its forced divisor differs from D_(h,h) by `product_(r=1)^h(x+r)`.
This is an algebraic diagnosis of the missing carry, not permission to
remove it from an actual moment. The checker confirms the distinction
already at h=2: the actual determinant has orders `(0,1)` at x=-1,-2,
whereas the no-carry determinant is divisible by `(x+1)(x+2)^2`.

## 4. First consumer: the actual endpoint-39 family

Now set h=6 and let `x=g-19>=1` be an integer. For nonzero complex
coefficients alpha,beta,gamma put

```text
f(u)=alpha*u^-39+beta*u^(2g-39)+gamma*u^(3g-39),
tau=alpha*gamma^2/beta^3,       X=alpha^x beta^18 gamma,
K=(2g)_26>0.
```

The complete mass-g and mass-2g count fibers are respectively

```text
(x+j,18-3j,1+2j),          j=0,...,6,
(2x+e,36-3e,2+2e),         e=-1,...,12.                     (15)
```

Their multinomial expansions give the exact coefficient identities

```text
CT(f^g)=X*binom(g,13)*p_(6,x)(tau),
CT(f^(2g))=X^2*K*q_(6,x)(tau).                              (16)
```

Thus all seven first channels and fourteen doubled channels, including
its lower carry, are retained. With `gcd(g,39)=1`, the charge equation
modulo g forces every positive admissible mass to be divisible by g.
The j=0 channel shows g is admissible. The gcd condition cannot be dropped:
at g=21, support `(-39,3,24)` has first admissible mass seven.

By THM-4436 with `A=2,B=3,h=6,r=0,z=1`, the first row has six
distinct strictly negative phase roots. To detect the doubled row at
all of them, the producer divides its six characteristic coefficients
by the structural factors (5). The residual degrees are exactly

```text
11, 21, 30, 38, 45, 57.
```

Every one of their **208 monomial coefficients in x** is strictly positive.
The full rational certificate, including positive content, ascending integer
coefficient lists and every removed factor exponent, is frozen in six
`DEFLATED_CERTIFICATE` lines in the [output](third_20260906_laurent.out).
These lists are part of the proof; the sign is not inferred from samples.
No translation of x is needed.

Each D_(6,k) is positive at x>=1. Hence every c_(6,k)(x)>0 there, and
`C_(6,x)(z)>0` for every real z>=0. Evaluation at a real first root rho
makes `q_(6,x)(rho)` a real eigenvalue of T_(6,x), so it must be strictly
negative. Equation (16) proves the doubled moment is nonzero at every
first cancellation. Its real sign is the normalized sign
`CT(f^(2g))/(X^2 K)<0`; no sign is assigned to a general complex raw moment.

Outside the six first roots, the g-th moment is nonzero. Taking all three
coefficients equal to one attains this first alternative; taking
`beta=gamma=1` and alpha equal to any first root attains first detection
at 2g. Both are below the endpoint width `3g`.

The factor theorem has a concrete audit consumer: interpolation of the
six residual polynomials needs only 58 parameter nodes, rather than the
73 needed for the largest unreduced degree. It does not make the next
height automatic; positivity of its residuals remains a separate obligation.

## 5. Why coefficientwise positive matrix pieces do not explain the result

A natural proposed structural explanation is to make the multiplication
operator symmetric and show each parameter coefficient matrix is positive
semidefinite. The canonical Bezout realization already obstructs that
particular approach at h=2.

Let `S_x(t)=(-p'_x(t)R_x(t)) mod p_x(t)`, and define B_x by

```text
[p_x(t)S_x(v)-p_x(v)S_x(t)]/(t-v)
 = sum_(i,j=0)^1 (B_x)_(i,j) t^i v^j.
```

At an actual first root rho its diagonal evaluation is
`-p'_x(rho)^2 R_x(rho)>0`; the known endpoint15 sign therefore makes
B_x positive definite for x>=1. Retain the natural phase scale t=x*y,
and set w=1/x. The normalized symmetric polynomial matrix is

```text
Btilde_(i,j)(w)=w^(6-i-j) (B_(1/w))_(i,j).
```

For w>0 this is a positive scalar times a diagonal congruence of B_(1/w).
Nevertheless its coefficient at w^4 is a positive rational multiple of

```text
[ 49138290   357162621  ]
[ 357162621 2402157940  ],
```

whose determinant is `-9527204358067041`. This coefficient matrix is
indefinite. The first failed implication is passing from positivity of
the whole parameter-dependent quadratic form to positivity of every
coefficient matrix. The strongest survivor is the entire symmetric form;
parameter weights or another basis may still matter. No obstruction to
every possible Gram representation is claimed.

This failure and the forced factors revise the board differently: the
former rejects one direct matrix-positivity explanation, while the latter
identifies a reproducible degeneration law that every successful all-height
proof must respect. The remaining all-height target is positivity of b_(h,k),
not their already-explained negative-integer factors.

## 6. Exact controls and reproduction

The [standalone producer](../../04-computation/third_20260906_laurent.py)
imports no repository implementation. For h=2,...,6 it reconstructs the
rows and polynomial multiplication operator, checks the carry cancellation,
coefficient degree bounds, Cayley-Hamilton identity, and every formal
boundary prerequisite. It verifies that the forced factor orders are exact
in this finite height range; only divisibility is proved for all heights.

For h=6, x=1,...,24, literal multinomials independently verify every first
and doubled coefficient, all first charge equations, and the lower carry.
The g=21 early-admissibility hostile and the exact h=2 Bezout and omitted-
carry controls are retained. These controls challenge the symbolic proof;
they are not the unbounded justification of the endpoint family.

```bash
python3 -B 04-computation/third_20260906_laurent.py
python3 -B -O 04-computation/third_20260906_laurent.py
```

Normal and optimized producer outputs agree byte for byte, with 668
always-active exact gates. Raw LF SHA256 values:

```text
source 29e6fe582f0301e91b9fe27023aaf210269975acc00c7b1c0251db3bfe587103
output 4020a649ea09986521314e5a7d36c1b8df6f9e8883e881ebb23f89487c635221
semantic 412cb08f795638eaf61bbcc8751afbfd8f8beb73f8e63d68641a2b33884fdc71
```

Independent written and coefficient reconstruction audits remain pending.
No external priority or Lean verification is claimed.
