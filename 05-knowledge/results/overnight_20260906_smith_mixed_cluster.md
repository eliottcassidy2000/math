# The first odd-prime mixed two-jet cluster

**Status: PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**
The claim concerns one complete residue layer with exactly one doubled
residue, including arbitrary further depth and arbitrary integer lifts.
The exact audit is supporting evidence; the proof below carries all-prime,
all-depth quantifiers. No external novelty or priority claim is made.

The root reviewer and the separate moments lane both accepted the complete
proof. The latter's
[independent analytic audit](overnight_20260906_moments_smith_audit.md)
checks the arbitrary-lift divided row, saturation/extension step, every
weighted-minor alternative, p=3 and p=2 boundaries, and the integral CRT
corollary. The exact independent algorithm and all-minor checks are recorded
in Section 7.

## Inheritance and active concepts

The closest proved mechanism is
[THM-4419 / twojet-prime-wall-precision-and-dyadic-triple-smith-law](../../01-canon/theorems/THM-4419-twojet-prime-wall-precision-and-dyadic-triple-smith-law.md),
including the concurrent
[complete-residue wall proof](synthesis_20260905_wildcard_smith_boundary.md).
The canonical hostile is THM-4419's dyadic adjoining example: adding a node
can change existing Smith factors. The corrected near miss is
[THM-4010 / confluent-consecutive-hasse-observer-kernel-index-and-smith-firewall](../../01-canon/theorems/THM-4010-confluent-consecutive-hasse-observer-kernel-index-and-smith-firewall.md)'s
false repeated-factorial diagonal. The least-used sidecar is the competing
cost of saturating a duplicated derivative row versus retaining one fewer
derivative row in a minor.

The board consists of: geometric cluster depth e; duplicate depth d; the
Frobenius derivative kernel; weighted determinantal divisors; and the
distinction between real and modular collinearity. The first decisive test
was e=2,d=1. It changes a preexisting invariant, so the full mixed observer
must be analyzed instead of appending two factors to an old partition.

## 1. Statement

Let p be an **odd** prime, e>=1, and let u_0,...,u_p be integers. Their
residues cover F_p, with exactly one residue occurring twice and every other
residue once. If u_0,u_p denote the duplicated pair, assume

```text
delta=u_p-u_0 != 0,       d=v_p(delta)>=1.
```

All remaining pairwise u-differences are p-adic units. Put

```text
x_i=a+p^e u_i,          a in Z,
m=min(e,d).
```

On Z[X] of degree below 2p+2, form the square observer

```text
J(P)=(P(x_i),P'(x_i))_(0<=i<=p),
```

where P'=D^[1]P. Its increasing p-Smith exponent list is

```text
0,0,
e,2e,...,(p-2)e,
(p-1)e+1,
pe+m,
(p+2)e,(p+3)e,...,(2p-2)e,
(2p-1)e-1,
2pe+d-m,
(2p+1)e+3d.                                           (1)
```

The middle range `(p+2)e,...,(2p-2)e` is **empty when p=3**. In that case
the complete eight-entry list is

```text
(0,0,e,2e+1,3e+m,5e-1,6e+d-m,7e+3d).                 (2)
```

No ranges overlap. The displayed list is already sorted. In particular
`pe+m<=(p+1)e`, and the last three gaps are positive except for the inherited
equal pair at e=1 immediately before the two new terminal factors. At
p=3,e=1, that equal pair is explicitly the two entries 4,4 in (2).

## 2. Determinantal-divisor form

Let nu_h=v_p(Delta_h), with nu_0=0. The equivalent exact formulas are

```text
nu_h = e*(h-1)*(h-2)/2,                     1<=h<=p,
nu_(p+1) = e*p*(p-1)/2+1,
nu_h = e*(h*(h-1)/2-(p+1))+1+m,           p+2<=h<=2p-1,
nu_(2p) = e*(2p^2-2p-1)+m,
nu_(2p+1) = e*(2p^2-1)+d,
nu_(2p+2) = 2ep(p+1)+4d.                             (3)
```

The range in the third line contains one rank, h=5, at p=3. Consecutive
differences of (3) yield (1), including the empty exponent range in (2).
The final equality agrees with the full confluent Vandermonde determinant.

## 3. Residual derivative lattices

We work over the DVR R=Z_(p). Translation in X is unimodular, so set a=0.
Let q_1<...<q_h be a selected set of monomial degrees, and let r be the
number of selected derivative rows. The exact valuation of a nonzero minor
is

```text
e*(sum_j q_j-r)+v_p(residual minor),                  (4)
```

with residual rows `u_i^q` and `q*u_i^(q-1)`.

Several elementary rank facts determine the residual cost. In every claim
below, the columns are consecutive degrees 0,...,h-1.

**A. One representative of each residue, p derivative rows.** For
`p<h<2p`, their row lattice has saturation index p. Its reduction has rank
p-1: derivative monomials q=1,...,p-1 give functions of degrees 0,...,p-2;
q=p vanishes; and q=p+1,...,2p-2 reduces to the same function span on F_p.
The columns q=1,...,p have determinant
`p!*product_(i<j)(u_j-u_i)`, of exact valuation one. At h=2p, column
q=2p-1 supplies the missing value function X^(p-1), and the rank is p.

**B. All p+1 derivative rows.** For `p+2<=h<2p`, their saturation index
is p^(d+1). At h=2p their saturation index is p^d.

Here is a proof that keeps the duplicate depth. Subtract the derivative row
at u_0 from the row at u_p and divide by delta. The new row is integral,
and modulo p it acts as P'' at the common residue. After this division the
rows modulo p are values of P' at every residue together with P'' at the
duplicated residue.

In the lower range, the P' value rows have rank p-1 by A. The additional
row is independent: the polynomial

```text
P_*(X)=X^(p+1)/(p+1)-X^2/2                            (5)
```

belongs to R[X] of degree p+1<h, and its derivative modulo p is X^p-X.
That derivative vanishes at every residue, while its derivative is -1.
Thus the divided row raises the rank to p, leaving exactly one missing
direction. All maximal minors of the divided derivative matrix are
divisible by p. The original minor on columns q=1,...,p+1 equals

```text
(p+1)! * product_(i<j)(u_j-u_i),                      (6)
```

of exact valuation d+1. Therefore the total saturation cost is exactly d+1.
At h=2p, A gives p independent value rows and (5) again adds the independent
extra row. The divided matrix has full rank p+1 modulo p, so the total
saturation cost is exactly d. This proof uses only the residue data and
v_p(delta); arbitrary lifts and the unit part of delta are retained.

**C. Extension by value rows.** For every h<=2p, the full residual observer
has rank h modulo p. One complete set of p residues already suffices: a
polynomial in its kernel has every residue as a double root, so is divisible
by `(X^p-X)^2`, of degree 2p. If a chosen derivative bank is saturated using
unimodular operations on its own rows followed by divisions by its invariant
p-powers, its new span modulo p contains its old span. Consequently it can
be extended to rank h by selecting the needed number of value rows.
Undoing the saturation gives a minor with exactly the derivative bank's
saturation valuation. This is a row-lattice argument; there is no assumption
that a chosen Smith basis survives adjoining a node.

## 4. Weighted lower bounds and attainment

For h<=p, the lower bound and witness from
[THM-4080 / confluent-two-jet-single-scale-smith-partition](../../01-canon/theorems/THM-4080-confluent-two-jet-single-scale-smith-partition.md)
apply unchanged: one value row and h-1 derivative rows at distinct residues,
on degrees 0,...,h-1, have a unit residual determinant. The factorial
factor is `(h-1)!`, still a unit at h=p. This proves the first line of (3).

At h=p+1 the least possible scale cost is e*p*(p-1)/2. To attain that
cost one must use one value and p derivative rows on degrees 0,...,p.
Their residual rank drops modulo p, so the minor costs at least one further
power of p. Any other viable row/degree choice costs at least one additional
geometric unit e>=1. Choosing the p derivatives at distinct residues and
one value row yields the p!-Vandermonde witness of exact residual valuation
one. This proves the second line of (3). An all-derivative minor has zero
column zero and necessarily has strictly greater scale cost.

Now let p+2<=h<=2p-1, and put `L=h*(h-1)/2-(p+1)`.

- With all p+1 derivative rows, a least-cost minor uses degrees 0,...,h-1
  and costs at least eL+d+1 by B. A nonminimal degree set costs at least
  eL+e+d, because the duplicated derivative rows always contribute delta.
- With p derivative rows, the least degree set costs at least eL+e+1:
  its derivative rows have rank at most p-1 modulo p, even if the duplicated
  residue is retained. Any more expensive degree set costs at least eL+2e.
- With at most p-1 derivative rows, the cost is at least eL+2e.

Since e,d>=1, all of these bounds are at least
`eL+1+min(e,d)`. By B and C, the bank of all p+1 derivative rows attains
eL+d+1; by A and C, p derivatives at distinct residues attain eL+e+1.
The smaller of these two witnesses proves the third line of (3).

At h=2p put `L=2p^2-2p-1`. A minor with p+1 derivative rows contains
the duplicated pair and hence costs at least eL+d. With at most p
derivative rows it costs at least eL+e. B and C attain the first bound;
using exactly p derivatives and p values at distinct residues attains the
second, because that full residual Hermite matrix is invertible modulo p.
Their minimum is eL+min(e,d), proving the fourth line.

At h=2p+1, any selected row set contains at least one equal-residue pair
of the same Hasse order: there are two such pairs among all 2p+2 rows, and
only one row can be omitted. Subtracting those two rows shows that every
residual minor is divisible by delta, whatever its column set. The scale
cost is at least e*(2p^2-1).

To attain both bounds, take all p+1 derivative rows and p value rows, one
value at each distinct residue, on degrees 0,...,2p. Divide the duplicated
derivative difference by delta. Modulo p the resulting matrix consists of
values and first derivatives at all p residues, plus second derivative at
the duplicated residue. Since 2 is invertible, its kernel has a triple
root at that residue and double roots at every other residue. Its degree
would be at least 2p+1. Thus the square matrix on degree below 2p+1 is
invertible, giving exact valuation d. This proves the fifth line of (3).

Finally the full determinant is
`product_(i<j)(x_j-x_i)^4`. There are p(p+1)/2 pairs of basic depth e,
and the exceptional pair contributes an additional 4d. This proves the
last line of (3) and completes the theorem.

## 5. Ordering, the adjoining boundary, and consequences

Delete either member of the duplicated pair. The remaining p-node observer
has the THM-4419 profile. Comparing (1) with that profile gives an exact
repair of the tempting unchanged-prefix statement:

```text
The first 2p Smith exponents survive adjoining
if and only if d>=e.                                  (7)
```

If d<e, precisely exponent p+2 decreases, from `(p+1)e` to `pe+d`.
The terminal new factors absorb the remaining determinant mass. Thus
adjoining is not generally a direct sum of the old observer with a tail.

At d=e, m=e: the first 2p exponents agree with the old profile, followed
by `2pe` and `(2p+4)e`. At d>e, the prefix still agrees, while the final
two exponents are `(2p-1)e+d` and `(2p+1)e+3d`, in strictly increasing
order. For d<e the penultimate exponent is exactly 2pe. These formulas
also verify the ordering across the transition d=e without guessing it.

For the smallest odd-prime hostile,

```text
p=3, e=2, d=1:
old nodes (0,9,18):    (0,0,2,5,8,9),
new nodes (0,9,18,27): (0,0,2,5,7,9,12,17).           (8)
```

The formula does not extend to p=2. At e=d=1 the actual nodes (0,2,4)
have THM-4419 profile `(0,0,2,2,5,7)`. In particular the h=5 valuation
is 9, whereas the odd-prime argument would predict 8. The failed step is
precisely the last rank argument: the ordinary second derivative is twice
the second Hasse derivative and is zero modulo 2. Formula (5) also uses
the nonunit denominator 2 there. The dyadic theorem remains separate.

**Consecutive-node extension.** For every odd prime p and
`1<=n<=p(p+1)`, the complete p-primary two-jet partition on n consecutive
nodes is now explicit. Split by residues modulo p using THM-4010's local
CRT. Each nonempty class has at most p+1 nodes. For sizes below p use
THM-4080; for size p use THM-4419; for size p+1 use (1) at e=d=1.
The latter class is `c+p*(0,1,...,p)`, so it has exactly the declared
duplicated residue and no unrecorded collision. At n=p^2+r, 0<=r<=p,
there are r classes of size p+1 and p-r of size p. The full partition is
their sorted multiset union. This extends the previously proved quadratic
range by the entire next p-node band.

**Sharp precision.** The worst coefficient loss in this mixed cluster is

```text
L=(2p+1)e+3d.                                        (9)
```

Observations modulo p^(N+L) determine all source coefficients modulo p^N,
and the uniform loss is sharp, by the Smith-coordinate argument of
THM-4419. Finite-level kernel counts are obtained from (1) as
`p^(sum_i min(N,alpha_i))`. These concern this fixed integral degree box.

## 6. Connection to the no-three-in-line seed, with its limit

There is a concrete determinant map, but it does not transfer extremal
grid occupancy. Send an integer node t to the parabola point (t,t^2).
For three distinct nodes, the oriented area determinant is the Vandermonde
product of their three pair differences. Thus every triple is noncollinear
over the rationals and reals, whereas reduction modulo p can make its area
zero at arbitrarily high p-adic depth. If V is the Vandermonde product of
all N nodes and A_ijk are these oriented area determinants, then

```text
product_(i<j<k) A_ijk = V^(N-2),
det(J_twojet)=V^4,
det(J_twojet)^(N-2)=product_(i<j<k) A_ijk^4.           (10)
```

The map preserves exact determinants and their p-adic valuations. It loses
the prescribed rectangular grid size, slope incidence on general point
sets, and any distribution of extremal grid configurations. Its necessary
sidecar is the actual point coordinates and real determinant predicate;
modular collinearity is not real collinearity. The present mixed-cluster
theorem refines the aggregate determinant valuation in (10) into every
interpolation precision layer. It proves no larger no-three-in-line set.

Primary-source reconnaissance checked the published records for Bhargava's
[*On P-orderings, rings of integer-valued polynomials, and ultrametric
analysis*](https://doi.org/10.1090/S0894-0347-09-00638-9) and Grinberg--Petrov's
[*A Greedoid and a Matroid Inspired by Bhargava's p-Orderings*](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v28i3p6).
They supply relevant background on arithmetic interpolation and ultrametric
greedy structures. No result from either paper is imported into this proof,
and no claim that they exclude a prior form of (1) is made.

## 7. Reproduction and remaining frontier

The companion `04-computation/overnight_20260906_smith_mixed_cluster.py`
records its exact universe, independent modular and rational reductions,
small all-minor audits, arbitrary residue lifts, p=2 and unchanged-prefix
hostiles, and direct consecutive-node checks. Its matching transcript is
`05-knowledge/results/overnight_20260906_smith_mixed_cluster.out`.

The completed bank contains 144 full mixed matrices over p=3,5,7,11,
e=1,2,3, d=1,2,3,4, and three distinct lift systems; 132 divided-derivative
rank checks; four exhaustive rank-eight all-minor cases; and eight direct
consecutive matrices through n=30. Rational DVR and independently compiled
modular DVR reductions agree on every full partition. The modular precision
exceeds the independently known total determinant valuation, so it does not
assume the predicted largest invariant. There are 1,821,396 exact gates.
Normal and optimized replays produce identical LF transcripts.

```text
source_sha256:
3f2c90a6109177b66faba9a9378f7ac8f1ae80962966e173e2e6f67aaf86f7de
output_sha256:
4a7c801c9578f954395fa9bbc696f55133f31e131d821a917f807171735c6f10
semantic_sha256:
90f614c1c30dda130488af1659295122a08bec217ffccb651904d683682ab54c
```

```text
python3 -B 04-computation/overnight_20260906_smith_mixed_cluster.py
python3 -B -O 04-computation/overnight_20260906_smith_mixed_cluster.py
```

On Windows aliases that do not execute script paths, use
`python3 -c 'import runpy; runpy.run_path("04-computation/overnight_20260906_smith_mixed_cluster.py",run_name="__main__")'`,
and add `-O` for the optimized replay.

Two doubled residue classes, a triple residue, arbitrary larger clusters,
higher jets, and a general cluster-tree Smith recursion remain OPEN. The
first new consecutive case outside this result is n=p(p+1)+1. The
all-depth result here includes every choice of duplicate depth and lift
within its stated single-doubled-residue shape.
