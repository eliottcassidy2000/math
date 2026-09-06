# A one-bit septenary Smith law and an elliptic residue obstruction for higher jets

**PROVED BY ANALYTIC IDENTITIES + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The [independent proof and computation audit](overnight8_20260906_jets_independent_audit.md)
passes without requested repairs. The full p=7 equilateral four-jet partition, the odd-prime
four-jet largest-loss formula, and the general prime-order ceiling test
below are separate statements with separate scopes. The last is not a
full precision or full Smith formula.

## 1. Inheritance, observer, and the first failed implication

The observer takes all Hasse derivatives of orders `0,...,m-1` at each
of n distinct integer nodes, on integer polynomials of degree below `nm`.
Smith exponents and valuations are local at a specified prime p. The
largest exponent L is the sharp extra coefficient-recovery precision.
The inherited arbitrary-jet reciprocal-inverse equality is
[THM-4443, arbitrary-jet precision and dyadic unit boundary](../../01-canon/theorems/THM-4443-arbitrary-jet-precision-and-dyadic-unit-boundary.md).

The closest proved mechanism is the
[seventh-round full odd-prime three-node three-jet law](overnight7_20260906_oddjets.md),
combined with the incoming
[dyadic three-node three-jet law](three-node-threejet-dyadic-smith-overnight-hexagon-sep05.md).
The canonical hostiles are THM-4443's dyadic three-jet twins and the
already proved
[ternary four-node two-jet residue family](overnight3_20260906_smith_triple_single.md).
Thus this report does not re-claim discovery that odd-prime metric-only
Smith laws can fail somewhere. The corrected near miss is MISTAKE-547's
warning that changing the observed derivative orders changes the
cancellation budget. The least-used sidecars are the *gap between adjacent
minor weights* and the *simultaneous residue of top reciprocal coefficients*.

The concept board is: exact inverse jets; full determinantal ideals;
companion-minor attainment; dilation weights; unit-residue packets; and
finite-field coefficient extraction. The cheap probe increased uniform
three-node multiplicity from three to four. At p=7 it immediately separated
the labeled-isometric nodes

```
(0,7,14): (0,0,0,0,1,2,4,6,7,8,10,10),
(0,7,21): (0,0,0,0,1,2,4,6,7,8, 9,11).            (1)
```

All three pair valuations equal one in both examples. Their determinant
valuation is 48, but L is 10 versus 11. Modulo `7^10`, their kernel orders
are `7^48` versus `7^47`. Thus the difference is visible in the actual
finite-precision observer, not only in an auxiliary polynomial.

Four is the least uniform jet count at which a three-node odd-prime
counterexample can occur: the seventh-round law covers three jets,
[THM-4429](../../01-canon/theorems/THM-4429-arbitrary-three-node-two-jet-smith-form-and-metric-precision.md)
covers two, and the three-node values-only list is `(0,e,2e+d)` at
depths `(e,e,e+d)`. No global integer-diameter minimality or smallest-prime
claim for every full four-jet partition is made.

Targeted canon/results/correction searches found no prior higher-jet
Deuring coefficient consumer or four-jet p=7 law in the truth surfaces
read for this task. No external priority claim is made.

## 2. Exact largest loss for three nodes and four jets at every odd prime

Let the three pairwise depths be `(e,e,e+d)`, where e,d are nonnegative
integers and p is odd. Translation and a p-adic unit change of variable
are unimodular on the complete Hasse observer, so normalized nodes are

```
0, p^e, p^e a.
```

For d>0 choose zero at a closest-pair endpoint, so `v_p(a)=d`; for d=0,
both a and a-1 are units. The normalized a may lie anywhere in Z_p, not
only in an integer sample.

**Close pair, d>=1:**

```
L=11e+7d-[p=5].                                    (2)
```

**Equilateral, d=0:** L=0 if e=0. For e>=1,

```
L=11e-beta,
beta=[p=3]+[p=5]+[p=7 and a mod7 is in {2,4,6}].     (3)
```

The indicators in (3) are mutually exclusive. The septenary condition is
intrinsic: after dividing out the common depth and reducing modulo seven,
one of the three node residues is the average of the other two. This is
unchanged by permutation and by a common affine unit transformation.

### Reciprocal proof and the exact one-digit cap

At a node with differences u,v to the others, the coefficients of
`((u+T)(v+T))^(-4)` through order three are

```
a0=1/(uv)^4,
a1=-4(u+v)/(uv)^5,
a2=(10u^2+16uv+10v^2)/(uv)^6,
a3=-20(u+v)(u^2+uv+v^2)/(uv)^7.                     (4)
```

THM-4443 gives `L=max(-v_p(a_l))` over all nodes and observed orders,
including attainment in the inverse observer.

For the closest pair in (2), u,v have depths e and e+d. The cubic in
the numerator of a3 has a unique lowest monomial after its common factor
20 is removed. Thus its valuation is `3e+[p=5]`, giving (2). Its
lower-order candidates are bounded by

```
8e+4d, 9e+5d, 10e+6d.
```

The last bound deliberately allows possible quinary cancellation in a2.
Each is at most (2), since `e+d-[p=5]>=0`. The outsider contributes at
most 11e, strictly less than (2). This proves the formula even at
`p=5,e=0,d=1`, where a crude order-two upper bound ties the top order.

For the equilateral case, after removing unit denominators and the common
factor 20, the three cubic numerators are, up to sign,

```
C0(a)=(a+1)(a^2+a+1),
C1(a)=(2-a)(a^2-3a+3),
Ca(a)=(2a-1)(3a^2-3a+1).                            (5)
```

They satisfy exact integer identities

```
C0+C1=7(a^2-a+1),
Ca-6C0=-7(3a^2+a+1),
C0=(a+3)(a^2-a+1)+2(2a-1).                         (6)
```

For an odd prime other than seven, simultaneous vanishing modulo p
forces `a^2-a+1=0` and `2a-1=0`, hence p=3. At p=3 the admissible
equilateral class is a=2 modulo three. The quadratic factors accompanying
a+1 and 2-a in (5) are units. Their two linear factors add to three, so
their minimum valuation is exactly one, for every higher lift.

At p=7 all three vanish together precisely at
`a=2,4,6`: modulo seven, `C0=(a+1)(a-2)(a-4)`, and (6) makes the other
two cubics its scalar multiples. At these three residues, `a^2-a+1`
is a unit. The first identity in (6) therefore prevents both C0 and C1
being divisible by 49. Again the simultaneous cap is exactly one.

Every other admissible prime/residue has a unit among (5). Restoring
20 adds one exactly at p=5. This proves that the maximum top-order
candidate is `11e-beta`. For e>=1 every lower order is bounded by 10e,
so the top candidate wins or ties. For e=0 all node differences are units,
and the whole matrix has a unit determinant; all exponents are zero.
This proves (3), including its shallow mask.

## 3. Full p=7 equilateral Smith law from nine minimal-weight minors

For the remainder of this section p=7, the three distances all have depth
e, and `a,a-1` are units. Set

```
kappa=[a mod7 in {2,4,6}].
```

At e=0 all twelve exponents are zero. At every integer e>=1 the complete
nondecreasing list is

```
(0,0,0,0,e,2e,4e,5e+1,7e,8e,10e-1+kappa,11e-kappa). (7)
```

Thus the metric plus **one bit** restores the full partition in this
entire equilateral family. Both bit values occur, so the bit cannot be
discarded. No p=7 close-pair full-partition formula is asserted here.

Clear the four Hasse rows at zero. The residual matrix has columns
`j=4,...,11` and rows `r=0,...,3` at each of the normalized nodes 1,a:

```
R_((x,r),j)=7^((j-r)e) binom(j,r) x^(j-r).
```

Rows 0,...,3 belong to node1 and rows4,...,7 to node a. A minor with
row set I and degree set J has the form `7^(We) P(a)`, where P is integral
and `W=sum J-sum_(i in I) r_i`. For residual ranks one through six the
least possible weights are

```
w=(1,3,7,12,19,27).                                (8)
```

Equality requires the first r columns and the r largest available
derivative orders. There are only nine equality cases in total:

| rank | rows I | W | exact P(a) |
|---|---|---:|---|
| 1 | 3 | 1 | `4` |
| 1 | 7 | 1 | `4a` |
| 2 | 3,7 | 3 | `40a(a-1)` |
| 3 | 2,3,7 | 7 | `200a(a-1)(2a-1)` |
| 3 | 3,6,7 | 7 | `200a^4(a-1)(a-2)` |
| 4 | 2,3,6,7 | 12 | `2100a^4(a-1)^4` |
| 5 | 1,2,3,6,7 | 19 | `560a^4(a-1)^4(7a^2-7a+2)` |
| 5 | 2,3,5,6,7 | 19 | `560a^9(a-1)^4(2a^2-7a+7)` |
| 6 | 1,2,3,5,6,7 | 27 | `1680a^9(a-1)^9` |

The degrees are always `J=(4,...,3+r)`. These are exact determinant
identities; the companion source verifies each one symbolically.

Every minor of ranks one through three has valuation at least `w_r e`.
At ranks four through six, every equality-weight polynomial in the table
is divisible by seven. Every other minor has `W>=w_r+1`, hence valuation
at least `w_r e+e>=w_r e+1`. This proves the universal lower bounds
without expanding those other minors. The requirement e>=1 is essential
to this use of the weight gap.

The table attains all bounds. At rank three the two remaining linear
factors cannot both be nonunits because `(2a-1)-2(a-2)=3` is a unit at
seven. At rank five the first quadratic is always 2 modulo seven, and
the second is `2a^2`; both are units. The constants 2100,560,1680 each
have seven-adic valuation exactly one. Therefore, writing D_j for the
sum of the first j full Smith exponents,

```
(D5,D6,D7,D8,D9,D10)
       =(e,3e,7e,12e+1,19e+1,27e+1).               (9)
```

The confluent determinant gives D12=48e. Equation (3) gives
`D11=D12-L=37e+kappa`. Taking successive differences proves (7).
They are sorted for all e>=1; the last two tie exactly when e=1,kappa=1.

The certificate enumerates the weights of all 12,804 residual minors of
ranks one through six, counts `(64,784,3136,4900,3136,784)`, and verifies
that the nine cases in the table are exactly the equality-weight cases.
It expands only those nine polynomials. The mathematical reduction is
the weight gap, not extrapolation from a small depth bank.

## 4. A short coefficient packet detects ceiling saturation at prime-order jets

This statement concerns a different, broader target. Let p be any odd
prime, `k=(p-1)/2`, `m=k+1`, and let n>=2 integer nodes have every
pairwise depth equal to e>=1. Normalize them as

```
x_i=c+p^e z_i,
```

where the residues `t_i=z_i mod p` are distinct in F_p, so n<=p.
Put

```
F(X)=product_i (X-t_i),
F(X)^k=sum_j c_j X^j,
q=floor((nk+1)/p),
P(F)=(c_(p-1),c_(2p-1),...,c_(qp-1)).               (10)
```

The exact assertion is

```
L <= (nm-1)e,
L = (nm-1)e  iff  P(F) is nonzero.                  (11)
```

If the vector is zero, L is at most `(nm-1)e-1`. This does **not** give
the size of a further drop or the full Smith list. The packet depends
only on one residue digit of the normalized nodes. Its vanishing is
intrinsic under affine unit changes, though the coordinate vector itself
depends on the chosen polynomial coordinate. No universal minimal-encoding
claim is made for that vector.

### Proof

All normalized reciprocal coefficients are integral, so the inherited
inverse law bounds the candidate of observed order l by
`((n-1)m+l)e`. At e>=1 only top order k can reach `(nm-1)e`.

Fix i and write `Q_i(X)=F(X)/(X-t_i)`. In F_p[[T]], because `m=p-k`,

```
Q_i(t_i+T)^(-m)
   =Q_i(t_i+T)^k/Q_i(t_i+T)^p
   =Q_i(t_i)^(-p) Q_i(t_i+T)^k modulo T^(k+1).       (12)
```

The last step uses Frobenius: the denominator is its nonzero constant
modulo T^p, and k<p. Since `F(t_i+T)^k=T^k Q_i(t_i+T)^k`, the top
coefficient is the unit `Q_i(t_i)^(-p)` times
`[T^(p-1)]F(t_i+T)^k`.

For an integer j, the binomial coefficient `binom(j,p-1)` is zero
modulo p unless `j=sp-1`, when it equals one. One elementary proof writes
j=ap+r, 0<=r<p, and uses
`(1+Z)^j=(1+Z^p)^a(1+Z)^r` modulo p. It follows that

```
[T^(p-1)]F(t_i+T)^k
  =sum_(s=1)^q c_(sp-1) t_i^((s-1)p)
  =sum_(s=1)^q c_(sp-1) t_i^(s-1).                 (13)
```

The last polynomial has degree at most q-1<n, because k<p. Its values at
all n distinct residues vanish if and only if every coefficient vanishes.
Thus a unit occurs among the top reciprocal coefficients exactly when
the packet in (10) is nonzero. This proves (11), with no elliptic-curve
theorem used.

For n=2 the packet is the nonzero leading coefficient of F^k, recovering
ceiling saturation. For n=3 or4 it has one coordinate; for more residues
it generally requires a vector. As a hostile to automatic saturation,
the complete residue set has `F=X^p-X`. Its kth power has exponents
`k+(p-1)j`, 0<=j<=k, none congruent to p-1 modulo p. Its packet is zero.
This example does not give its exact loss without retaining lower jets.

## 5. The three-node packet is the Deuring polynomial

For n=3 normalize the residues as `(0,1,a)`, with a different from zero
and one. Let

```
H_p(a)=sum_(j=0)^k binom(k,j)^2 a^j,
k=(p-1)/2.                                        (14)
```

Expanding `X^k(X-1)^k(X-a)^k` gives directly

```
c_(p-1)=(-1)^k H_p(a).                             (15)
```

Thus at every e>=1, for uniform multiplicity `(p+1)/2`,

```
L=(3m-1)e  iff H_p(a) is nonzero modulo p.           (16)
```

The zero set is invariant under the full three-node permutation group.
The binomial coefficients give `H_p(a)=a^k H_p(1/a)` directly.
Reflection of the cubic gives
`H_p(1-a)=(-1)^k H_p(a)`: coefficient degree p-1 is invariant under
translation of a polynomial of degree below2p-1, by the binomial argument
in (13), and `F_a(1-X)=-F_(1-a)(X)`. These two transformations generate
the six cross-ratio transforms. No group action on an unmarked ordered
parameter is silently discarded.

At p=7, (14) becomes
`1+2a+2a^2+a^3=(a+1)(a-2)(a-4)`, recovering exactly the arithmetic-
progression bit in (3) and (7). At p=3 its sole root2 is the only allowed
parameter; at p=5 it has no root in F5. Those facts explain why the
seventh-round metric-only three-jet conclusion does not contradict (16).

There is a uniform family of odd-prime metric failures: if p>=7 and
p=3 modulo four, k is odd. The reciprocal symmetry pairs terms to give
`H_p(-1)=0`. Its degree k is less than the p-2 admissible parameters,
so at least one admissible parameter has nonzero value. Their equally
spaced-depth observers have different precision by (16). This conclusion
needs neither a classification of all roots nor a bound on their higher
p-adic lifts. The number of observed jets varies with p.

**CITED interpretation, not an algebraic proof dependency:** Roland Auer
and Jaap Top, *Legendre elliptic curves over finite fields*,
[arXiv:math/0106273v1](https://arxiv.org/pdf/math/0106273), Section3,
printed pages5-6 (paragraph preceding Proposition3.2), identifies the
roots of this polynomial, up to a harmless global sign convention, with
supersingular parameters of `Y^2=X(X-1)(X-a)`. The primary PDF was read
after checking CORE-PAPERS. Only this three-node identification is
imported; the source does not establish any Smith law here.

The map therefore sends the normalized three-node observer to its
Legendre residue curve and preserves **failure to attain the top
precision ceiling** as supersingularity. It loses common depth, higher
unit digits, lower reciprocal orders, and the intermediate determinantal
ideals. For p7/four jets the separate cubic cap and minimal-weight proof
restore the exact loss and full partition. For general p/m those stronger
conclusions remain outside (16).

## 6. Computation, controls, and remaining question

The companion
[source](../../04-computation/overnight8_20260906_jets_residue.py)
uses no repository imports. Its mathematical identities are proved above;
finite controls provide independent representations rather than universal
sample extrapolation.

- Exact cubic identities, resultant `4116=2^2*3*7^3`, and the full joint
  cap on complete small residue-lift ranges at p3,5,7,11,13,17,19.
- 1,065 literal rational inverse-jet loss checks, spanning close pairs,
  equilateral configurations, several outer/inner depths, negative units,
  and all admissible small residue lifts.
- 51 literal12x12 integer Smith controls for the odd-prime precision law
  and the hostile twins; signed affine and node-reversal controls included.
- The weights of all12,804 relevant p7 residual minors; all nine exact
  equality-weight polynomial identities; complete joint-attainment controls;
  and175 additional full Smith matrices across35 unit lifts and five depths.
- 208 general coefficient-packet checks, using every residue subset at
  p<=7 and named deterministic configurations thereafter. Each packet is
  checked against both direct translation and separately multiplied
  reciprocal jets; 158 are saturated and50 have simultaneous cancellation.
- Complete Deuring root sets for p3,5,7,11,13,19,23,31, checked against
  all endpoint top jets, both generating permutation identities, and an
  independent finite Legendre point count. These finite point counts are
  controls, not the source of the cited general interpretation.

The companion prints **20,319 optimization-live gates PASS**. Reproduce:

```
python -B 04-computation/overnight8_20260906_jets_residue.py
python -B -O 04-computation/overnight8_20260906_jets_residue.py
```

Normal and optimized output equality and source/output hashes are recorded
in the freeze note below. Dyadic use of the odd-prime predictor is rejected.
The minimal-weight argument explicitly retains the e>=1 boundary; the
e=0 observer is handled by its unit determinant.

The strongest next question is now specific: at a zero of the coefficient
packet, which higher unit digits and lower reciprocal orders determine the
actual precision drop, and which additional minimal-weight minor residues
restore the full partition? A top-jet scalar alone cannot be assumed to
answer those questions. The
[sixth-round normalized-jet recursion](overnight6_20260906_jets_sidecar.md)
already provides an exact adaptive state for the largest loss; the present
packet is a much smaller exact observer only for ceiling saturation.

All files remain outside the repository for parent-managed integration.
No Git mutation, shared navigation edit, theorem reservation, or external
priority claim was made by this lane.

### Freeze note

Normal and optimized runs are byte-identical, with LF-only output:

```
source 102225c5476d59b3aec02c5b479093d105298000b9411fdb4526f9423efd3b65
output f77f84da627c5c1f33468326baf24cc93c1bcf222028a935a956728f1d878061
```

The independent referee reconstructs all12,804 minor weights and the nine
polynomials using Bareiss determinants and degree-bounded interpolation.
Its separate standard-library p-local elimination checks91 full Smith
matrices, and independently recursive inverse series verify155 coefficient
packets and33 complete-observer ceiling tests. It also verifies the complete
S3/Deuring residue controls at eight primes and reads the cited primary
Auer--Top passage directly. The proof audit accepts all four stated scopes,
including the source normalization, odd-prime lower-order competition,
weight-gap boundary, simultaneous packet test, and limited interpretation.
Normal and optimized audit runs are identical and pass26,520 live gates:

```
audit source 4e3d99d9ef917ca5026a4348e71058471e54ee20ff8af91a19a5ffb026197ca1
audit output a32847a7def96f074895538bbc4cffc6af908eb975481871c9a41955f9e60efa
```

**Filing:** root integrated these audited artifacts in the eighth checkpoint;
reproduction paths are relative to the repository root. Earlier outside-worktree
notes describe author provenance, not the present file location.
