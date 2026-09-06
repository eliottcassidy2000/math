# Four-node metric blindness and the surviving exact Hermite precision law

**PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.** The
all-scale counterexample is proved by a complete exact minor certificate;
the all-node precision formula is proved by explicit inverse columns.
Independent full proof reviews and evidence manifests are in Section 5.

## 1. Inheritance and the question that failed

[THM-4429, arbitrary three-node Smith form](../../01-canon/theorems/THM-4429-arbitrary-three-node-two-jet-smith-form-and-metric-precision.md)
proves that the full p-Smith partition at three nodes depends only on the
pairwise valuation tree. The closest incoming extension is the independently
audited [odd-prime mixed cluster](overnight_20260906_smith_mixed_cluster.md),
which keeps competing saturation costs rather than appending old factors.
Neither theorem asserts metric sufficiency for arbitrary four-node clusters.

The canonical hostile is THM-4419's change of old factors on adjoining a
node. The corrected near miss is the determinant/Smith confusion in
THM-4010. The least-used sidecar is now the **coefficient valuation of an
intermediate minor after its homogeneous scale has been removed**.

The board is: ultrametric tree; full Smith list; intermediate minor ideal;
global dilation; cardinal Hermite polynomials; largest precision loss.
The map from integer nodes to their p-adic distance tree preserves every
pairwise valuation and the determinant valuation. It destroys unit residues.
The cheapest decisive test is an isometric pair with different intermediate
determinantal divisors, not merely a matching determinant or largest factor.

An initial four-node census through diameter32 retained labelled distance
matrices and found no collision. The decisive audit canonicalizes the entire
distance matrix under all24 node permutations. It still has no counterexample
through32, but the first dyadic diameter is40. Changing the search quotient
and changing the depth are separate obligations.

## 2. An isometric all-scale counterexample family

For every integer e>=0, put

```text
A_e=2^e*(0,1,2,5),       B_e=2^e*(0,1,3,4).
```

The labelled map from `(0,1,2,5)` to `(1,0,3,4)` preserves every dyadic
distance. The six valuations in either set are four copies of e and one
each of e+1,e+2. Thus every corresponding three-node restriction also has
the same dyadic Smith list by THM-4429. Nevertheless at e=3,

```text
nodes (0,8,16,40): p-Smith exponents (0,0,4,7,12,16,19,26),
nodes (0,8,24,32): p-Smith exponents (0,0,4,7,12,17,18,26).
```

Their determinant valuations both equal84 and their largest exponents both
equal26. These invariants and all three-node restrictions miss the difference.
The distance-only claim for complete four-node two-jet Smith forms is
therefore **REFUTED**. No theorem at three nodes is contradicted.

Here is the complete all-scale law. Remove the two unit factors from the
zero-node rows. Let delta_j be the valuation of the jth determinantal
divisor of the residual six-by-six matrix, delta_0=0. For A_e and B_e,
respectively, all e>=0 give

```text
delta1=min(2e,e+1),
delta2=min(4e,3e+2),
delta3=min(7e+2,6e+6),
delta4(A)=min(12e+4,11e+6),
delta4(B)=min(12e+4,11e+7),
delta5=17e+7,
delta6=24e+12.
```

Successive differences, preceded by0,0, give the full partition. The lists
coincide exactly when e<=2, and differ exactly by transferring one unit of
valuation from the seventh factor to the sixth for every e>=3. In both
families the largest exponent is7e+5. The depth threshold is exact.

## 3. Native witness and all-scale proof certificate

For nonzero base nodes a,b,c, clearing the zero-node unit block leaves rows
`(x^j)_(j=2..7)` and `(j*x^(j-1))_(j=2..7)` at x=a,b,c. Under scaling
x->2^e*x, a minor on degrees J and r derivative rows is its integer base
determinant times `2^(e*(sum(J)-r))`. Thus every minor contributes one affine
cost `(sum(J)-r)*e+v_2(base determinant)`.

The companion independently evaluates **all923 nonempty minors of each
base matrix** by fraction-free elimination and literal signed permutations.
Its complete nondominated cost lists are exactly the two-term/single-term
minima above. Domination is componentwise in slope and intercept, so this
finite certificate proves the formulas for every real e>=0, in particular
every integer scale. It is not extrapolation from a finite depth sample.

The distinguishing least-degree four-minor uses value(a), derivative(a),
derivative(b), derivative(c) and columns2,3,4,5. Its exact factorization is

```text
-2*a^4*b*c*(a-b)*(a-c)*(b-c)
  *[3a^2-5a(b+c)+10bc].
```

At (a,b,c)=(1,2,5) it is16320=2^6*255; at(1,3,4) it is
12672=2^7*99. The common pair-distance contribution is unchanged, but the
remaining quadratic has a different valuation. The alternative degree12
minor has valuation4 in both families. At shallow depth that competing
minor hides the unit discrepancy; at e>=3 the degree11 witness exposes it.
This is the first failed implication in a metric-only recursion.

The strongest survivors are the exact three-node theorem, determinant
valuation, and (for these families) the largest precision loss. The repaired
full-partition question must retain the valuations of non-Vandermonde minor
coefficients, or an equivalent saturation sidecar. Their minimal sufficient
description for general cluster trees remains OPEN.

## 4. A positive replacement: largest precision for arbitrary integer nodes

Although the full partition is subtler, the largest factor has an exact
formula for **any number n>=1 of distinct integer nodes**. Let

```text
F(X)=product_(i=1)^n (X-x_i),    d_i=F'(x_i)!=0.
```

For values and first Hasse derivatives on integer polynomials of degree<2n,
the largest integer Smith factor is

```text
s_max = lcm_i [ |d_i|^3 / gcd(|d_i|,F''(x_i)) ].                 (H1)
```

Here F'' is the ordinary second derivative, including its factor two.
Consequently the sharp uniform p-adic loss is

```text
L_p=max_i { 2*v_p(d_i), 3*v_p(d_i)-v_p(F''(x_i)) },              (H2)
```

with v_p(0)=infinity and N>=1. Observation modulo p^(N+L_p) recovers every source
coefficient modulo p^N, and one less digit does not suffice uniformly.
This concerns the fixed degree box and full observer, not moving modules.
It is an elementary consequence of the classical Hermite cardinal basis;
no external priority or invention of that basis is claimed.

### Proof from explicit inverse columns

Put Q_i(X)=F(X)/(X-x_i), so Q_i(x_i)=d_i, and
tau_i=Q_i'(x_i)/d_i=F''(x_i)/(2d_i). The cardinal polynomials are

```text
K_i(X)=(X-x_i)*Q_i(X)^2/d_i^2,
H_i(X)=[1-2*tau_i*(X-x_i)]*Q_i(X)^2/d_i^2.
```

They have derivative/value data `(0,delta_ij)` and `(delta_ij,0)` at x_j.
Their coefficient columns are therefore the inverse of the square two-jet
observer. All coefficients are integral combinations of1/d_i^2 and
F''(x_i)/d_i^3. Both denominators are actually attained: the leading
coefficients of K_i and H_i are respectively1/d_i^2 and-F''(x_i)/d_i^3.
When the latter is zero, only the first denominator matters.

The least integer clearing this pair of denominator types is

```text
lcm(d_i^2, |d_i|^3/gcd(|d_i|^3,F''(x_i)))
 = |d_i|^3/gcd(|d_i|,F''(x_i)).
```

Taking the lcm over i gives the least denominator of the whole inverse,
which is the largest Smith factor. This proves(H1); valuations give(H2).
Sharp precision follows in Smith coordinates, as in THM-4419: a primitive
last coordinate multiplied by p^(N-1) defeats a proposed loss L_p-1.

### What this recovers and what it does not

Writing S_i=sum_(j!=i) v_p(x_i-x_j), the same formula is
`L_p=max_i[2S_i+max(0,-v_p(2*sum_(j!=i)1/(x_i-x_j)))]`.
If a node has a unique nearest neighbour with depth f_i, its reciprocal sum
has valuation-f_i. Its contribution is exactly
`2S_i+max(0,f_i-v_p(2))`. This explains the common7e+5 loss of the hostile
families even when intermediate factors differ.

The incoming odd-prime mixed-cluster loss `(2p+1)e+3d` follows immediately
by applying this unique-neighbour rule to the doubled pair and bounding
the other nodes. This recovers its precision, not its complete partition;
the all-minor saturation proof remains indispensable for those extra layers.
No metric-only statement for(H2) at general tied-neighbour configurations
is assumed. The evaluated F'' term is deliberately retained.

## 5. Exact reproduction and independent review

```bash
python3 04-computation/four_node_metric_hostile_overnight_hexagon_sep05.py --height 40
python3 04-computation/hermite_precision_overnight_hexagon_sep05.py
```

Repeat each with`python3 -O`. The first compares all1,846 residual minors
by two independent determinant definitions, full integer Smith forms at
depths0..20, all48 permutations at depth3, and the four corresponding
three-node restrictions. Its complete translated-to-zero dyadic head has
9,880 four-node sets of diameter<=40 and54 exact isometry types, with two
oriented-height witnesses at the first diameter40. The second independently
compiles full matrices and rational inverses for(H1), with original cardinal
value/derivative tests and no import from a minor-profile program.

Both normal/optimized replays have byte-identical output. The metric checker
passes175,798 explicit gates. The precision checker compares1,283 complete
integer observers of sizes n=1..8,71 literal rational inverses and cardinal
value/derivative controls, passing2,761 gates. Its universe retains arbitrary
translations, nonunit scales, n=1 and vanishing-F'' boundaries.

One independent agent reproduced the literal eight-by-eight Smith lists at
depths0..6 and checked the labelled isometry. A second fully audited both
proofs, separately checked all1,846 residual minors by SymPy domain Gaussian
elimination and the symbolic quadratic factor, and independently checked51
fresh full Smith/inverse examples for n=1..5 at p=2,3,5. All passed.

Raw LF-byte manifests:

```text
metric source 26967fd51486c47b70a98a1b00c6f6969b21c1f8a3240fe8fe1fa36068b0b9b8
metric output 5ce62e8f29b739d0527c6b8906039f8e896d23ca3daf7430a8c229d80aba53f5
metric semantic 48b5a19755729118aab7aa587ae14ea37d0328c187caa1ab224f21fe1a62903b
precision source 020f973c3a8ea90a6617a16df63cde4132bdeb9b2571d5ef4b59b93302a38606
precision output 2b34c5bc5751d7987ef72d6af82cdadf1e2e7ef068875204f990b15e415ed391
precision semantic ea42a165fd27a96e64deec663dee77c21ce32301edb9fe1f5098468c0e59a5b3
```

The incoming ternary two-plus-two theorem
[overnight2 Smith double-pair](overnight2_20260906_smith_double_pair.md)
is metric-only at p=3 under its stated two-cluster hypotheses. Our p=2
counterexample does not contradict it. The prime, depth and intermediate
minor coefficient are separate coordinates; the largest-factor formula
is a complementary all-node invariant, not a replacement full partition.
