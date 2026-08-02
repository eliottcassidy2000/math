# The transport deficit is one minimal fifth-order Schur circuit

Date: 2026-08-02  
Scope: exact scout and stopping boundary; **not a theorem and not a proved
dependency**

## 1. Inheritance and question

`THM-3110-arbitrary-anchored-product-gamma-dominant-tail-and-low-histogram-reduction`
constructs the two signed `24/25` response banks for an arbitrary anchored
support `{0,a,b}`, proves their common-root quotient, and finds a uniform
first-order row-majorization transport of mass `28`.  The unmatched positive
and negative masses are respectively `9` for `I1` and `11` for `I2`.

The natural next guesses were:

1. pair each unmatched debt edge with one forward edge;
2. split the aggregate into smaller row-slide circuits; or
3. regard the result as an ordinary totally-positive divided-difference or
   Peano-kernel matrix.

The exact answer is sharper.  Each residual core has a complete transport in
the *reverse* dominance direction.  Forward and reverse edges together form
one full-support circuit which annihilates every Schur flag through degree
four and exits the kernel positively in degree five.  There is no proper
subcircuit in any of the four chamber universes.  Ordinary coefficientwise
total positivity nevertheless fails already in a raw `2x2` minor.

## 2. The two chamber-independent reverse cores

Write a row as the multiset of interval lengths whose root alphabet is the
union of `[0,L)`.  The signed residual rows are:

```text
I1
+2 (2a,2a,a,b,b,b)
+2 (2a+b,a,a,a,b,b)
+5 (3a,a+b,a,b,b)
-1 (2a,a,a,a,3b)
-3 (2a,2a,a+2b,b)
-3 (2a+b,2a,a,2b)
-2 (3a,a+b,a+b,b)

I2
+3 (2a,a+b,a,a,b,b,b)
+3 (3a,a,a,2b,b,b)
+4 (3a,a+b,a+b,b,b)
+1 (3a,2a,b,b,b,b)
-3 (2a,a+2b,a,a,2b)
-2 (2a,2a,a,3b,b)
-1 (2a+b,2a,a+b,b,b)
-4 (3a,a+b,a,2b,b)
-1 (3a,2a,2b,b,b).
```

In both ratio chambers every negative supply can be transported to a
positive demand which it row-majorizes.  The reverse capacities are exactly
`9/9` and `11/11`.  Thus the residual functional itself has the opposite
Schur sign.  For example,

```text
Phi_res,1(s_(1))=(a+b)(2a-9b)<0,
Phi_res,2(s_(1))=3a^2-2ab-14b^2<0                     (1)
```

for `0<a<b`.  The reverse core is a genuine debt subtracted from the
28-mass forward surplus, not a second positive transport.

## 3. Explicit edge order

Within each original bank, index positive rows `P0,P1,...` and negative rows
`N0,N1,...` in the exact lexicographic order returned by the THM-3110
companion.  `Pi>Nj x c` denotes `c` copies of the forward majorization edge;
`Nj>Pi x c` denotes a reverse-core edge.  One deterministic set of flows is:

```text
I1 tight F:
P0>N0x1 P1>N3x4 P2>N3x2 P3>N0x1 P3>N1x3 P3>N5x2
P5>N2x1 P6>N8x1 P7>N5x3 P8>N5x1 P8>N6x3 P8>N9x2
P9>N11x1 P10>N9x1 P10>N10x1 P11>N10x1
I1 tight R:
N2>P4x1 N4>P4x1 N4>P6x2 N7>P9x3 N10>P9x2

I1 wide F:
P0>N0x1 P1>N3x4 P2>N3x2 P3>N0x1 P3>N1x3 P3>N5x2
P5>N2x1 P6>N8x1 P7>N5x3 P8>N5x1 P8>N6x3 P8>N10x2
P9>N11x1 P10>N9x2 P11>N9x1
I1 wide R:
N2>P4x1 N4>P4x1 N4>P6x2 N7>P9x3 N10>P9x2

I2 tight F:
P0>N2x4 P1>N4x1 P1>N8x2 P2>N4x2 P3>N6x1 P4>N1x2
P5>N9x3 P6>N3x2 P6>N10x4 P7>N8x4 P7>N9x2 P9>N0x1
I2 tight R:
N5>P3x3 N7>P8x2 N9>P8x1 N11>P10x4 N12>P11x1

I2 wide F:
P0>N2x2 P0>N3x2 P1>N0x1 P1>N4x1 P1>N8x1 P2>N4x2
P3>N6x1 P4>N1x2 P5>N9x3 P6>N2x2 P6>N10x4 P7>N8x5
P7>N9x1 P9>N9x1
I2 wide R:
N5>P3x3 N7>P8x2 N9>P10x1 N11>P8x1 N11>P10x3 N12>P11x1.
```

Expanding forward differences positively and reverse differences negatively
recovers the original signed `24/25` bank exactly, row by row.

## 4. The exact flagged matrices

For an edge `e=(H>L)` put

```text
d_e(lambda)=s_lambda(S_H)-s_lambda(S_L).               (2)
```

For each partition `lambda` of sizes one through four, expand `(2)` exactly
in the binomial basis

```text
binom(a,i) binom(b,j).                                  (3)
```

A matrix row is one nonzero coefficient flag
`(degree,lambda,i,j)`; a column is one edge type in the order above.  Since
`s_lambda` has parameter degree at most `2|lambda|`, the rectangular forward
difference extraction through degree eight is exact, not interpolation
evidence.  The matrices are:

| circuit | flagged matrix | rank | nullity | forward/reverse columns |
|---|---:|---:|---:|---:|
| I1 tight | `330 x 21` | 20 | 1 | 16 / 5 |
| I1 wide | `330 x 20` | 19 | 1 | 15 / 5 |
| I2 tight | `328 x 17` | 16 | 1 | 12 / 5 |
| I2 wide | `328 x 20` | 19 | 1 | 14 / 6 |

The kernel vector is the edge-capacity vector, positive on forward columns
and negative on reverse columns.  Every coordinate is nonzero.  Consequently
the aggregate is the unique universal relation among these edge currents
through degree four, and its support is minimal.  A smaller edge-slide
circuit or a common-intermediate decomposition inside these flows cannot
exist: it would give a second kernel vector or a kernel vector with proper
support.

Deleting the common root alphabet is an invertible triangular change on
Schur flags up to any fixed degree.  Hence the undeleted coefficient-matrix
kernel and the common-root-quotient kernel are the same.

## 5. Equality through four and the degree-five exit

The circuit annihilates `s_lambda` for every `|lambda|<=4`.  At degree five
it leaves the kernel.  With

```text
L1=3a+5b,       base1=a^4 b^2 (b-a)^3,
L2=3a+4b,       base2=a^3 b^2 (b-a)^4,
kappa_lambda=sum_i lambda_i(lambda_i-2i+1),             (4)
```

the direct circuit convention is

```text
Phi_j(s_lambda)
 =base_j f^lambda/2 (L_j+kappa_lambda/4),  |lambda|=5. (5)
```

This is the conjugate-index form of THM-3110's displayed formula for
`A_(j,mu)=Phi_j(s_(mu'))`; conjugation reverses `kappa`.  Since
`kappa_lambda>=-20` and `0<a<b`, `(5)` is strictly positive.  Thus the
object really is an order-five Schur circuit, not merely a numerical
nullspace.

This order-five exit is a promising rooted-flag invariant.  It is not by
itself an all-degree B-spline proof.

## 6. The smallest total-positivity hostile

Order columns canonically at the representative tight support `(a,b)=(2,3)`
by the high Ferrers partition, then the low Ferrers partition, both in
dominance/lexicographic order.  Three columns in that order have their
degree-one edge currents

```text
c0: (2a,2a,a,3b) > (2a,a,a,a,3b),
    d_c0(s_(1))=a^2;

c2: (a,a,a,a,a,3b) > (a+b,a,a,a,a,2b),
    d_c2(s_(1))=2b^2-ab;

c4: (2a,2a,a+2b,b) > (2a,2a,a,b,b,b),
    d_c4(s_(1))=b^2+2ab.                              (6)
```

Take the raw monomial coefficient rows `[a^2]` and `[ab]`.  The three
columns are respectively `(1,0)`, `(0,-1)`, and `(0,2)`.  Therefore

```text
det(c0,c2)=-1,             det(c0,c4)=+2.              (7)
```

The same coefficient rows and canonical column order already contain both
minor signs.  In the binomial-Newton flags the values are `-2,+4`, so the
hostile is not a monomial-basis scaling artifact.

This kills an ordinary totally-positive matrix, positroid, extended complete
Chebyshev system, or coefficientwise nonnegative Peano kernel before any
`3x3` or maximal-minor test.  Notice that every edge function in `(6)` is
positive on the chamber; the hostile comes from the chamber-adapted
combination `2b^2-ab`.  Any successful positivity mechanism must therefore
remember the chamber Newton cone or the rooted rank-four flag.  Raw
coefficientwise total positivity forgets exactly that sidecar.

## 7. Reproduction and scope

Run

```text
python 04-computation/gmc_product_gamma_order5_transport_circuit.py
python -O 04-computation/gmc_product_gamma_order5_transport_circuit.py
```

and compare with

```text
05-knowledge/results/gmc_product_gamma_order5_transport_circuit.out.
```

The companion derives both transports from the live THM-3110 banks, checks
the row-current identity, builds every exact coefficient flag, verifies the
four ranks and full-support kernels, checks the degree-five exit on
degree-bound-determining chamber grids, and derives `(6)--(7)` symbolically.

This scout neither extends the proved degree range nor proves arbitrary
anchored width-three product-Gamma goodness, GMC(2), NC2, LRC(14), JC(2), or
DC(2).  Its positive content is the minimal order-five circuit.  Its stopping
content is the exact TP2 sign conflict.
