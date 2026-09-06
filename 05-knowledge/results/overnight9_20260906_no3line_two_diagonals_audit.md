# Independent audit of the two-direction diagonal covariance theorem

**Verdict: PASS — analytic proof accepted, with independent exact controls.**
No mathematical correction to the producer report was needed. This audit
accepts the covariance and uniform random-density conclusions, not an
extremal nonexistence theorem or optimality of the probability constant.

## Scope, inheritance, and independence

The audited report is
`overnight9_20260906_no3line_two_diagonals.md` in this directory. Its model
is a fixed simple bipartite 2-regular graph on n+n vertices with independent
uniform shore labels. The eighth diagonal-density theorem supplies the
uniform matching and colored-overlap estimates and the one-direction
variance. The new information is a grid cell selected by both factorial
tuples, together with its integer crossing parity.

The producer source was neither read nor imported. I derived the overlap
orders and geometric integrals, implemented a separate verifier, and only
then read the producer proof. The producer enumerates shore labelings;
this referee enumerates distinct boards by unordered column pairs in each
row, then classifies components by graph traversal. Within one cycle type
the two measures agree because every board has the same number of shore
labelings. Thus automorphism multiplicity cannot alter their comparison.

The concept comparison is: the LRC anchor retains endpoint owners lost by
scalar mass; the Laurent niche retains a carry coordinate; this wildcard
retains a shared cell lost by an edge-disjoint covariance projection. Here
the decisive cheap test is the exact two-edge moment with one common cell.
The inherited C8 versus 2C4 hostile prevents inferring a finite-size sign
or cycle-independent law from a uniform leading asymptotic.

## Analytic review

For fixed a,b>=1 put k=a+b. A pair of partial matchings T,U has lengths
L,M, common-row/column counts R,C, and at most one shared cell Z. A cell
selected once in each tuple must be an isolated edge of the selected
incidence graph: neither matching can contain another edge at either of
its endpoints. There is at most one such repeated cell.

Before excluding other intersections its tuple count is exactly
`ab Z (L-1)_(a-1)(M-1)_(b-1)`. For k>=3 the fully disjoint remainder has
count `ab Z L^(a-1)M^(b-1)+O(n^(k-3))`; for k=2 it is exactly Z.
There are only k-1 distinct cells, so multiplying by their matching
probability creates the positive term
`ab 2^(k-1) Z (L/n)^(a-1)(M/n)^(b-1)/n`.
This is a first-order contribution and cannot be discarded.

For a selected colored pattern with v vertices and c components, target
placements are O(n^c), source embeddings O(n^c), and the containment
denominator has order n^v. A prescribed repeated cell removes one target
choice, giving O(n^(2c-v-1)). Unless all components are isolated edges,
this is O(n^-2). Without a duplicate the only first-order nonmatching
pattern is one two-edge path with isolated edges. Replacing its R+C count
by R+C-2Z changes only the remainder. These constants depend on a,b,
not on G or the lengths; zero or bounded lengths cause no division issue.

Opposite-direction line unions can contain a 4-cycle. For example,
`{(0,0),(1,1),(2,2)}` and `{(0,2),(1,1),(2,0)}` contain both a 4-cycle and
an isolated common cell. Such patterns obey the above remainder bound;
the parallel-line absence of cycles must not be inherited. The producer
explicitly retains them and is correct.

The resulting covariance is therefore

```
ab 2^(a+b-1)/n * [theta^a eta^b
 + theta^(a-1)eta^(b-1)(Z-(R+C)/n)] + O_(a,b)(n^-2).
```

As an independent normalization check, the exact a=b=1 formula is

```
E(Y_T Y_U) = 2Z/n + 2(R+C-2Z)/(n(n-1))
             + (LM-R-C+Z)(4n-6)/(n(n-1)^2).
```

The referee checks this on each cycle type at n=3,...,6 using two
antidiagonals, including both shared-cell parity possibilities.

For triple counts divide the a=b=3 formula by 36. Summing O(n^2) line
pairs leaves an O(1) error uniformly in the carrier. The three limiting
geometric terms are independently recovered as follows:

* Length cubes give `(integral_(-1)^1 (1-|t|)^3 dt)^2=1/4`.
* Splitting length quadrants gives row/column overlap
  `4*(1/14+71/1260)=23/45`.
* Direct integration over the triangle `0<=x<=1/2, x<=y<=1-x`, with
  integrand `4(1-y+x)^2(x+y)^2`, gives the common-cell term `19/90`.

The crossing term uses actual grid cells, not a smooth crossing indicator
with an assumed parity density. Consequently
`8*(1/4-23/45+19/90)=-2/5`, each inherited variance is `40n/9+O(1)`,
and the total variance is `364n/45+O(1)`. The exact total mean is
`(4n-10)/3`. The squared mean correction in the L2 statement is of order
n^-2. No triple has both slopes, so S_++S_-<=X without overcounting.
Chebyshev yields the uniform limsup constant `(364/45)/(16/9)=91/20`.
Arbitrary mixtures of cycle types have the same mean, so no between-type
mean variance is missing.

## Additional exact geometry certificate

Set L_d=n-|d| and M_e=n-|e-(n-1)|. Write

```
A_n = sum_d L_d^3,
O_n = sum_(d,e) L_d^2 M_e^2 (R_de+C_de),
Z_n = sum_(r,c) (n-|c-r|)^2 (n-|c+r-(n-1)|)^2.
```

The independent rational arithmetic certificate gives, for n>=1,

```
A_n = n^2(n^2+1)/2,
O_n = (46n^7+85n^5+49n^3)/90,
Z_n = (19n^6+25n^4-44n^2)/90 + (n mod 2)n^2.
```

Here `R+C=min(L,M)+max(L+M-n,0)`. Splitting those two summands at
L=M and L+M=n produces polynomial summation regions with integer affine
boundaries, so O_n is a polynomial of degree at most seven. For Z_n,
split the square along r=c and r+c=n-1. The only moving half-integer
boundary has period two; polynomial weight degree four and two sums give
a quasipolynomial of degree at most six on each parity. These degree bounds
justify exact rational interpolation. Direct sums on n=1,...,50 check
both the line-pair parity definition and the separate grid-cell definition,
including held-out values beyond the interpolation nodes.

Thus the normalized leading geometric kernel is exactly

```
8*(A_n^2/n^8-O_n/n^7+Z_n/n^6)
 = -2/5 - 4/(3n^2) + (8*(n mod 2)-94/15)/n^4.
```

This identity concerns the leading geometric kernel; it does not improve
the actual covariance's O(1) remainder, which still contains carrier data.

## Finite controls and failure boundaries

The independent census contains all 70,086 distinct simple degree-two
boards for n=3,...,6: respectively 6,90,2040,67950. It computes both
diagonal statistics for every board and literal integer determinants for
n<=5. Full all-slope counts at n=6 are deliberately outside this census.

| n | Cycle type | Boards | Cross covariance | P(both zero) | P(full zero) |
|---|---|---:|---:|---:|---:|
| 3 | C6 | 6 | -1/9 | 1/3 | 1/3 |
| 4 | 2C4 | 18 | 1 | 1/2 | 1/2 |
| 4 | C8 | 72 | -2/3 | 1/36 | 1/36 |
| 5 | C4+C6 | 600 | -34/225 | 1/75 | 0 |
| 5 | C10 | 1440 | -203/180 | 7/120 | 1/45 |

All four n=6 cycle types are retained in the output. In particular 3C4
has positive covariance 109/135, another check against an illicit finite
negative-sign claim. The independent five-by-five hostile

```
{(0,0),(0,1),(1,2),(1,3),(2,0),(2,4),(3,1),(3,2),(4,3),(4,4)}
```

has both diagonal counts zero and exactly two other collinear triples.
Thus the one-sided zero-event implication cannot be reversed.

## Reproduction and frozen identities

```
python -B 04-computation/overnight9_20260906_no3line_two_diagonals_audit.py
python -B -O 04-computation/overnight9_20260906_no3line_two_diagonals_audit.py
```

Both runs pass **72,622 always-active gates**, with byte-identical LF output.

* Referee source SHA256:
  `2428e10f6df7b7dca45f043b97188796ead06246c55e2cff48b2a3140b698ceb`.
* Referee output SHA256:
  `c4d9ff3b16234dbd66e097cfbc5cbadcf5cfff01cd5c0c6e0ed063184e85f3ec`.
* Semantic SHA256:
  `05c03390a43ff6b872a5ce6862c75492ce409bea088124898d0cb7cf5a458dd9`.
* Producer source identity supplied for this review, not imported:
  `b9d6dba0402216ac055c37658a5d0c4ed85ad2975c948e862355fcb7fe7398a8`.
* Producer output identity:
  `a74650e096fcc5881d8750c4c7c55326339e7b90b2011e64d175688cd9817d32`.

The accepted all-size statements rest on the uniform pattern proof and
geometric identities; the finite census supplies independent adversarial
and normalization checks rather than an extrapolation argument.

**Filing:** root integrated these audited artifacts in the ninth checkpoint;
reproduction paths are relative to the repository root. Earlier outside-worktree
notes describe author provenance, not the present file location.
