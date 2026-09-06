# Crossings strengthen the random saturated-board diagonal obstruction

**Status: PROVED ANALYTICALLY / INDEPENDENTLY AUDITED.** The
[independent audit](overnight9_20260906_no3line_two_diagonals_audit.md)
accepts the proof and independently checks 70,086 distinct boards. The exact
controls below are not a substitute for the uniform all-size proof.
This is a random-density result, not an extremal nonexistence theorem.

## Inheritance and statement

The closest proved mechanism is the uniform colored-overlap lemma in the
[eighth diagonal-density report](overnight8_20260906_no3line_diagonal_density.md).
The inherited hostile is the finite C8 versus 2C4 difference: cycle types
retain different covariances and zero probabilities at finite size.
The corrected near miss is using the edge-disjoint covariance formula for
two intersecting lines. The least-used sidecar is a selected grid cell
appearing once in each of two factorial tuples. It contributes at the same
order as the row and column exclusion terms and cannot be dropped.

The concept board is: injective shore labels; factorial tuples; repeated
cells; row/column exclusions; orientation; and the one-sided zero-event
gate. The map sends a labeled saturated board to its two diagonal triple
counts. It preserves the implication that a positive count is a collinear
triple, loses all other slopes, and requires shared-cell and row/column
incidences as sidecars. A concrete five-by-five board below has both counts
zero but another collinear triple, preventing reversal of the implication.

Let G be any simple bipartite 2-regular graph on n+n vertices. Independently
label each shore uniformly with 0,...,n-1; B is its set of 2n grid cells.
Define

```
D_d={(r,c):c-r=d},       E_e={(r,c):r+c=e},
S_+=sum_d binom(|B intersect D_d|,3),
S_-=sum_e binom(|B intersect E_e|,3),       S=S_++S_-.
```

Let X count all nonaxis integer-collinear triples. No triple has both
distinct slopes, so `0<=S<=X` with no multiplicity factor.

**Theorem.** For every G and n>=3,

```
E_G S_+=E_G S_-=(2n-5)/3,       E_G S=(4n-10)/3.       (1)
```

Uniformly over G as n tends to infinity,

```
Cov_G(S_+,S_-) = -(2/5)n + O(1),                      (2)
Var_G S = (364/45)n + O(1),                            (3)
E_G[(S/n-4/3)^2] = 364/(45n)+O(n^-2),                 (4)
limsup_n n sup_G P_G(X=0) <= 91/20.                   (5)
```

The absolute remainders allow every cycle type, including a linear number
of C4 components. Every mixture of cycle types inherits the probability
bound; in particular it holds for the uniform distribution on distinct
saturated boards. The constant improves the eighth variance bound 10.
Neither constant is claimed optimal, and neither rules out rare successful
labelings. The asymptotic extremal no-three-in-line problem remains OPEN.

## The shared-cell covariance lemma

Take any two grid partial matchings T,U, with lengths L,M and at most one
common cell. Write

```
Y_T=|B intersect T|,  Y_U=|B intersect U|,
R=|rows(T) intersect rows(U)|, C=|cols(T) intersect cols(U)|,
Z=|T intersect U|<=1,
theta=L/n, eta=M/n, omega=(R+C)/n.
```

For fixed a,b>=1 and k=a+b, the strengthened uniform formula is

```
Cov((Y_T)_a,(Y_U)_b)
 = ab*2^(k-1)/n *
   [theta^a eta^b + theta^(a-1)eta^(b-1)(Z-omega)]
   + O_(a,b)(n^-2).                                  (6)
```

The exponent zero convention for empty matchings is the ordinary polynomial
one. Uniformity includes bounded or zero lengths; no division by L or M is
used in the proof.

Expand into ordered distinct cells within each color. If a cell occurs in
both tuples, it is isolated from every other selected cell in their
underlying incidence graph: T and U are each matchings, so no other cell
of either color can use its row or column. There can be only one duplicate.
The number of tuples with this duplicate and all other chosen cells
pairwise vertex-disjoint is

```
ab Z L^(a-1)M^(b-1)+O_(a,b)(n^(k-3)),                 (7)
```

for k>=3; at k=2 the count is exactly Z. Before requiring all the other
cells to be disjoint, the exact count is
`ab Z (L-1)_(a-1)(M-1)_(b-1)`. At Z=1 no other selected cell can meet the
duplicate; an additional intersection among the remaining cells costs
another free target choice.

There are k-1 distinct cells in (7). The inherited matching probability is

```
p_j=(2/n)^j[1+j(j-1)/(4n)+O_j(n^-2)].                 (8)
```

Multiplying (7) by p_(k-1) gives exactly the positive Z term in (6), with
O(n^-2) remainder. Removing these O(n^(k-2)) tuples from an all-k-matching
count affects its contribution by only O(n^-2).

For completeness, every more complicated pattern containing the duplicate
also lies in the remainder. The fixed repeated cell removes the n choices
for its target component. For an underlying graph with v vertices and c
components, there are O_k(n^(c-1)) target placements and O_k(n^c) injective
source embeddings, with containment denominator of order n^v. Its order
is at most n^(2c-v-1). Unless every component is an isolated edge, this is
O(n^-2). Patterns without the duplicate use the inherited n^(2c-v) bound.
These bounds are uniform over every carrier G.

Among tuples without a duplicate, the only other first-order pattern is
one two-edge path plus isolated edges. Replacing R+C by R+C-2Z changes its
leading target count only by O(n^(k-2)); thus it has the same first-order
term as the edge-disjoint lemma. Its count is
`ab(R+C)L^(a-1)M^(b-1)+O(n^(k-2))`, and its containment probability is
`2^(k-1)n^-k+O(n^(-k-1))`. All remaining patterns, including 4-cycles,
are O(n^-2). Opposite-direction line unions can contain such 4-cycles;
the parallel-line absence of cycles is not inherited here.

Combining these counts with (8) and subtracting the separate factorial
means proves (6). At a=b=1 there is also an exact diagnostic:

```
E(Y_T Y_U)= Z*(2/n) + (R+C-2Z)*2/[n(n-1)]
          +(LM-R-C+Z)*(4n-6)/[n(n-1)^2].              (9)
```

Every summand is a distinct geometric type: one cell, a two-edge path,
or two disjoint edges. The central diagonals at odd n give Z=1 and show
directly why the edge-disjoint formula has the wrong leading covariance.

## Sum over the two orientations

For a=b=3 divide (6) by 36. With L_d=n-|d| and M_e=n-|e-(n-1)|,

```
Cov(Z_d,Z_e)=(8/n)[theta_d^3 eta_e^3
                    +theta_d^2 eta_e^2(Z_de-omega_de)]
                    +O(n^-2).                       (10)
```

There are (2n-1)^2 pairs. Uniformity of the error gives O(1) after summing.
The three leading geometric sums divided by n^2 converge to the following
integrals, each with O(1/n) error.

First, the product of length-cube sums tends to

```
[integral_(-1)^1 (1-|u|)^3 du]^2 = 1/4.             (11)
```

For row and column overlaps use Fubini. The squared-length mass of the
slope-one lines through a normalized row r is

```
f(r)=integral_(-r)^(1-r) (1-|u|)^2 du
    =1/3+r-r^2.
```

Reflection gives the same function for the other direction and for the
column contribution. Therefore the overlap sum tends to

```
2 integral_0^1 f(r)^2 dr = 23/45.                    (12)
```

Finally, do not replace a crossing by a smooth indicator over line pairs:
half of geometric crossings may miss the integer grid because of parity.
Instead there is an exact, parity-preserving identity

```
sum_(d,e) L_d^2 M_e^2 Z_de
 = sum_(r,c in grid) (n-|c-r|)^2 (n-|c+r-(n-1)|)^2.  (13)
```

After division by n^6 this is a square Riemann sum, tending to

```
J=integral_[0,1]^2 (1-|c-r|)^2(1-|c+r-1|)^2 dr dc
 =2 integral_(u>=0,v>=0,u+v<=1) (1-u)^2(1-v)^2 du dv
 =19/90.                                            (14)
```

The factor 2 is four quadrants times the Jacobian 1/2 of
`u=c-r, v=c+r-1`. It already incorporates the exact cell parametrization,
so no additional parity factor is allowed. The piecewise polynomial
integrands in (11)-(14) are bounded Lipschitz functions on fixed compact
regions; ordinary mesh estimates give O(1/n) errors. The discrete overlap
sum can equivalently be written `2 sum_r A_r^2`, where
`A_r=sum_(d=-r)^(n-1-r) (n-|d|)^2`, and is divided by n^7. This also
establishes the claimed error without a line-pair boundary convention.

Thus the cross covariance coefficient is

```
8*(1/4-23/45+19/90)=-2/5.
```

The eighth theorem gives each single-orientation variance `40n/9+O(1)`;
column reflection transfers it exactly to the other orientation. Hence

```
Var(S_++S_-)=2*(40/9)n+2*(-2/5)n+O(1)
           =(364/45)n+O(1).
```

The exact means follow by the same reflection. Equation (4) includes the
squared mean correction 100/(9n^2). Since X=0 implies S=0, Chebyshev gives
`P(X=0)<=Var(S)/(E S)^2` for n>=3, proving (5) uniformly. Equal means
across G also ensure that arbitrary mixtures introduce no between-type
variance in this argument.

## Controls, failure boundary, and reproduction

The standalone verifier independently constructs literal shore labelings
for C6,C8,2C4,C10,C4+C6: 29,988 labelings are checked individually in the
transcript totals. For every board it computes all integer determinants
and checks S<=X. It verifies the exact means and both covariance terms,
the shared-cell two-edge formula, the symbolic first-order expansion for
all 1<=a,b<=3, and the three exact rational integrals. Exact line-pair
geometry through n15 is checked against row and cell Fubini sums, including
crossing parity; larger sums through n400 are finite controls only.

For C10 at n5 it finds `P(S=0)=7/120` while `P(X=0)=1/45`.
One actual board with S=0 but X>0 is

```
{(0,0),(0,2),(1,2),(1,4),(2,1),(2,4),(3,1),(3,3),(4,0),(4,3)}.
```

For C8 and 2C4 the cross covariances are respectively -2/3 and 1. Thus
the negative limiting coefficient is not a pointwise covariance sign
theorem or an assertion that finite cycle types are interchangeable.

The first failed implication in the discarded formula is 'all factorial
positions name distinct cells.' The strongest survivor is the inherited
edge-disjoint lemma. The repaired form is (6), retaining Z, and the new
question is which other fixed directions or bounded-change zero-event
statistics give stronger density bounds. No growing-direction or central
limit theorem follows merely from this two-direction calculation.

```
python -B 04-computation/overnight9_20260906_no3line_two_diagonals.py
python -B -O 04-computation/overnight9_20260906_no3line_two_diagonals.py
```

Both runs pass 34,607 always-active gates with identical LF output.
Source does not import repository code.

```
source b9d6dba0402216ac055c37658a5d0c4ed85ad2975c948e862355fcb7fe7398a8
output a74650e096fcc5881d8750c4c7c55326339e7b90b2011e64d175688cd9817d32
```

The independent verifier passes 72,622 gates with normal/optimized agreement.
It uses distinct boards through n6, instead of this producer's labeled carriers,
and independently derives the parity-sensitive finite geometry polynomials.
No correction to the analytic proof was needed.

**Filing:** root integrated these audited artifacts in the ninth checkpoint;
reproduction paths are relative to the repository root. Earlier outside-worktree
notes describe author provenance, not the present file location.
