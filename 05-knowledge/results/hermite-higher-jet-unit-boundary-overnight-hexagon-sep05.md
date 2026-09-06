# Higher jets: exact inverse precision and the first unit-sensitive boundary

**PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.** The full
written proof passed an independent review; literal inverse/Smith controls
and both normal/optimized pairs passed. The canon route is
[THM-4443](../../01-canon/theorems/THM-4443-arbitrary-jet-precision-and-dyadic-unit-boundary.md).
No all-order metric-only assertion or literature-priority claim is made.

## 1. Inheritance and the operation that found the boundary

[THM-4439, terminal-cluster two-jet precision](../../01-canon/theorems/THM-4439-all-node-twojet-metric-precision-by-terminal-clusters.md)
proves that the worst precision loss is metric-only when every node supplies
exactly two Hasse jets. Its closest mechanism is the inverse denominator
from [THM-4435, universal Hermite precision](../../01-canon/theorems/THM-4435-four-node-metric-blindness-and-universal-hermite-precision.md).
The canonical hostile is the isometric four-node full-Smith pair in that
theorem: intermediate factors differ, while the largest factor agrees.
The corrected near miss is to transfer the full-partition failure to the
two-jet precision, or conversely to extend the precision theorem to a
different observer. The least-used operation here is **change jet
multiplicity while retaining exactly the same node metric**.

[THM-4010, confluent Hasse observer and Smith firewall](../../01-canon/theorems/THM-4010-confluent-consecutive-hasse-observer-kernel-index-and-smith-firewall.md)
already distinguishes Hasse from ordinary derivatives and the determinant
from the Smith partition. Its consecutive-node higher-jet observer does
not settle arbitrary node sets or heterogeneous multiplicities.

The live concepts are: observer multiplicity; inverse Taylor coefficients;
terminal distance clusters; simultaneous reciprocal cancellation; unit
residue; sharp loss. The map from a multiplicity-labelled node set to its
labelled valuation tree preserves neither the general inverse denominator
nor the examples below. The repaired sidecar is explicit: a reciprocal
Taylor jet in general, and a four-valued unit class in the first uniform
dyadic three-jet case.

## 2. Exact sharp precision for arbitrary multiplicities

Let x_1,...,x_n be distinct integers, let m_i>=1, and put M=sum_i m_i.
Observe D^[r]P(x_i) for 0<=r<m_i on integer polynomials of degree<M.
Define the monic integer polynomials and rational local coefficients

```text
Q_i(X)=product_(j!=i)(X-x_j)^(m_j),
a_(i,l)=[T^l] Q_i(x_i+T)^(-1),       0<=l<m_i.           (1)
```

The largest integer Smith factor is exactly the least common multiple of
the reduced denominators of **all** a_(i,l). Equivalently, for any prime p,

```text
L_p=max_(i,l<m_i) {-v_p(a_(i,l))},  v_p(0)=infinity.    (2)
```

The maximum is nonnegative because the l=0 terms are reciprocals of
nonzero integers. For every N>=1, data modulo p^(N+L_p) determine all
coefficients modulo p^N, and using one less digit fails uniformly when
L_p>0. When L_p=0 no extra digits are necessary.

For proof, the cardinal polynomial for the (i,r) datum is

```text
Phi_(i,r)(X)=(X-x_i)^r Q_i(X)
             *sum_(l=0)^(m_i-r-1) a_(i,l)(X-x_i)^l.    (3)
```

It has degree at most M-1, its prescribed Hasse jet is one, and every
other observed jet is zero. Thus these are precisely the inverse matrix
columns. Formula (3) bounds their denominators by those in (1). Conversely
the coefficient of X^(M-1) in Phi_(i,r) is exactly a_(i,m_i-r-1): every
factor is monic and all other summands have smaller degree. Every proposed
denominator is therefore attained as an actual inverse entry. This proves
(2) and the integer least-common-multiple statement. The sharp precision
claim also follows directly by integral unimodular Smith changes on source
and target. A coordinate with exponent L_p gives the one-digit hostile.

The independent referee pointed out this direct leading-coefficient
attainment; triangular recovery from a single cardinal column is unnecessary.
For one node, all inverse columns are integral translated monomials, so L=0.

For uniform multiplicity k, set d_i=product_(j!=i)(x_i-x_j),
S_i=sum_(j!=i)v_p(x_i-x_j), and let h_l be the complete homogeneous
symmetric function of the multiset containing each reciprocal
1/(x_i-x_j) exactly k times. Then

```text
a_(i,l)=(-1)^l d_i^(-k) h_l,
L_p=max_(i,0<=l<k) {k S_i-v_p(h_l)}.                  (4)
```

Under dilation x_i -> p^e x_i, this particular candidate increases by
[k(n-1)+l]e. The full precision is the upper envelope, not a chosen
largest-degree coefficient. For k=1 this is max_i S_i, hence metric-only;
k=2 recovers THM-4435 and THM-4439. For two nodes at distance d and any
uniform k, (4) becomes the metric-only formula

```text
L_p=max_(0<=l<k) [(k+l)v_p(d)-v_p binom(k+l-1,l)].     (5)
```

## 3. Complete dyadic three-node three-jet precision law

Now n=k=3 and p=2. Every three-node dyadic metric has pair depths
(e,e,e+d), with e>=0 and d>=1. Let x,y be the closest pair and z the
remaining node. If d>=2, then

```text
L_2=8e+5d-1.                                        (6)
```

If d=1, the unit

```text
t=2(z-x)/(y-x) in Z_2^*
```

is well-defined modulo16. Set

```text
Gamma(t)=3  if t=3 mod4;
         2  if t=5 mod8;
         1  if t=1 mod16;
         0  if t=9 mod16.
```

These disjoint classes exhaust all odd residues. The exact law is

```text
L_2=max(7e+4, 8e+Gamma(t)).                           (7)
```

The residue class is invariant under unit-affine changes of coordinates.
Exchanging x and y replaces t by 2-t and preserves Gamma. Thus the
four-valued sidecar is intrinsic to this metric shape, not an ordering
artifact. It is sufficient for the loss at every depth, and all four
eventual branches occur, at the triples 2^e*(0,2,t) with t=3,5,1,9.
At small e the maximum can hide a branch difference; no claim of necessity
at a fixed shallow scale is made.

### Proof of the unit law

At a node with differences u,v to the other nodes, (1) gives

```text
a_0=1/(uv)^3,
a_1=-3(u+v)/(uv)^4,
a_2=3(2u^2+3uv+2v^2)/(uv)^5.                        (8)
```

At a closest-pair node write u=2^e U, v=2^(e+d)V with U,V odd
(rational two-adic units suffice after normalization). The l=0,1
candidates are 6e+3d and 7e+4d. If d>=2, the three terms in the
quadratic numerator have valuations 2e+1, 2e+d, 2e+2d+1. The first
is uniquely minimal, so the l=2 candidate is 8e+5d-1. It exceeds
the first two. At the outsider all differences have valuation e;
its three candidates are at most 6e,7e,8e, so it cannot dominate.
This proves (6).

For d=1, put t=U/V at one closest-pair node. This is the same t as
defined above (the signs in the two ratios cancel). The two relevant
quadratic numerators after removing 2^(2e+1)V^2 are

```text
P(t)=t^2+3t+4,
P(2-t)=t^2-7t+14.                                   (9)
```

Let nu=min(v_2(P(t)),v_2(P(2-t))). Direct residue arithmetic gives
nu=1,2,3,4 in the respective four classes displayed above. For detail:
if t=3 mod4 then P(t)=2 mod4. Otherwise write t=1+4a. The two
values are 4(2+5a+4a^2) and 4(2-5a+4a^2). Odd a gives nu=2.
For a=2b they become 8(1+5b+8b^2),8(1-5b+8b^2). Even b gives
nu=3; odd b makes both inner values even, while their difference
10b has valuation one, so their minimum valuation is exactly one
and nu=4. This argument works for arbitrary two-adic units t, not
just integer examples. Consequently the maximal closest-pair l=2
candidate is 8e+4-nu=8e+Gamma(t).

The outsider's l=2 candidate is exactly 8e: its normalized differences
are odd, making 2u^2+3uv+2v^2 odd after the common square is removed.
Since Gamma>=0 it does not exceed the pair maximum. The l=1 pair
candidate 7e+4 exceeds every l=0 candidate and every outsider l=1
candidate. Formula (7) follows.

## 4. Minimal uniform boundary and the nonuniform warning

Consider the isometric node sets

```text
A_e=2^e*(0,1,2),    B_e=2^e*(0,1,3), e>=0.
```

An explicit labelled isometry sends A's ordered points (0,1,2) to
B's (1,0,3). For three jets their losses are

```text
L(A_e)=max(7e+4,8e+1),
L(B_e)=max(7e+4,8e+3).                               (10)
```

They agree for e=0,1 and differ at every e>=2. At e=2 the full dyadic
Smith exponent lists, independently computed from literal 9x9 integer
matrices, are

```text
A_2: (0,0,0,2,6,8,11,18,18),
B_2: (0,0,0,2,6,8,11,17,19).
```

Thus k=3 and n=3 are the smallest uniform jet multiplicity and node
count allowing a metric-only precision failure: k=1 and k=2 are metric-only
for every node count, while (5) excludes every two-node uniform example.
This is not a global smallest-diameter assertion.

Uniformity is itself load-bearing. Give the closest two nodes multiplicity
two, and the outsider multiplicity one. For the ordered weighted sets

```text
A'_e=2^e*(0,2,1),  B'_e=2^e*(1,3,0),  m=(2,2,1),
```

the weighted metrics agree, but (1) gives

```text
L(A'_e)=max(3e+2,4e+1),
L(B'_e)=max(3e+2,4e).                                (11)
```

For A', the closest-pair reciprocal weighted sum is plus/minus2 before
dilation; for B' the two sums are 0 and4/3. Their denominator candidates
are therefore 4e+1 for A' and at most4e for B'. The order-zero pair
candidate is3e+2, and the outsider's order-zero candidate is4e, proving
(11). At e=2 the losses are9 and8. This does not contradict the complete
uniform two-jet theorem: one derivative was deliberately not observed.

The first failed implication in both hostiles is that a whole terminal
cluster has the same one-digit simultaneous cancellation budget after the
observer is changed. The surviving statement is the exact inverse law
(1)-(4), with the observer multiplicity retained. The uniform three-jet
dyadic law (6)-(7) identifies precisely which first unit class repairs the
first failure. Higher node counts, other primes and full partitions are
not inferred from this boundary calculation.

## 5. Reproduction and audit record

The initial exploratory script
[higher-jet probe](../../04-computation/hermite_higher_jet_probe_overnight_hexagon_sep05.py)
compares all nonnegative integer dilation envelopes, not only the undilated
loss. Its stated head found the three-jet dyadic hostile immediately;
passing heads at other multiplicities are not unbounded statements.
The [primary verification](../../04-computation/hermite_higher_jet_precision_overnight_hexagon_sep05.py)
and [frozen output](hermite_higher_jet_precision_overnight_hexagon_sep05.out)
cover every zero-translated three-node set of diameter<=80 (3,160 rows),
all four unit classes at17 scales,1,000 signed/unit-affine pairs,49 literal
full inverse/Smith matrices and112 two-node controls:5,947 explicit gates.

The [independent referee](../../04-computation/hermite_higher_jet_referee_overnight_hexagon_sep05.py)
imports no primary routine. It separately constructs126 literal integer
Hasse matrices and reciprocal series, comparing their **largest integer
factors**, not just one prime valuation. Its1,262 gates include89 named
dyadic matrices, five singleton multiplicities and32 seeded signed mixed
observers. The [stored output](hermite_higher_jet_referee_overnight_hexagon_sep05.out)
also retains the complete9x9 and5x5 hostile partitions at six scales.
Independent full written proof review: PASS, including denominator
attainment, all four residue classes, endpoint exchange, shallow masking,
minimal-uniform scope and the nonuniform warning. The root replayed the
referee. Primary and referee normal/optimized outputs each agree exactly.

```bash
python3 04-computation/hermite_higher_jet_precision_overnight_hexagon_sep05.py
python3 -O 04-computation/hermite_higher_jet_precision_overnight_hexagon_sep05.py
python3 04-computation/hermite_higher_jet_referee_overnight_hexagon_sep05.py
python3 -O 04-computation/hermite_higher_jet_referee_overnight_hexagon_sep05.py
```

Frozen raw-LF SHA256:

```text
primary 6263c2b9b5fc65036e7b5bd3bf2d751d529defb63635f07aa8beba9ee1a52a97
output 0fa167e65476acbe873812aae7e7a52dd48a9451cd49d9a7ec0d3d5475cf11b6
referee f14c3ac44f34cbbe275408895de781eb6363d7586fd03e4fffae6c86a095f697
referee output c0c03411b4c78f6101d660930eb2a8e5bc3c7467abc34864d913a675c3fed1e7
```
