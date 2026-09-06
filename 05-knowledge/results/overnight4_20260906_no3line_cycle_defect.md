# Cycle defects, signed squares, and the failure of component-additive cumulants

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The [independent referee](overnight4_20260906_no3line_cycle_defect_audit.md)
uses graph BFS and an Euler-log recurrence, checks the formal all-length
proof, and reproduces all retained geometric profiles and both obstructions.
Its 278,884 optimization-live checks pass normally and with optimization. The previous third-moment coefficients are inherited, not new.
The new mechanism explains their place in a connected-cluster expansion
before geometric averaging. It also supplies a negative homogeneous
weight-four square supported entirely on realizable weight-eight forest
types, and a nonzero component-additivity contrast for third cumulants.

## 1. Inheritance and the target of the reframe

The closest proved mechanism is the
[short-cycle copy theorem](overnight_20260906_moments_pairprofile_theorem.md).
The [third factorial-moment calculation](overnight2_20260906_no3line_third.md)
and [independent audit](overnight2_20260906_no3line_third_audit.md) supply the
complete n=8 event-union multiplicities. They prove

```text
M3(G)=E[(X3)_3]=A+B*c4+D*c6+E*c8+F*c4^2,
E=172483/529200,  F=11881/50400.
```

The canonical hostile is `C16` versus two `C8` components: the first two
moments agree, but the third does not. The corrected near miss is that a
cycle monomial permitted by the edge budget need not survive geometric
cancellation. The least-used sidecar is the multiplicative generating
function of complete non-induced edge-subset types, before assigning the
whole-event geometric weights.

The live board is: cycle-component counts; path-component variables;
connected clusters; whole-event multiplicity; global label denominators;
and ordinary versus factorial cumulants. The question is whether a positive
or component-additive cumulant interpretation survives the averaging map.

## 2. A universal formal identity before geometric averaging

Work over the rational ring completed by positive edge weight. Its variables
are `p_(2j+1)` for an odd-edge path, `p_(2j)^L,p_(2j)^R` for an even-edge
path with the indicated shore larger, and `z_L` for an L-cycle (`L>=4`
even). Each variable has weight equal to its number of edges. A monomial
records a disjoint union of components, with repetition retained. In
particular `p1` is a single edge; `p2^L` has two left and one right vertex.

For a simple bipartite 2-regular graph G with n vertices per shore, define

```text
Z_G = sum_(S subset E(G)) monomial(type(S)),
P_L = Z_(C_L).
```

The constant term is one. Its monomial coefficient is exactly the number
of non-induced edge subsets of that type. Disjoint skeleton components give

```text
Z_G=product_L P_L^(c_L),
log Z_G=sum_L c_L log P_L.                             (1)
```

There is a universal path-only formal series b such that, in each weight
m, `[log P_L]_m=L[b]_m` whenever `L>m`. Consequently, with
`Delta_L=log P_L-Lb`,

```text
Delta_L starts in weight L,
Z_G=exp(2n b)*exp(sum_L c_L Delta_L).                  (2)
```

**Proof of the locality assertion.** Regard each nonempty connected edge
subset of a cycle as a polymer, with its component-type variable as
activity. Two polymers are incompatible exactly when they share a vertex.
Compatible families are precisely the connected components of an edge
subset, so their partition function is P_L. The formal logarithm contains
only connected incompatibility clusters. This follows by expanding the
pairwise compatibility factors and applying the exponential formula to
connected components of the resulting incompatibility graph. Repetitions
are allowed in the logarithm and a polymer is incompatible with itself.

A cluster of total edge weight m has connected union containing at most m
distinct edges. If `L>m`, it cannot contain the full cycle; its union is a
proper path arc. It therefore lifts uniquely to the infinite alternating
path after marking its initial edge. Every local pattern has `L/2` possible
placements for each shore orientation. Its coefficient is thus L times a
quantity independent of L. This proves locality coefficient by coefficient;
no convergence or positivity assumption is needed. Equation (2) follows
from (1) and `sum L*c_L=2n`.

The formal connected-cluster identity is the classical Mayer expansion;
a checked primary source is Fernández--Procacci, *Cluster Expansion for
Abstract Polymer Models*, equations (2.3)--(2.4),
[published paper](https://webspace.science.uu.nl/~ferna107/papers/ferpro07.pdf),
[DOI 10.1007/s00220-007-0279-2](https://doi.org/10.1007/s00220-007-0279-2).
Only its formal connectedness identity is used. No analytic convergence
criterion or sign conclusion from that paper is imported into this model.

## 3. Explicit four-cycle defect and its next correction

Let `D=[Delta_4]_4` and write `u=p2^L`, `v=p2^R`,
`r=p4^L`, `s=p4^R`. Direct cyclic-run enumeration, or expansion of (1), gives

```text
D = z4+p1^4-2p1^2(u+v)+4p1*p3+u^2+v^2-2(r+s),
T = p1^2*p3-p1*u*v+p1*z4-p1(r+s)+(u+v)*p3-p5,
[Delta_4]_5 = -8p1*D+4T.                             (3)
```

All monomials in D have weight four; those in T have weight five. Since
`[b]_1=p1`, equation (2) determines the full quadratic four-cycle
coefficient through weight nine for every fixed n:

```text
[c4^2] Z_G = D^2/2 + (n-8)p1*D^2 + 4D*T  modulo weight>=10. (4)
```

The first term has weight eight and the other two have weight nine. In
particular at n=8 the middle term vanishes. This is a coefficient of the
universal cycle-count expression, not an assertion that a logarithm of a
geometrically averaged statistic has the same coefficients.

For the eight-cycle coefficient, let `d8=[Delta_8]_8` and
`d9=[Delta_8]_9`. The analogous complete expression is

```text
[c8] Z_G = d8+d9+2n*p1*d8 modulo weight>=10.           (5)
```

The verifier computes d8 and d9 directly from `log P8-8b` (45 and 74
monomials respectively). Both (4) and (5) agree as whole polynomials with
the old skeleton contrasts, including types absent from the geometric
census. Thus this is an exact structural expression for both coefficients,
not an extra fitted census.

## 4. The averaging map and its signed-square obstruction

At n=8 let `W8(H)` be the number of unordered sets of three distinct nonaxis
collinear grid triples with complete union type H, divided by

```text
N8(H)=(8)_(left vertices)*(8)_(right vertices)/aut_shore(H).
```

Extend W8 linearly to the component-type ring, assigning zero to impossible
types. This normalization retains unordered event multiplicity, non-induced
copies, and every shore-preserving automorphism. The actual moment is
`M3=6 W8(Z_G)`; the factor six restores all ordered event triples. For
(4)--(5), only the inherited complete eight- and nine-edge bank is needed.

Applying this map to (4) recovers

```text
F = 3 W8(D^2)+24 W8(D*T)
  = 456371/2116800 + 42631/2116800
  = 11881/50400.                                    (6)
```

This looks like a quadratic expression. Positivity on the entire unital
ring is already trivially impossible because `W8(1)=0` while W8 is nonzero.
The useful question is narrower: can W8 be positive on homogeneous
weight-four squares in the retained weight-eight sector, where `D^2`
lives? It cannot. Take the homogeneous `Q=u*(p1^2-u)`. Every term of its
square is a realizable weight-eight forest type. Three inherited counts give:

| Monomial | Event triples | Grid copies | W8 value |
|---|---:|---:|---:|
| `p1^4*u^2` | 66988 | 4233600 | `16747/1058400` |
| `p1^2*u^3` | 32988 | 2822400 | `2749/235200` |
| `u^4` | 1200 | 176400 | `1/147` |

Therefore

```text
W8(Q^2)=16747/1058400-2*(2749/235200)+1/147
        =-397/529200<0.                              (7)
```

Every individual count and averaging weight here is nonnegative. Q is a
formal combination of union types, not a random variable whose square is
being averaged. Multiplication in this ring means disjoint union of
component types; W8 counts whole geometric event triples and divides by
global injective-label denominators. It is not an algebra homomorphism or
a positive moment functional on this formal multiplication. Thus (7) does
not assert a negative probability or expected square of an actual variable.

The positive value in (6) requires its actual arithmetic geometry. Merely
calling D a cycle defect and observing D squared cannot prove positivity.
The same loss explains why acyclic union types contribute to short-cycle
coefficients: their embedding-exclusion corrections detect ambient cycles.

## 5. A precise obstruction to component-additive scalar cumulants

Write `mu=E[X3]`, `M2=E[(X3)_2]`, and `M3=E[(X3)_3]`. The third factorial
and ordinary cumulants are respectively

```text
K3=M3-3mu*M2+2mu^3,
kappa3=M3+3(1-mu)M2+mu-3mu^2+2mu^3.                  (8)
```

At fixed n the mean is cycle independent and M2 is affine in c4,c6.
Hence neither subtraction can change the c8 or c4-squared coefficient.
Both cumulants inherit E and F above. A fixed-n coefficient alone can be
misleading if the cycle basis is redundant; here the five basis functions
`1,c4,c6,c8,c4^2` have rank five on valid n=8 skeletons (a five-row minor
using `C16,C4+C12,C6+C10,2C8,4C4` has determinant 24).

There is also a direct contrast annihilating *every* component-additive
function, without any chosen polynomial basis. The cycle multisets satisfy
`(2C8) disjoint-union (4C4) = 2*(2C4+C8)`. Thus an additive function of
individual cycle components would have zero contrast. Instead, for M3
and for either cumulant in (8), the exact contrast is

```text
stat(2C8)-2stat(2C4+C8)+stat(4C4)
    =8F=11881/6300>0.                                (9)
```

All three skeletons have n=8. This refutes a decomposition of these scalar
cumulants into a sum of contributions from individual skeleton cycles,
even allowing arbitrary n-dependent contributions for each cycle length.
It does not refute a connected expansion whose clusters may involve more
than one skeleton component. Global label exclusion and line events
spanning components are precisely the information such an expansion needs.

## 6. Exact verification and scope

[Standalone verifier](../../04-computation/overnight4_20260906_no3line_cycle_defect.py)
and [frozen output](overnight4_20260906_no3line_cycle_defect.out) construct
all cyclic edge masks for lengths 4,6,8,10,12,14,16, decompose proper masks
into runs, and perform exact rational formal logarithms through weight nine.
The independent long-cycle controls agree coefficientwise. The finite
checks support the locality proof; they are not its all-length justification.

Every one of the 150 retained geometric profiles agrees with the old copy
coefficient and its weighted contribution. No third-event geometry census
is repeated or enlarged. The inherited
[certificate](overnight2_20260906_no3line_third_certificates.json) has SHA-256
`f2c566ac1b2bcedb530af72f3a290479841e7db87844b331558bed68c93ba727` and
was already independently audited. Collapsing every component variable to
the corresponding power of a single edge variable kills all cycle defects,
as it must because `Z_G` then becomes `(1+t)^(2n)`.

```text
python -B 04-computation/overnight4_20260906_no3line_cycle_defect.py
python -B -O 04-computation/overnight4_20260906_no3line_cycle_defect.py
```

All **71,420** optimization-live gates pass, and normal and optimized
outputs are byte-identical. SHA-256 values refer to LF bytes:

```text
source c57585de9a836abbc40387d8c7fc23c300921ba3707daa3b16d07f16c5427c4d
output c0682d15f8ec8bb488f86f6827e7ef9ba5e8a69d8112441d5cce98171822c14c
formal semantic 50a0f39da9447a60a63fb9397ab8c717ec2dc90063fa3a94baa01cf0c7dc1d9d
```

The general formal identity and the bounded counterexamples are separate
claims. No all-size sign theorem,
zero-defect probability estimate, convergence radius, or independence of
geometric line events is inferred.
