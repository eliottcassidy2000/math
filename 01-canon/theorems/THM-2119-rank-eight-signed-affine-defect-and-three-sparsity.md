---
id: THM-2119
title: "Signed affine defect and three-sparsity gates for rank-eight torus covers"
status: >
  PROVED from settled LRC below the target rank and exact one-circle measure.
  In THM-2098's pure-transverse rank-eight coefficient plane, every complete
  signed affine-line family has a strict mixed-safe point, including a line
  parallel to the guard and a line through the origin. More generally, a
  guard-transverse signed affine core of at least five terminals escapes when
  every outlier is transverse to its direction; a guard-parallel core of at
  least four escapes; and a maximal rank-one core of at least four escapes.
  Consequently every projective terminal direction and every absolute
  quotient coefficient has multiplicity at most three in a cover. Combined
  with THM-2104, the quotient multiset is three-sparse and crosses every
  2-, 3-, and 5-adic valuation wall. This does not close rank eight or LRC(14).
source: codex-2026-07-22-LRC-signed-affine-defect
depends_on:
  - THM-2098
  - THM-2099
related:
  - THM-2103
  - THM-2104
  - THM-2105
  - THM-2115
  - THM-2117
  - THM-2120
---

# THM-2119 -- signed affine defect and three-sparsity

## 1. Setup and affine cores

Retain THM-2098's pure-transverse rank-eight plane. Thus

```text
g,c_1,...,c_8 in Z^2,                                  (1)
```

each `c_i` is rationally independent of the guard `g`, and some integer
direction `d` makes `g.d` odd positive and the eight `c_i.d` distinct
positive integers. The putative terminal bands cover the guard-safe region:

```text
{X:||g.X||>1/7}
  subset union_i {X:||c_i.X||<1/14}       almost everywhere. (2)
```

For a sign vector `sigma_i in {+1,-1}`, a subset `I` is a **signed affine
core** if there are characters `p,q`, with `p` primitive and nonzero, and
integers `k_i` such that

```text
sigma_i c_i=q+k_i p                         for i in I. (3)
```

There is no loss in the integral normalization: if all signed differences
lie in one rational line, its primitive lattice generator is `p`, and
`Qp intersect Z^2=Zp`.

## 2. Every complete signed affine line escapes

Assume first that `I={1,...,8}`.

### Nondegenerate line, direction transverse to the guard

Suppose `p,q` are independent and `p,g` are independent. Since `p` is
primitive, the fiber

```text
K={X:p.X=0}                                             (4)
```

is a connected circle. The restrictions of `q` and `g` are nonzero integer
characters on `K`, hence Haar-measure preserving. On `K`, all terminal
distances in (3) equal `||q.X||`. The closed bad sets have measures

```text
measure_K{||g.X||<=1/7}=2/7,
measure_K{||q.X||<=1/14}=1/7.                          (5)
```

Their union has measure at most `3/7<1`. A point outside it is strictly safe
for the guard and every terminal.

### Nondegenerate line parallel to the guard

Suppose `p,q` are independent but `p` is rationally proportional to `g`.
Primitivity gives

```text
g=m p,                         m in Z\{0}.              (6)
```

Because `g.d=m(p.d)` is odd, `m` is odd. The full-rank torus map
`X |->(p.X,q.X)` is surjective. Choose

```text
p.X=1/2,                       q.X=1/4.                 (7)
```

Then

```text
||g.X||=1/2,
||c_i.X||=||sigma_i c_i.X||
          =||1/4+k_i/2||=1/4.                         (8)
```

Thus the guard-parallel exception left open after THM-2099 equation (21) is
also always safe. The odd specialized guard is the load-bearing parity input.

### Line through the origin

Finally suppose `q in Qp`. Then every `c_i` is an integral multiple of the
primitive character `p`. Orient `p` so that

```text
c_i=a_i p,                 a_i distinct positive integers. (9)
```

Positivity and distinctness follow from the specialization direction `d`.
Transversality makes `p,g` independent. Pad the eight integers `a_i` by four
new positive integers and apply settled LRC for twelve nonzero speeds. There
is a `t` with

```text
||a_i t||>=1/13>1/14                 for every i.      (10)
```

Surjectivity of `(p,g)` realizes `p.X=t`, `g.X=1/2`, proving a strict mixed
escape. A collapsed direction `p=0` is incompatible with eight distinct
positive specializations.

The three cases prove unconditionally that a rank-eight cover must satisfy

```text
rank_Q{sigma_i c_i-sigma_1 c_1:2<=i<=8}=2
                              for every sign vector sigma. (11)
```

## 3. Exact near-pencil defect bounds

The same fiber argument tolerates outliers.

### Guard-transverse affine core

Let `p,q` and `p,g` be independent, and suppose every outlier `c_j`, `j notin
I`, is also independent of `p`. Restrict to (4). The guard costs `2/7`, the
entire core costs one common `1/7` arc, and each of the

```text
r=8-|I|                                                  (12)
```

outliers costs `1/7`. Therefore the closed bad union has measure at most

```text
(3+r)/7.                                                (13)
```

If `r<=3`, equivalently `|I|>=5`, (13) is strictly below one and a mixed-safe
point exists.

### Guard-parallel affine core

If `p,q` are independent and `g=mp`, then `m` is odd as above. Work on the
fiber `p.X=1/2`. The guard is safe. The core phases

```text
q.X+k_i/2
```

occupy at most two parity arcs, of total measure at most `2/7`. Every outlier
is transverse to `p` because the pure-transverse setup has `g=mp`, and costs
`1/7`.
The bad mass is at most

```text
(2+r)/7<1                         whenever r<=4.        (14)
```

Hence every maximal guard-parallel affine core of size at least four escapes.

### Maximal rank-one core

Let a maximal projective class be `c_i=a_i p`, `i in I`. Settled LRC at the
at most eight core speeds gives a `t` at which all core distances are strictly
greater than `1/14`. On the connected fiber `p.X=t`, the guard and each
outlier are nonzero characters. They cost at most

```text
2/7+r/7<1                         whenever r<=4.        (15)
```

Thus every maximal rank-one core of size at least four also escapes.

All inequalities use closed bad arcs. A strict total below one leaves positive
Haar measure outside their union, so equality endpoints cannot spoil the
strict conclusion.

## 4. Three-sparsity and valuation-wall consequences

The rank-one result immediately gives

```text
every projective terminal direction has multiplicity <=3. (16)
```

There is a second, basis-aware consequence. Write a possibly imprimitive odd
guard and transverse terminals as

```text
g=m p,                    m odd, p primitive,
c_i=a_i p+n_i q,          n_i!=0, p,q independent.     (17)
```

If four terminals have the same absolute quotient coefficient `|n_i|=N`,
choose `sigma_i=sign(n_i)`. Then

```text
sigma_i c_i=Nq+(sigma_i a_i)p                           (18)
```

is a guard-parallel signed affine core of size four. The other four terminals
are allowed to be arbitrary outliers, so (14) excludes a cover. Hence

```text
#{i:|n_i|=N}<=3                         for every N>0.  (19)
```

THM-2104 independently says a cover must occupy at least two `ell`-adic
valuation shells for each `ell in {2,3,5}`. Together, (19) and THM-2104 say
the quotient multiset is exact-value three-sparse and crosses all three
small-prime walls. THM-2105's affine clocks, THM-2115's Fourier--Toeplitz
moments, and THM-2117's separation gate retain further residue information;
none is implied by multiplicity alone.

As a direct compatibility check with the incoming frontier, both ledgers in
THM-2120 equation (17) violate (16): the all-line case has eight projectively
parallel terminals, while the exceptional case has the blocker and six
residuals on the same projective line. Thus once THM-2120 derives its global
phase-rigidity dichotomy, three-sparsity supplies a shorter final contradiction.
THM-2120 remains load-bearing for deriving that dichotomy.

## 5. Circuit-potential sidecar

The affine condition has a useful exact relation coordinate. For an
independent terminal pair, write its signed primitive circuit as

```text
A_ij g=B_ij c_i+C_ij c_j,                 A_ij>0.       (20)
```

Fix a sign gauge. Along any connected spanning tree of such pairs, the signed
columns lie on one affine line not through zero if and only if there is one
rational potential `R` satisfying

```text
A_ij R=B_ij sigma_i+C_ij sigma_j          on every edge. (21)
```

For the forward direction, choose an affine functional `ell` with
`ell(sigma_i c_i)=1` and set `R=ell(g)`. Conversely, prescribe
`ell(g)=R`, `ell(c_root)=sigma_root`; the circuit equations propagate
`ell(c_i)=sigma_i` along the tree. Non-tree edges test the cycle holonomy.
The value `R=0` is exactly the guard-parallel affine line; `R!=0` is the
guard-transverse case. A through-origin rank-one class is the separate
degenerate branch of Section 2.

This turns affine defect into incompatible edge potentials rather than an
unlabelled maximum-tree statistic. It is the lawful sidecar suggested by
THM-2065 and THM-1226 after THM-2103 refuted the original tree-or-pencil
dichotomy.

## 6. Assumption challenge and Tournament Analysis

The pairwise observable is the circuit potential in (21), and the switch is
the terminal sign gauge. Equality of edge potentials, not an orientation, is
the predicate that detects an affine core. Orienting edges by potential and
breaking ties by labels produces a tournament, but its score histogram,
cycles, SCCs, edge flips, and Hamiltonian paths depend on that arbitrary
tie-break and do not preserve (21). The faithful object is the edge-labelled
graphic matroid with cycle holonomy.

Candidate vertices challenged here were terminal columns, projective
directions, absolute quotient fibers, valuation walls, and proof obligations.
Projective quotienting preserves (16) but destroys lift parity; absolute
quotienting preserves (19) but destroys the sign and affine offset. Those
losses explain why the clock and Toeplitz sidecars remain necessary.

This theorem removes complete signed affine lines and the stated near-pencils,
and gives two unconditional multiplicity laws. It does not remove a guard-
transverse near-pencil with an outlier parallel to its direction, classify the
remaining coefficient planes, force a Toeplitz violation, discharge the
vertical branch, or prove LRC(14). QED.

