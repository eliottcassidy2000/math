# Affine column orbits: an exact carry decoder and a pairwise-law obstruction

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The new results concern a specified operation on actual integer grid boards:
replace column c by ac+b modulo an odd prime p, then take the representative
in {0,...,p-1}. They give an exact decoder for successful affine images, a
minimal-prime obstruction to replacing a uniform column permutation by a
pairwise-uniform one, and a modular-line capacity certificate. They do not
settle no-three-in-line extremal asymptotics, furnish a general successful
construction, or assert independence of collinearity events.

[Independent proof and exact audit](continuing4_20260906_wildcard_affine_columns_audit.md) passes.

## 1. Inheritance and the concrete operation

The primary historical source remains
[Guy--Kelly, The No-Three-In-Line Problem, DOI10.4153/CMB-1968-062-3](https://doi.org/10.4153/CMB-1968-062-3).
The maintained source analysis
`05-knowledge/results/overnight_20260906_no3line.md` separates its grid-triple
count from the independence heuristic and records the later correction to
the heuristic entropy constant. Neither heuristic is used here.

The closest proved repository mechanisms are the exact conditional matching
law in `overnight11_20260906_no3line_rowfreeze.md` and its consumer
`overnight12_20260906_no3line_count_restart.md`, both under
`05-knowledge/results/`. Their random-column conclusions retain a conditionally
uniform full permutation. Purposeful changes within a board are explicitly
outside that restart theorem. The canonical hostile is the cycle-blind mean:
at n=4 the two possible cycle profiles have the same triple mean and different
zero probabilities. The present operation exposes an even smaller lost
coordinate: identical one- and two-column laws need not preserve success.

Anchor: find a faithful within-board operation. Niche: retain the integer
carry after a finite-field reduction. Wildcard: recover midpoint closure and
short lattice directions from a complete modular line. The live board is:

| Object | Preserved information | Information requiring a sidecar |
|---|---|---|
| Fixed-row degree-two board | Every row/column degree and incidence skeleton | Euclidean collinearity after a new column order |
| Affine column group | Every one-/two-column uniform cylinder | The affine ratio of three source columns |
| Modular triple | Determinant zero modulo p | Primitive row span and signed lifted slope |
| Forbidden-translation intervals | Exactly the actual bad affine images | Their union, not only total length |
| Full modular line | Midpoint closure for equal-parity endpoints | Euclidean line layers within that modular line |
| Tangent lattice | Congruence class and exact line direction | Integer width of the selected grid box |

The operation has a concrete failure and repair: the fixed 5-by-5 cycle board
below has no successful affine column image, while a single non-affine swap
followed by an affine map produces a successful board. This is a bounded
research result, not another density or random-greedy proposal.

A primary-source cross-check also read the proof-ideas and preliminary
geometric statements in Grebennikov--Kwan's author-hosted
[No-(k+1)-in-line problem for large constant k](https://mkwn.github.io/NPLC.pdf).
Its spread distributions and regular-subgraph extraction concern large k;
they are not invoked for k=2 here. The distinction reinforces why preserving
only pairwise columns supplies no transfer. No external novelty claim is made.

## 2. Exact affine images of a single triple

Let p be an odd prime, and write [p]={0,...,p-1}. For a nonzero a modulo p
and any b modulo p, let

```text
Phi_(a,b)(x,y)=(x, [ay+b]_p).
```

This is a column permutation. Let T consist of three points with distinct
rows and distinct columns, ordered by rows as (x0,y0),(x1,y1),(x2,y2).
Integer determinant zero is the target predicate throughout.

If det(T) is nonzero modulo p, every Phi_(a,b)(T) has nonzero integer
determinant: its determinant modulo p is a det(T). Now suppose det(T)=0
modulo p. There is a unique nonzero modular slope u such that
yi-y0=u(xi-x0) modulo p. Define

```text
g=gcd(x1-x0,x2-x0),
m=(x1-x0)/g, n=(x2-x0)/g, 0<m<n, gcd(m,n)=1,
K=floor((p-1)/n).
```

**Theorem 1 (exact event).** Fix a. The transformed triple is integer
collinear if and only if there is an integer k satisfying

```text
0<|k|<=K,             k=a u g modulo p,
max(0,-kn) <= [a y0+b]_p <= min(p-1,p-1-kn).             (1)
```

There is at most one such k for each a, since n>=2 implies 2K<=p-1.
When k exists, its forbidden b values form one cyclic interval of exactly
p-|k|n residues; otherwise this triple forbids no translation for that a.

Proof: any integer line through these rows has column differences km and
kn for an integer k. Indeed n times the first difference equals m times
the second, and gcd(m,n)=1. Distinct columns force k nonzero. Their total
column span is |k|n, which must be at most p-1. The congruence for either
difference gives k=a u g. Conversely that congruence and the last interval
in (1) force the unique representatives of all three column residues to
be [ay0+b]_p, [ay0+b]_p+km, [ay0+b]_p+kn. They are integer collinear.
This proves both directions and retains weak endpoints of the interval.

An immediate exact probability consequence, for uniformly chosen (a,b), is

```text
Pr(det Phi_(a,b)(T)=0)
 = 2/[p(p-1)] sum_(k=1)^K (p-kn)
 = [2Kp-nK(K+1)]/[p(p-1)].                            (2)
```

The expression depends on the primitive row span n, not merely on the
fact that the determinant vanishes modulo p. At p=5 the row triples
(0,1,3) and (0,2,4), with source columns equal to their rows, give respectively
1/5 and 2/5. The first triple under a=2,b=0 becomes
(0,0),(1,2),(3,1), whose integer determinant is -5, despite modular zero.

## 3. An iff decoder for the full affine orbit

Let B be any subset of [p]^2. If a row or column contains three points,
every affine image is already bad. Otherwise any potentially collinear
triple has distinct rows and columns. For each such modularly collinear
triple, enumerate the 2K signed integers in (1), calculate

```text
a = k (u g)^(-1) modulo p,
b in [max(0,-kn),min(p-1,p-1-kn)] - a y0 modulo p,      (3)
```

and add this cyclic interval to a forbidden-residue union U_a.

**Theorem 2 (full orbit).** Phi_(a,b)(B) has no three integer-collinear
points if and only if b is outside U_a. In particular

```text
number of successful affine parameter pairs
 = sum_(a=1)^(p-1) [p-|U_a|].                          (4)
```

This is the union of the exact triple events, so it requires no independence
or inclusion-exclusion truncation. Summing (2) over the board's triples
also gives its exact expected integer triple count, but that scalar alone
does not determine (4). The union sidecar is essential.

The supplied implementation uses p-bit masks for U_a. It generates only
the 2K possible signed slopes per modular triple, instead of trying every
multiplier. Its universe is the p(p-1) affine parameter pairs, not all p!
column permutations. Consequently an empty decoded orbit is not a proof
that its incidence skeleton admits no successful column ordering.

## 4. Pairwise-uniform columns can erase every success

The affine group is sharply two-transitive: prescribed images of two
distinct source columns determine exactly one (a,b). Consequently uniform
affine columns and a uniform member of S_p agree on every cylinder
prescribing the images of at most two columns. Equivalently, they agree on
expectations of all degree-at-most-two polynomials in board-cell indicators.
This includes the first two factorial occupancy moments used by the
conditional matching representation.

They disagree already at three columns. The affine ratio
(y2-y0)/(y1-y0) is preserved modulo p; a uniform unrestricted permutation
does not preserve it. A compatible prescribed triple has probability
1/[p(p-1)] under the affine law, versus 1/[p(p-1)(p-2)] under the full law;
an incompatible prescribed triple has affine probability zero.

The exact success hostile is the saturated p=5 board

```text
B={(i,i),(i,i+1 modulo 5): 0<=i<5},
```

whose incidence graph is one ten-cycle. Exhaustion of its 120 column
permutations gives four successful boards. All 20 affine maps fail, with
the following complete integer-triple histogram:

| Number of triples | Affine parameter pairs |
|---:|---:|
| 2 | 8 |
| 3 | 4 |
| 5 | 4 |
| 14 | 4 |

Thus the full-success probability changes from 1/30 to zero, and the mean
integer triple count changes from 64/15 to 26/5, although every one-/two-column
cylinder agrees. The six affine orbits on the full column-permutation
space, represented by fixing the first two images to 0 and 1, have success
counts 0,2,0,2,0,0 in lexicographic representative order. This checks the
full lost coordinate rather than a single bad image.

A minimal operation repairs the example. Swap source columns 3 and 4,
so tau=(0,1,2,4,3), then apply c->2c+1 modulo 5. The resulting column
permutation is (1,3,0,4,2), and its board is successful by both literal
determinants and repeated primitive-direction tests. A purely affine move
could never do so. At p=3, AGL(1,3)=S_3, so p=5 is the smallest odd prime
at which this affine/full-column distinction can occur.

The first failed implication is "pairwise column uniformity preserves the
actual geometric success event." The strongest survivor is the exact
agreement of linear/quadratic observables. The missing coordinate is the
three-column affine ratio, restored for this operation by (1)--(4).
No conclusion about a general search algorithm's running time is made.

## 5. Complete modular lines have an exact midpoint obstruction

Let L_(u,v)={(x,[ux+v]_p):x in [p]}, where u is nonzero. It is the complete
set of representatives of a non-axis affine line over F_p. Whenever two
points of L have the same two parity coordinates, their midpoint is an
integer point in the box and is also on L: 2 is invertible modulo p.
Every such unordered endpoint pair yields one distinct integer three-term
arithmetic progression; conversely every such progression has one endpoint
pair of equal parity. Therefore, if n_epsilon are the four parity counts,

```text
number of integer three-term AP triples in L
 = sum_(epsilon in {0,1}^2) binom(n_epsilon,2).          (5)
```

This counts equal-step triples, a subset of all integer-collinear triples.
It is not asserted to count unequal-step triples.

Write p=4q+r, r in {1,3}. Balancing four integer class sizes minimizes (5),
so every such line has at least

```text
A(p)=r binom(q+1,2)+(4-r)binom(q,2)                    (6)
```

AP triples. This lower bound is sharp over every affine column orbit of a
non-axis modular line. Indeed an affine column map can send its slope and
intercept to (2,0). For y=[2x]_p, the two no-wrap/wrap row intervals give
four parity classes of sizes q or q+1, with exactly r larger classes.

For p>=5, A(p)>0. Hence every affine image of a board containing a complete
modular line fails the actual no-three-in-line predicate. The saturated
two-parallel-line board {(i,i),(i,i+1 modulo p)} contains two disjoint such
lines and has at least 2A(p) genuine collinear triples under every affine
column map. This is an unbounded family of empty affine orbits. It does
not claim that unrestricted column permutations succeed for every p.
The p=5 instance supplies the exact positive unrestricted control above.

## 6. A separate capacity certificate for modular-line multiplicity

Let ell be any affine F_p line and Lambda its tangent preimage in Z^2.
This is an index-p lattice. Let w=(u,v) be a shortest nonzero vector for
the L1 norm, and put r=|u|+|v|. Then

```text
r <= floor(sqrt(2p)) < p, and gcd(u,v)=1.              (7)
```

An elementary pigeonhole proof avoids importing a geometry-of-numbers
bound. Apply (x,y)->(x+y,x-y) to Lambda. Its image Lambda' is an index-2p
subgroup of Z^2. If h=floor(sqrt(2p)), the (h+1)^2 points of {0,...,h}^2
occupy only 2p cosets, so two have a nonzero difference in Lambda'. The
inverse transform is integral on Lambda', and the L1 norm of its inverse
is the maximum absolute coordinate of that difference, at most h. Finally,
if a shortest vector had a nontrivial coordinate gcd d, then d<p and is
invertible modulo p; dividing by d would leave a shorter tangent vector.

Let S be any no-three-in-line subset of {0,...,N-1}^2; N need not equal p.
Points of S reducing into ell all have the same residue of vx-uy modulo p.
The integer functional vx-uy varies over an interval of width (N-1)r.
Each exact value is one Euclidean line, containing at most two points of S.
Consequently

```text
|{P in S: P modulo p lies in ell}|
 <= 2 [floor((N-1)r/p)+1].                            (8)
```

In particular, for N=p this is at most 2r<=2 floor(sqrt(2p)). The finer
directional certificate can use the exact occupied layer lengths of the
complete modular line, giving sum_j min(2,length_j). These are necessary
capacity bounds, not a sufficient characterization of no-three-in-line.

The short-direction bound itself is attained whenever
p=k^2+(k+1)^2 is prime. The orthogonal vectors (k,k+1) and (-(k+1),k)
generate an index-p lattice. The four unit basis choices have L1 length
2k+1; every other nonzero basis combination has Euclidean length at least
sqrt(2p)>2k+1. Thus r=2k+1=floor(sqrt(2p)). The verifier includes
p=5,13,41,181. No infinitude of primes of this form, or sharpness of the
resulting no-three-in-line cardinality bound, is asserted.

## 7. Exact verification, boundaries, and next question

The standalone verifier imports no repository implementation. Its universe
contains all 7956 triples with distinct rows and columns at p=3,5,7;
all 2100 row triples at p=11,13,17,19 on a modular-line representative;
the full 120-element S5 orbit and all its affine suborbits; every one-/two-column
cylinder; and 8253 slope/intercept tangent-lattice controls through p=43.
It checks 8244 midpoint hostiles, exact AP counts against literal triples
through p=13, and all safe modular-line subsets at p=3,5. Literal determinant
tests and an independent primitive-direction test agree on every S5 board.
The continuous/unbounded claims are proved above, not inferred from these
finite banks. All **117652 always-active gates** pass normally and with
optimization; reproduction and frozen hashes are recorded below.

```text
python3 -B 04-computation/continuing4_20260906_wildcard_affine_columns.py
python3 -B -O 04-computation/continuing4_20260906_wildcard_affine_columns.py
```

Source and outputs share this report's stem, with `.py` and `.out` extensions; the normal and optimized
runs reproduce the same retained `.out`. Frozen raw-LF SHA-256 hashes:

- source: `4c3317fb0aa7989927887bc7eced1cbd2a045975c9fb589f765ab063799927d0`;
- each output: `68290e80e06968f81ca4c4619e0c0ee751563fc752d06dc3c087b06f3c8a0844`.

The successful research move is to restrict the column operation while
retaining its exact lost ratio and carry, then test the actual zero event.
It uses the maintained "existence is a maximum or tail" and "respect
symmetries" cards: affine maps preserve the residue predicate, while the
integer lift is not equivariant without its interval sidecar. No new general
meta-pattern is claimed on this one application.

The surviving constructive question is whether useful non-line modular
boards have uncovered translation residues in (4), or whether a bounded
number of non-affine swaps can be certified to create one. The p=5 repair
shows what information such a method must change; it is not evidence for
an all-prime repair theorem. No generic density, random-greedy, Smith, or
finite-length stability problem is reopened here.
