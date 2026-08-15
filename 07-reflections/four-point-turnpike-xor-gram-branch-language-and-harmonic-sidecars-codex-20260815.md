# Four-point turnpikes: XOR incidence, Gram structure, branch languages, and harmonic sidecars

**Research synthesis, 2026-08-15.**  The truth-bearing atlas is the proved,
exact, and independently hostile-audited
[THM-3457](../01-canon/theorems/THM-3457-four-point-line-metric-turnpike-preorder-atlas.md).
During the session the
incoming [THM-3456](../01-canon/theorems/THM-3456-left-permutive-trace-bijection-and-rule30-seed-boundary.md)
Rule 30 package was independently audited and promoted; its trace/seed split
is therefore used below as a proved comparison with a different target.

## Outcome

For four ordered line points there are three positive gaps and six distances:

```text
(a,b,c) |-> (a,b,c,a+b,b+c,a+b+c).
```

The five possible equality walls cut the positive projective gap triangle
into exactly 25 preorder strata:

```text
10 strict + 10 one-tie + 4 two-tie + 1 maximal-tie.
```

The stable Fibonacci window is only one wall stratum.  Its equality
`c=a+b` is visible to the six-distance preorder, but its marked-origin seam
and Pell unit are not.  This explains both the positive connection and the
failure of the converse without appealing to numerical coincidence.

## Inheritance pass

| role | object |
|---|---|
| closest proved mechanism | THM-3454's same-`U`-ray distance isometry and marked Fibonacci seam |
| canonical hostile | gaps `(1,3,4)` at origin index `1`: Fibonacci tie order, failed origin seam |
| corrected near miss | MISTAKE-402: spine index `t` has actual rooted depth `t-1` |
| least-used sidecar | XOR weight on the six edge vectors, which remembers adjacent versus opposite `K4` edges |
| wildcard input | the `4 x 4` Gram organization in Rybin--Zhang--Luo's RXTX algorithm |
| incoming proved comparison | THM-3456's free-trace tree versus distinguished-seed boundary |

## Concept board

1. **Four ordered points.**  Their geometric order is a transitive `T4`.
2. **Six interval sums.**  Comparing the six `K4` edge distances gives a
   total preorder on six vertices, not another tournament on four vertices.
3. **Boolean edge layer.**  The six edges are the weight-two vectors of
   `F_2^4`; XOR detects their incidence orbit.
4. **Projective arrangement.**  The atlas is the positive-cone sign
   stratification of five rational linear forms.
5. **Rank-one Gram locus.**  Squared line distances are linear functionals
   of a centred `4 x 4` rank-one PSD Gram matrix.
6. **Sparse recurrence selector.**  Fibonacci chooses spine indices on one
   Berggren branch; it does not define a second child operation.
7. **Boolean/harmonic branch carrier.**  A subset of the naturals is a bit
   branch and a reciprocal subseries, but a scalar sum forgets its labels.

The board changed after the exact atlas was derived.  XOR stopped being a
candidate orientation and became the exact incidence sidecar.  The Gram view
turned the paper's size-four block structure into a typed map rather than a
cardinality analogy.  The recurrence wall exposed the missing origin.  The
Rule 30 input then supplied a second proved branch-versus-language instance
with a different missing coordinate.

## What “XOR is essentially size four” can mean exactly

Encode edge `ij` by

```text
v_ij=e_i+e_j in F_2^4.
```

For two distinct edges, `wt(v XOR w)` is two when they meet and four when
they are disjoint.  Complementation by `1111` pairs the three opposite-edge
pairs.  Thus the weight-two layer is the octahedral graph `L(K4)`, and `S4`
acts by permuting the four Boolean coordinates.

This gives a clean division of labour:

```text
XOR weight     -> adjacent/opposite incidence,
distance value -> orientation or tie,
line order     -> positional gauge.
```

XOR is symmetric, so it cannot orient a pair by itself.  The atlas adds the
strong constraint that a lone tie is always between adjacent edges.  The
only opposite-edge equality is paired:

```text
01=23 iff a=c iff 02=13.
```

That forced companion is the first non-cosmetic payoff of combining the XOR
and turnpike views.

## Size four, size six, and missing/both-way edges

If all six distances differ, increasing value gives a transitive `T6`.  A
tie deletes one strict comparison.  If weak comparisons are drawn instead,
the same tie becomes a two-way edge.  Therefore “tournaments with missing
edges” and “semicomplete digraphs with both-way edges” are two presentations
of the same total preorder only when the tie blocks are retained.

The four atlas types have respectively `15,14,13,11` strict arcs and
`1,2,4,12` linear extensions.  The last row has tie blocks of sizes three,
two, and one; saying “four ties” means four tied unordered pairs, not four
blocks.

The original four points also carry a transitive `T4` by ancestry/order, but
that object forgets every gap magnitude and every equality among the six
edge distances.  Neither tournament determines the other without sidecars.

## The hyperplane-arrangement mechanism

After scaling to `a+b+c=1`, the only equality forms are

```text
a-b,      b-c,      a-c,      a-b-c,      c-a-b.
```

Their feasible signs in `{-1,0,+1}^5` are the atlas strata.  Exact rational
feasibility, rather than an integer box scan, is the right companion because
every stratum is a relatively open rational polytope.  This representation
suggests a generalization: for `m+1` ordered points, compare all contiguous
gap sums and classify the positive part of the resulting interval-sum
arrangement.

Recurrence laws then become selected walls.  At four points, the Fibonacci
gap law is `c=a+b`; at longer windows, the complete recurrence produces a
chain of interval-sum equalities plus one origin alignment.  The arrangement
sees the homogeneous walls but not the absolute alignment.

## Reduced fractions, the Berggren ternary tree, and Fibonacci

THM-3454 identifies

```text
t/(t+1)  <->  P_t=U^(t-1)(3,4,5)
```

on one fixed-cusp Farey leaf chain and one parabolic branch of Berggren's
ternary tree.  Their common distance is `|s-t|`.  Selecting
`t=F_n` therefore transfers Fibonacci index gaps into this line-metric atlas.

The typing matters:

- the Berggren tree is ternary under its three child matrices;
- the selected `U`-spine is one unary branch;
- Fibonacci is a sparse subset of the spine indices;
- actual rooted depth is `F_n-1`.

So the “tree of Fibonacci primitive triples” is a sparse selection on a
distinguished Berggren branch, not a new ternary ancestry law.  Farey
cross-determinants and Lorentz chords preserve the index differences, while
the six-distance preorder further forgets their magnitudes.

## The rank-one Gram and RXTX connection

For centred coordinates `z=Hx`, the matrix

```text
G=zz^t
```

is a rank-one `4 x 4` PSD Gram matrix and

```text
d_ij^2=G_ii+G_jj-2G_ij.
```

Thus the ten upper-triangular Gram coordinates split naturally into four
diagonal vertex terms and six off-diagonal `K4` edge terms.  The line-distance
preorder is the pullback of six linear functionals on the rank-one Gram cone.

Rybin--Zhang--Luo's
[*XX^t Can Be Faster*](https://arxiv.org/abs/2505.09814v2) gives a different
use of the same size-four symmetric organization: its RXTX scheme uses eight
recursive symmetric calls and 26 general products.  The repository imports
its exact span-enumeration and minimum-cover search architecture, not a
tournament claim.  RXTX works for general Gram products; it neither imposes
rank one nor orients the six edge functionals.  The common `4+6` indexing is
real, and the target predicates are different.

A useful next experiment is to feed the rank-one Gram identities into an
RXTX-style exact span-and-cover engine and ask for the smallest certificate
of each preorder wall.  The hostile gate is mandatory: a shorter algebraic
span that loses positivity or the line chart is not a turnpike proof.

## Every subset of the naturals as a harmonic branch

For `A subset N_(>=1)`, the following three objects are faithful before
aggregation:

```text
characteristic branch:  chi_A in {0,1}^N,
harmonic subseries:      (chi_A(n)/n)_n,
Dirichlet germ:          D_A(s)=sum chi_A(n)n^(-s), s>1.
```

The phrase “a subset of the harmonic series” is exact at the sequence level:
choose the terms `1/n` with `n in A`.  It is not exact after replacing that
labelled subseries by one scalar.  At `s=1`, some choices diverge and finite
collisions such as

```text
1/2=1/3+1/6
```

destroy injectivity.  The full Dirichlet germ is faithful by the
least-differing-index limit as `s` tends to infinity.

THM-3455 gives a nontrivial periodic example.  Its cap-seven spine subset has
spine-index harmonic coefficient `76/187`; Fibonacci pullback has labelled
index coefficient `43/90`, while reciprocal Fibonacci values, positive
rooted depths, and nonlinear branch labels yield convergent subseries.  The
Boolean set is the same only after the carrier map is declared.

## Proved Rule 30 branch-language comparison

THM-3456's left-permutive trace bijection says every binary branch occurs as a
free-input centre trace, while the distinguished single-seed trace is selected
only after restoring an inverse-boundary word.

It gives a second exact carrier for every subset of the naturals:

```text
A -> chi_A -> some free-input trace,
```

but not a realization by the named Rule 30 seed.  The missing coordinates in
the live examples would then be:

| quotient | forgotten coordinate |
|---|---|
| six line distances | one origin coordinate |
| distance preorder | magnitudes and origin |
| one harmonic scalar | all term labels |
| free trace language | inverse boundary selecting the named seed |

This is a shared controlled-forgetting pattern, not an isomorphism of the
underlying trees.  A binary trace-cylinder tree, a ternary Berggren ancestry
tree, and the unary `U`-spine have different branch operations and different
targets.

## Typed connection ledger

| source | target | map | preserved | lost | sidecar | cheapest hostile |
|---|---|---|---|---|---|---|
| four line points | six-edge preorder | compare interval sums | all ties/order | scale, magnitudes, translation | exact gaps and origin | `(1,3,4)` |
| `K4` edges | Boolean four-cube | `ij -> e_i+e_j` | adjacency/opposition | distance order | numerical weight | XOR is symmetric |
| line points | rank-one Gram | `x -> (Hx)(Hx)^t` | squared distances | translation, reflection sign | chart/orientation | arbitrary higher-rank Gram |
| Fibonacci indices | atlas wall | successive gaps | `c=a+b`, `a<b` | origin seam and Pell unit | `x_0` and Cassini norm | gaps `(1,3,4)` |
| subset `A` | harmonic carrier | `n -> chi_A(n)/n` | labelled membership | labels after summation | coefficient sequence/germ | `{2}` versus `{3,6}` |
| free trace language | named seed trace | THM-3456 inverse | prefix compatibility | seed boundary | inverse-boundary word | Rule 60 comparison |

## New frontiers

### 1. The `m`-point interval-sum arrangement

For positive gaps `g_0,...,g_(m-1)`, all line distances are contiguous sums.
Classify realizable weak orders, their tie-incidence graphs, and reversal orbits
without enumerating a large integer box.  The first decisive test is five
points; determine which apparently isolated opposite-edge ties are forced to
bring companions.

### 2. Recurrence-wall recognition

Given a linear recurrence, identify the smallest interval-sum wall packet it
forces and the minimal absolute sidecar needed for a converse.  Test Fibonacci,
Lucas, and one genuinely non-order-two recurrence.  A wall alone must face a
translated or wrong-origin hostile.

### 3. RXTX-style exact certificate search

Use rank-one Gram identities as candidate generators, enumerate exact linear
relations, and solve a minimum-cover problem for the 25 wall certificates.
Retain PSD rank, line chart, and strict positivity as target constraints.

### 4. LRC speed-window transplant

Sorting four speeds gives line gaps, but multiplying by time and reducing
modulo one changes the cut and can reorder the points.  Any transport must
retain winding integers, the circular cut, owners, and common time.  The first
hostile is a wall crossing under modular wrap; a static preorder is not a
physical LRC cover.

### 5. Full-germ branch observers

Compare finite jets and full Dirichlet/Fourier germs on the periodic THM-3455
rank word and on free Rule 30 traces.  Ask which smallest
observer separates named branches while retaining the operation that selects
them.

## Boundary

The atlas is a precise answer to the four-point/six-edge/XOR puzzle.  It does
not move the LRC(14) ledger, realize the
`7 x 13` bispectrum, construct the D5 `H^1` map, or classify Keller maps.
Its export is a reusable warning with a positive core: edge incidence may be
Boolean, orientation may be metric, and recurrence may be a tie wall, but
these are three different coordinates and none should be silently substituted
for another.
