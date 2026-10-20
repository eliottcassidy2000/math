# A one-chord extension of the tribonacci partner

**Status.** PROVED: the two-tail Hamiltonian-path formula, rational transfer,
rigidity, self-converse classification for order at least six, the central
blue-count law, and the exact leaf-to-partner line family below.
FINITE-EXACT: the stated independent controls and small exceptions.
OPEN: a grammar for unrestricted bounded-range defects, or a class-counting
recurrence for all tournaments. This note does not claim either.

## Inheritance and the new object

The closest proved mechanism is the March 5 S11 run decomposition of the full
tiling tournament, in
`03-artifacts/drafts/tiling-class-structure.md` and
`04-computation/verify_tribonacci_proof.py`. Its current recovery is
**THM-4731**, `01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md`,
statement/proof 6. The tribonacci sequence itself is inherited, not new here.
The canonical hostile is the previously corrected guess
`H(T*_n)=2^(n-2)+1`, first failing at order seven. The relevant corrected near
miss is the lost clone-allocation coordinate in
`05-knowledge/results/weighted_core_transport_20261020/allocation.md`.
The less-used sidecar is the March decomposition's endpoint condition: a
singleton run cannot border another singleton run.

Use vertices `1,...,n` and a fixed Hamiltonian path in increasing order.
Let `R_n` have the forward arcs `i -> i+1` and every longer arc backward.
Write `tau` for complementing the free, nonpath arcs and `sigma` for
converse followed by reversing the labels. These are the tiling operations
of **THM-643**, `01-canon/theorems/THM-643-blue-black-line-parity-structure.md`.
Thus `R_n=tau(TT_n)`. The path-fibre dictionary is **THM-1430**,
`01-canon/theorems/THM-1430-the-tiling-class-metagraph-dictionary-and-which-tricks-pay.md`:
ordinary fibre size is `H(T)/|Aut(T)|`; a blue tiling is sigma-fixed.
All uses here also follow the setting of THM-4731.

For integers `m,r>=0`, put

    n=m+r+3,  j=m+1,
    S_(m,r) = R_n with the one arc (j+2 -> j) reversed.

The defect is one extra forward chord of length two. The parameters are the
numbers of vertices below and above its central triple. This is genuinely
beyond changing only adjacent arcs in a backward transitive reference:
deleting an adjacent forward arc creates a uniform cut and merely factors
the tournament into an ordered join of smaller `R` blocks. The length-two
defect crosses that cut and creates a new path grammar.

The live coordinates are the defect location, ordered run boundaries,
endpoint permissions, converse symmetry, and the distinction between paths
and isomorphism classes. The connection developed here maps a simple
`TT_m -> C3 -> TT_r` core through tau to a growing tribonacci-like fibre.
It preserves a particular line and its reflection; it does not identify the
other lines incident at that fibre.

## 1. Exact run grammar and Hamiltonian-path formula

Define three tribonacci sequences, including their zero terms:

    h: 1,1,1,3,5,9,17,31,...
    a: 1,0,1,2,3,6,11,20,...
    F: 0,1,2,3,6,11,20,37,...

Each obeys `v_t=v_(t-1)+v_(t-2)+v_(t-3)` for `t>=3`.
Here `h_t=H(R_t)`, with `h_0=1`. For `t>0`, `a_t` counts paths of `R_t`
not beginning at the top vertex, equivalently not ending at the bottom
vertex. The value `a_0=1` means that an absent boundary imposes no restriction;
it is not a statement about an endpoint of an empty path. Also
`F_t=sum_(u=0)^(t-1) h_u`.
For `t>=1`, the first descending interval contains the top vertex; avoiding
that top vertex as the path's start is exactly the condition that this first
interval has length at least two. Removing it leaves an arbitrary `R` path,
so `a_t=F_(t-1)`. This verifies the boundary sequence and its stated initials.

**Proposition 1 (PROVED).** For every `m,r>=0`,

    s_(m,r) := H(S_(m,r))
      = h_(m+r+3)
        + F_m F_(r+1) + F_(m+1) F_r
        - h_r a_m - a_r h_m.                         (1)

**Proof.** In `R_n` every maximal ascending run is an increasing consecutive
interval. The disjoint intervals appear in decreasing order of their
positions. At the boundary between two intervals the required downward arc
exists unless both intervals are singletons. Thus the old paths are exactly
the descending lists of consecutive intervals with no two adjacent
singleton blocks. This is the inherited March grammar.

Partition the change in paths according to whether they use the changed arc.
An old path that uses `j+2 -> j` must have the two successive runs

    [j+1,j+2], [j]       or       [j+2], [j,j+1].

The vertex `j+1` must occur in one of these two intervals because the
intervals are consecutive in their descending list. In the first case the
upper tail is unrestricted, while the lower tail may not start with its top
singleton. This gives `h_r a_m`. In the second case the upper tail may not
end with its bottom singleton, giving `a_r h_m`. The cases are disjoint.

A new path using `j -> j+2` has a unique ascending run with one hole.
Put `c=j+1`. Its run has the form

    [c-L,...,c-1], [c+1,...,c+R],       L,R>=1.

Every other ascending run is a consecutive interval, so the skipped vertex
`c` is a singleton run. All intervals above this special run precede both
special objects, and all intervals below it follow both. There are precisely
two possible middle arrangements:

* `[c]`, then the run, allowed exactly when `L>=2`;
* the run, then `[c]`, allowed exactly when `R>=2`.

Each external junction automatically has downward length at least two.
The lower and upper remaining tails are arbitrary `R` paths. Summing the
first arrangement over `2<=L<=m+1`, `1<=R<=r+1` gives `F_m F_(r+1)`;
summing the second gives `F_(m+1) F_r`. The two arrangements are disjoint.
Unchanged paths avoid the changed arc in both directions. Subtracting the
lost paths and adding the gained ones proves (1). ∎

This proof is an exact path-set grammar, including the zero-tail boundaries.
It is stronger than an empirically fitted recurrence. The script separately
constructs the old interval paths and the new hole paths, then compares the
resulting vertex-path sets with a direct graph search.

For example, the rows indexed by `m=0,...,n-3` are

| order | `H(S_(m,n-m-3))` |
|---|---|
| 3 | 1 |
| 4 | 5, 5 |
| 5 | 9, 13, 9 |
| 6 | 15, 23, 23, 15 |
| 7 | 29, 41, 41, 41, 29 |
| 8 | 53, 77, 73, 73, 77, 53 |
| 9 | 97, 141, 137, 129, 137, 141, 97 |

## 2. Two exact tribonacci directions and a minimal linear transfer

Put `D(z)=1-z-z^2-z^3`. The three one-variable generating functions are

    H(z)=(1-z^2)/D(z),
    A(z)=(1-z)/D(z),
    F(z)=z(1+z)/D(z).

**Proposition 2 (PROVED).** The bivariate series is

    sum_(m,r>=0) s_(m,r) x^m y^r = P(x,y)/(D(x)D(y)),          (2)

    P = 1+4x+4y+3x^2+4xy+3y^2+2x^2y+2xy^2+x^2y^2.

Consequently `s_(m,r)` satisfies the tribonacci recurrence separately in
either tail once that index is at least three. The complete initial block is

              [1   5   9]
    B       = [5  13  23],       det B = -4.                  (3)
              [9  23  41]

**Proof.** Formula (1) is a linear combination of products and translates of
tribonacci sequences, proving both recurrences. For an explicit algebraic
check of (2), set `G(z)=sum_(t>=0) h_(t+3) z^t=(3+2z+z^2)/D(z)`.
Then the `h_(m+r+3)` term has generating function
`(xG(x)-yG(y))/(x-y)`; the other four terms have generating functions
`F(x)F(y)/y`, `F(x)F(y)/x`, `-A(x)H(y)`, and `-H(x)A(y)`.
The divisions by x and y are valid since `F(0)=0`. Their sum simplifies to
(2). Equivalently, multiply the independently proved tail recurrences and
the 3-by-3 initial block by `D(x)D(y)`. ∎

Let `u_t` be the vector-valued tribonacci solution with
`u_0=e_0,u_1=e_1,u_2=e_2`. Then

    s_(m,r) = u_m^T B u_r.                                  (4)

The three coordinates on each side therefore suffice for this entire
infinite family. Since the first 3-by-3 minor has nonzero determinant, no
representation of this bivariate array as a sum of two separable terms
`f_i(m)g_i(r)` is possible over a field of characteristic zero. This is a
rank-three linear/bilinear lower bound, not a lower bound on arbitrary
nonlinear encodings.

The positive numerator also gives a positive weighted counting grammar:
choose a seed `(i,j)` with weight `P_[i,j]`, then independently append words
of steps `1,2,3` totalling `m-i` and `r-j`. This counts (1) exactly.
The direct correspondence to individual Hamiltonian paths is the run
construction in Proposition 1, not an asserted extra natural bijection
between those paths and these weighted seed words.

## 3. Rigidity, reflection, and the blue coordinate

**Proposition 3 (PROVED).** Every `S_(m,r)` is rigid. For `n=m+r+3>=6`,

    S_(m,r) is self-converse  iff  m=r.                       (5)

Thus the ordinary fibre size of this family is exactly (1), rather than a
path count still awaiting an automorphism correction.

**Proof.** An automorphism preserves outdegrees. For an internal chord
`2<=j<=n-3`, the vertices `j,j+1,j+2` form a transitive equal-degree triple,
and all other degree cells are singletons or the two transitive endpoint
pairs inherited from `R_n`. For an endpoint chord, the altered degree cells
are transitive pairs. Orders three and four have the same property by direct
inspection. Each degree cell is therefore a rigid transitive tournament;
every automorphism fixes every vertex.

Label reversal gives

    sigma(S_(m,r)) = S_(r,m).                               (6)

If `m=r` it is an anti-automorphism. For an internal chord, the degree
histogram of `R_n` changes by
`-X^(j-1)+2X^j-X^(j+1)`. The original histogram is palindromic under
`d -> n-1-d`. The perturbation is palindromic exactly when
`j=(n-1)/2`, equivalently `m=r`.

For the left endpoint chord `j=1`, the degree histogram remains palindromic,
so this necessary condition alone is insufficient. The degree-one pair is
`2 -> 3`, and the degree-`n-2` pair is `n-1 -> n`. An anti-automorphism would
have to send `3` to `n-1`. For `n>=6`, vertex 1 is the unique degree-two
vertex and `n-2` the unique degree-`n-3` vertex, so it must send `1` to `n-2`.
The new arc `1 -> 3` would then require `n-1 -> n-2`, contradicting the
fixed path arc. The right endpoint case is its converse. ∎

**Proposition 4 (PROVED).** The central family's blue fibre multiplicity is

    b(S_(m,m)) = h_(m+2),       m>=0.                         (7)

For order at least six, every noncentral member has blue multiplicity zero.

**Proof.** Rigidity makes label reversal the unique anti-automorphism of the
central tournament. A blue presentation is therefore a Hamiltonian path
whose reversed vertex list is obtained by applying this involution.
In an odd-length list, two involution-paired vertices occupy positions with
even separation. They can never be consecutive. The only flipped chord
joins just such a pair, so no compatible path uses either orientation of
that chord. The compatible path set is unchanged from `R_(2m+3)`.
THM-4731 statement/proof 6 gives its size `h_(m+2)`; `m=0` is the direct
three-vertex case. Noncentral members of order at least six are not
self-converse by (5), so have no blue presentations. ∎

The path-set argument proves a more general elementary survivor: for a fixed
reversing involution on an odd-order tournament, flipping any selection of
arcs joining its paired vertices preserves all Hamiltonian paths compatible
with that involution. The involution remains arc-reversing after those flips.
This does not by itself preserve the entire blue fibre, because other
anti-automorphisms or the automorphism group might change.

FINITE-EXACT small exceptions: every member of orders three through five is
self-converse, with blue rows respectively `(1)`, `(3,3)`, and `(1,3,1)`.
In particular an equal Hamiltonian count does not establish equality of
blue counts or isomorphism classes.

## 4. Exact metagraph leaves and their growing partners

**Proposition 5 (PROVED).** In the indicated fixed-path presentations,

    tau(S_(m,r)) = J_(m,r) := TT_m -> C3 -> TT_r.             (8)

Every `J_(m,r)` has `H=3`, automorphism order 3, and ordinary fibre size 1.
It is self-converse exactly when `m=r`.

**Proof.** Complementing the free arcs restores the transitive tournament
except for `j+2 -> j`, which makes the three consecutive central vertices a
directed cycle. Every earlier vertex dominates this cycle and every later
vertex, and the cycle dominates every later vertex. A Hamiltonian path must
visit these ordered components consecutively; the transitive pieces have
one path each and the cycle has three. Automorphisms preserve the ordered
strong components, fix each transitive singleton, and may rotate the cycle.
Converse reverses the strong-component order, giving exactly the condition
`m=r`. ∎

For `n>=6`, the `n-2` tournaments `S_(m,n-m-3)` are pairwise nonisomorphic.
Distinct internal positions have distinct degree histograms; the endpoint
histograms are equal but those two tournaments are nonisomorphic, since
they are converse and neither is self-converse by (5). Endpoint histograms
are different from every internal one. Hence this family contains
`ceil((n-2)/2)` distinct converse-merged partner nodes. The same number of
merged leaves `J` follows from the position of their unique cyclic strong
component.

Equation (8) determines their full outgoing incidence:

* If `m=r`, the `J` node has one tiling, fixed by sigma. Its single line is
  blue and joins it to the central `S` node.
* If `m!=r`, the merged `J` node consists of two converse classes, each with
  one tiling. These two tilings give two black lines into the same merged
  `S` node. They are distinct sigma twins.

For these orders the leaf and partner are distinct, since `J` has three
automorphisms and `S` is rigid. These are exact incidences, not just degree
parities. The rest of the incident lines at the much larger `S` fibre
remain uncomputed by this formula.

This is the useful small-object/large-fibre connection: a single cyclic
component with two transitive tails produces an explicitly calculable
partner under tau, with exponential path growth and a rank-three
two-sided transfer. It generates a family of metagraph incidences without
enumerating its full tiling cubes. It supplies only linearly many nodes at
each size; it is not a recurrence for A000568.

## 5. Boundaries of compression and operation closure

Two exact hostile cases show where a generalization must keep more state.

First, start with `R_6` and reverse two gap-two arcs. The forward-chord sets
`{(1,3),(3,5)}` and `{(2,4),(3,5)}` have the same defect count and the same
degree multiset `(1,2,2,3,3,4)`, but Hamiltonian counts 25 and 29.
Thus defect number plus degree distribution does not determine even `H`.
The missing coordinate is the relative position and overlap of the chords.

Second, local boundary degrees do not suffice to update Hamiltonian
connectivity. Use vertices `c=1,a=2,u=3,v=4,b=5,d=6`. The two partial covers

    a -> u -> b, c -> v -> d
    a -> u -> d, c -> v -> b

have identical incoming and outgoing degree at every named vertex. Both
occur in the same tournament obtained from `R_6` by forwarding
`(1,4),(3,5),(3,6),(4,6)`, all of range at most three.
Adding the same edge `b -> c` completes a Hamiltonian path in the first
cover and closes a directed three-cycle in the second. An exact
path-cover transfer must distinguish these boundary connectivity states
or retain other information that determines the same distinction. This
is an obstruction to a degree-only state, not a proof against every
bounded-state construction.

The precise operation claims are:

* Increasing `m` or `r` adds a prescribed endpoint vertex in `S`, and adds a
  transitive singleton on that side in `J`. The complete two-tail family
  is closed under these two additions. Equations (2)-(4) give its exact
  counting transfer.
* Reflection exchanges `m,r`, preserves (1), and gives the blue condition
  (5)-(7).
* Tau exchanges these two families at the distinguished tilings in (8).
  It is not asserted to map every presentation in an `S` fibre to the
  same `J` class; indeed the sizes already forbid that.
* Uniform transitive doubling takes the family outside its one-chord
  normal form. Each changed coarse arc is replicated across four clone
  pairs, and Hamiltonian paths can interleave the clones. The exact
  dilation sidecar is the allocation word and per-core shuffle in
  `weighted_core_transport_20261020/allocation.md`.
  The fixed-path correction from
  `metagraph_growth_20261020/growth.md` remains valid:
  `tau_(kn) D_k(t)=D_k(tau_n t) xor E_(n,k)`, where
  `|E_(n,k)|=(k-1)(3nk-2k-2)/2`. Ordinary tau therefore does not commute
  with dilation; no claim of one-chord closure under doubling is made.

The strongest next target suggested by this packet is a grammar for several
interacting hole-runs, retaining their endpoint connectivity and the order
of the remaining holes. A fixed finite motif with arbitrarily long ordinary
tails is a more controlled extension than guessing a recurrence solely
from defect counts. Its state size and its automorphism correction must be
proved independently; the rank-three result here does not supply them for
free.

## Reproduction and exact scope

Run from the repository root:

    python 05-knowledge/results/defect_boundary_20261020/grammar_check.py
    python -O 05-knowledge/results/defect_boundary_20261020/grammar_check.py

The code uses explicit exceptions rather than assertions. Its inherited
bitset subset DP is independent of the new interval/hole-run enumerator.
The recorded output covers:

* all 55 chord positions at orders 3 through 12 against subset DP;
* all 21 chord positions at orders 3 through 8 as complete vertex-path sets;
* a 31-by-31 coefficient grid for the generating function and both recurrences;
* all 171 positions at orders 3 through 20 for rigidity, self-converseness,
  literal tau/reflection, and class distinctions;
* 36 independent leaf path counts, and central blue counts through order 15;
* both hostile controls above, including their common ambient tournament.

Finite controls support the proofs and catch boundary mistakes. They are
not the proof of the infinite claims.
