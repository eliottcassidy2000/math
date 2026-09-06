# Independent audit of saturated-board counting and uniform-column restarts

**Status: PASS, proof-level independent audit.** Equations (1)--(5) of the
[producer report](overnight12_20260906_no3line_count_restart.md)
are valid with its stated measure and algorithmic scope. No mathematical
repair is requested. The companion supplies a different exact counting path,
literal configuration fibers, and a conditioning hostile. Its finite results
are FINITE-EXACT; the all-n statements use the proofs and inherited audited
probability theorems, not extrapolation.

Producer source reviewed: SHA256
`047e5cc5c9d374db311ddb41cce741be791de950302a4244c2a687bfa7556540`.
Producer output reviewed: SHA256
`a3837d9b8539426ef467e9bd9a72c84008a0c3a8d40c19fd113bad4f957fd6ce`.
Neither producer file was changed or imported by the independent program.

## 1. Inheritance, types, and the actual consumer

The relevant inherited statements were checked directly in the
[eleventh conditional theorem](overnight11_20260906_no3line_rowfreeze.md),
especially its equations (1)--(5), and the
[tenth two-permutation theorem](overnight10_20260906_no3line_defect.md),
equations (1)--(3). Both are already independently audited. Their different
randomness hypotheses remain different in the new corollary.

For a simple bipartite two-regular skeleton on n+n vertices and every fixed
row labeling rho, the eleventh result gives, for n>=4,

```
P(no collinear triple | rho) <= exp(-E(n)),
E(n)=b(n)^2/[8(n-1)],
b(n)=max(0, alpha*n-21, 2(n-9)/15+2/[n(n-1)]),
alpha=1-5*exp(-2).
```

The tenth result supplies the additional exponent n^2/[900(n-1)] when both
shore labelings are independent and uniform. A uniform bound passes to a
mixture of its eligible measures. A bound for one measure does not pass just
because the samples have the same support. This is the exact reason E_all
is valid for uniform uncolored boards but cannot replace E in arbitrary
row-frozen restarts.

The source is a probability estimate on full integer no-three-line success.
The map to counting multiplies by the cardinality of the correct uniform
orbit. The map to restarts uses the success hazard conditional on the entire
previous history and the current skeleton/row choice. These maps preserve the
zero event; they discard constructive strategies for choosing columns. Their
required sidecars are respectively orbit stabilizers and conditional laws.

N_n means the number of saturated boards that themselves have no collinear
triple. The phrase "count their subsets" in the producer's first definition
should be read in that sense; changing it to "count those boards" would remove
a harmless wording ambiguity. It cannot mean counting all point subsets of
each board. Every 2n-point no-three-line set on an n-by-n grid is saturated:
each of n rows contains at most two points, so all contain exactly two, and
the same argument applies to columns. There is no dihedral quotient.

## 2. Three distinct multiplicities and an independent representation

Let c be the number of connected cycles of the skeleton and c4 its number
of four-cycles. The following factors are all correct and are not
interchangeable.

* The full two-shore orbit has size
  `(n!)^2 / product_r[(2r)^c_(2r) c_(2r)!]`.
* With rows fixed pointwise, the column orbit has size `n!/2^c4`.
* A board has `2^c` ordered decompositions into two disjoint perfect matchings.

For a cycle of length 2r, its bipartition-preserving automorphism group has
order 2r. Permutations of equally sized components supply the factorials in
the first formula. This includes C4, whose shore-preserving group has order
four. A column permutation with rows held pointwise fixed can only permute
columns with identical neighborhoods. A pair of such columns consumes the
two available edges of both of its row neighbors and is exactly an isolated
C4. Three identical columns violate row degree two. Thus each C4 supplies
one independent swap and nothing else contributes. Finally an alternating
coloring of each even cycle has two choices, proving the third formula.
C6 has c=1 and c4=0, so it is the smallest separating control.

The independent program collapses each column to an unordered edge between
its two row neighbors. This gives a loopless two-regular multigraph on the
labeled rows. A doubled edge is allowed and is exactly a C4 of the original
board. This multigraph is precisely a fixed-row column orbit. Conversely,
assigning its edge multiset to the n labeled columns recovers all the boards
in that orbit, each once. Hence its number of boards is n!/2^c4.

Orienting each row cycle gives a derangement. Each row cycle of length at
least three has two orientations, whereas a doubled edge has just one
derangement orientation. Thus a multigraph with parameters (c,c4) has
exactly 2^(c-c4) preimages under this derangement map. This independently
explains both identities

```
sum_over_row_multigraphs 2^(c-c4) = !n,
sum_over_row_multigraphs (n!/2^c4)*2^c = n! * !n.
```

In particular the ordered-permutation model is a biased uncolored-board
measure. One may still apply a probability estimate uniform over all cycle
types to that mixture, but its weighted count is not the uncolored count.

The exact B_n generating function also follows from connected components:
on r specified row and r specified column labels the cycle count is
`r!(r-1)!/2`, for r>=2. The bivariate labeled component series is therefore
`sum_(r>=2) z^r/(2r)`, whose exponential is
`exp(-z/2)(1-z)^(-1/2)`. The independent program uses the corresponding
distinguished-row recurrence, rather than the producer's recurrence or
generating-function convolution:

```
B_0=1, B_1=0,
B_n=sum_(r=2)^n binom(n-1,r-1) binom(n,r)
                    * r!(r-1)!/2 * B_(n-r).
```

These arguments prove the all-n orbit and count formulas directly. Taking
the actual orbit-size mixture then proves producer equations (1) and (3).

## 3. Configuration fibers and the entropy term

Give each row and each column two distinguishable stubs. A simple board
has two choices for assigning its incident edges to the row stubs at every
row, and independently two choices at every column. These choices determine
a unique stub bijection. Thus every simple board has exactly 4^n preimages
among the (2n)! bijections. Bijections with repeated row-column edges are
discarded, establishing B_n <= (2n)!/4^n with the stated normalization.

For large n the term alpha*n-21 is the largest term of b(n), because
alpha>0.323>2/15. With c=alpha^2/8, exact polynomial division gives

```
E(n)=c*n+(alpha^2-42*alpha)/8+(alpha-21)^2/[8(n-1)].
```

Also c>1/900, so E_all(n)=E(n) eventually. Both equal c*n+O(1).
Consequently log(N_n/B_n)<=-E_all(n) proves the claimed limsup, including
the convention log0=-infinity. Stirling's factorial estimate yields

```
log((2n)!/4^n)=2n*log(n)-2n+O(log n),
log N_n <= 2n*log(n)-(2+c)n+O(log n).
```

This is an entropy reduction among saturated labeled boards. The resulting
upper bound still tends to infinity and does not imply nonexistence or an
extremal no-three-line constant.

## 4. Adaptive hazards, including unsuccessful termination

Let H_(j-1) contain all earlier information available to the algorithm.
The current skeleton and row order may be arbitrary measurable or randomized
choices from that information. Condition once more on those choices. The
uniformity hypothesis on the new column permutation then gives success
hazard at most q=exp(-E(n)) on every still-alive history. Integrating this
inequality over the choices retains the same bound conditional on H_(j-1).

Writing s_m=P(T>m), the tower property gives
`s_m >= (1-q)s_(m-1)`. Induction gives the survival and success inequalities
in producer equation (5); summing s_m for m>=0 gives E T>=1/q. This argument
requires no independence among success events and no fixed skeleton across
attempts. Constant hazard q is an equality control for the abstract hazard
lemma. For E=0, q=1 and the statements are valid but trivial.

If an algorithm stops unsuccessfully, defining T=infinity means its
first-success index is infinite. The theorem concerns that variable, not
the finite number of attempts actually executed before giving up. The
producer explicitly makes this convention. No claim about arbitrary
algorithmic runtime follows.

A concrete hostile isolates the first failed implication when conditional
uniformity is dropped. At n=4 take fixed column neighborhoods

```
({0,1}, {0,2}, {1,3}, {2,3}).
```

This is one C8. Among its 24 column permutations exactly two are successful;
the selected-slope defect mean is 5/6. List the 24 permutations in
lexicographic order, pick a uniform starting offset, then visit the list
cyclically. Each individual attempt has a uniform marginal distribution.
Nevertheless success occurs by attempt 21 for every offset, and the mean
first-success index is 79/8. The constant marginal success probability is
1/12, so the incorrectly imported geometric mean bound 12 fails. Every
positive geometric survival lower bound also fails by attempt 21.

The exact-mean conditional concentration exponent for this orbit is
`(5/6)^2/[8*3]=25/864>0`; importing that exponent from marginal uniformity
would likewise assert a positive tail after 21 attempts. The producer does
not make that import. Its coarse b(4) is zero, so this hostile is not a
counterexample to its literal small-n bound. It identifies exactly why the
conditional hypothesis is load-bearing.

## 5. Independent exact controls and freeze

The companion uses no producer import and checks 28,336 always-active gates.
It constructs row-neighborhood multigraphs from derangements through n=8.
The orbit totals are:

| n | Row multigraphs | Labeled boards B_n |
|---:|---:|---:|
| 2 | 1 | 1 |
| 3 | 1 | 6 |
| 4 | 6 | 90 |
| 5 | 22 | 2,040 |
| 6 | 130 | 67,950 |
| 7 | 822 | 3,110,940 |
| 8 | 6,202 | 187,530,840 |

All board counts, profile weights, orientation multiplicities and stabilizer
formulas are checked. This does not enumerate 187 million boards: the larger
rows are exact orbit sums. Literal boards are reconstructed only for n=2..5
(2,137 boards total), using unique edge-multiset orderings. A primitive
affine-line dictionary, separate from the producer's determinant and
direction routines, checks all integer collinearity. It recovers the
producer's success and selected-zero counts for this universe.

For each of those 30 fixed-row orbits the source also checks a stronger
rational small-orbit inequality: if x=mu^2/[8(n-1)], then
`good/orbit_size <= 1-x <= exp(-x)`; all tested x lie in [0,1).
The comparison uses the coefficientwise bound exp(x)<=1/(1-x), without
floating-point arithmetic or the producer's Taylor routine.

Literal stub bijections through n=4 reconstruct exactly the same simple
boards and have fibers 4^n. Their numbers of simple bijections are 16, 384,
and 23,040 among 24, 720, and 40,320 total bijections. The distinguished-row
component recurrence checks the producer's second-order count identity and
configuration bound through n=60. These finite checks supplement the
all-n combinatorial proofs above.

An independent binary observation tree through depth ten uses hazards
depending on the entire history, bounded by 3/7. It checks survival, the
union bound, and truncated tail-sum expectations with rational arithmetic.
Constant-hazard controls check equality. The actual C8 hostile above checks
the weaker-hypothesis failure. A fresh 20-term rational series encloses
alpha and verifies the eventual asymptotic branches.

Reproduce:

```
python -B 04-computation/overnight12_20260906_no3line_count_restart_audit.py
python -B -O 04-computation/overnight12_20260906_no3line_count_restart_audit.py
```

Both modes have identical LF output. Independent source SHA256:
`9c65c292dbd9dfd9a35928febd8af29404e577fbf8f31512ba4ace77fe0f860f`.
Independent output SHA256:
`900f1691974ab628dd2c48e76e46719d069e39ab9e056660efc64605ef9ce6fa`.

This audit adds a different representation and a sharp conditioning
firewall, not a new probability exponent or a lower bound for purposeful
within-board search. The remaining rate question must act on the actual
one-permutation law or preserve a sufficient conditional replacement for it.

**Filing:** root integrated these independently audited artifacts in the twelfth
checkpoint. Reproduction commands are relative to the repository root.
