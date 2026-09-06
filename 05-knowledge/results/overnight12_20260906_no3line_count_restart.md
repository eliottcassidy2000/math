# Counting successful saturated boards and the limit of uniform-column restarts

**Status: PROVED elementary consequences of the audited eleventh conditional
theorem; independently audited.** The finite counts below are FINITE-EXACT.
These are enumerative and algorithmic consumers, not a new extremal bound on
the maximum number of points and not a lower bound for arbitrary algorithms.

## Inheritance and the retained measure

The [initial source analysis](overnight_20260906_no3line.md)
proved the exact saturated-board count and the nonuniformity of ordered
permutation-pair sampling. The
[eleventh theorem](overnight11_20260906_no3line_rowfreeze.md)
proves its probability bound for every fixed row ordering, retaining the
induced row-subgraph matching polynomials. The least-used sidecar here is
the column stabilizer with the rows held pointwise fixed. It differs from
the number of alternating decompositions of a whole cycle skeleton.

The board is: uncolored boards, cycle profiles, fixed-row column orbits,
conditional probability, adaptive histories, and counting entropy. The map
from the conditional theorem to a search lower bound retains a fresh uniform
column permutation after every chosen history. It loses all information about
an algorithm that improves the current column order. The latter is outside
the theorem, not covered by calling it another restart.

Put alpha=1-5/e^2 and, for n>=4, define

```
b(n)=max(0,alpha*n-21,2(n-9)/15+2/[n(n-1)]),
E(n)=b(n)^2/[8(n-1)],
E_all(n)=max(E(n),n^2/[900(n-1)]),
c=alpha^2/8=0.013067267481528453... .
```

Here E is valid with arbitrary rows fixed in advance. E_all combines it
with the tenth theorem only in the original two-permutation model or in
uniform uncolored-board sampling. The distinction is load-bearing.

## 1. Exact orbit normalization and the counting bound

Let B_n count simple n-by-n zero-one boards with every row and column sum
two, and let N_n count those boards with no collinear triple over the
integers. These are labeled grid boards, without a dihedral quotient.
Every 2n-point no-three-in-line set is one of these saturated boards.
The inherited exact identity is

```
B_n=(n!)^2 [z^n] exp(-z/2)(1-z)^(-1/2).
```

If a skeleton has c_(2r) components of length 2r, its row-and-column orbit
has size `(n!)^2/product_r[(2r)^c_(2r) c_(2r)!]`. Uniform independent
label permutations are uniform on that orbit. Decomposing the uniform
uncolored measure into those orbits and applying the uniform eleventh and
tenth bounds therefore proves

```
N_n <= B_n exp(-E_all(n)),                         (1)
limsup (1/n) log(N_n/B_n) <= -c,                  (2)
```

with log0=-infinity. This mixture argument weights each orbit by its actual
size. Uniform ordered pairs of disjoint permutation matrices would instead
weight each board by 2^(number of cycle components).

There is also a useful fixed-row statement. Hold a skeleton and every row
label fixed. Its stabilizer under permutations of the column vertices has
size exactly `2^c4`, where c4 counts its four-cycles. Indeed a stabilizer can
only permute columns with identical two-row neighborhoods. Any identical
pair consumes both edges of those two rows and forms one C4 component;
three identical columns are impossible because row degrees are two. Every
C4 gives one independent swap and there are no other stabilizers. Thus the
column orbit has size `n!/2^c4`, uniformly sampled by a random column
permutation. In that orbit the number of successful distinct boards is at
most

```
(n!/2^c4) exp(-E(n)).                             (3)
```

The two powers of two above count different objects. A C6 has two ordered
alternating decompositions but no nontrivial column stabilizer when its
rows are fixed. It is the smallest cycle separating them.

For a direct coarse entropy consequence, label the two half-edges at each
row and column. There are (2n)! bijections of these half-edges. Every simple
board has exactly `2^(2n)=4^n` preimages, while some bijections create a
double edge. Hence B_n<=(2n)!/4^n. Combining with (1) and the elementary
factorial asymptotic gives

```
log N_n <= 2n log n -(2+c)n + O(log n).             (4)
```

The exponent E_all(n)=cn+O(1). The upper bound still grows like
exp(2n log n). Consequently an exponentially small success fraction does
not imply N_n=0, nor the Guy--Kelly conjectured extremal constant. Equation
(4) changes a linear entropy term and leaves the leading n log n scale.

## 2. Arbitrarily adaptive restarts within the stated sampling rule

Consider any randomized search at fixed n. Before attempt j it may use the
entire past to select a simple degree-two skeleton G_j and a row ordering
rho_j. Conditional on that past and that choice, its current column
permutation must be uniform. Let T be the first successful attempt, allowing
T=infinity. The eleventh theorem gives, on every history with no earlier
success, conditional success probability at most q=exp(-E(n)). Therefore

```
P(T>m) >= (1-q)^m,
P(T<=m) <= 1-(1-q)^m <= m q,
E T >= 1/q = exp(E(n)).                            (5)
```

Proof: multiply the conditional failure lower bound 1-q inductively,
then sum P(T>m), m>=0. The proof uses a conditional hazard bound, not
independence of the successive skeleton choices or success events. An
algorithm that stops without success can be assigned T=infinity, and the
same inequalities hold. In particular the expected attempt count is at
least exp(cn+O(1)). Time per attempt is not included.

The quantifier is essential. Choosing columns after inspecting the current
board need not leave them conditionally uniform, even if an unconditional
permutation somewhere in the procedure was uniform. The bound cannot be
imported into local search, backtracking, SAT solving, or a construction.
Similarly the extra tenth exponent in E_all is not automatically valid
for a search that fixes arbitrary rows before each attempt.

## 3. Finite controls and reproduction

The standalone source constructs every board by choosing an unordered pair
of column neighbors at each row, pruning only impossible remaining column
degrees. It imports no previous enumerator. Every board is checked by two
geometric definitions: literal integer determinants on triples and repeated
primitive unoriented directions from a point. The complete universe is
70,087 boards over n=2,...,6:

| n | B_n | N_n | Selected-slope zero | Fixed-row column orbits |
|---:|---:|---:|---:|---:|
| 2 | 1 | 1 | 1 | 1 |
| 3 | 6 | 2 | 4 | 1 |
| 4 | 90 | 11 | 35 | 6 |
| 5 | 2,040 | 32 | 545 | 22 |
| 6 | 67,950 | 50 | 12,155 | 130 |

The 160 column orbits check (3)'s stabilizer and the exact-mean conditional
tail inequality with a rational Taylor upper bound on its exponential.
Every cycle-profile size is independently checked by the orbit formula.
Weighting each board by its literal number 2^components of ordered
decompositions gives n! times the derangement count, preserving the hostile
sampling bias. The exact B_n recurrence

```
B_0=1, B_1=0,
B_n=n(n-1)B_(n-1)+n(n-1)^2 B_(n-2)/2
```

is checked through n40 against a separate formal convolution of the two
generating functions. Rational Taylor arithmetic verifies alpha>0.323 and
prints conservative finite exponents at n64 through4096. A finite branching
bank checks the survival induction for genuinely history-dependent hazards.
Those arithmetic examples supplement the general proofs; no asymptotic or
algorithmic theorem is extrapolated from a finite universe.

```
python -B 04-computation/overnight12_20260906_no3line_count_restart.py
python -B -O 04-computation/overnight12_20260906_no3line_count_restart.py
```

Both modes must pass 238,047 always-active gates with identical LF output.
The [independent referee](overnight12_20260906_no3line_count_restart_audit.md)
accepts every proof and scope restriction, with 28,336 exact normal/optimized
gates from a separate row-multigraph enumerator and literal stub bijections.
It also supplies a fixed-C8 restart strategy whose every marginal permutation
is uniform but whose mean success time is79/8, below the reciprocal orbit
success probability12. This records why marginal uniformity cannot replace
the conditional hypothesis in (5).

Frozen source SHA256:
`047e5cc5c9d374db311ddb41cce741be791de950302a4244c2a687bfa7556540`.
Frozen output SHA256:
`a3837d9b8539426ef467e9bd9a72c84008a0c3a8d40c19fd113bad4f957fd6ce`.

The next substantive probability question is a stronger rate or a bound
covering purposeful within-board changes. Repeatedly recounting the same
small grid sizes is not needed for that question.

**Filing:** root integrated these independently audited artifacts in the twelfth
checkpoint. Reproduction commands are relative to the repository root.
