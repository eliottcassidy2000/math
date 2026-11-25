# A positive weighted obstruction to finishing the 18-branch seed

**FINITE-EXACT.** The stored cyclic seed on 18 branches of 56 sites has 2,479 legal normalized exterior neighbor tuples. Every allowed difference between every pair of branches occurs in some tuple. Despite this local coverage, the seed admits no degree-57 Moore completion, including a completion that completely breaks its cyclic symmetry outside the seed. A 39-entry integer table certifies the contradiction. A sharper certificate bounds the number of compatible exterior sites by 1,604, enough for at most 28 additional full branches.

## What an exterior site must supply

An exterior site here is a site in a further branch of the same rooted chart. It has exactly one neighbor in each of the 18 old branches. Writing their cyclic coordinates as \(x_i\), normalize its tuple by \(t_i=x_i-x_0\pmod{56}\). This normalization uses the symmetry only to enumerate possibilities on the **old** branches. Let a cell \(c=(i,j,d)\) prescribe one pairwise difference \(d=t_j-t_i\). Of the 56 differences for each \(i<j\), 17 are already reached by the direct matching or a two-step route through one of the other 16 old branches. The 39 remaining differences each describe 56 ordered physical site pairs, one for each absolute starting slot. Their midpoints must be distinct for a fixed site pair.

Let \(m_t\) count all exterior sites using tuple \(t\). Each site occupies one cell per old branch pair. In any partial extension, occupancy of a cell is **at most 56**: every site charged to that cell uses a distinct one of its 56 physical pairs, by girth. In a full completion every cell has occupancy exactly 56 and there are \(39\cdot56=2184\) exterior sites. Different absolute shifts of one normalized tuple may be used. The argument imposes no probabilistic, integrality, or equivariance assumption on the nonnegative weights \(m_t\).

## A certificate of 39 small integers

Enumerate all normalized tuples from the seed. For each of the \(153\cdot39=5967\) cells, let \(f_c\) denote the number of legal tuples that occupy it. The frequencies range from 14 to 152. Give each cell the following integer weight according to \(\lfloor f_c/4\rfloor\), indexed from zero:

```
-17 -17 -17 -15 -17 -17 -3 -3 -2 -2 -3 -1 0 3 2 3 3 7 4 8
  9  10   2  11  17   6  6  5 16  3 10  1 2 -7 10 17 10 17 4
```

The sum over all cells is 19,598. The sum over the 153 cells of any legal tuple is at least 503 and at most 956. Thus a full completion would require a total score of \(56\cdot19598=1097488\), while its 2,184 sites must score at least \(2184\cdot503=1098552\). The strict excess is 1,064. Equivalently the exactly uniform target is outside the convex hull of the 2,479 tuple incidence vectors. This explains why checking every difference separately misses the obstruction.

## At most 1,604 exterior sites

The sparse [certificate data](seed18_balance_weights.txt) give a nonnegative weight to each cell. The sum of these weights is 2,865,315. Every legal tuple costs at least 99,994 under those weights. Since each cell can be used at most 56 times, the total number of exterior sites is at most
\[
 \left\lfloor\frac{56\cdot2865315}{99994}\right\rfloor=1604.
\]
This leaves at least 580 of the required 2,184 sites unfilled. Each further complete branch needs 56 sites, so **at most 28 complete branches can be added**, and at most 46 including the seed. This upper bound still allows repeated normalized tuples and arbitrary noncyclic matchings.

The same file certifies that 1,604 is the integer part of the optimum of this fractional relaxation. Give the 973 listed tuples their listed integer multiplicities; every cell then has load at most 9,999,991, and the total multiplicity is 286,532,033. Scaling by \(56/9999991\) gives a feasible fractional packing of size strictly above 1,604. The two integer inequalities bracket the optimum, without assuming that the rounded weights themselves are optimal or that an integer packing of that size exists. Floating-point linear programming suggested the weights; all claimed bounds use only the displayed integral comparisons.

The [integer verifier](seed_balance.py) reads all normalized words from the seed enumerator. The upper weights in each branch-pair row use the 39 allowed differences in increasing order. The lower certificate indexes the lexicographically sorted words starting at zero. Run in this directory:

```sh
c++ -O3 -std=c++17 seed_transversals.cpp -o /tmp/moore57-transversals
/tmp/moore57-transversals --words < seed18.txt | python3 -B seed_balance.py
```

[Stored seed](seed18.txt) · [Exhaustive legal-word enumerator](seed_transversals.cpp)
