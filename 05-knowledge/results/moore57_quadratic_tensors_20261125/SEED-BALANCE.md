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

## The cyclic ceiling of the same seed is exactly 20

If the new branches also use translation matchings, each new branch is the full orbit of one normalized word, labeling its sites so that its old branch-0 contact is its own slot. Two words \(u,v\) can coexist only if the 18 differences \(u_i-v_i\) are all distinct: repeated differences give a quadrangle between the two new branches and the two implicated old branches. There are 43,592 unordered word pairs passing this test.

For such a pair the shift \(z\) from the first new branch to the second must avoid every \(u_i+s_{ij}-v_j\), where \(s_{ij}\) is the old matching shift and \(s_{ii}=0\). At \(i=j\) this forbids a triangle. At distinct \(i,j\) it forbids a square using the new edge. These are all remaining short cycles involving both new branches, since each word is individually legal. Exactly 642 word pairs remain, with 648 choices of shift in total; each is an actual 20-branch chart. Their compatibility graph has no triangle, so three cyclic branches cannot be added simultaneously. This is an exact ceiling for this seed under the translation condition. Append `--cyclic` to the verifier command to reproduce it.

## Changing two branches yields 21

**FINITE-EXACT positive construction.** The distinct stored [20-branch seed](seed20.txt) can grow in order when two branches are released and replaced by three. Delete old indices 11 and 15 (zero-based), retaining the other 18 branches in order. Enumerate their normalized words; take a triangle of the compatible-pair graph just described and select its three interbranch shifts. All direct and two-step routes for the assembled chart must be tested again, including routes involving all three new branches. The resulting full upper-triangular array is [seed21.txt](seed21.txt). Its first 18 branches are literally the retained subarray of seed20. For every one of its 210 unordered branch pairs the direct route and 19 intermediate routes have 20 distinct translations. Thus the 1,176 sites form a 20-regular matching cover of girth at least five. Adjoining their 21 heads and a root realizes this as a partial degree-57 Moore chart on 1,198 vertices.

This 21-branch seed has **zero** possible normalized words for any further exterior site, including one with unrestricted new matchings. Its positive size and its nonextendibility are both certified by `/tmp/moore57-transversals < seed21.txt`. The 18-branch core retained from seed20 here is different from the original [seed18.txt](seed18.txt), whose exact cyclic ceiling is 20.

The 21-branch chart also has a short spectral obstruction, unlike the strictly bounded spectra of the three earlier stored seeds. Set \(B_{ij}=(-1)^{s_{ij}}\) for \(i\ne j\) and \(B_{ii}=0\). For the integer vector
\[
 v=(2,1,1,0,1,-1,1,2,1,0,0,1,-2,1,2,2,3,0,2,1,-1)
\]
the array gives \(v^tBv=302>7v^tv=301\). The site vector taking value \((-1)^x v_i\) at slot \(x\) of branch \(i\), zero on other branches, is centered on every fiber; its Rayleigh quotient is \(302/43\). An actual completion would impose the upper bound 7 on this subspace. Thus the 21st branch increases the girth-five subsystem size while already failing a necessary Moore test.
