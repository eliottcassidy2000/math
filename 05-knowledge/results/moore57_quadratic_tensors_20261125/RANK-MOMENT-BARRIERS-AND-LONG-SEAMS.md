# Sharp Selmer-moment limits, quantitative twist coverage, and a maximal signed seam

2026-11-26. **PROVED** reductions and height estimates; **FINITE-EXACT**
on the marked graph. The all-height rank split and the full marked
neighborhood remain open. This extends
[Half-rank averages and signed physical lifts](HEIGHT-RANK-AND-SIGNED-LIFTS.md).

## 1. Even perfect algebraic parity does not make the known moments decisive

Use the usual Bhargava–Shankar height for this section. Their first
moments at [2](https://arxiv.org/abs/1006.1002),
[3](https://arxiv.org/abs/1007.0052), [4](https://arxiv.org/abs/1312.7333)
and [5](https://arxiv.org/abs/1312.7859) are respectively 3, 4, 7 and 6.
For this test add the hypothesis that half the **algebraic** ranks are
odd. This is an extra hypothesis, not an inference from the moments
or from the root number alone.

Write b=r mod 2. The rank contribution to the 5-Selmer group gives
E[5^r] at most 6 in the limit supremum, while E[5^b] tends to 3.
For every rank,

    5^r - 5^b >= 24 * 1_{r>=2},
    5^r - 5^b >= 12 * (r-b).

For the second inequality increase r from b in steps of two. Every
step increases the left side by at least 24. Thus the tails of rank 2
and higher have limsup at most **1/8** and the mean rank has limsup at
most **3/4**. These bounds are sharp for the stated data, even if all
four average sizes are imposed simultaneously.

Here is an abstract model on 32 equiprobable points. Take no rational
torsion, the listed free rank and the listed finite group S.

| Rank | Finite S | Number of points |
|---:|---|---:|
| 0 | trivial | 4 |
| 0 | (Z/2)^2 | 4 |
| 0 | (Z/3)^2 | 4 |
| 1 | trivial | 14 |
| 1 | (Z/2)^2 | 2 |
| 2 | trivial | 3 |
| 2 | (Z/2)^2 | 1 |

Use the split group (Z/n)^r + S[n] as the n-Selmer model. Each nontrivial
S carries a perfect alternating pairing; for m divisible by n the
model lift map is reduction on the free-rank part and multiplication
by m/n on S[m]. The four size totals are exactly (96,128,224,192).
There are 16 odd ranks and four ranks equal to two; their arithmetic
mean is 24/32. This is a countermodel to an inference from **these
four means and parity**, not a proposed distribution of actual curves.
It does not incorporate every theorem about Selmer groups. For
example, an extra persistence estimate on the lifts can distinguish
it from a minimalist distribution.

The 5-Selmer bound also improves the quantitative tail transfer:
if δ is the probability of rank at least two, conditional Jensen
gives E[r; r>=2] <= δ log_5((6+o(1))/δ). This is an unconditional
uniform-integrability estimate in the usual height, but not a proof
that δ tends to zero.

## 2. How much height mass do collections of small bases capture?

Use the unweighted height max(|A|^3,|B|^2) here, and the canonical signed
squarefree parameter d with base height h from the previous note.
Let M(Y) count twist-minimal bases up to height Y. Expanding the two
floors in Möbius inversion shows

    M(Y) = (2/ζ(5)) Y^(5/6) + O(Y^(1/2)).

Indeed the floor terms sum to O(Y^(1/2) Σk^-3 + Y^(1/3) Σk^-2
+ Y^(1/6)); truncating the main Möbius series costs O(Y^(1/6)).
Removing singular and special-j models does not affect that error.
For positive squarefree n <= X^(1/6) the exact unsigned count at
|d|=n is 2M(X/n^6), with uniform error O(X^(1/2)n^-3).
Summing over n and including the limiting mass above X^(1/6)
proves total-variation convergence of the twist-label law at rate
**O(X^-1/3)**.

The same argument gives a joint refinement. Write H for the resulting
curve height, not the base height. In the limit |d| and H/X are independent,
with H/X having distribution function u^(5/6) on [0,1]. In detail,
the normalized joint count of |d|=n and H/X<=u is obtained from
2M(uX/n^6). The sum over n of absolute errors from the product law is
O(X^-1/3) uniformly in u. This is a CDF statement for H/X, not total
variation against a continuous height law. In particular the limiting
CDF of h/X is

    (ζ(10)/ζ(5)) Σ_{squarefree n>=1} min(t,n^-6)^(5/6).

More sharply, let B tend to infinity with B=o(X). Take every signed
squarefree twist of **every** base of height at most B and keep curves
of height at most X. The fraction of all X-curves covered is

    [5ζ(10)/(4ζ(2)ζ(5))] (B/X)^(2/3)
        * [1+O(B^-1/3+(B/X)^(1/12))].

The coefficient is 0.7335753721. For a base at height h the signed
squarefree count is 2(X/h)^(1/6)/ζ(2) + O((X/h)^(1/12)),
by Möbius inversion for squarefree integers. Partial summation gives
Σ_{h<=B} h^-1/6 = 5 B^(2/3)/(2ζ(5)) + O(B^(1/3))
and Σ_{h<=B} h^-1/12 = O(B^(3/4)). Divide the sum of twist counts
by (4/ζ(10))X^(5/6)+O(X^(1/2)). In particular even uniform results
over all base heights <=X^(1−ε) affect only order X^(-2ε/3) of
the full height ensemble.

There is an equivalent **budget of bases** version. Each base's
twist count is nonincreasing in its height, so among K bases the
optimal coverage takes the K smallest. If K tends to infinity and
K=o(X^(5/6)), the fraction covered by the best choice of K bases is
asymptotic to

    [5ζ(10)/(4ζ(2)ζ(5))] [ζ(5)/2]^(4/5) (K/X^(5/6))^(4/5).

Ties at the cutoff do not change the asymptotic since M(Y) is
asymptotic to a positive power of Y. This makes the fixed-base
coverage limit quantitative even for an optimally chosen collection.

## 3. The longest chordless two-phase contradiction has length 36

On the marked board, use a vertex for each two-phase option (an
exterior branch and a sheet residue mod 28). For this test a pair of
options on different branches may be adjacent if they overlap in
exactly one center residue. The edge carries the XOR demand needed
to hit opposite physical phases. A nonedge in an induced cycle must
have no overlap and the branches must all be distinct. An odd sum of
edge demands around the cycle prevents a physical lift. Deleting any
vertex breaks the cycle and leaves a path, so the obstruction is
minimal within that set. The paired options have no unary constraints.

The exact four-cycle count is **488,133 inconsistent** and
**487,750 consistent**. Count each cycle by its least vertex, its
nonadjacent diagonal vertex, and the unordered two common neighbors.
Require zero overlap across the diagonals, different exterior branches
and exactly one overlap along each edge.

An induced cycle of length L>=4 uses exactly 2L distinct residues:
one on each of the L cycle edges, used twice, and one further residue
per vertex, used once. If an edge residue were shared with any
further vertex, it would make a chord. There are just 73 double-capacity
residues accessible to the two-phase options. Therefore L<=36.
Starting from a checked length-35 cycle, freeing a 13-vertex segment
and searching for an induced 14-vertex replacement produced this
length-36 witness (branch,residue) in cyclic order:

```text
(38,16) (24,2) (48,24) (43,19) (15,19) (4,8)
(12,8) (50,18) (8,18) (33,4) (37,15) (9,15)
(23,22) (29,22) (42,11) (6,17) (18,21) (30,5)
(7,9) (32,9) (39,12) (46,13) (21,13) (44,23)
(52,23) (3,3) (20,3) (31,6) (17,6) (19,14)
(55,27) (41,27) (34,20) (26,20) (10,26) (13,16)
```

The witness has zero repeated branches, no chords, and odd signed
cycle sum. It uses (23,24,25) distinct residues in the three center
families. The replay independently reconstructs the actual 3,250-vertex
partial graph and inserts the 35 physical contacts, one by one, for
each of the 36 one-deletion cases. All keep the marker and reach source
degree 40 without introducing triangles or squares. So checking every
35-choice subset of these 36 contacts is not enough to certify a lift.
This proves **36 is the exact maximum** for this induced two-phase
obstruction class, not for obstructions using unary constraints.

There is an important limitation. The entire 36-choice witness cannot
be extended to a complete residue profile even before the phase test;
the residual capacity search exhausts after 42 nodes. Moreover every
one of the 72 oriented one-deletion lifts already leaves some outside
branch with no physical contact. The shorter 35-cycle has the same
limitation. Thus these are exact tests for local-to-global inference
on partial choices, not examples of full residue profiles that fail
to lift. For exact-cover search the on-the-fly signed cycle check
remains preferable to preloading all long cycle clauses.

## Replay

From this directory, with Node.js:

```sh
node height_moment_barriers.js
node signed_cycle_extremum.js period_fourteen_site.txt
```
