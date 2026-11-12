# Coupled path-edit correction: keep correlation before taking absolute values

Status: PROVED algebraic identities; FINITE-EXACT controls. No global coverage claim.

## Inheritance and advance

THM-116 (`01-canon/theorems/THM-116-arc-flip-contraction-reduction.md`)
already gives the single-arc contraction difference. THM-517
(`01-canon/theorems/THM-517-ocf-mayer-cluster-expansion.md`) already gives
formal hard-core cumulants and explicitly warns against evaluating a divergent
Taylor series at fugacity two. THM-4757
(`01-canon/theorems/THM-4757-linear-forest-loop-expansion-of-hamiltonian-paths.md`)
already compiles signed forced-arc insertions through directed linear forests;
the contraction mechanism below is inherited/classical, not a new forest identity. We do not rename either as new. Here we use a
finite, terminating coupled inclusion-exclusion ledger for simultaneous physical
arc changes. Its purpose is a correction certificate that cancels shared error
before estimating it. Inclusion-exclusion/Bonferroni themselves are standard.

## Definitions and exact correction theorem

Let T be a finite tournament, E a set of its arcs to reverse, T' the result,
and U the digraph containing every arc in T or T'. In particular both directions
are present on edited pairs; U is not claimed to be a tournament. A directed
Hamilton path p in U uses a(p) original edited arcs and b(p) reversed edited
arcs. Put N_ab=#{p:a(p)=a,b(p)=b}. Necessarily a+b<=n-1. Set

    Q(x,y)=sum N_ab x^a y^b.
    H(T)=Q(1,0), H(T')=Q(0,1).

Thus the desired adjustment is supported only on the two axes: all paths using
both kinds of edited arcs disappear from the exact difference. The joint state
is not replaceable by its total mass.

Define M_j^a=sum N_ab binom(a,j), M_j^b=sum N_ab binom(b,j) and

    D_r=sum_{j=0}^r (-1)^j (M_j^a-M_j^b).
    w_r(k)=binom(k-1,r) if k>r, and 0 otherwise.

Then, for every r>=0,

    D_r - (H(T')-H(T))
      =(-1)^r sum_ab N_ab (w_r(a)-w_r(b)).                 (1)

Proof: for k>0, sum_{j=0}^r(-1)^j binom(k,j)=(-1)^r
binom(k-1,r), with zero binomials past k; for k=0 the sum is1.
Apply to each actual path. This is an all-order finite identity, not a
convergence assertion. At r=n-1 it terminates exactly.

The correlation-preserving certified error bound is

    |D_r-Delta| <= sum_ab N_ab |w_r(a)-w_r(b)|,             (2)

which is at most the separate-tail bound
sum N_ab(w_r(a)+w_r(b)). The exact saving in that triangle bound is
2 sum N_ab min(w_r(a),w_r(b)). Every diagonal cell a=b cancels at every
truncation order. Bound (2) is sharp whenever the nonzero cell differences all
have the same sign. It is not asserted sharp for every tournament.

More economical partial information also gives certificates: any certified
upper/lower bounds on the signed cell sums in (1) propagate directly to Delta.
With only marginals, ordinary alternating Bonferroni bounds give an interval
for each H, and subtraction gives a valid, generally wider interval for Delta.
The distinction is exact: correlated error can be smaller than either error.

## Contracting the moments without enumerating Hamilton paths

M_j^a is the sum, over j-subsets S of original edited arcs, of the number C(S)
of U-Hamilton paths containing S. Likewise for b. If S has an in-degree or
out-degree exceeding1, or a directed cycle, C(S)=0. Otherwise its nontrivial
components are disjoint directed paths. Contract each to a block, retaining
its ordered first/last vertices. Give a block-to-block arc precisely when U
has an arc from the last vertex of the first block to the first of the second.
Then C(S) is the Hamilton-path count of this contracted digraph.

Proof: paths containing the prescribed arcs are exactly concatenations of
these blocks. The endpoint interface is essential; arbitrary graph contraction
is incorrect. Parallel or bidirectional arcs do not create multiplicity: the
underlying adjacency remains binary. Singletons are retained. A contracted
one-vertex graph has one Hamilton path. This supplies an exact compiler for
low-order corrections, whose complexity depends on edit budget and residual
contraction sizes; no polynomial-time claim is made.

## A cancellation example

For the four-row tournament 0001 / 1010 / 1001 / 0100, reverse
(1,2),(2,0),(1,0),(3,1), with indices starting at0. Counts change5 to1,
so Delta=-4. The joint nonzero cells are

    (a,b,N): (0,3,1),(1,1,4),(1,2,1),(2,0,4),(2,1,1),(3,0,1).

D_1=-8: first-order correction overshoots by-4. D_2=-4 is already exact,
although each of the two separate order-two Bonferroni approximations still
has error1. The two cubic tails cancel. This is the smallest possible vertex
order allowing a nonzero cubic shared tail. We claim no minimality among all
other measures or examples. This cancellation is invisible to independent
absolute-error estimates and visible immediately in the paired ledger.

## Consequence and limitation for Chan–Pak

Given a physical source and an actual edit menu, the ledger certifies the
absolute correction, not just a residue or asymptotic mean. It lets a search
prioritize menus whose diagonal/shared tails cancel and compute only surviving
mixed-interface contractions. But obtaining an edit menu that realizes every
prescribed target remains a separate construction problem. Neither the
ability to compute Delta nor the existence of some pair for every even gap
anchors an arbitrary odd target. This is a reusable local correction compiler,
not a proof of Chan–Pak.

## Reproduction

`python3 05-knowledge/results/session278_signed_corrections_20261112/edits/experiment.py`

The script uses seed278, 255 tournament/edit inputs of orders3..7, independently
compares direct U-path inventories with subset dynamic programming on T,T',
and verifies (1) at every truncation order. OUTPUT.json records eight concrete
nontrivial corrections. All proof identities are unrestricted in order; only
these computational controls have the stated finite scope.
