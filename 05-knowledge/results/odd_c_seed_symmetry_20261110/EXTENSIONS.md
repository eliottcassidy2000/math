# From qualifying paths to actual Hamiltonian counts

**PROVED IN NOTE, elementary refinement of inherited path-cover laws.
The target-existence step remains OPEN.**

## The exact extension compiler

Fix a labelled tournament T on n>=1 old vertices and k>=1 named new
vertices. There are r=kn+binom(k,2) free orientation bits. For each pair
i<j involving a new vertex, x_ij=1 means i->j.

Call a permutation p of all n+k vertices *qualifying* if every consecutive
old-old pair is already an arc of T. Its remaining consecutive pairs
prescribe a partial assignment of the r bits. Let A_p be the bits demanded
to equal one and B_p those demanded to equal zero. Then

    F_T,k(x) = sum_(qualifying p)
                  product_(e in A_p) x_e product_(e in B_p) (1-x_e)

satisfies F_T,k(x)=H(T_x) at every Boolean orientation x. Each actual path
is counted once and each incompatible qualifying path contributes zero.
This is an exact compiler; no averaging, sign cancellation, or division by
a presumed common number of extensions is needed.

**Degree bound.** This multilinear polynomial has degree at most 2k,
independently of n. Delete the new vertices from p. The old vertices form
j nonempty runs, 1<=j<=k+1. Exactly j+k-1 adjacency requirements involve
a new vertex, hence each cylinder has degree j+k-1<=2k. Equivalently,
each of the k new vertices touches at most two path adjacencies.
The two-degree one-birth law is inherited from THM-4114; the bound here
states the same mechanism for every k. Exact degree four already occurs
for k=2,n=3 in the checked families. No literature priority is claimed.

## Why an aggregate qualifying count is insufficient

Let P_j(T) count ordered j-tuples of disjoint nonempty directed paths
covering the old vertices. Choose their j positions among the k+1 gaps
of an ordering of the new vertices. The number of qualifying permutations
with j old runs is

    k! binom(k+1,j) P_j(T).

Each such permutation is a Hamiltonian path in exactly
2^(r-j-k+1) extensions. Thus summing the cylinder formula recovers the
**KNOWN** extension mean of
[THM-5126 C](../../../01-canon/theorems/THM-5126-unit-offset-rays-two-adic-lines-isolation-of-7-and-21-extension-mean-law.md):

    mean H = k! 2^(1-k) sum_j binom(k+1,j) 2^(-j) P_j(T).

The multiplier depends on j. More fundamentally, different qualifying
permutations demand overlapping and sometimes contradictory orientations.
An arbitrary desired subcollection of N qualifying permutations cannot
be declared to be the whole Hamiltonian-path family of one tournament.
Even compatible demands can force additional paths.

A three-vertex example makes the latter defect explicit. With old vertex
0 and new vertices 1,2, select the qualifying words 012 and 120. Together
they force 0->1->2->0, which also forces word 201. The selected count two
becomes the actual count three. Positivity or compatibility alone does not
solve the exact-total problem.

The faithful existence object is the support of

    Z_T,k(z) = sum_(x in {0,1}^r) z^(F_T,k(x)).

The coefficient [z^N] is positive iff a k-vertex extension has H=N.
The mean retains only Z'_T,k(1); it loses this support. Constructing a
Boolean assignment with F=N, while controlling every forced extra path,
is precisely the missing lifting step in this approach.

## Ordered covers need endpoint labels

The supplied Chen-Mogielnicki-Ono paper, Theorem 1.1, keeps the canonical
order of record blocks. Its N(1,2)=20 and N(2,1)=10 are different. This
suggests retaining ordered boundary states, but simply replacing path
partitions by compositions is insufficient here.

If P_alpha(T) counts unrestricted ordered covering paths with lengths
alpha=(a1,...,aj), permuting the paths is a bijection, so P_alpha depends
only on the multiset of the a_i. Such a noncommutative encoding would
still have exchangeable coefficients. The asymmetry arises only after
retaining attachment roles or endpoint restrictions.

For the transitive chain 0->1->2, unrestricted two-path covers of lengths
(1,2) and (2,1) each number three. Require the second path to start at
vertex0: the counts become two and one. `verify.py` enumerates these
covers and all alternating permutations of order six independently.
For a general extension, retain the ordered new-vertex word, the old path
runs, and their endpoints. These determine every cylinder requirement.
Block sizes alone do not.

## A concrete sufficient Chan-Pak induction

The incoming, audited
[decomposition16 induction J4](../decomposition16_odd_spectrum_20261110/induction/PROOF.md)
already gives a particularly sharp target:

> For every realised odd h>=95, find a one-vertex extension from some
> realiser of h to 3h+2, and another from some realiser of h to 3h+4.

This is **OPEN**. Together with the known initial values through 287 and
ordinal-join multiplication by3, it proves the full conjecture by strong
induction: for odd N>=289, choose delta in {0,2,4} so that
h=(N-delta)/3 is odd; then 95<=h<N. The two source realisers may differ.

The incoming [55/105 source classification](../chanpak_55_105_golden_20261110/README.md)
makes this qualification essential: at h=105 no single source supports
both 317 and 319 at any order. The weaker rule allowing different sources
survives. Its concrete reset 105->19->105->9->105 changes the source while
retaining six original vertices; it does not establish a universal reset.

The missing source pair (93,281) prevents lowering this proposed uniform
threshold below 95. The same inherited packet proves 13->27 impossible
for a single induced vertex birth at **every order**, although both
13 and 27 are Hamiltonian counts.

In the compiler above, this becomes the demand for actual Boolean
solutions of F_T,1=3h+2 and F_T',1=3h+4, allowing T,T' to range over
realisers of h. Larger k gives more attachment freedom but requires a
different proved descent rule; a larger mean does not discharge it.

## Consequence for seed searches

A selected base (B,O) is a unit seed exactly when

    H(B)-W_B(O)=1,       U_B(O)=c*2^j.

These two integer equations retain the offset and the slope. This
session's search solves them on literal one-vertex extensions of known
deletion cores. More generally, fix an old base and let k added vertices
and a new root have arbitrary attachments; fix the auxiliary transitive
tails. If F1 and F2 are the resulting literal path counts with tails of
length one and two, respectively, then

    offset = 2 F1-F2,       slope = (F2-F1)/2.

As polynomials in the attachment bits of the k added vertices and the
root, F1,F2 have degree at most 2(k+1): every variable edge is incident to
one of those k+1 vertices. Thus seed construction has a bounded-degree
Boolean feasibility formulation with exact endpoint content. It has no
proved solution for every odd c. Even universal seed existence would
still leave the points below the rays.

## Checks and loss ledger

The external target is Conjecture 6.2 in
[Chan and Pak, section 6.5](https://www.math.ucla.edu/~pak/papers/DominoPermanent16.pdf),
which attributes the 2016 conjecture to bof and Gordon Royle. This primary
statement was checked directly; none of the extension-existence claims
is imported from that conjecture.

`verify.py` compares the full cylinder polynomial to independent subset
path counts for 2,130 extensions in 84 base families, including every
three-vertex base with two new vertices. Boolean finite differences
verify the degree bound, including a degree-four positive control.
It recovers the inherited exact mean using rational arithmetic.

Source: actual ordered paths and attachments. Target: Boolean cylinder
polynomial. Preserved: every labelled path and every exact total.
Lost on averaging: shared-edge incidence and target support. Restoration:
the complete bit-demand sets, not only P_j or their first moment.
