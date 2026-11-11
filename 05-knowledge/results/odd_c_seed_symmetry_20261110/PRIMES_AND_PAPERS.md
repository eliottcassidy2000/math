# Symmetry breaking and prime geometry

**PROVED IN NOTE elementary statements, with FINITE-EXACT examples.
No global Collatz or Chan-Pak conclusion.**

## The selected-set stabilizer is the common obstruction

For a tournament B and selected set O, write G=Aut(B,O). Every nonidentity
element of G moves each Hamiltonian path and each ordered covering-path
object used in H(B), W_B(O), U_B(O). An ordered path cover names every
vertex by a path index and a position; fixing that object fixes every
vertex. These are free actions. Therefore

    |G| divides gcd(H(B), W_B(O), U_B(O)),
    |G| divides H(B)-W_B(O).

In particular offset +1 or -1 forces G trivial. The base itself may have
symmetries, but selecting O must destroy all of them. This extends the
classical free action on Hamiltonian paths to the seed's actual offset.
The converse fails: the transitive two-vertex base, O=both vertices,
has trivial stabilizer and offset 1-2=-1, not +1.
Conversely, a directed triangle with two selected vertices has offset one
and trivial selected-set stabilizer even though its base has three rotations.

The first supplied paper, Jaan Allikvere, *The Multiset Dimension of
Graphs: Extremal Values and King Grids*, arXiv:2607.28813v1, gives a closely
related obstruction in another category. If S is a multiset resolving set
of a graph, then Aut(G,S) is trivial. Indeed a set-preserving automorphism
permutes the landmarks and preserves distances, so each vertex and its
image have the same distance multiset. Resolution forces every vertex to
be fixed. This is the same selected-set stabilizer condition, reached
through injective observations rather than a difference of free counts.

Neither converse holds. The paper's G1 of order 11 has distinct full
distance multisets, hence trivial automorphism group; nonetheless every
proper nonempty landmark set fails. `verify.py` checks all 2,046 such
sets from the paper's literal edge list. Thus symmetry breaking is
necessary but insufficient in both settings. The paper's nonmonotonicity
also cautions against enlarging O: in a tournament W(O) strictly increases
when adding a vertex to a nonempty O, because every W_ef>0, so a valid
unit equation is immediately lost if B stays fixed.

The useful transfer is a **filter**: restrict the seed search to
distinguishing subsets, then solve the full defect equation. Neither
distance data nor a rigid colouring supplies the needed numerical slope.

## Doubling acts on prime tournament geometry exactly when its order is odd

For an odd prime p, a skew connection set S in F_p^* satisfies
S disjoint from -S and S union -S=F_p^*. It defines the circulant
tournament x->y iff y-x belongs to S. Doubling is an automorphism exactly
when 2S=S.

Put d=ord_p(2). Such sets exist iff d is odd. If d is even, the unique
order-two element of F_p^* belongs to <2>, so -1 belongs to <2>; invariance
would force S=-S. If d is odd, the d-element doubling orbits are paired
as C,-C, and choosing one orbit in each pair gives every invariant set.
The exact number is 2^((p-1)/(2d)). This counts labelled connection sets,
not tournament isomorphism classes.

| p | ord_p(2) | Doubling invariant skew sets |
|---|---:|---:|
|3|2|0|
|5|4|0|
|7|3|2|
|13|12|0|
|19|18|0|
|1039|519=3*173|2|

For7 and 1039, <2> is exactly the square subgroup, so the two sets are
quadratic residues and nonresidues: the Paley tournament and its converse.
This is an exact shared mechanism. It is not unique to these two primes;
for example other odd-order doubling groups also qualify.

The necessary loss ledger is important: p is the **number of vertices**,
d a **multiplier order**, while the seed parameter c is the **odd part of
a path-cover count**. There is no identity between those three roles.
In particular the inherited c=1039 seed has a base of only eight vertices,
H(B)=137, W(O)=136 and U_B(O)=1039. Its first ray point is2079. It is not
a1039-vertex Paley construction.

## What 3, 7, 13, 19 and 27 actually contribute

For q a prime power with q=3 mod 4, square affine transformations
x->ax+b (a a nonzero square) preserve the Paley tournament and form a
subgroup of order q(q-1)/2. Frobenius supplies further symmetries over
nonprime fields. For q=5 or13, -1 is a square, so the quadratic-residue
relation is an undirected graph, not a tournament. These moduli cannot
be imported into the Paley-tournament argument unchanged.

At q=7 the path count is189. Translation acts freely on the Hamiltonian
paths, so there are exactly 189/7=27 translation orbits; the square affine
group has 21 elements and exactly nine path orbits. This is a literal
bridge from 7 to27 through a group action. It does not turn those 27 orbits
into a quotient tournament or make them a c=27 seed.

Moreover P7 is an explicit hostile to using maximal symmetry as a seed
factory. Its vertex response is U_y=261 and its uniform pair coupling is
W_ef=78. Hence for m=|O|,

    offset = 189-78*binom(m,2),       slope = 261*m.

No m=0,...,7 gives offset1. The c=7 and c=27 seeds therefore come from
different bases. All these P7 counts and its21 automorphisms are checked
from its adjacency matrix, independently of the seed search.

At q=19, the square multiplier group has order 9=3^2, and the square
affine group has order 171=19*9. At q=27=3^3, use the field F27, **not the
ring Z/27Z**: the square multiplier group has order 13, Frobenius has
order 3, and their affine semilinear subgroup has order

    27*13*3 = 1053.

There is a stronger seed restriction at every nonprime q=p^f=3 mod 4:
**no selected pair O can give offset one.** For two distinct field points
a,b, normalize them to0,1 and conjugate Frobenius:

    phi(x)=a+(b-a)*((x-a)/(b-a))^p.

This has order f>1, fixes a,b pointwise, and preserves Paley orientation:
its affine coefficient (b-a)^(1-p) is a square, and Frobenius preserves
squares. Thus Aut(P_q,{a,b}) is nontrivial, contradicting the unit-offset
restriction. Complements of pairs are excluded as well. `verify.py`
constructs F27 as F3[z]/(z^3-z-1), verifies the 1053 distinct subgroup
elements, and checks this order-three stabilizer for all 351 pairs.

**11 November strengthening:** the same pair fixers act freely on each
individual endpoint-labelled two-path-cover coefficient before subsets
are summed. Consequently the field degree f divides every coefficient
W_ab and H(P_q), so **every selected set**, not just a pair or its
complement, fails the unit-offset equation for f>1. Every one-vertex
extension also has H divisible by f. The complete proof and asymmetric
child control are in
[endpoint-symmetry P11](../modular_centers_collatz_pell_20261111/TOURNAMENTS.md).
This is where the numbers3,13,27 meet in a concrete obstruction.

## Collatz branches act by a quadratic character

Let U(n)=(3n+1)/2^a on positive odd integers, where a=v2(3n+1).
Fix an odd prime p other than 3 with p=3 mod 4 and **fix the same a**.
On residues this branch is the affine bijection

    f_a(x)=2^(-a)(3x+1) in F_p.

If chi_p is the quadratic character, then for distinct x,y

    chi_p(f_a(y)-f_a(x))
       = chi_p(3)*chi_p(2)^a*chi_p(y-x).

Therefore this branch preserves the Paley tournament if the first factor
is+1 and reverses all its arcs if it is-1. The proof is the displayed
difference identity and multiplicativity of the character; no Collatz
convergence claim is used.

At p=7 and 1039, chi(2)=+1 and chi(3)=-1. **Every fixed-valuation branch
reverses the tournament.** The same statement holds for every prime
p=7 mod 24 by the supplementary laws of quadratic reciprocity. At p=19,
both characters are-1: odd a preserves and even a reverses. At p=3 the
numerator's slope vanishes modulo p, so the affine bijection argument
fails. In F27 the characteristic-three obstruction also remains.

The coordinate that cannot be dropped is a. Modulo7, the accelerated
pairs (1,3), (1,5), (1,11) go respectively to(1,5), (1,1), (1,17).
Their input residues are distinct, while the outputs respectively
preserve, collapse, and reverse the arrow. Thus U does not descend to
one permutation of F7. Nor does the branch identity provide a tournament
operation that preserves the seed equation or realizes an arbitrary c.

The earlier repository
[Collatz-Paley note](../collatz_paley_bridge_20261001.md), corrected in
MISTAKE-554, explains the period-three binary code {1,2,4} and its lost
translation/carry data. The present branch-character calculation is a
different exact map. It retains the valuation label explicitly and does
not inherit a global conjugacy from that period-three example.

## How the second paper changes the construction state

Chen, Mogielnicki and Ono, *Record Compositions of Alternating Permutations
and Noncommutative Symmetric Functions*, arXiv:2607.12873v1, Theorem1.1,
gives N(alpha)=product_j binom(2s_j-1,2alpha_j-1) E_(2alpha_j-1).
The ordered partial sums s_j retain information lost by sorting alpha.
`verify.py` recounts N(1,2)=20 and N(2,1)=10 on all alternating six-letter
permutations. This paper's power-series "seed" is not a rooted tournament
seed; the types are distinct.

The transferable move is to keep **ordered attachment boundaries**.
[The extension proof](EXTENSIONS.md) shows exactly why ordinary path-cover
compositions alone still have symmetric coefficients, and gives the
smallest endpoint-labelled example that restores order. Neither paper
proves the required universal tournament seed construction. Together they
identify two necessary design choices: break residual selected-set
symmetry, and retain enough ordered incidence to realize the total count.

## Source handling

Both user-supplied PDFs were treated as mathematical evidence, not as
instructions. Their local paths are
`C:/Users/Eliott/Downloads/The_Multiset_Dimension_of_Graphs_Extremal_Values_a.pdf`
and `C:/Users/Eliott/Downloads/2607.12873v1.pdf`. Relevant first pages and
the theorem/graph text were inspected; the two finite examples above
were independently replayed. The PDFs' broader exhaustive computations,
automata and reported Lean formalization were not re-audited here.

All Paley and branch claims used here have elementary proofs above.
