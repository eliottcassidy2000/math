# Prime valuations and composition of tournament symmetries

**Status:** PROVED by the elementary arguments below; independently reviewed
with a full checker rerun in [the fibre agent audit](../fibre/AUDIT.md).
Exact controls are FINITE-EXACT and are distinguished from the all-order proofs.
This is a classification of the locus attaining one prime's maximum symmetry,
with an exact arithmetic formula for symmetry lost under composition. It is
not a classification of all prime tournaments or all tournament automorphism
groups.

For every odd prime p and every n >= 1,

    max over n-vertex tournaments T of v_p(|Aut(T)|)
      = v_p(n!) = (n - s_p(n))/(p-1),                         (1)

where s_p is the base-p digit sum. There is a complete structural description
of equality. For strongly connected tournaments and n >= 3 the sharper answer is

    max over strong n-vertex T of v_p(|Aut(T)|)
      = v_p(n!) - 1_{s_p(n)=2}.                             (2)

In particular the largest 3-part for any 12-vertex tournament is 243, whereas
the largest 3-part for a strong 12-vertex tournament is 81. A strong tournament
with automorphism order 81 is therefore extremal for this prime, while one
with order 27 has one further missing factor of 3.

## Inheritance and scope

The closest proved mechanisms are the canonical strong-module tree and its
prime/linear quotients in
[session56 hierarchy, Sections 1 and 2](../../session56_rmq_lca_connections_20261025/tournament-hierarchy/REPORT.md),
the cycle-index substitution in
[A000568 families](../../a000568_intransitivity_families_20261026/FAMILIES.md),
and the ambient-Sylow normalizer count in
[THM-4721, Theorem S](../../../../01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md).
The published decomposition input was checked in McConnell and de Montgolfier,
[Linear-time modular decomposition of directed graphs, Sections 2 and 4](https://www.cs.colostate.edu/pubserv/pubs/McConnell-rmm-digraphDecomp.pdf).
The proofs below derive the arithmetic and Sylow fixed-family descriptions;
they do not attribute these extensions to that paper.

The canonical hostile is loss of the child-type stabilizer when a quotient
rotation no longer preserves the attached children. The corrected near miss is
simultaneous prime saturation: the maxima in (1) are separate maxima, and
[the current uniform equivariance obstruction](../../uniform_equivariance_obstruction_20261020/UNIFORM_EQUIVARIANCE_OBSTRUCTION.md)
proves that, for n >= 5, no tournament carries both full ambient Sylow-3 and
Sylow-5 subgroups. The least-used relevant sidecar is the actual permutation
action and cycle index, distinguished from its order by
[THM-4764, modular cancellation depth](../../../../01-canon/theorems/THM-4764-modular-cancellation-depth-of-odd-regular-gap-actions-and-the-fifty-vertex-realization.md).

The anchor is exact modular composition, the niche is prime-specific symmetry
extremality, and the wildcard is base-p carrying as an obstruction to strong
connectivity. The live concepts are the strong-module tree, child isomorphism
colors, ambient Sylow action, digit carries, and action-sensitive observers.
The carry calculation connects the first four; the equal-order hostile below
shows why the fifth cannot be discarded.

## The exact local automorphism formula

Let the children of a canonical strong-module node be nonempty tournaments
T_1,...,T_m of sizes n_1,...,n_m, and let Q be their quotient tournament.
Write

    H = {g in Aut(Q): T_i is isomorphic to T_{g(i)} for every i}.

The children are canonical, so every automorphism permutes them. Restriction
to their interiors gives the exact sequence

    1 -> product_i Aut(T_i) -> Aut(T) -> H -> 1.             (3)

The right arrow is onto: choose a common model and identifications for every
child isomorphism class, and lift every color-preserving quotient
automorphism through those identifications. These choices also give a split
extension, although the splitting is not canonical. In particular,

    |Aut(T)| = |H| product_i |Aut(T_i)|.                    (4)

At a linear node H is trivial. At a prime node it is the automorphism group
of the quotient **with the complete child isomorphism type as a color**.
Equal child size or equal child group order is not enough. Formula (3) is
asserted for the canonical partition, or another partition proved invariant
under every automorphism, not for an arbitrary chosen module partition.

Define the nonnegative p-deficit

    delta_p(T) = v_p(|V(T)|!) - v_p(|Aut(T)|).

Equation (4) gives the exact recursion

    delta_p(T) = sum_i delta_p(T_i) + c_p(Q; T_1,...,T_m),

    c_p = v_p(n! / product_i n_i!) - v_p(|H|) >= 0.         (5)

For the inequality, H acts freely on the ordered partitions of an n-element
set into blocks of the indicated sizes, by reindexing blocks. The action is
well-defined because H preserves sizes. A fixed partition would require two
different nonempty disjoint blocks to be equal, so only the identity fixes
one. Thus |H| divides the displayed multinomial coefficient. This also proves
the inequality without assuming a formula for the prime factorization of H.

Consequently a canonical tree is p-saturated (delta_p = 0) if and only if
every local c_p is zero. More generally delta_p(T) is the sum of its local
c_p values over all internal nodes. This is an exact accounting identity;
computing H still requires the child types and quotient action.

## The carry terms

The elementary factorial identity

    v_p(n!) = sum_{j>=1} floor(n/p^j) = (n-s_p(n))/(p-1)

follows by counting the multiples of each power of p among 1,...,n. It yields

    v_p(n! / product_i n_i!)
      = (sum_i s_p(n_i) - s_p(n))/(p-1).                   (6)

For two summands this is the number of base-p carries when adding their
sizes. For several summands it is the total number of carry units; a column
carrying more than one unit contributes its full multiplicity.

For an ordered sum A -> B, there is no quotient symmetry, so

    delta_p(A -> B) = delta_p(A) + delta_p(B)
                     + v_p(binomial(|A|+|B|, |A|)).        (7)

For a circulant tournament Q on p vertices and any nonempty T of size m,

    delta_p(Q[T,...,T]) = p delta_p(T) + s_p(m) - 1.       (8)

To justify the full group formula here, every vertex-transitive tournament
of prime order is prime: its canonical maximal proper strong modules form
an invariant partition into equal sizes dividing p, and its quotient cannot
be a nontrivial transitive tournament while the action is transitive.
Thus its maximal proper strong modules are singletons. Substitution into a
prime quotient preserves the substituted blocks as the maximal proper strong
modules, by the inherited decomposition lemma. Hence (4) applies. The
circulant p-cycle gives v_p(|Aut(Q)|)=1, with the reverse bound from Aut(Q)
being a subgroup of S_p. Substitute (6) into (5) to obtain (8).

It follows that saturated children and a saturated p-vertex quotient produce
a saturated homogeneous substitution **exactly when m is a power of p**.
For example C3[TT2,TT2,TT2] has six vertices and automorphism order 3.
Its children and quotient are 3-saturated, but its total deficit is 1,
because s_3(2)-1 = 1. This is the first failed implication of the tempting
claim that saturation composes without a size condition.

## Classification of the p-saturated locus

Write n = sum_k a_k p^k with 0 <= a_k < p, and let r = sum_k a_k = s_p(n).
Define W_0 = 1 and, on p^k points,

    W_k = W_{k-1}^p semidirect C_p,

where C_p cyclically permutes p blocks. Its order is p^{(p^k-1)/(p-1)}.
On n points take the direct product P_n of a_k copies of W_k, with disjoint
supports. Its order is p^{(n-s_p(n))/(p-1)}, so it is a Sylow p-subgroup
of S_n.

**Theorem.** An n-vertex tournament is p-saturated if and only if, after a
relabelling, it has the following form:

1. There are a_k distinguished blocks of size p^k for every k, each a module.
2. The quotient on these r blocks is an arbitrary tournament.
3. A block of size p^k is a lexicographic tower of k circulant p-vertex
   tournaments, one possibly different circulant tournament at each level.
   The same lower tower is attached to all p children at any level.

This describes a Sylow-adapted presentation. It need not be the canonical
strong-module tree and need not be unique up to changing the Sylow subgroup.

**Proof.** A p-saturated automorphism group contains an ambient Sylow
p-subgroup. By Sylow conjugacy in S_n, relabel to make it contain P_n.
Between two different P_n-orbits, the independent transitive direct factors
force every cross arc to be uniform. Each orbit is therefore a module.
Within a p^k orbit the independent W_{k-1} factors similarly force uniform
arcs between the p children. The top p-cycle forces their induced
tournaments to be identical in the standard coordinates and their quotient
to be circulant. Apply induction inside one child. Conversely every
displayed tournament is fixed by P_n, hence its automorphism group has
the full allowable p-part. QED.

For a **fixed standard labelled subgroup** P_n the number of fixed tournaments
is therefore

    |Fix_T(P_n)| = 2^{binom(s_p(n),2)
                       + (p-1)/2 * sum_k k a_k}.           (9)

There is one bit per pair of distinct Sylow orbit blocks and (p-1)/2 bits
per tower level in each orbit. These are independent choices. This is not
an isomorphism-type count: the normalizer orbit quotient, including its
decoration action, is still required. THM-4721 supplies that inherited
normalizer reduction; (9) supplies the explicit fixed family.

Equation (1) now follows: Aut(T) <= S_n gives the upper bound, and the theorem
constructs an attaining tournament for every n and odd p. At p=2 all tournament
automorphism groups have valuation zero: any involution swaps some pair and
would reverse that pair's arc; Cauchy's theorem excludes even group order.

## The sharp strong-connectivity refinement

Assume n >= 3. If s_p(n)=1, a single tower is strong. If s_p(n)>=3, choose a
strong quotient on the digit blocks and any tower inside each block. The
result is strong: the quotient paths transport between any two blocks, and
a quotient directed cycle through a block supplies a path between any two
of its vertices if an internal path is missing. Strong tournaments exist
on every r>=3; orient a Hamiltonian cycle and orient the remaining pairs
arbitrarily. These choices attain the full bound (1).

If s_p(n)=2, the classification gives two nonempty modules with uniform arcs
in one direction. Every saturated tournament is consequently non-strong.
The integer p-valuation of a strong tournament is at most v_p(n!)-1.

To attain this smaller bound, start with the two digit blocks. Since n>=3,
one of them has size p^k with k>=1. Replace that block by p blocks of size
p^{k-1}, retaining the other digit block. There are now p+1 blocks. Choose
any strong quotient on them and put standard towers inside them. The
independent block actions form a subgroup whose valuation is

    sum over new blocks (|block|-1)/(p-1) = v_p(n!)-1,

since replacing W_k by W_{k-1}^p loses exactly one factor p. The preceding
upper bound makes this exact even if the final tournament has additional
automorphisms. This proves (2).

At n=12, the base-3 expansion is 110, so a saturated tournament has a
9-vertex tower and a 3-vertex tower with one directed edge between them.
Its automorphism order is 81*3=243 and it is non-strong. The explicit strong
tournament C3[C3 -> 1, C3 -> 1, C3 -> 1] has order 3*3^3=81 and realizes
the strong bound. No claim is made here that every strong extremizer has
this latter presentation.

## What arithmetic does and does not classify

Two 12-vertex constructions both have automorphism order 81:

    A = C3[C3,C3,C3] -> 1 -> 1 -> 1,
    B = C3 -> C3 -> C3 -> C3.

Aut(A) is the nonabelian group C3 wr C3, whereas Aut(B) is C3^4. Thus order,
size and every prime valuation together still fail to determine the group,
its permutation action, or the tournament's composition tree. This is
independent of the emblematic stronger failure in which even the entire
labelled cyclic-triple support is the same while the orders differ.

Likewise, being a prime tournament does not imply having a simple
automorphism group. The prime Paley tournament on seven vertices, with
arcs i -> j when j-i lies in {1,2,4} modulo 7, has automorphism group of
order 21. Its translation C7 is normal and proper: direct conjugation by
the multipliers 1,2,4 exhibits C7 semidirect C3, and the exact arc census
below rules out additional automorphisms. The appropriate composition
analogy retains prime tournament quotients and their actions, not merely
abstract simple-group factors.

The random historical reflection
[Sylvester Mersenne lossless](../../../../07-reflections/sylvester-mersenne-lossless.md)
suggests a useful comparison without importing its numerology. Distinct
pairwise coprime integer bases give an injective exponent code by unique
prime factorization. In contrast, (4) adds valuations of many different
nodes using the same primes: the map from a colored decomposition tree to
its valuation vector deliberately forgets where the factors came from.
The A/B example is the cheapest exact hostile to injectivity. The sidecar
that repairs reconstruction is the complete quotient-labeled tree, its
child types, and the induced actions; the valuation vector alone cannot.

THM-4764 supplies a distinct observed consequence: regular cyclic and
elementary-abelian actions of order 25 produce different subset-orbit
polynomials and modular cancellation depths. The transfer here is from a
decomposition node to its actual finite permutation group H and then to
the cycle-index/subset-orbit observer. Equations (4)-(9) preserve order
and its prime valuations; they destroy cycle structure. Retaining the
permutation action restores the input required by that observer. No
modular-depth formula is asserted for an unspecified action.

## Reproduction and exact controls

Run from the repository root:

```text
python3 -B 05-knowledge/results/tournament_composition_classification_20261026/arithmetic/check_arithmetic.py --out 05-knowledge/results/tournament_composition_classification_20261026/arithmetic/controls-normal.json
python3 -B -O 05-knowledge/results/tournament_composition_classification_20261026/arithmetic/check_arithmetic.py --out 05-knowledge/results/tournament_composition_classification_20261026/arithmetic/controls-optimized.json
```

The stdlib verifier has no inherited census filters. It checks signed pair
orbits directly from generators, independently of the digit formula, for
every n=1,...,40 and p in {3,5,7,11}. It enumerates all fixed tournaments
when n<=12 and the fixed-cube dimension is at most 12, totaling 8,930
objects; larger cases receive three specified controls. Every enumerated
fixed tournament is checked against the tower/module description.

It also enumerates every labelled tournament through five vertices and
compares a score-and-arc automorphism search with an independent enumeration
of all vertex permutations, checking the ordinary and strong maxima. The
strong attaining construction is checked through n=40 for the four primes,
including actual reachability and each claimed subgroup generator's action
on every arc. Eight explicit larger matrices have their full automorphism
groups enumerated, including the order-81 group distinction, the carry
failure, the decoration failure, the strong extremizer, and Paley7.

All tests use runtime checks that remain active under optimization. The
JSON files record matrices, group element-order histograms, finite universes
and exact outputs. The all-order classification and bounds follow from
the proofs, not from extrapolating those controls.
