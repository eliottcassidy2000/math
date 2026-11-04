# A verified mixed-prime defect and forced symmetry loss

Session 190, 2026-11-04. **FINITE-EXACT independent audit at degree ten; PROVED consequences for equivariant maps and doubly stochastic correspondences.** This is a continuation of the separate `session190-bijection/NOTE.md`, not a replacement or a claim of a new global enumeration theorem.

## 1. Convention, objects, and exact source scope

For an automorphism g of a simple graph X, its orientation character is

    epsilon_X(g) = (-1)^{#{edges {u<v} of X with g(u)>g(v)}}.

The graph is even when this character is +1 for every automorphism. This is **not** the sign of the permutation of the edge set. The orientation-character convention and global equality between unlabelled tournament types and even-graph types are from [Royle, Praeger, Glasby, Freedman, and Devillers](https://doi.org/10.1007/s10801-022-01197-0), [arXiv:2204.01947](https://arxiv.org/abs/2204.01947).

Let T_n be the set of labelled tournaments and E_n the set of labelled even graphs, with S_n acting by vertex relabelling. For H <= S_n let

    A_H = {tournament types [T]: Aut(T) contains an S_n-conjugate of H},
    C_H = {even-graph types [X]: Aut(X) contains an S_n-conjugate of H}.

Write tau_T(H)=|A_H| and tau_E(H)=|C_H|. These are numbers of **unlabelled types having a subgroup**, not numbers of labelled H-fixed objects, nor counts by exact automorphism group.

The incoming pre-audit sources were read only at pinned commit `7363371836804e94f3189ea56d4a9a09fb4cf894`:

* [THM-4871](https://github.com/eliottcassidy2000/math/blob/7363371836804e94f3189ea56d4a9a09fb4cf894/01-canon/theorems/THM-4871-cpq-defect-family-of-hyp-9340-normaliser-orbit-parametrisation-disjoint-orbit-lemma-for-odd-order-groups-first-defect-at-the-prime-distance-and-the-census-refuting-hyp-9344-2.md), especially the definitions, statements 1, 3, 4, and the (3,5) census row.
* [HYP-9463](https://github.com/eliottcassidy2000/math/blob/7363371836804e94f3189ea56d4a9a09fb4cf894/05-knowledge/hypotheses/HYP-9463-cpq-defect-laws-first-defect-at-the-prime-distance-the-b-k-minus-d-formula-the-sign-and-first-order-circulant-rules.md), statements 1 and 2.
* [Producer census source](https://github.com/eliottcassidy2000/math/blob/7363371836804e94f3189ea56d4a9a09fb4cf894/04-computation/experiments/cpq_defect_20261104/cpq_defect_census.py), inspected for conventions, not executed or imported.

The metadata for commits `7363371836804e94f3189ea56d4a9a09fb4cf894` and its child `d6244c82968115e4509c0916708322165bdb2a70` was read through the Git-data commit endpoint, not a generic file-content diff. The later commit's additional computations have not been independently reproduced here. The incoming result already attributes the degree-ten obstruction to earlier repository work, so **no priority is claimed for the counterexample**.

## 2. A concrete witness and independently checked counts

On vertices 0,...,9 put

    h = (0 1 2)(3 4 5 6 7),     H = <h> ~= C15,
    P = {0,1,2},   Q = {3,4,5,6,7},   F = {8,9}.

The normalizer is

    N = Hol(C3) x Hol(C5) x S2,     |N|=240.

To verify this description, a normalizing permutation must preserve the H-orbits P,Q,F because their sizes differ. Its restrictions to P and Q normalize their regular cyclic groups and hence are affine permutations. Conversely independent affine multipliers combine by the Chinese remainder theorem to conjugate h to a power of h.

There are nine H-orbits on unordered vertex pairs: one inside P, two inside Q, one between P and Q, four between F and the two nontrivial orbits, and one inside F. None is orientation-reversing because H has odd order. Thus there are exactly 512 labelled H-fixed tournaments and 512 labelled H-fixed graphs.

The independent exhaustive audit finds:

    normalizer orbits of H-fixed tournaments                 32
    normalizer orbits of N-even H-fixed graphs               32
    unlabelled H-symmetric tournament types tau_T(H)          32
    unlabelled H-symmetric globally even graph types tau_E(H) 31

Here N-even means epsilon_X is trivial on N intersect Aut(X). It is weaker than being globally even.

The unique N-even but globally odd graph type among those 32 is K_{5,5}, with parts

    A = P union F = {0,1,2,8,9},       B = Q = {3,4,5,6,7}.

Orient every edge from A to B. Every element of N preserves both parts and hence reverses no edge: the graph is N-even. But

    sigma = (0 3)(1 4)(2 5)(6 8)(7 9)

exchanges the two parts and reverses all 25 edges. Therefore epsilon_X(sigma)=-1, so K_{5,5} is not a globally even graph.

As a further sign control, with the numeric reference orientation u->v for u<v, the same sigma reverses13 edges. Its permutation on unordered edges fixes5 edges and transposes the remaining20 in10 pairs, so its **edge-permutation sign is +1**, although its orientation character is -1. Both reference orientations give the same required parity.

All 32 N-orbit representatives on each side are pairwise non-isomorphic under the **full S10 action**. Consequently there are no multiplicity excesses at this order, and the deficit is exactly

    tau_T(H)-tau_E(H)=1.

This conclusion does not rely on the incoming paper's all-order disjoint-support lemma, its claimed general multiplicity formula, or its nauty output.

## 3. What the defect actually forbids

**Proposition.** There is no S10-equivariant map f:T_10 -> E_10 whose induced map on unlabelled types is injective. In particular there is no such map inducing a bijection of unlabelled types, or an orbit-bijective equivariant surjection.

*Proof.* For every labelled tournament T, equivariance gives

    Aut(T) <= Aut(f(T)),

because gT=T implies g f(T)=f(gT)=f(T). Thus the induced type map sends A_H into C_H. An injective type map would inject 32 objects into 31 objects, impossible.

The precise scope matters:

* This obstructs even a **one-way** preservation rule: an injective type correspondence cannot preserve the presence of this conjugacy class of C15 subgroups from every tournament to its matched graph.
* It does not forbid S10-equivariant maps which collapse types, such as sending every tournament to the empty graph.
* It does not refute global equinumerosity or an unrestricted bijection between unlabelled types. Every such bijection must send at least one of these 32 tournament types outside the 31 graph types carrying H.
* It does not, by itself, settle the broader question of a natural unlabelled correspondence. In particular canonical rank choices or the supported matchings constructed from the Walsh matrix remain logically possible.

## 4. Forced symmetry loss for the stochastic correspondence

The earlier session190 note constructs the canonical doubly stochastic matrix B indexed by tournament types T and even-graph types G:

    B_(T,G) = pi_T Phi_G(T)^2 / r_G.

The following deduction uses only its nonnegativity and unit row and column sums, so it applies to every doubly stochastic type correspondence, including every permutation matrix of a type bijection.

For subsets A of rows and C of columns write B(A,C) for the sum of entries in that rectangle. Then

    B(A,C-complement) - B(A-complement,C) = |A|-|C|.     (1)

Indeed the row sums give |A|=B(A,C)+B(A,C-complement), while the column sums give |C|=B(A,C)+B(A-complement,C). Subtraction proves (1).

Take A=A_H and C=C_H for the explicit C15 above. The verified defect gives the exact balance law

    B(A_H,C_H-complement) - B(A_H-complement,C_H) = 1.   (2)

In particular B(A_H,C_H-complement)>=1. If a tournament type is chosen **uniformly among these 32 H-symmetric types**, then the probability that its B-sampled graph does not have an S10-conjugate of H in its automorphism group is at least

    1/32.

Some H-symmetric tournament row therefore loses H symmetry with probability at least 1/32. This is a lower bound, not a claim that equality holds for the Walsh matrix, and not a statement about uniformly random labelled H-fixed tournaments. No order-ten Walsh matrix was computed.

More generally any verified positive fixed-subgroup defect delta supplies the universal bound delta/tau_T(H). Thus fixed-subgroup defects quantify an unavoidable symmetry trade-off for the probabilistic correspondence, instead of contradicting its two stochastic marginals.

For a deterministic type bijection, (2) states that the number of pairs crossing out of the H-symmetric strata exceeds the number crossing in by exactly one. Every supported perfect matching in the earlier note must obey this rule. The support graph itself must contain an edge from A_H to C_H-complement; otherwise its total available column capacity would be at most31 for32 rows.

## 5. Pre-audit formula issues: do not import them wholesale

These are scope checks on the pinned incoming text, separate from the independently verified degree-ten theorem.

1. THM-4871 statement 3(c) writes E_T(p)=4 t_p t_q. Its own parametrization leaves independent circulant choices on P and on the p fixed points F. Under its stated two-class multiplicity argument, the count is instead

       2 cross orientations x t_p(P) x t_q(Q) x 4 colours x t_p(F)
           = 8 t_p^2 t_q normalizer orbits,

   hence E_T(p)=4 t_p^2 t_q. The missing factor is invisible for p=3,5, where t_p=1. This is a correction to the count **conditional on that general multiplicity classification**, not an independent audit of all group-theoretic inputs.

2. HYP-9463 statement 1 appends a unique-witness/zero-excess assertion at k=|p-q|. It conflicts with its own cited THM-4871 (3,7), k=4 table, which lists W=17, E_T=64, E_even=48. Those numbers have not been independently rerun here; the point is that the proposed strengthening is internally inconsistent even with the packet's declared evidence. The proposed exact defect-one law must not be conflated with uniqueness of its contributing witnesses.

3. If j=k-|p-q|, all tournament coincidences at level j are the stated E_T^F(j), and R(j)=2 t_p t_q b_j(4), then the conjecture

       tau_T(k)-tau_E(k)=b_j(4)-E_T(j)

   algebraically rewrites as

       tau_E(k)=tau_T(k)-tau_T(j)+(2 t_p t_q-1)b_j(4),

   not with coefficient1 on b_j(4) except when t_p t_q=1. This correction is purely algebraic; it does not prove the conjecture.

4. The incoming statement E_T(k)=0 for k<p, by itself, gives only nonnegativity of the defect through

       tau_T(H)-tau_E(H)=W+E_even+E_wit-E_T.

   It does not prove vanishing below |p-q| without separately eliminating all three nonnegative terms. The packet's scope/status language should not be read as such a proof.

5. The pinned displayed table contains53 rows across its11 families (7+7+7+5+4+4+5+4+4+5+1), while its prose says64. Some displayed families do not reach k=|p-q|. Consequently those rows support an observed first defect at that distance only **in families where the distance was reached**. Additional source calculations may exist, but they are not supplied by that displayed table and were not audited here.

The incoming C15 values84 rather than88 at k=5 and672 at k=6 remain unreviewed source claims in this packet. They are not needed for (1), (2), or the degree-ten obstruction.

## 6. Reproduction and independent checks

Run from the enclosing workspace:

    python3 session190-cyclic/check.py > session190-cyclic/results.json
    python3 -O session190-cyclic/check.py > session190-cyclic/results_optimized.json
    cmp session190-cyclic/results.json session190-cyclic/results_optimized.json

The script uses only the standard library and imports no producer code. It explicitly generates all512 H-fixed objects on each side and all240 normalizer elements. Normalizer representatives are compared under unrestricted vertex isomorphism by exhaustive backtracking. Refinement is used only for safe rejection/pruning; a found isomorphism is checked against every adjacency entry.

Global graph evenness is checked through a full stabilizer-chain transversal: at level i, with vertices0,...,i-1 fixed, find one automorphism for every possible image of vertex i. These transversals over all levels generate the full automorphism group. The orientation character is checked on each;547 successful nonidentity transversals are tested before the first negative witness where applicable. This is not a random or sampled automorphism check.

A separate exact graph check merges vertices with identical open neighbourhoods into false-twin independent sets. The32 N-even graph representatives have at most4 such classes. All weighted quotient permutations are exhausted to classify graph types and test their full orientation characters independently. Both graph methods give32 full types,31 even. Both normal and optimized Python runs agree byte-for-byte, with explicit runtime checks retained under optimization.

No unbounded census, full audit of THM-4871, exact order-ten Walsh probabilities, or global natural bijection is claimed.
