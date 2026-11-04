# Independent audit of Session 190 group-action packet

2026-11-04. **PASS. No mathematical correction required.**

Reviewed target: `session190-groups/PROOF.md`, `verify.py`, `results.json`, `SOURCES.md`, and README. Exact reviewed hashes are in `REVIEWED_SOURCE_SHA256SUMS`. The producer packet was not modified. The audit has no repository imports and did not run a startup helper or read the forbidden ledger.

## Executive findings

1. The planar/spherical distinction is correct. The five families concern finite SO(3) subgroups; planar finite point groups are cyclic or dihedral. The convention |D_n|=2n is explicit, and the order-two degeneracy is handled.
2. The Paley7 group and its signed extension are exactly C7 semidirect C3 and AGL(1,7), orders 21 and 42. The orbit sizes 7,21,42 genuinely come from point, unordered-pair and ordered-pair actions. They are not identified with forbidden Hamiltonian-path counts.
3. Ordinary path-complement edge stabilizers are trivial. Signed stabilizers are 1 or C2, and the negative element cannot exchange the endpoints. Both the proof and finite checks are sound.
4. The all-order unsigned/signed edge-orbit formulas have the correct factors of two and the stated n>=3 boundary. The n=8 producer row is explicitly formula-only.
5. The real minimum faithful dimension six and binary minimum dimension three are correct. Maschke is used with its characteristic restriction, and the warning that reducing skew signs modulo two loses orientation is essential and correctly present.
6. The Fermat-prime assertion keeps the prime-degree affine theorem as an explicit inherited input. Trivial multiplier stabilizer alone would not prove the full automorphism statement without that input.

No claim of new literature priority or a five-to-five correspondence is justified or made.

## Proof review

### A. Classification scope

A finite plane-isometry group fixes the average of an orbit, so conjugation by a translation reduces it to O(2). The embedding O(2)->SO(3) used by the packet is A->diag(A,det A), so an odd-order group with a faithful planar action is covered by the SO(3) obstruction.

Every nonidentity involution of a finite vertex set exchanges some pair and reverses that tournament arc. Thus tournament automorphism groups have odd order. In the finite SO(3) classification, only odd cyclic groups have odd order. For a signed symmetry group with an odd-order kernel of index at most two, the only possible faithful SO(3) group types are cyclic of odd or twice-odd order, or dihedral of order twice an odd number. The exceptional A4,S4,A5 fail this condition. These are necessary realization restrictions, not claims that every permitted action is realized by the specified tournament.

The five Platonic solids are the five ordered (p,q) solutions of the stated curvature inequality; duality has three orbits. Consequently their rotation groups provide three exceptional groups, not five distinct group families.

### B. Exact edge action

For T_P obtained by reversing all nonpath arcs, the common directed arcs are exactly P. Any permutation preserving the unordered pair {T,T_P} must preserve that directed path and therefore fixes every vertex. This proves ordinary edge freeness, not just common-endpoint-stabilizer freeness.

A negative signed symmetry must reverse the common directed path, hence its underlying permutation is the unique path reversal r_P. The pair of endpoints of P is an off-path pair for n>=3. Applying r_P and global converse leaves its original arc unchanged, whereas T_P reverses it. The candidate cannot interchange T and T_P. It therefore exists exactly when it fixes T, equivalently r_P is an anti-automorphism. Since reversing twice is the identity, the stabilizer is 1 or C2.

After normalizing P, q off-path bits remain. Passing from an ordered incidence to an unordered edge quotients by the fixed-point-free bit complement, giving 2^(q-1) ordinary edge orbits. The converse involution acts by coordinate reflection, with f fixed coordinates and c=(q+f)/2 coordinate cycles. It fixes 2^c bit words, hence 2^(c-1) complement-pairs. A reflection-to-complement word is impossible because the apex coordinate is fixed. Burnside then gives precisely

    (2^(q-1)+2^(c-1))/2

signed edge orbits. The formula counts the actual carrier described in the proof. It does not count simple edges after collapsing all parallel incidences, black support, or arbitrary flip adjacency.

### C. Paley and field-sensitive representations

The 7/21/42 diagram is orbit-stabilizer for an actual affine action. The order-two stabilizer of an unordered pair is the affine reflection z->x+y-z. It exchanges the two ordered lifts, explaining why the full affine group cannot preserve one tournament orientation. Restricting multipliers to the squares removes that obstruction.

For the minimum real dimension, any faithful real representation is nontrivial on normal C7. One nontrivial character forces its order-three multiplier orbit and, by reality, the disjoint complex-conjugate orbit. Their sizes are three plus three. This is a valid lower bound for every faithful real representation, not merely orthogonal permutation models. The explicit monomial complex matrices give a faithful SU(3) model, attaining real dimension six.

For the binary model, x^7-1 splits into coprime factors of degrees 1,3,3 over F2. The multiplier two preserves both cubic character orbits and negation exchanges them. The faithful model GammaL(1,8) has order 21 inside GL(3,2), and the lower bound follows from the absence of order-seven elements in GL(1,2) and GL(2,2).

The bad-characteristic examples establish existence of nonsplit modules. They do not imply that every module is nonsplit. Over F7, an invariant functional for the seven-cycle has constant coefficients and evaluates to zero on the all-ones vector. Thus no invariant retraction to that constant line exists. The analogous cyclic modules over F2 and F3 justify the other stated controls.

## Independent numerical verification

`independent_check.py` imports neither the producer verifier nor its result file. It uses directed arc sets, Hamiltonian-path depth-first search, generator-based orbit traversal, enumeration of all binary matrix entries, modular row reduction, and symbolic monomial matrices.

Independent checks passed:

* All 2+8+64 normalized path-complement words for n=3,4,5, testing all 6,24,120 vertex permutations respectively. Ordinary stabilizers are trivial; negative stabilizers are exactly the permitted reversal; ordinary edge orbits are 1,4,32 and signed edge orbits 1,3,20.
* Every one of the 5,040 Paley7 vertex permutations, verifying the full ordinary and anti-automorphism sets.
* All 189 Paley7 Hamilton paths found by DFS; the endpoint matrix is 1 on arcs and 8 on reversed arcs. There are 189 distinct labelled path-complement neighbors and no surviving original nonidentity automorphism.
* Nine ordinary path orbits of size 21; signed path orbits 21,21,21,42,42,42, found from generators rather than applying the producer's complete group list.
* Seven anti-automorphism involutions, each fixing nine reversed paths.
* Exact integer multiplication S^2=J-7I.
* All 512 binary 3-by-3 matrices, retaining the 168 invertible ones. A separately constructed companion matrix and Frobenius matrix generate 21 elements and the complete normalizer of the companion's C7.
* All 128 binary seven-coordinate vectors for the 1+3+3 decomposition and negation exchange.
* Exact symbolic SU(3) monomial group of order 21, including PDP^-1=D^2 and determinant one.
* Modular linear-system inconsistency for invariant retractions onto constants in dimensions/characteristics 2,3,7.

Normal and optimized independent runs are byte-identical. The producer's normal and optimized `--check` runs also passed, with results SHA256 `7f4fcd05b9d1324f70716a70cbd4ac861ded89d8c5e2476cff6d698a11bbd199`. Its original SHA256 manifest passed before review.

## Primary-source verification

The source URLs and theorem numbers were independently opened and checked:

* [Etingof, Groups Around Us, Theorem 10.2](https://math.mit.edu/~etingof/groups.pdf): the SO(3) classification. Etingof writes D_(2n); the packet explicitly uses D_n of order 2n, so there is no mathematical discrepancy.
* [Arapura, Algebra, Theorem 14.1](https://www.math.purdue.edu/~arapura/algebra/algebra14.pdf): SO(2) cyclic and O(2) cyclic/dihedral, including the degenerate dihedral conventions.
* [Etingof et al., Introduction to Representation Theory, Theorem 4.1.1 and Proposition 4.1.2](https://math.mit.edu/~etingof/representationtheorybook.pdf): Maschke and its converse. The semisimplicity statement used here is valid over arbitrary fields with the stated characteristic condition; no algebraic-closedness assumption is needed for that portion.
* [Akin, Rock, Paper, Scissors, Etc., Theorems 3.13–3.14](https://arxiv.org/html/1806.11241v4): the multiplier criterion and construction from odd automorphism subgroups. These do not by themselves replace the separately acknowledged prime-degree affine theorem.

## Audit limits

The classification theorems are cited classical inputs, not reproved computationally. The all-order edge formulas and dimension bounds were checked by proof, with finite controls; finite census is not their justification. The inherited raw Hamiltonian 7/21 exclusion and its original proof were not re-audited in this assignment. No new conclusion about the full raw Hamiltonian spectrum follows.

No producer file was changed. A fresh audit is needed if the mathematical claims, verifier, or sources are materially edited after the reviewed hashes.
