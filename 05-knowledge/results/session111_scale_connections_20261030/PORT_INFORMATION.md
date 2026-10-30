# Session111 checkpoint: 25 cycles, 50 free bits, 100 and 200 constraints

**PROVED conditional on the inherited socle-shift construction; FINITE-EXACT literal control.** No literature-priority claim.

## Inherited map

The concurrent [ramified symmetry theorem S3b](https://github.com/eliottcassidy2000/math/blob/5a6d7be9b8b1e7939bb6a016f93b0fb07620ed13/05-knowledge/results/ramified_symmetry_spectral_20261030/PROOFS.md) proves that every simply ramified clock-invariant bank on R_k=O/(p^k) has automorphisms
sigma_h(x)=x+eta h(rho(x)),
where eta=g^(p^(k−1))−1, rho:R_k→Fp, and h:Fp→Fp is arbitrary.
The independent fibre shifts form C_p^p. This source result concerns the literal arcs, not a symmetry inferred from numerical spectra.

A nonzero shift on one residue fibre has p^(2k−1) moved vertices, grouped in L=p^(2k−2) cycles of length p. Its permutation index is D=(p−1)L.

## Linear constraint theorem

For r binary named ports, each equality between a vertex word and its shifted word is a linear constraint over F2. The full constraint operator is (I−sigma)⊗I_r. Its rank is rD, because each p-cycle contributes p−1 independent equality constraints per bit.

At p=5,k=2, the active fibre has125 vertices and 25 cycles. One port has125 input bits and 25 free bits:100 constraints. Two independent ports have250 input bits and 50 free bits:200 constraints. Here25,50,100,200 refer to dimensions of explicitly related spaces. They are not interchangeable with vertex counts, group orders or the conservative constant200 in Session110's spectral theorem.

The fixed coloring space on a cycle is the repetition code. The quotient by it is the difference-information space. Thus the apparent missing information is an exact linear rank, not a heuristic entropy analogy.

## Joint law, not separate marginal laws

Suppose vertex words are IID from an arbitrary finite joint distribution pi on the r-bit palette. Write c_p=sum_a pi(a)^p. A fibre shift fixes all words with probability
u=c_p^L.
The p fibres use disjoint input vertices, so
E |Stab_(C_p^p)|=(1+(p−1)u)^p,
P(Stab_(C_p^p)=1)=(1−u)^p.
More generally for every real s,
E |Stab_(C_p^p)|^s=(1+(p^s−1)u)^p.

For uniform independent binary ports, u=2^(−rD). For two identical fair bits, the support has size 2 rather than 4 and u=2^(−D). The second named port adds no independent constraint.

A sharper hostile: three independent fair bits and the uniform even-parity distribution on {000,011,101,110} have identical one- and two-port marginals. Nevertheless, at p=5,k=2 their fixation probabilities are2^(−300) and 2^(−200), respectively. Pairwise summaries miss a full100-constraint distinction.

These formulas count stabilizers in the specified fibre subgroup. They do not automatically equal full graph-stabilizer probabilities. S3b identifies the full group for one pinned625 source, and its larger group includes additional permutations.

## General code-valued ports

If the joint word is uniform on a binary linear code C of dimension d, then c_p=2^{−d(p−1)}, hence u=2^{−dD}. The length r of the named word is irrelevant to this probability beyond the dimension and distribution actually supplied.

This can be strengthened as a hostile family: the full r-bit cube and its even-parity subspace have identical every proper-coordinate marginal, but their exponents differ by D. Thus bounded-order port observations can miss information that a single legal symmetry test detects. This is a concrete instance of the repository's joint-port information boundary.

The independence assumption is across vertices. If words at different vertices are correlated, c_p^L need not hold; an everywhere-constant word is a compulsory counterexample.

## Exact verification

[ports.js](ports.js) builds R=(Z/25)[mu]/(mu²−mu−1), g=mu⁴=(2,3), eta=g⁵−1=(5,15). For all five fibre-shift generators it verifies the ordered difference stays in the same literal g-orbit for every unordered pair:975,000 checks. This proves those generators preserve every choice of bank orientation without enumerating all bank states.

It independently computes binary matrix ranks100,200,300 and compares the triple laws' pairwise marginals. It does not enumerate all 625-vertex graph automorphisms.

Reproduce with a JavaScript runtime supporting BigInt: node ports.js.
