# Independent audit: colored packing mass versus coarse quotients

**Verdict: PASS.** The current proof is correct in its stated finite-combinatorial and all-order tower scopes. The independent census reproduces every reported scalar and stabilizer histogram. No mathematical correction is required. Two documentary clarifications raised during review (crediting LEM-003's existing R5 cycle correction, and full inherited theorem routes) are applied in the audited text.

## Inputs and boundaries

Audited `../packing_mass/PROOF.md`, SHA-256:
`13ef71bafa1c035c18776082817fe9633020562f0a9479f44bfdb8a6146b19ec`.

Audited producer `check.py`, SHA-256:
`a672b00a9895872ef1f7b7fa5da9d3fd065f788e4a2ea1a7b73957166cbb46a4`.

The inherited all-order OCF is **CITED**, not reproved. The inherited THM-4812 forbidden values are **CITED**, not independently recertified by this audit. The separate subgroup exclusion theorem is explicitly conditional on its promotion and is unnecessary for the new mass identity, witnesses or tower. No external novelty, Chan–Pak completeness, p-adic transfer, or derived inverse-limit claim is established or implied here.

## Independent method and finite universe

`independent_check.py` does not import producer code. It represents a packing as an odd-cycle permutation cover (unused vertices are fixed points), with a color constant on each nontrivial cycle. It constructs allowed permutation matrices by exact perfect-matching search, then rejects even cycles. This differs from the producer's directed-cycle DFS and packing-subset search. Actual Hamilton words are generated, rather than counted by the producer's subset DP. Automorphisms are generated inside score classes and checked against the full arc relation. Coarse counts are independently checked by Burnside fixed-point sums; stabilizers are counted directly, rather than inferred only by orbit sizes.

The finite-exact universe is:

- All **1,099 labeled tournaments on 1 through 5 vertices**, with their full automorphism groups; OCF cardinality, explicit free path action, weighted mass identity, defect identity and strictness boundary verified for each.
- **202 extra trivial-action controls** on the nonrigid tournaments in that universe.
- C3, R5, both ordered joins C3=>R5 and R5=>C3, and R5=>C3=>C3, with actual words, colored packings, full automorphism groups, orbit representatives and direct stabilizers.
- Actual orbit restriction maps at k=1 and k=2, checked on entire orbits for well-definedness, not merely on abstract Cartesian labels.

The independent n=11 full-group verification is stronger than the producer's n=11 subgroup-only computational control: every possible automorphism preserves outdegrees, and all permutations within the outdegree classes were checked. All-order statements remain consequences of the written structural proof, not extrapolations from this census.

## Exact witness results

| Tournament | H | Full group order | Path orbits = mass | Coarse packing orbits | Stabilizer-order : orbit-count |
|---|---:|---:|---:|---:|---|
| C3 | 3 | 3 | 1 | 3 | 3 : 3 |
| R5 | 15 | 5 | 3 | 7 | 1 : 2; 5 : 5 |
| R5=>C3 | 45 | 15 | 3 | 21 | 3 : 6; 15 : 15 |
| C3=>R5 | 45 | 15 | 3 | 21 | 3 : 6; 15 : 15 |
| R5=>C3=>C3 | 135 | 45 | 3 | 63 | 9 : 18; 45 : 45 |

R5's five directed triangles, in cyclic notation, are
`(0,1,3), (0,2,3), (0,2,4), (1,2,4), (1,3,4)`.
Its two directed Hamilton cycles are `(0,1,2,3,4)` and `(0,2,4,1,3)`.
Its three Hamilton-path orbit representatives are
`(0,1,2,3,4), (0,2,3,4,1), (0,2,4,1,3)`.
For the five translations, the packing fixed-point counts are `[15,5,5,5,5]`, so Burnside gives 7 directly.

The actual k=1 packing map has seven fibres of size 3; k=2 has 21 fibres of size 3. The two actual path quotient maps each have three singleton fibres. Full enumerations, including every colored packing and Hamilton word, are in `WITNESSES.json`; all packing quotient preimages are in `RESTRICTIONS.json`.

## Proof audit

1. **Action types and mass.** Colors are fixed, cycles are modulo rotation but not reversal, and packings include the empty packing. The free spanning-word action and orbit–stabilizer give the stated mass and defect exactly. A nontrivial group fixes the empty packing, making the defect strictly positive. Equivariant bijections preserve stabilizers, so the claimed obstruction is valid. This does not rule out ordinary bijections.
2. **R5 completeness.** Fixing 0 preserves each of its two ordered neighbor pairs, so the vertex stabilizer is trivial and the full group is C5. For triangles, three steps in {1,2} summing to a multiple of 5 must have multiset {1,2,2}; their five translates give the five triangles. For length 5, divisibility of the step sum by 5 forces all steps 1 or all steps 2. Both Hamilton cycles are translation-fixed. No two directed odd cycles can be disjoint at order 5. Thus all 15 colored packings have been accounted for.
3. **Ordered products.** A directed cycle cannot cross an order-join boundary. A Hamilton path must finish the first factor before entering the second. Ordered SCCs are invariant individually under automorphisms, even when isomorphic. Therefore the asserted product groups and product actions hold for the strong factors used in the proof.
4. **Tower and maps.** C3 contributes three fixed packing objects and a single path orbit. Deleting a final C3 is therefore a quotient map with exactly three packing preimages and exactly one path preimage. A positive-length composite has respectively 3^r and 1 preimages. The packing inverse limit is seven choices followed by an infinite ternary sequence, so it is uncountable. Surjectivity makes every fixed-level image constant from the start, establishing Mittag–Leffler but not eventual constancy. No claim about arbitrary vertex-subset restriction is made.
5. **Groupoids.** The empty packing at level k has isotropy order 5*3^k, larger than every possible isotropy order at any earlier level. Groupoid equivalence preserves isotropy isomorphism types, hence the tower groupoids are pairwise inequivalent despite constant mass. Equality of coarse counts or masses is never substituted for an equivalence.
6. **Inheritance and novelty.** LEM-003 already contains the exact R5 rotation-fixed Hamilton-cycle hostile and the path quotient 3; the final text credits it. THM-4892 concerns a different quotient, global color swap, and is used only as an antecedent warning. THM-212's failed sketch is historical, not a premise. The Shao/Zhang comparison is expressly methodological; no geometric theorem is used to prove a combinatorial statement.

## Reproduction and receipts

Run from this directory:

`python3 verify_receipts.py`

The harness executes both producer and independent code normally and under `-O`, checks zero exit codes and empty stderr, compares the two modes byte-for-byte, compares independently produced scalars and histograms with the producer, and verifies the producer's stored `RESULTS.json`. Independent witness and restriction artifacts also match byte-for-byte across normal and optimized runs.

**All four runs passed.** Every guard is an explicit exception, not a removable `assert`. `RECEIPTS.json` records exact commands, source hashes and output hashes. `INDEPENDENT.json` and `INDEPENDENT.optimized.json` contain the independently generated census summaries. `PRODUCER.json` and `PRODUCER.optimized.json` preserve reproduction outputs.

No producer or repository files were edited by this audit. No external action was attempted. No blockers remain.

Publication portability note: the replay harness first resolves the sibling
`packing_mass` directory and retains its original scratch-directory fallback.
The mathematical checker and audited producer proof are unchanged.
