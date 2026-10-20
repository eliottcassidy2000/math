# Source audit and count conventions

Checked 2026-10-20. This packet's new arguments are independent of any unproved classification conjecture.

## Published inputs

1. John D. Dixon, *The Maximum Order of the Group of a Tournament*, Canadian Mathematical Bulletin 10(4) (1967), 503-505.
   - DOI: https://doi.org/10.4153/CMB-1967-048-9
   - Primary PDF: https://www.cambridge.org/core/services/aop-cambridge-core/content/view/682531E83B0A3AA042CFC5852474E914/S0008439500067461a.pdf/div-class-title-the-maximum-order-of-the-group-of-a-tournament-div.pdf
   - Re-read the primary PDF, including Theorem 1 on page 503 and its proof. Imported statement: |Aut(T)|<=3^((n-1)/2) for every n-vertex tournament. No transitivity or regularity assumption. The proof invokes Feit-Thompson. This is a cited published input, not a newly proved group theorem.

2. Carsten Thomassen, *Hamilton Circuits in Regular Tournaments*, in *Cycles in Graphs*, North-Holland Mathematics Studies 115 / Annals of Discrete Mathematics 27 (1985), 159-162.
   - DOI: https://doi.org/10.1016/S0304-0208(08)73005-6
   - Primary publisher page: https://www.sciencedirect.com/science/article/pii/S0304020808730056
   - Rechecked the indexed primary publisher text. It explicitly gives a lower bound 4^k k! for strong tournaments with minimum outdegree at least 3k+3. The indexed exponent appears as uppercase K, but the statement has just one parameter k. The earlier packet independently corroborated lowercase k in Theorem 14.3 of Bang-Jensen and Gutin's survey.
   - Direct page/PDF fetch returned HTTP 403; indexed publisher content remained readable. The full four-page proof was not retrieved during this source check. The result is therefore used as a CITED published theorem, with the precise statement and hypotheses verified in its primary publisher abstract, not as an independently reproved or full-text-audited theorem.
   - The relevant bound is about distinct Hamiltonian circuits/cycles, not Hamiltonian paths or only regular tournaments. Our HC convention is cycles modulo cyclic starting point. Even a rooted-cycle convention would change the comparison only by a factor n and not the finiteness conclusion, but we retain the ordinary unrooted cycle convention used in the inherited audited packet.

3. Camion's theorem, inherited in both input packets: every nontrivial strong tournament has a directed Hamiltonian cycle. It is used only for existence of one cycle in the degree argument. There is no use of a Hamiltonian-cycle counting strengthening of Camion.

## Local mathematical inputs

- `../paired_black_leaf_reduction_20261020/CYCLE_CUT_NEAR_PURITY_OBSTRUCTION.md` (staged input: `session14-black-leaf-reduction/CYCLE_CUT_NEAR_PURITY_OBSTRUCTION.md`): the exact q-orbit / gcd(2,q) blue-cut calculation, reproved here for arbitrary defect and without a self-converse restriction.
- `../paired_black_leaf_reduction_20261020/BLACK_QUOTIENT_LEAF_CLASSIFICATION.md` (staged input: `session14-black-leaf-reduction/BLACK_QUOTIENT_LEAF_CLASSIFICATION.md`): exact paired-black degree identities, rigid-wrapper leaf classification, unknown center series, and the complete deficit-two exclusion through order nine. The finiteness conclusion in the new packet strictly upgrades that note's old infinite-residual boundary; it does not alter its structural proofs or turn its finite census into a classification.
- `../merged_metagraph_finite_pure_atoms_20261020/FINITE_ATOMS_AND_RIGID_CLOSURE.md` (staged input: `session14-pure-blue-finiteness/FINITE_ATOMS_AND_RIGID_CLOSURE.md`): count conventions, prior finiteness mechanism, cited inputs, and the proved strong rigid alphabet {1,C3}. Its later superseded pure-blue classification frontier is not inherited as current open status.

Input hashes are recorded in `PUBLICATION_MANIFEST.json`. The source hashes pin the precise inherited mathematical statements.
