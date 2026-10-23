# Independent audit: Session43 restricted coupling

Verdict: PASS. No theorem correction requested. This audit does not claim a native-metagraph spectral theorem.

## Reproduction and independent arithmetic

The supplied checker was run normally and with Python optimization, with bytecode writes disabled. Both outputs match the supplied receipt byte for byte. All source manifest hashes and the two inherited source-provenance hashes match.

The independent `audit.py` does not import the supplied checker or its orbit code. It constructs graph and directed-tournament adjacency arrays, canonicalizes by actual vertex permutations, counts Hamilton paths by subset dynamic programming (rather than permutation path enumeration), reconstructs every incidence table, and computes exact characteristic polynomials from the opposite, tournament-side Gram product. The AB/BA nonzero-spectrum identity supplies the graph-side zero multiplicity. All 26 supplied report records, including marked actions and soft weights, match. No numerical eigenvalues enter this audit.

Independent normal and optimized receipts are byte-identical.

## Proof review

1. Observer conditioning: a transitive tournament has exactly one directed Hamiltonian order. Inside the fixed-path cube its only representative is the all-forward tournament, so the empty graph and transitive tournament give an isolated incidence block. Other words exist for n>=3. Complement pairing cannot add the absent complete graph; converse preserves the transitive class. The all-size zero-gap claim is correct.
2. At n=4 the nontransitive component has seven words, graph weights (3,2,1,1), and eigenvalues 1,8/15,0,0. Thus its gap is exactly 7/15, not an approximate numerical estimate. Both paired and unpaired checks agree.
3. Hamiltonian counts modulo 8: the independent enumeration covers every word for n=2,3,4,5. At n=5, residue 7 has mass64 and nontrivial eigenvalue89/315; residue5 has mass360 and nontrivial eigenvalue143903/1801800. No residue violation occurs earlier. These are conditioning results, not coarse-observation results.
4. Coarse observations: projection to a smaller tournament-statistic subspace can only reduce the projected norm. Together with the inherited full-cube bound, this proves the claimed at-least8/9 gap. Conditioning changes the underlying inner product and marginals, so the projection argument does not transfer.
5. Marked actions: fixing two individual vertices makes their edge bit common information in the unpaired model. With three fixed vertices, a product of two edge characters is complement-even and supplies the paired obstruction. The note correctly does not infer the paired obstruction from the two-vertex witness alone.
6. Exact PSD criterion: for graph-class f, the Dirichlet form is f^T(D_r-C D_c^-1 C^T)f and weighted variance is f^T(D_r-r r^T/|A|)f. Their comparison for all real f is exactly the displayed PSD condition, assuming only the retained positive-mass rows and columns.
7. Positive density: minimizing E_w(f-h)^2 over tournament-class h bounds the weighted conditional variance below by a times the original minimum. Minimizing E_w(f-c)^2 over constants bounds weighted variance above by b times the original variance. Dividing gives8a/(9b). The normalization of w cancels from a/b. Hard conditioning sets the lower density to zero and is outside this nonzero guarantee.
8. Observer fibre masses H(T)/|Aut(T)| follow from the free action on ordered spanning paths and agree with the separate native observer packet. Normalized converse reflects coordinates without complementing the free mask.

## Retained boundaries

The full-cube Walsh theorem is an inherited premise, whose source hashes were verified; its all-size proof was not re-audited here. The statements concern quotient-incidence coupling, not tournament metagraph adjacency. A zero global gap does not determine the gaps within all nontrivial components. The PSD test is a complete finite certificate, not an efficient all-size structural classification. Exact small-order results do not establish a universal mod8-filter law.

Commands from the workspace root:

    python -B session43-audit-coupling/audit.py
    python -B -O session43-audit-coupling/audit.py

Source reruns and independent receipts are retained beside this audit. No source files or Git state were modified.

## Publication portability

The independent checker now locates its coupling sibling relative to its own file, supporting both workspace and published names. Inherited provenance is read from the pinned publication Git object, with an optional MATH_REPOSITORY location and verified extracted-source fallback. A temporary published sibling layout was run from an unrelated working directory in both normal and optimized modes; both match the original independent receipt exactly. No mathematical source was changed.
