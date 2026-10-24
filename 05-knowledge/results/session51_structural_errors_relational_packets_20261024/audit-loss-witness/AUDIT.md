# Independent audit of the symmetry-decoder loss witness

Status: PASS, finite-exact. The final producer manifest is
`dfb2fadd452327a7681c3338b456e9caf335a9a8a91f436b2abf61e3d5abdaf8`.
All six payloads were verified and physically copied before execution.

The two supplied ordinary tournaments are strong, have the exact same actual
stable partition into two consecutive cyclic five-cells, and have degree-two
crossings. Independent arc-preserving backtracking finds only the identity
in each full automorphism group. Each has ten fixed vertices in five-cells,
so the defining refinement-error sum is40.

More strongly than agreement of abstract group orders, the entire JSON output
of the frozen cyclic-cell decoder is identical on the two tournaments:
stable cells, recovered cyclic coordinates, crossing gain/pin graph, frozen
cells, free components, generators, vertex-orbit sizes, group order and error.
There is exactly one rigid crossing edge and no free component.

An independently written DP distinguished by initial vertex retains the full
10-by-10 endpoint-count matrix. A separate literal depth-first traversal of
all directed simple path prefixes counts every Hamiltonian path, with no
memoization or orbit quotient. Both algorithms agree on total and terminal-
vertex counts. They give

    H(A)=15465, H(B)=15549, difference84.

The endpoint matrices and terminal counts are preserved in normal.json.
Converse and random relabelling controls preserve the DP totals. The full
ordinary automorphism groups are trivial, so these are also the ordinary
Hamiltonian-path fibre masses H/|Aut|. The Hamiltonian path action of Aut is
free: an automorphism fixing the ordered path fixes every vertex.

The differing invariant proves the tournaments are nonisomorphic and that
this symmetry decoder is not a lossless tournament or Hamiltonian-path-fibre
decoder. It also cannot reconstruct a full metagraph object that retains
those masses. It does not alone prove that two particular bare unweighted
adjacency graphs differ. The producer's final wording explicitly respects
this distinction after the audit raised it. No minimal-order or novelty
claim is made, and no forbidden-count conclusion is inferred from84.

Normal and optimized producer outputs match their pinned receipts; normal
and optimized independent outputs match byte for byte. Reproduce with

    python3 -B audit.py > normal.json
    python3 -B -O audit.py > optimized.json
    python3 -B replay.py

The final replay verifies all audit hashes, copies precisely the manifest
payloads to a fresh temporary directory, and runs both independent and
producer programs in both modes against the pinned outputs.
