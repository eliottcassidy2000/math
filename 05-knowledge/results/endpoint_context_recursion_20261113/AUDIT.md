# Validity audit and reproduction

**Self-audit with independent computational routes; no external review or
Lean formalization is claimed.** General results are proved in the notes;
the enumerations have the explicit finite scopes listed below.

## Proof checks

- PG1 counts maximal nonroot-input runs, not arbitrary walks. The shared
  root can be used once. The second-bin-empty case is included. The
  cancellation of the temporary end-root coordinate has an explicit
  insertion/deletion bijection. The theorem is scoped to nontrivial
  operands; singleton identity contexts are separate.
- A context transfer has an actual entry and exit vertex. Every claimed
  word is a literal tournament construction. Root changes invalidate
  the fixed-root sign-cone hypothesis; this is the mechanism of escape,
  not a contradiction of that theorem.
- The H5 seed pair is checked by permutations and is strong. Seven U/W
  words through length2 are constructed and recounted, with strongness
  checked. The larger inherited H605 pair undergoes two literal returns
  at orders14 and17. The three selected stranded-pair escapes are
  recounted at orders16,17,18.
- The four-coordinate quotient to(h,E) is used only on F+s=h+q and under
  matrices that preserve the nonzero half-slope repair. It is sufficient
  for this language's future observations and entry, but not for unit
  readiness or ordinary isomorphism classes.
- The inverse selector is an iff statement about a defined finite
  alphabet and seed. Positive numerical inversion alone is never called
  a general tournament-realization theorem.
- The sparsity estimate counts construction words, so collisions only
  strengthen the upper bound. It allows a finite initial seed set but
  does not cover arbitrary larger contexts or interleaved preparations.
- The Collatz inverse parity guards, positive-domain guard and graph
  directions are explicit. Finite-box misses, greedy-strategy tests,
  and universal reachability are three separate assertions.
- Ordinary adjacency powers and squarefree Hamiltonian extraction are
  distinguished. The oriented7/21 witnesses lose one pairwise edge and
  therefore are not tournament counterexamples.

## Exact universes and independent routes

| Engine | Universe or check |
|---|---|
|context_head|Complete inherited ordinary types through8; every root canonically quotiented; fixed-root response determinant|
|typed_glue|7,952 physical labelled small gluings;112 permutation controls;216 associativity controls;56,030 distinct response vectors through8|
|feedback_probe|All64 common one-vertex masks for the aligned six-core pair; specified one-exception masks on four inherited stranded pairs; every output root|
|reroot_probe|All roots of eleven literal pairs; four selected later contexts physically counted|
|port_compiler|Universal finite run-word expansion,5,920 independent physical controls, selected order8 matrix certificate|
|port_search|Complete ordinary context head through8, four specified input slopes;10,139 eligible entry/exit occurrences; every first rule independently compiled|
|renewal_verify|Three selected escapes; two actual repeated returns; exact characteristic polynomial and decreasing inverse identities|
|renewal_catalog|All half-slope plane-preserving return matrices through8; literal U/W context factorization|
|language|Complete forward closure throughH=10^8; independent backward membership on every1,346 positive state and1,000 controls|
|collatz_graphs|Three induced graphs on1..10^6;270,000 reverse-edge guard controls; separate greedy certificates for starts through10^5|
|walk_controls|All4,096 loopless digraphs on four vertices: permutations, subset DP, squarefree propagation; selected one-arc witnesses|
|audit|Eleven strict sign-cone checks; seven literal strong words;6,087 primitive gcd-gate controls; cross-receipt checks|

Largest literal subset DP order is18, and18!<2^63. All its path, marked
two-bin, and deletion two-bin counts fit signed64-bit integers. Larger
symbolic words use Python arbitrary-precision integers. No huge tournament
was silently recounted using overflowing fixed-width arithmetic.

The ordinary bank is inherited, not re-enumerated from all2^28 labelled
tournaments. Its type counts1,1,2,4,12,56,456,6880 are checked. Completeness
of its all-root expansion uses the fact that every tournament has at least
one Hamiltonian ending vertex. Its frozen receipt and all imported helper
files actually used by this packet are hashed in REPLAY.json.

## Reproduction

From the repository root:

    python3 -X utf8 05-knowledge/results/endpoint_context_recursion_20261113/replay.py
    python3 -X utf8 05-knowledge/results/endpoint_context_recursion_20261113/replay.py --optimized

The driver verifies pinned inputs before running. The optimized run also
requires byte-identical normal/optimized compressed receipts, including
the response catalogue. Gate checks use explicit exceptions, not Python
assert statements. Final completion status is recorded below after the
run finishes.

**Completion: PASS.** All twelve engines completed in normal and optimized
Python modes, and all thirteen compressed receipt pairs were byte-identical.
The replay manifest pins26 executable/input files. The final optimized
replay was run after the substantive code changes above.

Corrections caught during development: the inherited list has eleven
pairs, not ten; the length of a first saved witness word is not the
longest word after quotient deduplication. The latter was replaced with
an exact dynamic program on strictly increasing target H (maximum8).
Neither correction changes a mathematical proof or the1,346-state count.
