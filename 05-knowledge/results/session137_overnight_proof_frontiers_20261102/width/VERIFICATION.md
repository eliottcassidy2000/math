# Verification and source boundary

## Reproduction

```
python3 -B context_residue.py > normal.json
python3 -B -O context_residue.py > optimized.json
cmp normal.json optimized.json
sha256sum -c SHA256SUMS
```

The main program uses only the standard library and no producer imports.
All load-bearing checks use explicit `require`, not removable assertions.
Recorded outputs are byte-identical and contain:

- 1,099 labelled contexts and 1,099 independent literal cover checks;
- coefficient gcds1,2,12;
- sharp detector rows through eleven vertices;
- 320 direct actual-graph substitution checks of the modular transition;
- all positive representatives for the exact class counts5,3,14,11;
- the exact seven-state restricted grammar and H11/H55 escapes;
- 5,421 repeated-support divisor checks, sharp balanced examples, and the
  three-repeat cyclic hostile;
- the L3/L4 modulus8 old-decoder failure and sharp-detector repair.

The independent reviewer `/root/audit_overnight_proof_frontiers` independently
accepted the all-order sharp-divisor, contextual-equivalence, transition,
low-modulus classification and branching-obstruction proofs, and separately
accepted the late inactive-repeat-support refinement. Its separate
literal x-word checker is
`/tmp/session137_research/independent/audit_width.py`, with 1,099 quotient
checks and sharp detectors through c4. It does not import the producer code.
The parent session owns publication and placement of those audit files.

## Inherited mathematics versus new consequences

The Session136 substitution kernel, Rédei parity, OCF, reciprocity mechanism,
and known spectral widths are inherited. This packet's advances are the
sharp factorial-return coefficient ideal, the exact minimal modular context
quotient and explicit finite algebra, the small universal context-class
counts, the quotient-aware repeated-support divisibility, and the infinite
same-width branching-resource separation. The low-modulus upper bounds are
all-order proofs, not extrapolation from labelled tournaments through5.

The sparse baseline contains the complete Session136 packet but omits some
older canon files named inside it. The exact input statements are present
in the inherited packet; those missing older files have not been pretended
read. A targeted search of the available packet and canon subset found no
prior exact modular contextual quotient or sharp excursion divisor. This
is not an exhaustive external literature search or priority claim.

## Discarded or stopped exploratory routes

- A larger unary-orbit sweep was terminated with process exit137 before its
  requested modulus range finished. No large-modulus counts from it are
  load-bearing or included in the final receipts. The bounded verifier above
  finishes normally in both execution modes.
- A speculative half-weighted polynomial shift formula failed already on
  the cubic monomial. The full path-cover coefficient law(9), proved by
  counting joining edges and independently checked against actual graphs,
  is used instead. No failed formula remains in the proof graph.
- Fixed-width density and globally bounded/unbounded spectral width were
  not resolved. The exact modular theory and infinite constructor hierarchy
  are substitutes for neither of those open problems.

No git fetch, mutation, commit or push was performed by this child. All
intentional files are confined to its assigned directory, below0.1MB.
