# Exact scope of positive physical context operators

Two independently reviewed packets distinguish what the universal context algebra can and cannot do.

- `PERIODIC_CONTEXT_DICHOTOMY.md`: repeating any fixed whole-module context word from any nonempty seed is exactly geometric when every root is passive, and superexponential when a root is in a directed triangle. Hence no nonconstant affine clock with a nonzero offset is realized by a fixed whole-module period. Exposed-port births and adaptive rebuilding are outside this statement.
- `PHYSICAL_NON_TN_CONTEXT.md`: the actual context C3[input,Tr3,Tr2] has a negative3-by-3 central-factorial minor. Adding one final transitive vertex makes it parity-preserving and keeps the negative minor. Entrywise positivity therefore cannot be promoted to total nonnegativity, and the padded operator cannot equal a parity-returning D/X word. No negative output-Hankel claim follows.

Run the two top-level verifiers normally and with `-O`. Their named executable dependencies already reside in this directory. Verbatim independent PASS reports, standalone checkers and receipts are in `periodic-context-review/` and `non-tn-review/`. For the periodic review, compile `check_physical.cpp`, then run `check_transitions.py` in that review directory; its comparison receipts are included. Source replay receipts are retained there, while the unchanged executable source is the corresponding top-level verifier and dependencies.

`CONTEXT_SCOPE_DEPENDENCIES.json` binds those previously published dependencies. `CONTEXT_SCOPE_MANIFEST.json` binds all new payload files. The later adaptive active-depth supplement and paired three-cover construction are separate and not included here.
