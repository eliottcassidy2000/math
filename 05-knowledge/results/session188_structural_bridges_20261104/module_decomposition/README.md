# Session 188: strong-module reduction

`PROOF.md` proves hereditary Hamilton/path-cover monotonicity under arbitrary
simultaneous XOR of modules of the original tournament, using the classical
Gallai decomposition and the positive occupied-run substitution identity.

It proves the following are equivalent for such an output U:

* H(U)=H(T).
* Every originally transitive strong-tree quotient remains transitive.
* Every induced-subset Hamilton count and every ordered cover count agrees.
* The physical directed-triangle sets agree.
* U lies in the legal sequential module-inversion class of T.

The final equivalence is consistent with the classical C3-structure theorem;
an explicit legal sequence is constructed and computationally replayed. The
scalar-preserving class is a **residual full-response problem**, not a claimed
collision or a complete reconstruction theorem.

Run `python verify.py > normal.json` and
`python -O verify.py > optimized.json`, then compare the files. The checker
enumerates every named tournament and every original-module XOR family through
order five, with modules built directly from their definition. It checks every
induced-subset ordered-cover coordinate and verifies legality at every step
of the constructed current-module inversion sequences.

Both modes pass with identical receipts. There are 132,041 original-module
families: 21,629 scalar-equal and 110,412 strictly increasing. The same-session
`interval_modules` packet preserves full cubic/quartic contrast controls,
including a nested scalar-equal example that is not a full named collision.
