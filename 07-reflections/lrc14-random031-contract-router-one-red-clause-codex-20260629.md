# LRC14 random031 proof-contract one-red-clause audit

HYP-3528 executes the one-red-clause audit for the contract-router idea after
the namespace was cleaned up: HYP-3525 remains the guarded-emission atlas,
HYP-3527 carries the prefix-free proof-contract router, and HYP-3528 supplies
the independent one-red-clause executable synthesis.

The main result is a change of interface.  random031 should not be thought of
as a row-to-route table at this point.  It is a sidecar-to-safe-emission
contract:

```text
terminal clause + required sidecars -> emitted proof token + remaining tail
```

The executable contract imports the current spigot/carry stack and finds:

```text
79 component events stream 77 terminal certificates
doublet predigits have rank gaps (1,1)
owner carry opens at rank 45 and reduces at rank 46
final owner carry is (45,173)
I/Q prove private status but do not reconstruct route
R is still required for terminal dispatch
vertical halfturn glues 69 components with 7 mixed components
```

The useful compression is severe but honest.  The proof can stop keeping most
of the component table alive, but it cannot forget:

```text
R
residual_pair_(45,173)
no_hidden_tail_guard
vertical_halfturn_guard
```

The single open clause is now explicit:

```text
residual_pair_close_tail
```

That clause must show that after transport `(23,93,113)` and branch-boundary
lift `(147,169)`, the residual pair `(45,173)` cannot hide in a downstream
quotient when route sidecar `R` is retained.  Residue buckets, sliced-box
volumes, owner counts, and raw counts are not substitutes for that lemma.

Creative reframe: the forbidden seam is not a wall.  It is a carry boundary in
a mirror-punctured cylinder.  Phase flow moves in the complement, emits
ordinary/free-hole/bypass certificates, and leaves a typed residual tail.  The
proof ABI is the list of sidecars each clause consumes before it is allowed to
print a theorem-facing token.
