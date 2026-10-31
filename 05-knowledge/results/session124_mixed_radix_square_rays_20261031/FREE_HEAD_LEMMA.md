# Independent semantic audit: separated heads need independent shared colors

**The relaxation is valid.** It does not weaken any already frozen theorem: the earlier within-word linkage was a stronger sufficient condition, not a necessary invariant.

Let R=q² be odd and h=(R-1)/2. In a compatible input pair, a common source label x has the same one-based input-position parity in both members. Thus its representative `R*x+(-1)^j*c` is the same whichever member is selected for offset c. The bands `[R*x-h,R*x+h]` are disjoint, and within each band the offset is unique. The low singleton labels are also disjoint from those bands.

When a block is cut after source label1, the head is R-c, in band1. Every tail label has source x>=2 and lies strictly above band1. Hence a head cannot secretly coincide with a tail, a head of another offset, or a low singleton. The source-type decomposition of every common output label is unique.

For a tail beginning at output-position parity s and an input member of length L, the output-minus-input parity shift is s-2 for a forward traversal and s+L for a reverse traversal (modulo2). Whole-block shifts are s-1 forward and s+L reverse. These shifts are constant throughout their respective parts.

Consequently two output words are compatible if and only if:

1. Each low singleton has the same output parity in both words.
2. Each split head has the same output parity in both words.
3. Each tail or whole block has the same input-to-output parity shift in both words.

The converse uses the nonemptiness of each common source set; it holds automatically for the B inputs here, of size at least42. **No condition equating a head's within-word shift to its tail's shift follows or is required.** Later recursive lifts require only common-label positional parity, which the three conditions establish exactly.

All other obligations remain: interval ownership, forced member selection, complete traversal, symbolic length, endpoint squares, protected prefix and terminal rules. This lemma does not turn an arbitrary port matching into a Hamiltonian path.

## Real square-sum counterexample to necessity of the old linkage

The independently replayed scale25 candidate at input parity0, residue18 splits offset6. Its head is19 and occurs in position4, so its output parity is0. The tail shift is0, although the old linkage would require the head parity to be1+0=1.

Both output words agree on the head's own parity and on the tail shift separately. Direct independent expansions on B42 give compatible square-sum Hamiltonian paths on[1,1068] and[1,1069]. Thus this is an actual Hamiltonian-pair counterexample to the old linkage's necessity, not merely an abstract coloring example.

The changing-block rule still holds: the extra output label has parity N+1 and its source label n+1 has parity n+1, forcing the changing tail/block shift to be d modulo2. This constrains the tail's shift, not the detached head's color.

## Independent checker and hostile test

`B/audit_independent.py` imports no producer module. Its cross-word profile comparison treats head entries separately from tail shifts. It never imposes the obsolete within-word equality. A hostile semantic-certificate mutation flips exactly one head-color entry in one output profile, leaving every other profile entry unchanged. The same comparison function rejects it specifically as `shared head parity` failure.

The partial receipt initially checks42 available paired states and all50 terminal states, including direct full-vertex expansions from B42/B43:317,115 explicit checks, byte-identical normal/optimized receipts. This partial check deliberately does not assert an all-size four-prefix theorem. Complete-table existence and full constructor replay remain separate obligations.

The later frozen complete table has now passed the full independent audit: all50 paired and50 terminal states,51,028,631 checks normal/-O, plus103 streamed cycles. See `B/INDEPENDENT_AUDIT.md`. The earlier partial receipt remains a research-only intermediate semantic check, not a dependency of this package.
