# All greedy endpoint descent, via reverse-injective past guards

19 October 2026. **PROVED, scoped; independent audit PASS.**

[ALL_GREEDY_ENDPOINT_DESCENT.md](ALL_GREEDY_ENDPOINT_DESCENT.md) proves:

- Every greedy V/F first-crossing word W_m strictly descends on every
  positive lawful source, for all m>=1.
- More generally, fix any nonempty expansion template U over {K,V}. Append
  the maximal allowed number of F resets after each U, then one final F.
  Every resulting word with m>=2 copies strictly descends. The arbitrary-U,
  m=1 case is explicitly outside this corollary.
- Completed-iteration checkpoints are pairwise distinct: their terminal
  image types identify the last code branch, and injectivity lets an assumed
  return be decoded back to an impossible return to the original source.
- A finite minimal-left-occurrence proof gives an equivalent way to remove
  the cycle alternative from any sufficiently long repeated symbol factor.

The mechanical reset sequence has factor complexity at most k+1. At
k=ceil(m/100), finite pigeonhole provides the required repetition, and the
variable-length encoding preserves its long literal future. No infinite
rotation recurrence estimate or numerical episode cap supplies the proof.

This closes the endpoint-descent goal for the whole greedy family, strengthening
the preceding image-guard packet’s upper-convergent subsequence. The stricter
different-types-only finite occurrence hypothesis remains separate and open. The complete realized
S_d boundary and arbitrary bank entry remain open.

## Dependencies and controls

The proof uses the audited
[matched-block theorem](../extensions/MATCHED_BLOCK_EXTENSION.md), plus the
exact [bank image guards](../image-guards/IMAGE_GUARD_DESCENT.md).
The former inherits the primary explicit logarithmic-form input.

Run:

    python3 verify_all_greedy.py
    python3 -O verify_all_greedy.py

Both modes agree with the supplied receipts: 40,345 explicit checks.
Controls include all 2,046 nonempty binary words through length ten at every
factor length; selected greedy words through m=1,024 with actual lawful
macrobranch replay; and 21 fixed-template controls, including templates
requiring multiple resets per iteration. Direct literal Terras replay is
separately marked in the receipts. These finite controls audit the lemmas
and arithmetic; the written proof supplies every universal quantifier.

## Independent reproduction and scope overlay

Read the [independent proof review](review/ALL_GREEDY_ENDPOINT_DESCENT_REVIEW.md). The exact [scope-reviewed text](review/ALL_GREEDY_ENDPOINT_DESCENT.scope-reviewed.md), patch and overlay metadata are retained in `review/`. The main paper differs from that reviewed text only by its completed-audit status line.

Using Python 3 and its standard library, from this package directory:

    python3 -B review/independent_all_greedy_controls.py
    python3 -B -O review/independent_all_greedy_controls.py

Compare stdout with `review/independent-receipt.json` and `review/independent-optimized.json`. Both runs perform 156,750 explicit checks, including 3,138 finite binary words and 107 lawful whole-word controls for 34 fixed templates. The checks remain effective under optimization. `SOURCE_PROVENANCE.json` records the exact original, scope-reviewed and published identities, the package-relative manifest lookup, and the receipt provenance changes. No verifier mathematics was altered.
