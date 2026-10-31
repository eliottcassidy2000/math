# Exact integer-potential producer and scheduling controls

The inherited height-adaptive graph is unchanged. This experiment uses
integer excess 8*S_E-5*A, rather than the source's large rational numerator.
Every lifted edge, paid entry phase, and original positive-letter clock
remains the same. A potential is a proof only after a separate checker
verifies all its required inequalities and infinite tails.

The first FIFO E10 run was stopped before completion because it was slow.
A duplicate-bearing bucket queue was also stopped before completion when
its stale entries threatened unnecessary memory growth. Neither has a
certificate result. The final candidate uses an intrusive integer bucket
queue with at most one queued copy per state. The previous and next arrays
use eight bytes per state, the potential two, plus the height table.
The E10,R2 graph has245312500 states, so the main arrays require2453125000
bytes, plus9765625 height bytes. It runs under a3000000KiB virtual-memory
limit and reports peak resident memory during execution.

This is a label-correcting schedule, not Dijkstra's settled-once algorithm:
a processed state can be inserted again whenever its label improves.
Integer priorities merely determine update order. They do not discard
edges, alter the bound, or identify relaxed paths with actual words.

E5,R2, E6,R2 and E7,R2 final potential files match the inherited FIFO
scheduler byte for byte, even though the tested producer uses letters1..40
and the original controls used1..64. Exact hashes:

- E5:716d658df261477daa29bc81f1d8cd6e9c8b2212f83329b1f81815a00edeab52
- E6:fb1d060ebc0e73a8b4e0611cb6023abf991c57bb2bb1e65059ff9a1533b058fc
- E7:04dafa22053fc5ad88d8de3cacb18859e6bb22a3a7bd5d0047e8598e577f3088

The independent full-letter checker, rather than these regression matches,
will determine whether a final E10 output is valid. In particular its
all-letter tail checks must pass. A producer's current terminal label is
not an upper bound on its eventual label and cannot establish the desired
terminal fee before completion and checking.

## Completed E10 verification

The producer finished with277792217 pops, terminal-30, minimum-128 and
maximum-10. Its coarse40-letter entry tail prints false, because-125 is
not at most-128. The independent64-letter checker nevertheless passes all
18850000000 edge inequalities, all630 paid entries and the65+ tails.
Thus the output is a valid all-letter potential despite the producer's
insufficient coarse proposal-only tail test. The1,962,500,000-byte potential
is not committed. Its SHA256 is
`9afc878aaaf0ff11b075678493f6ef935c1735afe4f81bdcd9a2e1282c652c31`.

The independent near-tight graph and semantic multiplier exclusions, not
the raw terminal score alone, supply the stronger actual endpoint fee
needed for E14. See `E14_CLOSURE.md` and the companion independent audit.
