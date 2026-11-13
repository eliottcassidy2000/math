# Session 316: Marked square-sum hinges and exact signed tournament repairs

## The square-sum input

[Square-sum hinges](SQUARE_SUM_HINGE.md) identifies $Q_8$ exactly:
$4-5$ and $2-7$ are direct sum-9 edges, while $8-1-3-6$ is
the composite reflection $R_9R_4R_9$. The middle edge consumes real
vertices: replacing it by an unmarked virtual connection gives the
wrong spanning-path test. The whole marked segment first appears
in a Hamilton path at order 32. A new marked 25-fold glue word keeps
all three literal motifs in Hamilton cycles at infinitely many orders,
starting $887,22187,554687$. This uses older reflection-orbit,
Pell, marked-ear, and Gerbicz lift work.

[Squarefree port elimination](SQUAREFREE_PORTS.md) makes the analogy
exact. A new interior vertex flips the only spanning two-port pairing
from straight to crossed in both a square-sum graph and a tournament.
Every individual port-to-port reachability stays positive. Nilpotent
interior labels prevent using the same vertex twice. With all exterior
seams available, a $k$-path cover with $p$ odd paths has weighted
odd-only closure count $2p(k-1)!$, with an exact parity condition on
the exterior permutation.

## Exact tournament control

* [One disagreeing contact](NEAR_TWIN_GATE.md) costs at least twice
  the reservoir count. Exactly +2 occurs only for absolute pairs
  $(1,3),(3,5),(9,11),(13,15)$.
* [Two opposing contacts](OPPOSED_GATES.md) give strong absolute pairs
  $6\cdot2^L-1,\ 6\cdot2^L+1$ for every $L\ge1$, by a single
  arc reversal. The next two-slot response is never equal across the
  pair: the exact defect changes sign at module order 4.
* [The neutral mirror bank](NEUTRAL_BANK.md) gives
  $\binom{d}{k}$ distinct rigid strong self-converse tournaments
  with identical full path-cover polynomial and identical *named*
  apex interface. Local probes distinguish the position of a switch
  and can have opposite signed responses. [Sequential reentry](SEQUENTIAL_REENTRY.md)
  shows that the unprobed choices still admit exact neutral exchanges
  after any fixed pattern of revealed ports; probing right to left\n  keeps a strictly positive readout at each nontrivial branch.
* [Cyclic wrap and the 7-ray](CYCLIC_WRAP.md) shows that the first
  three entries $7,21,63$ are the only self-converse holes on
  $7\cdot3^j$. The seeds at 189 and 567 are strong. Ordinal
  propagation is not strong, and the [mixed-cycle cost](PARTITION_COST.md)
  proves it cannot be made so at the same count just by changing
  the cross arcs. The exact cyclic wrap operator gives a separate
  mod-9 obstruction to repeated refill from the 567 seed.

## Replay and scope

Run the Python replays in this folder: square-sum, ports, near twins,
opposed gates, neutral bank, and wrap. The larger literal 21/22 vertex
neutral-bank counts can also be checked by compiling
`replay_literal.cpp` with C++17 and feeding each line of
`literal_rows.txt` separately. The finite scoped extension scanner
takes input `n target`, the old outgoing rows, and the old involution
as a zero-based permutation.

This session gives reusable exact port control on explicitly anchored
families, not a way to choose the anchor for every odd target. In
particular the signed probe effects have not been shown to cancel at
arbitrary prescribed magnitudes. Chan–Pak remains open.
