# Seven-return bounds, native seam recursion and delayed activation

2026-10-29; inherited `201fa96604d5`, with incoming `f8846701cd` and the
relevant native-operator results of `2a81e9c29d` read and integrated.
**PROVED scoped elementary results / FINITE-EXACT / CONDITIONAL / OPEN.**
Same-agent proof review with independent algorithms; no independent-agent
proof audit or literature-wide novelty claim. No new scarce theorem ID is
reserved or promoted here.

* [COLLATZ](COLLATZ.md): a finite two-stage test of the uniform seven-return
  candidate, the all-original-source mass bound through horizon74, and a
  sharp reverse decoder on the 81-divisible layer, retaining the source and
  both integer clocks.
* [SEAM](SEAM.md): three-coordinate transfer for arbitrary directed blocks
  with marked entry/exit; composition and converse; a five-vertex hostile
  showing that even the entire transfer operator does not survive native
  complement as a faithful representation.
* [COUNTING](COUNTING.md): exact delayed activation and half-deficit packing
  for inherited A000568 engines, with matched-output benchmarks and the
  boundary between fewer visits and a measured runtime improvement.
* [RECOVERY](RECOVERY.md): inheritance, prior stronger results, concept board,
  and the precise connections and next obligations.
* [PRECISION_EXTENSION](PRECISION_EXTENSION.md): two proved bounds on future
  multiplicity refunds extend the incoming low-bit recurrence, discarding
  whole partial-partition branches exactly.

## Main results and their limits

The bounded Collatz search excludes all positive odd sources below
53456334770533302271 from 74 consecutive growing returns. The resulting
full original-source survivor weight is at most

    641060333941629/2^148 < 1.797e−30.

This checks a consequence of two proposed seven-return steps from the exact
horizon60 minimum. It does not establish the actual horizon67-to74 ratio,
because the exact horizon67 minimum remains unknown. The uniform inequality,
A_k→0, and repayment below the original source remain **OPEN**.

On positive m divisible by81, the valuation-adjusted residue
`2^v3(m) * (m/3^v3(m)) mod81` identifies the last nonzero return among four
types, even after a run of zero letters. Its values are 61,13,70,34.
The weaker 27-layer fails by an explicit endpoint collision; a separate
marked-source hostile shows why reverse natural-number descent is not a
forward convergence proof.

The tournament observer now applies to arbitrary block sizes and coincident
ports, with exact matrix composition and converse. It has minimal linear
dimension three for seam path counts. Two five-vertex marked blocks have
the same full matrix but complements with path counts11 and13. Thus the
observer does not reconstruct classes or close the unrestricted metagraph.

For A000568, parking inactive buckets avoids88915 state visits at N=200;
the full-prefix engine's three-run median improves from1.376 to1.306 seconds
on this machine. Half-deficit projective packing saves a bit per response
digit but gives essentially no measured speed gain. All outputs agree exactly.
The same identity-stage contraction arrived concurrently and is credited.
Building on that packet's stronger precision theorem gives a larger gain:
the exact eight-normalized-bit request at order401 falls from2.169 to0.345
seconds in matched medians, about6.3-fold. This is a partial count, not the
complete a_401; all4983 contributing tuples are retained.

## Reproduction and audit scope

From this directory:

    python3 -B seam_transfer.py
    python3 -O -B seam_transfer.py
    python3 -B reverse_reset.py
    python3 -O -B reverse_reset.py
    python3 -B check_counts.py --bench
    python3 -O -B check_counts.py
    python3 -B check_colours.py
    python3 -B check_precision.py
    python3 -O -B check_precision.py
    python3 -B check_refund.py
    python3 -B check_full_precision.py
    c++ -std=c++17 -O3 window74.cpp -o /tmp/seven_return_window74
    python3 -B check_engine.py /tmp/seven_return_window74
    /tmp/seven_return_window74 74 53456334770533302271

The last command is the large finite exclusion, visiting763878259 prefixes.
It was completed once; small native controls and a high positive control
are independently checked, but a full second horizon74 algorithm has not
been run. The width proof and exact source hash are recorded. The C++ code
also passes undefined- and unsigned-overflow sanitizers at horizon30 with
output identical to the optimized executable.

The seam universe comprises2769 marked blocks: every labelled tournament
through4 with every ordered pair of ports including equal pairs; every
loopless directed graph on3 vertices with every pair of ports; and every
fixed-spine tournament presentation of orders5 and6. One hundred actual
compositions are checked against literal path/cut enumeration. The full
matrix hostile is checked by literal permutations. Normal and `-O` receipts
are byte-identical.

The reverse-decoder controls include every legal nonzero return from a
positive81-multiple at most1000000 and all its legal zero suffixes, plus
independent endpoint enumeration through100000. Native odd steps verify
the branch words and clocks. Normal and `-O` receipts are byte-identical.

Counting controls include every N=0..30, both parking modes, both projective
packing modes, all three vertex-colour counts1,2,3, direct odd-partition
counts, and an independent grouped-core coloured kernel. Larger benchmark
outputs at100 and200 are compared exactly. No truth-bearing test depends
on Python `assert`. `manifest.json` pins the scripts, proofs and receipts;
`source_pins.json` pins the inherited baseline snapshots.

Precision checks compare four channels through30 against full rational
Burnside sums at five precisions. Refund inequalities are also checked
on every prefix of every odd partition through30, independently of the
modular totals. Complete counts verify the partial outputs at101 and201;
no full count at401 was calculated.

The next proof targets are a uniform source-preserving bound for all k,
a legal reset-payment rule, and a marked-path carrier closed under both
seams and native complement/re-marking. They are not claimed solved by
the finite data or by the smaller observers.
