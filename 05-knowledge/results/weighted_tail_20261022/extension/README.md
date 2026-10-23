# Sixteen-edge native extension with immutable source owners

October 23, 2026. **FINITE-EXACT:** bounded extensions of the four frozen
ordinary native states, grounding additional inherited source fibres and
improving their all-future weighted gates. **OPEN:** uniform 153N for every
N >= 4096, arbitrary-source frontier entry, and the Collatz conjecture.
The N7 counterexample remains a counterexample.

## Mechanism and retained information

The input is the ordinary synchronous state from
[future_tail_blocks_20261021b](../../future_tail_blocks_20261021b/README.md),
with the optimizer from [the weighted packet](../README.md). This is not an
extension of the optional-pulse relaxation. The latter is a different type
of object and cannot be silently interpreted as a native state.

Each original pending node at the input depth D owns one immutable fibre of
original primitive sources. We supply up to 16 additional native edges from
that node, stopping exactly at its first visit to the pinned frontier V.
Every original owner appears once. At intermediate cut D+r, its value is the
original four-vector multiplied by the sign of its actual extra integer
clock. A first-arrived owner contributes to its literal terminal once;
otherwise it contributes to its actual current node. Owners at the same
node may be aggregated for the energy calculation, but remain separate in
the ownership table. A four-coordinate zero aggregate may be omitted from
that calculation; this does not ground or erase any original source.

This representation preserves source identities, separate old/new channels,
mod15/mod35 registers, actual native depths, full integer clocks, and the
first-arrival terminal. It forgets only the distinction between owners when
forming their synchronous energy aggregate; the owner table is the sidecar
that recovers that distinction. Its cheapest decisive test is independent
owner transport versus synchronous aggregate transport at every one of the
17 cuts, including the original cut. These agree exactly.

The source identities are the inherited lossless `(N,D,s)` carrier, not an
unlabelled collection of counts. Their original prefix is reproducible by
the frozen native rule. For a member s of owner u, its full continued clock
is its own original clock plus the authenticated common suffix clock for u;
the original residues remain `s mod 15` and `s mod 35`. Each owner retains
its inherited multiplicity, minimum and maximum full clocks, and literal
source witness. The consumer independently replays every witness's complete
original prefix and charge, checks both residues, then appends the suffix.
The original finite source census and every fibre membership are inherited
from the pinned prefix receipt, not independently re-enumerated here. A
hash binds that established input; it does not replace its original proof.

The finite universe is precisely the previously pending nonzero owner bank.
Sources omitted earlier because their initial charge or synchronized linked
aggregate was zero are outside these new grounding counts. There is no claim
that those omitted sources enter V.

## Exact results

| N | Input D | Output D | Newly grounded owners | Original sources in those owners | Remaining owners | Remaining nonzero nodes | Final weighted bound | 153N |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 7 | 12 | 28 | 5 | 5 | 0 | 0 | 1097 | 1071 |
| 4096 | 105 | 121 | 146 | 302 | 49 | 49 | 441842 | 626688 |
| 8192 | 110 | 126 | 247 | 554 | 101 | 101 | 837316 | 1253376 |
| 32768 | 130 | 146 | 471 | 954 | 289 | 286 | 3173911 | 5013504 |

The three positive cases ground 864 original owner fibres containing 1810
original source addresses, counted **within their separate N-specific
carriers**. This total is not a claim of 1810 globally distinct integers
across the overlapping scales. Their final bounds improve the preceding
ordinary weighted bounds by 109327, 280974, and 959532 respectively.

The maximum actual increments among the new cuts are 423991, 795841, and
2844812. Every new cut is below its unchanged budget. The inherited receipt
covers earlier cuts, and the final weighted gate covers all later cuts,
including arbitrary nonarrival of the remaining owners. There is no
requirement that a computed bound decrease at every intervening cut.

The final optimizing parameters alpha are respectively 0, 63/100, and
173/285. The corresponding angular gates are 441898, 837433, and 3177849;
the absolute envelopes are 442192, 839113, and 3189619. These are separate
upper bounds on the same final state. Their improvements are not added to
one another, and no source or pair credit is recycled.

At N32768, the 289 remaining original owners occupy 286 synchronous nodes.
This is genuine coalescence at the same depth, with all their owner labels
retained. It is not inferred from equal vectors or equal residues.

The grounded frontier is bound to exactly 627 primitive integers. The
consumer checks its rank function independently: every nonroot point maps
inside V to a point of smaller nonnegative rank, and the root is 1. Thus
every newly arrived owner has a finite authenticated route into this
already grounded set. Reaching a frontier point and producing a weighted
tail certificate are distinct conclusions.

N7 supplies a particularly useful hostile: after the extension no owner is
unresolved, yet the actual final increment is 1097 > 1071. Grounding all
children does not by itself establish the proposed budget. This also shows
why a successful finite grounding computation cannot substitute for the
uniform interaction inequality.

## Verification boundary and reproduction

`produce.py` performs bounded path discovery using the original native
implementation. `check.py` only consumes those paths, with an independent
valuation-based edge implementation. It checks:

* literal integer types, exact case and owner schemas, and the complete
  one-use owner partition;
* every supplied transition, full integer clock, primitive-node condition,
  and first-arrival stopping rule;
* all original witness prefixes, charge formulas, clock ranges and residue
  guards;
* all 17 new synchronized cuts by two transport implementations;
* the pinned 627-point frontier and strictly decreasing grounding ranks;
* exact final weighted certificates and the retained N7 hostile.

Hostile controls corrupt clocks, duplicate or omit owners, truncate a path
before its unresolved cutoff, use floating-point/Boolean aliases, add a
spurious schema field, and continue a path past its first arrival. All are
rejected. Checks use explicit exceptions and remain active under `python -O`.

From the repository root:

```text
python -B 05-knowledge/results/weighted_tail_20261022/extension/produce.py
python -B 05-knowledge/results/weighted_tail_20261022/extension/check.py
python -B -O 05-knowledge/results/weighted_tail_20261022/extension/check.py
```

Both consumer modes reproduce the receipt with 676143 checks, including
the saved-receipt comparison. Receipt LF SHA256:
`16efb7ba793ddc9c8965315e7883614e803e58c6a6b2f431249cee858c8c767b`.

The next uniform obligation remains an estimate that guarantees a finite
passing cutoff for every intended N, together with its earlier cuts. The
present extension is a bounded exact improvement at three existing scales;
it does not fill that quantifier gap.
