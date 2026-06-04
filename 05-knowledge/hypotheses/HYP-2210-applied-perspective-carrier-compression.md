# HYP-2210: Perspective-Flip Carriers Give Practical Compression for LRC, Unit-Distance Cores, and Rooted Tournament Counts

**Status:** PARTLY PROVED / VERIFIED through S634.  The Burnside and
source/sink pairing statements are formal counting consequences; the
unit-distance statement is an exact graph audit of the five stored `n=21`
cores; the extension claims remain open.

## Claim

The S629/S630 perspective-flip carrier should be used as an algorithmic
compression principle, not only as a self-converse tournament curiosity.

It gives three immediate workhorses:

1. **LRC:** THM-381's observer-source target and the observer-sink obstruction
   are the same converse-paired obligation.  A source/sink ledger should carry
   one merged target family of size `U(n-1)`, not two independent lists.
2. **Tournament counts:** rooted tournament perspectives modulo all-edge
   reversal are counted by Burnside:

   ```text
   Q_n = (P_n + F_n)/2,
   ```

   where `P_n` is the number of rooted perspectives and `F_n` is the number
   of rooted perspectives fixed by the SC perspective flip.
3. **Unit distance:** the exact/Moser `n=21` row is not only numerically
   `57 = 20 + 37`; every one of the five stored `57`-edge graph6 cores is
   traceable, so the `20`-edge spine is graph-real and the remaining `37`
   edges form the `C_hex(3)` bulk carrier.

## Evidence

S634 computed the rooted/converse table through `n=6`:

| n | classes | rooted perspectives `P_n` | SC fixed roots `F_n` | rooted modulo converse `Q_n` | source targets | sink targets |
|---|--------:|--------------------------:|---------------------:|-----------------------------:|---------------:|-------------:|
| 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| 2 | 1 | 2 | 0 | 1 | 1 | 1 |
| 3 | 2 | 4 | 2 | 3 | 1 | 1 |
| 4 | 4 | 12 | 0 | 6 | 2 | 2 |
| 5 | 12 | 48 | 8 | 28 | 4 | 4 |
| 6 | 56 | 296 | 0 | 148 | 12 | 12 |

The source and sink target columns match `U(n-1)`, as expected from deleting
the unique source or sink.  Converse pairs the two marked target lists, and
THM-409 makes the pairing internal on self-converse classes by a perspective
involution.

S634 also audited the five exact `n=21` unit-distance cores from S614:

| core | edges | traceable | split | degree histogram | deletion-edge deck |
|------|------:|:---------:|-------|------------------|--------------------|
| 1 | 57 | yes | `20 + C_hex(3)` | `{5:18, 8:3}` | `{49:3, 52:18}` |
| 2 | 57 | yes | `20 + C_hex(3)` | `{3:1, 4:2, 5:8, 6:7, 7:3}` | `{50:3, 51:7, 52:8, 53:2, 54:1}` |
| 3 | 57 | yes | `20 + C_hex(3)` | `{3:1, 4:2, 5:9, 6:5, 7:4}` | `{50:4, 51:5, 52:9, 53:2, 54:1}` |
| 4 | 57 | yes | `20 + C_hex(3)` | `{3:1, 4:2, 5:8, 6:7, 7:3}` | `{50:3, 51:7, 52:8, 53:2, 54:1}` |
| 5 | 57 | yes | `20 + C_hex(3)` | `{3:1, 4:3, 5:6, 6:8, 7:3}` | `{50:3, 51:8, 52:6, 53:3, 54:1}` |

Thus HYP-2206's requested exact-core audit now passes for the stored cores.

## Why This Helps

The common operation is not raw equality of numbers.  It is retaining the
side channel that makes an apparent scalar collision harmless:

- LRC retains the marked observer/source threshold payload.
- Tournament counting retains rooted perspectives and the converse action.
- Unit distance retains the Hamiltonian spine plus centered-hex bulk split.

Once those channels are retained, several state spaces shrink before any hard
search begins.  The LRC source/sink target list halves under converse pairing;
rooted tournament counts get a direct quotient formula; and unit-distance
`n=22` extension work can focus on endpoint-compatible ears from traceable
`21`-cores rather than rediscovering the whole graph.

## Next Tests

1. Compute `P_7` and the rooted-converse quotient `Q_7`, preferably by a
   rooted canonical generator or by reusing the S630 SC generator only for the
   fixed-term `F_7`.
2. Upgrade the LRC quotient from unthresholded source/sink roots to
   threshold-decorated observer-source obligations: root, safe box, endpoint
   owners, denominator shields, and carry fibers.
3. For each exact `n=21` core, enumerate endpoint-compatible Hamiltonian
   spines and candidate one-vertex ears.  The `n=22` frontier asks whether a
   degree-4 extension of a `57`-core or a degree-5 extension of a `56`-core can
   be embedded faithfully.
4. Add `bulk = edges-(n-1)` and `bulk_shell` to the unit-distance beam ledgers,
   so centered-hex carrier events are tracked separately from raw edge counts.

## See Also

`04-computation/applied_perspective_carrier_s634.py`;
`05-knowledge/results/applied_perspective_carrier_s634.out`;
`07-reflections/applied-perspective-carriers-s634.md`;
THM-381; THM-409; HYP-2206; HYP-2205; HYP-2204; HYP-2203; HYP-2176.
