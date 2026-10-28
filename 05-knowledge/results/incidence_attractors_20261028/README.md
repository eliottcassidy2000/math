# From divisor-height cycles to guarded codes and attachment algebras

2026-10-28. **PROVED elementary maps and obstructions; FINITE-EXACT minimal
observers in an explicitly finite language; INHERITED results credited below;
PROPOSED larger research programme.** No claim of literature priority or of
an unrestricted tournament/arithmetic equivalence.

The recovered six labels are not a literal six-cycle. Session43's divisor
height map has the five-cycle `4->2->5->3->7->4`, with `13->2` an extra arrow.
Its common-forward-image fold identifies4 and13 only for the specified
forward observation. Their arithmetic attachments differ.

## Concrete results

1. **A literal graphic model.** On the atom bank `{4,2,5,3,7,13,11}`, the
   square classes of the divisor sums are the columns of a four-vertex
   multigraph's incidence matrix over F2. The triples `{2,5,7}` and
   `{4,7,13}` are circuits; `{2,3,11}` is a parallel pair plus a loop.
   Both66 and70 have sigma144, but their dependencies have different types.
2. **The realization guard is visible.** There are16 zero-syndrome subsets,
   of which12 are pairwise coprime. The symmetric difference of the two
   displayed triangle circuits contains both2 and4: its formal sigma is42²,
   whereas its actual product520 has sigma1260. The arithmetic subset family
   is not closed under the graph's cycle-space addition.
3. **An exact smallest finite carrier.** For guarded multiplication by
   those seven named atoms, starting from their96 pairwise-coprime subset
   products, and observing square sigma plus each legality response, the
   minimal deterministic quotient has **86 states**. Also observing the
   Boolean predicate F=U+S gives87; observing tau gives96. These are state
   counts, not linear ranks or memory-bit counts, and not an all-integer result.
4. **A Fano extension.** The prime bank `{2,5,7,19,29,59,89}` realizes all
   seven nonzero parity vectors on coordinates `{2,3,5}`. Its square-sigma
   subsets form the Hamming[7,4,3] code, with weight enumerator
   `1+7z^3+7z^4+z^7`. The seven three-element circuits are the Fano lines;
   their seven four-element complements are also circuits, for14 in total.
   This matroid is not graphic: ordinary graph cycle spaces cannot model
   every extension of the first graphic chart.
   Uniform laws on two code cosets agree on every at-most-three-prime
   selection marginal and on the marginal joint law of(mu,lambda), yet
   have square-sigma probabilities1 and0: a finite arithmetic analogue of
   the earlier Vitali alignment obstruction.
5. **Frobenius phase is explicit.** In F8=F2[alpha]/(alpha³+alpha+1), the
   chosen normal basis has arithmetic point labels `(5,89,29)`. Frobenius
   rotates them. It preserves the zero-syndrome predicate but sends source70
   to11837, changing sigma144 to14400. A chosen Singer order gives an actual
   Paley7 tournament with189 Hamilton paths and automorphism order21; the
   ordering and orientation are additional choices, not intrinsic prime order.
6. **A source connecting FUS, square sums and Berggren.** The inherited
   mixed balanced minimum364=4*7*13 has `(F,U,S)=(10,3,7)` and aliquot edge
   `364->420`, total28². Applying the inherited norm-one Berggren absorber
   `(5,5,7)` yields the literal square-sum path
   `19236--364--420--19180`, replacing the host edge19236--19180. This is a
   local insertion identity, not a proof of a spanning path in a whole host.

[THEORY.md](THEORY.md) gives maps, proofs, guards and the proposal.
[RECOVERY.md](RECOVERY.md) traces June, August, September and October work,
the live concept board, and the concurrent Session84 contributions.
[AUDIT.md](AUDIT.md) records scopes and reproducibility.

The resulting proposal is to organize the repository by **operations and
observable contexts**: use incidence relations for compatibility, retain the
transport action, and attach the actual arithmetic or geometric realization
guard. Classify the resulting response algebras, then study their quotient
maps. This is a research programme informed by recognizable graph algebras,
connection matrices and binary codes, not a newly proved universal theorem.

## Reproduce

```text
python -B 05-knowledge/results/incidence_attractors_20261028/verify.py --replay
```

The independent audit compares all128 atom subsets, all12288 finite
state/context pairs, and all128 Fano-bank subsets by a different divisor-sum
method. It also checks norm-one descent directly through c150. The producer
replays the full Petersen local-response comparison over1024 quotient edge
sets and the six K3,3 cycles, and checks364 branch addresses through depth5.
Normal and optimized outputs agree. Source pins are immutable Git blobs.

Inherited boundaries remain: the full unnamed tournament carrier and
strong-atom classification are open; the81/27 pair still defeats coarse
arithmetic summaries under path complement. Graceful-tree existence for
every prescribed shape, LRC(14), and global aliquot/Collatz behavior are not
settled by this packet. A finite prime bank supplies no prime-infinitude claim.
