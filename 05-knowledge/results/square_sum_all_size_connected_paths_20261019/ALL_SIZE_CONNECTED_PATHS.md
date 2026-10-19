# Simultaneous connected three-vertex paths in every sufficiently large interval

All98 weak-to-strong B states have been independently replayed, including all49 exact duals. The proof combines this finite certificate with the already complete weak B, protected A and final-cycle constructions.

## Theorem

Let S be the integer interval [−1200,1200] with the two intervals [−24,24] and [74,122] deleted. Thus |S|=2303. For every integer N≥6294221 there is a square-sum Hamiltonian cycle on [1,N] containing consecutive482,3362,6242 and simultaneously every path

    (2401+t, 7203−t, 14406+t),  t∈S.

These2303 paths are pairwise vertex-disjoint. Their adjacent sums are9604=98² and21609=147². The labels are fixed independently of N. No optimality of the sufficient threshold is claimed.

## 1. One-time weak-to-strong transition

The existing eventual theorem constructs weak B_m pairs for all m≥53, with prefix1,3, shared positional parity, and even-length terminal8 / odd-length terminal10. The new98-state table takes precisely these existing weak inputs. It has the same block/segment semantics and exact partitions as the old table, but every output word begins with singleton1, singleton3, singleton6. Thus the pair output has the stronger prefix1,3,6.

The universal checker explicitly validates all three prefix vertices and the common parity profiles. No new base pairs and no stronger input hypothesis are required. For every n≥2621, write n=49m+d with24≤d≤72. Then m≥53, so one new transition applied to the already available weak B_m gives a strong B_n pair.

The parity1 states are exact duals of the49 parity0 states: replace d by96−d, exchange the two output words, complement every selected input-member index, and retain offsets, orientations and singleton atoms. Therefore singleton prefix1,3,6 is preserved. The checker replays all98 dual-derived states individually rather than relying solely on the duality assertion.

## 2. Protected A lift

Take a strong B_m input. In its first three positions, P=(1,3,6) has alternating signs (−,+,−) under the affine lift. At each offset c, the copy is

    (49−c, 147+c, 294−c).

The old Aplus/Aminus tables use only whole blocks except their already certified prefix cuts at offsets0 or−2. Conservatively omit both offsets in every state. All47 remaining copies of P stay wholly inside uncut blocks, irrespective of their orientation or selected input member.

Set u=−c. The47 retained first-level parameters are u∈[−24,24] excluding0 and2. This sign is fixed by the odd first position of1 in the strong input; it is not an unspecified parity convention.

## 3. Final lift

Both A paths contain every retained copy with the same vertex labels and positional parities. All its vertices are at least25, so the final2–7 cut does not remove any of its edges. For each u, the next49 offsets produce parameters

    t=49u+σ_u c',  −24≤c'≤24,

where σ_u=±1 is determined by the common positional parity of the copy's starting vertex. Regardless of σ_u, these fill exactly[49u−24,49u+24]. Their union is [−1200,1200] minus the intervals centered at0 and98, namely[−24,24] and[74,122]. This is exactly S. Reversed blocks or final segments reverse only path traversal; the undirected prescribed paths remain.

The same final certificates give seed482–3362–6242. The prescribed path at t=961 is3362–6242–15367: it shares the seed edge3362–6242 and continues it. The shared vertex3362 is an endpoint of this prescribed path, so its other edge can remain482–3362. Hence the construction also forces the four-vertex seed482–3362–6242–15367. Pairwise disjointness is asserted among the2303 prescribed three-vertex paths, not between every such path and the original seed.

## 4. Threshold and disjointness

Writing N=49n+d2 and n=49m+d1 with both residues in24..72, the least m is floor((N−1200)/2401). Hence N≥2401·2621+1200=6294221 ensures m≥2621 and makes the strong pair available.

The three coordinate ranges are[1201,3601], [6003,8403], and[13206,15606]. They are disjoint, and each coordinate is injective in t, proving that the paths are pairwise vertex-disjoint.

## Verification and attribution

`verify_complete_strong_table.py` checks complete state coverage, exact duals, and invokes the independent typed checker on every state. The stronger prefix is explicitly checked. Both normal and optimized receipts are included. A streaming construction at N6294221 checks every cyclic edge, exact vertex coverage, all4606 prescribed path edges, and both original seed edges, without storing the large cycle. The four-vertex seed follows because the t961 prescribed path is verified together with the original seed.

This is a certificate-based extension of the previously credited square-sum construction, not a claim that the underlying lift, parity-pair method, or search heuristic is novel. No historical-priority claim is made for this prescribed-path statement without further literature work.

## Repository layout and commands

Place this folder next to `../seeded_eventual_cycle_theorem_20261019/`, `../seeded_final_split_table_20261019/`, and `../square_sum_B_base_interval_20261019/`. Run `python verify_complete_strong_table.py` and its `python -O` counterpart. Run `python verify_all_size_control.py` for the threshold cycle, or pass any integer N≥6294221 to construct and verify that interval's cycle in streaming form. No large cycle is saved. The finite table checker is self-contained; the streaming constructor uses the stated sibling dependencies.
