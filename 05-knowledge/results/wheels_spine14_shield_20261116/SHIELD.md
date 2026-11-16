# Factor-1 shields, +2 moves, and coordinated multi-bit edits in the synchronized binary reader

Session klein-2026-11-16-wheels2. The owner asked for "a factor-1 shield that still has a +2 move" and for "a
coordinated multi-bit edit or a different seam". This note fixes the terms from the repository's own packets, proves
what is immediate, reports an exact census of address-level moves, and states what a third coordinate would have to be.

## 1. The terms, as the repository uses them

The synchronized reader (`variable_reader_coverage_20261116/PROOF.md`, `live_exterior_reader_20261116`) carries a state
`(I, Z, g)`; an `n`-bit address `b` followed by one final ordinary sink has the exact count

`H(b) = I + Z·(2D − rev_n(b)) + g·(2b − J(b)),   D = 2^n,   J(b) = Σ_{i<j set} (2^(j−i−1) + [j = i + 1]),`

with `(I, Z, g) = (−11, 82, 9)`, `(−9, 82, 9)`, `(−7, 82, 9)` for the three synchronized seeds: the **+2 and +4 moves** are
the seed swaps, which change one or three old arcs and leave the address untouched. A **block** `B` appended to all three
readers acts by `(I, Z, g) ↦ (h_B I, q_B Z, h_B g)` (near-sink block) or multiplies `Z` by `q_B` and adds `(q_B − h_B)g`
(ordinary sink block), where `h_B = H(B)` and `q_B` counts ordered two-slot path covers of `B` (BLOCKS.md). The seed
monoid of THM-5302 says the same: the corner entry, which carries the step, is multiplied by the other factor's count.

**A factor-1 shield is a block with `h_B = 1`, i.e. a transitive block.** It is the only kind of common block that keeps the
unit step (BLOCKS.md: "to preserve that gap with a common block, one must have `h_B = 1`"), and for a transitive block
of `s` vertices `q_B = 2^s` (the `2^s − 1` nonempty subsets each give one ordered split of the unique Hamiltonian path,
plus the empty slot). So a factor-1 shield is exactly a run of `s` ordinary-sink bits: it preserves the +2 and +4 moves,
multiplies the slope by `2^s`, and changes the radix of nothing. This is the `m` extra zeroes of SWITCH.md (`Z = 82·2^m`).
Thus **a factor-1 shield that still has a +2 move exists trivially — every zero run is one — and, by BLOCKS.md, no
factor-1 shield can do more**: `h_B = 1` forces `q_B` to be a power of two, so no transitive block supplies a new odd
radix or an independent additive coordinate. The search the owner asked for terminates here with a proof of
non-existence for the uniform-append architecture (BLOCKS.md's barrier, restated).

## 2. Address-level exact moves: a census (FINITE-EXACT)

Where CONTACT.md studies single-bit derivatives (`Δ = 18W − ZV − 9(e(X) + e(U))`, small only at the one site
`2i = n + m + 2`), we asked the complementary question: for which addresses `b` does **another address** `b′` give
`H(b′) = H(b) + 2` exactly — an exact +2 move by editing the word rather than the seed? Enumeration of all `2^n` addresses
at `(I, Z, g) = (−11, 82, 9)` (`scripts/address_moves.py`, receipts `address_moves_*.out`):

| `n` | addresses | with a `+2` partner | closest partner at Hamming distance 1 | partners that avoid the top and bottom `n//5` bits | cells (top, bottom `n//5` bits) reaching all residues mod 9 |
|---|---|---|---|---|---|
| 12 | 4,096 | 120 (2.9 %) | 32 (bit 7) | 46 | 16 of 16 |
| 16 | 65,536 | 1,560 (2.4 %) | 0 | 222 | 64 of 64 |
| 18 | 262,144 | 5,929 (2.3 %) | 256 (bit 10) | 894 | 64 of 64 |
| 20 | 1,048,576 | 23,184 (2.2 %) | 0 (the critical bit 11 gives +4 for 512 addresses) | 1,331 | 256 of 256 |
| 22 | 4,194,304 | 92,294 (2.2 %) | 0 (the critical bit 12 gives +8 for 1,024 addresses) | 6,080 | 256 of 256 |

Reading. (i) Single-bit exact small partners occur only at the critical bit `i = (n + 2)/2` of CONTACT.md, and the gap
they realise depends on the width: `+2` at `n = 12, 18` (32 and 256 addresses, one in 128 and one in 1024), `+4` at `n = 20`
(512), `+8` at `n = 16, 22` (128, 1,024); the census confirms the contact bound and shows that the small derivative is
attained exactly for a positive proportion of addresses, but not always with the gap 2 (pre-registered P-A1's second
clause is refuted at `n = 20, 22`). (ii) Multi-bit +2
partners (Hamming distance 3–9) exist for about one address in forty, and about one address in three hundred has a +2
partner that leaves both exterior phase blocks untouched — the kind of edit the live-exterior construction would need.
(iii) The gaps +4, +6, +8 behave alike (about 2 % each). (iv) Every phase cell reaches every residue modulo 9, in line
with EXTERIORS.md's mod-9 strengthening. These are counts of addresses with partners, not a selection rule: an exact
difference of two counts, `ΔH = −Z·Δrev(b) + 9·(2Δb − ΔJ(b))`, must satisfy `Δrev(b) ≡ −ΔH·82^(−1) ≡ −ΔH (mod 9)`, so a
+2 partner needs `Δrev ≡ 7 (mod 9)`, and `ΔJ` depends on the whole word (every pair of set bits contributes
`2^(distance − 1)`), which is why the edits that work are address-specific.

## 3. What a "different seam" would have to be

The count is two registers reading one word in opposite directions: `82·(2D − rev b)` big-endian and `9·(2b − J(b))`
little-endian, which is why a split address `b = x·2^q + rev_q(y)` gives a product form `S²AB + K + E` (EXTERIORS.md) with
two free exteriors and an error `O(N^(2/5))`. An exact target needs an additive coordinate of range at least that error
and step 2 that is independent of the product. The seed swaps give range 4; the gap theorem (SWITCH.md) gives any even
gap up to `V` but only at anchors of scale `V²` that move with the gap; blocks multiply. HYP-9842 (a) records the
absence of a third independent digit in the seed monoid's `3 × 3` unipotent form `[[h, b, c], [0, G, a], [0, 0, h]]`,
whose single corner `c` is the only additive register (its update `h₂c₁ + b₁a₂ + h₁c₂` is the cross term above). A
different seam therefore means a reader whose transfer is at least `4 × 4` unipotent — two corners — so that a second
additive register reads the word independently; the repository's rank-four return language
(`chanpak_eight_braid_20261116/RANK4_RETURN.md`, `RANK4_STEERING.md`) is the existing candidate. This is a precise
reformulation, not a construction; we did not build such a reader.

Typed verdicts: factor-1 shield with +2 — EXISTS (zero runs) and CANNOT give a new coordinate (PROVED, by BLOCKS.md);
address-level exact ±2 moves — EXIST for about 2 % of addresses, one in a thousand single-bit at the critical site
(FINITE-EXACT to `n = 22`; interior +2 edits for 0.13–0.34 % of addresses, pre-registered P-A2 confirmed); a third independent coordinate — OPEN, needs rank ≥ 4.
