# The first 32 odd numbers: 3 | 1 | 6 | 1 | 21, the two solo separators, and what the macrocosm does with them

Session klein-2026-11-16-micro (2026-11-16). Owner: *"Consider microcosmically the first 32 odd numbers split into
3, 1, 6, 1, 21 and note that these are all triangular numbers and pick out 7 and 21 as solo separators. Odd numbers
greater than 64 are apparently unaffected by the gaps at {7, 21} arranged in tripling chains, or at least they are
affected in a more periodic macrocosmic way."*

## 1. The split ([receipts/microcosm.out](receipts/microcosm.out))

The odd numbers below 64 with their least orders (the census to order 9 of THM-5266):

| block | members | size | least orders `n(h)` |
|---|---|---|---|
| 1 | 1, 3, 5 | 3 = T₂ | 1, 3, 4 (the rigid fibres of THM-5268) |
| separator | 7 | 1 = T₁ | hole |
| 2 | 9, 11, 13, 15, 17, 19 | 6 = T₃ | 5, 5, 5, 5, 6, 6 |
| separator | 21 | 1 = T₁ | hole |
| 3 | 23, …, 63 | 21 = T₆ | 6–8; 63 is the last, born at 8 |

All five block sizes are triangular (M-1 confirmed); the block sums 9, 7, 84, 21, 903 are triangular only for the
separator 21 and the last block (903 = T₄₂); the 32 numbers sum to 32² = 1024 = 2¹⁰. The two separators are the
`k = 1` near-joins with borrowed responses of THM-5268 (`7 = 1 + 2·3`, `21 = 3 + 2·9`), the first block is the three
rigid fibres whose responses are pure dyadic shifts, the second block begins with 9 (the first flexible fibre) and holds
the four sterile-only or fertile-only hard counts 11, 13, 17, 19, and the third block ends at 63 = 2⁶ − 1 = 7·9 = 3·21,
the Zsigmondy exception and the atom through both holes (wave 55). Arithmetic, typed as such.

## 2. The tripling chains heal at 63

`7 → 21 → 63 → 189 → 567 → 1701`: the first two are the holes, every later member is a count (63 at order 8, 189 = `H(P_7)`
at order 7, 567 at order 8, 1701 by order 9). Tripling is the ordinal sum with `C_3` (THM-5264 (E)), so a hole can never
propagate upward once a multiple is a count; `63 = 7·9` is a count because 9 is a count and 63 has strong realisers of
its own (its two first creators). Among the multiples of 7 below 200, the only non-counts are 7 and 21. More strongly:
**every odd number in [65, 2881] is the count of a tournament of order ≤ 9** — the only odd numbers below the order-9
maximum 2881 that are missing are 7 and 21. So below 2881 the "macrocosm" is unaffected in the literal sense.

## 3. Is the macrocosm affected periodically? ([receipts/residues9.out](receipts/residues9.out))

The near-join family from bases ≤ 9 (THM-5265: complete to 398,129) misses 4.5 % of the odd numbers in `[2^19, 2^21)` and
43 % in `[2^21, 2^23)`. By residue class the miss rate is flat: relative rates 1.01 and 0.99 for the classes 7 and 21
modulo 64 (ranks 14 and 16 of 32), 1.00 against 0.98 for 7 against 1 modulo 8, within 1 % across classes modulo 7, and
the only signal is modulo 3, where the multiples of 3 are missed 5 % less often (the closure `3·NJ ⊆ NJ`). M-2 refuted:
the holes do not recur as residue classes in the near-join macrocosm. What *is* periodic above 64 is the concurrent
packet's `mod 18` phase for the ready `7p` cores (`ready_core_growth_and_mod18_20261117`), a property of one
construction's addresses, not of the counts; and THM-5267's seam-centre misses, all below 263 and none `≡ 2 (mod 8)`.
With bases ≤ 10 the family reaches 19,190,825 and its 1,197 misses below 2^25 all lie in `[2^24, 2^25)`, again flat by
residue: 33 misses in the class 7 and 44 in the class 21 modulo 64 against 37 expected, the top classes being 19, 45, 29
([receipts/residues10.out](receipts/residues10.out)). The macrocosm is unaffected by the separators except through the
lane `7p` (THM-4094), which starts at `p = 5`.
