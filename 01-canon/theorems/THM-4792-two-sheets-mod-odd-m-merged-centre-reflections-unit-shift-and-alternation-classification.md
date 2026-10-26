---
id: THM-4792
title: "The two sheets modulo an odd modulus, merged. (a) Centre: for odd m = 2k+1 the centre labels of the mod-m multiplication table are k = -1/2 and k+1 = 1/2, the fixed points of 3x+1 and 3x-1; the centre products k^2 = (k+1)^2 = 1/4 and k(k+1) = -1/4 form the 2-cycle 1/4 -> -1/4 -> 1/4 of the alternating pair (3x-1, 3x+1) (2026-10-05 Prop. 1). (b) Symmetries of the one-step maps floor(3x/2) = (3x-c)/2 and ceil(3x/2) = (3x+c)/2 (c = x mod 2): x -> 1-x is a symmetry of the first and x -> -1-x of the second (reflection about the centre labels); negation exchanges them; the unit translation x -> x-1, which is the root-27 section change z = x-1, conjugates the first into the second; it does not conjugate the accelerated Syracuse maps. (c) Classification: the halving-shift alternation of the modular edge tables, E^-_(k+t)(m) = {(-x, y) : (x, y) in E^+_k(m)} for all k (kind-pasteur Theorem A), exists iff -1 is a power of 2 mod m, iff every prime p | m has ord_p(2) even with one common 2-adic valuation, iff the two centre values +-1/4 lie in one coset of <2>. Primes p = 3 (mod 8) have valuation 1, p = 5 (mod 8) valuation 2, p = 7 (mod 8) valuation 0 (never alternate); 15 = 3*5 is the least odd modulus whose prime factors alternate separately but not jointly."
status: >
  PROVED (elementary, proofs below). (a) is 2026-10-05 Prop. 1 (opus) with the fixed-point reading added; (c)'s first
  and last equivalences restate kind-pasteur-2026-10-25-S4 Theorem A.3-A.4, whose open exercise (the classification of
  the alternating moduli) is the middle equivalence proved here; (b) records the identities behind the root-27 section
  change (root27_machine_bridges_20261026 anchor, T(z) = ceil(3z/2)) and THM-4474 (E).
  FINITE-EXACT controls: sheet_merge_20261026.py (identities (b) for |x| <= 1e5; (a) and all four conditions of (c),
  including the literal table shift for m <= 301, for the 1000 odd m in [3, 2001], no mismatch).
session: opus-2026-10-26-S57
source: 05-knowledge/results/s57_untouchables_root27_20261026/README.md
depends_on:
  - 01-canon/theorems/THM-4474-strategy-cube-provability-by-cycle-densities.md
related:
  - 05-knowledge/results/collatz_spectrum_equality_and_stitched_sheets_20261005.md
  - 05-knowledge/results/sheets_root27_rank_untouchables_20261025.md
  - 05-knowledge/results/collatz_two_sheet_receipts_20261005.md
  - 05-knowledge/results/root27_machine_bridges_20261026/anchor/README.md
scripts:
  - 04-computation/untouchables_s57/sheet_merge_20261026.py (+ sheet_merge.json)
  - 04-computation/experiments/sheet_tables_mod_odd_kps_20261025.py (kind-pasteur, Theorem A controls)
---

# THM-4792: the two sheets modulo an odd modulus, merged

The owner asked to merge three threads. The 2026-10-05 stitched-sheets note found that "tables mod odd `n`
alternate between `3n-1` and `3n+1` edges": the centre of every odd multiplication table is `+-4^-1`, swapped by the
two affine steps. kind-pasteur's Theorem A (2026-10-25) related the edge tables by negation and a halving shift. The
root-27 packet changes section by `z = x - 1`. This theorem puts them in one statement.

**(a) Centre.** Let `m = 2k+1`. Then `2k = -1`, so `k = -1/2` and `k+1 = 1/2 (mod m)`. These are the fixed points of
`x -> 3x+1` and `x -> 3x-1` respectively. Their products are `k^2 = (k+1)^2 = 1/4` and `k(k+1) = -1/4`. The
alternating pair maps `1/4 -> 3/4 - 1 = -1/4 -> -3/4 + 1 = 1/4`. Over `Q`, this is the rational 2-cycle of the
sign conjugacy `U_-(x) = -U_+(-x)` (THM-4474 (E)), and the table centre is its reduction mod `m`.

**(b) One-step maps.** For integers `x`,

    floor(3(1-x)/2) = 1 - floor(3x/2),      ceil(3(-1-x)/2) = -1 - ceil(3x/2),
    -floor(-3x/2) = ceil(3x/2),             ceil(3(x-1)/2) = floor(3x/2) - 1.

*Proof.* For an integer `y`, `ceil((y-1)/2) = floor(y/2)` and `floor((y+1)/2) = ceil(y/2)`. The first identity is
`floor((3-3x)/2) = 1 + floor((1-3x)/2) = 1 - ceil((3x-1)/2) = 1 - floor(3x/2)`. The second is the first conjugated
by negation, which is the third identity. The fourth is `floor((3z+3)/2) - 1 = floor((3z+1)/2) = ceil(3z/2)` with
`z = x - 1`. `QED`

So the 3x-1 sheet's one-step map `floor(3x/2) = (3x - c)/2` (Antihydra's macro on `x`) and the 3x+1 sheet's
`ceil(3z/2) = (3z + c)/2` (the root-27 packet's `T`) are conjugate by the unit translation. That translation is the
packet's change of section; the source `x = 8 -> 12 -> 18 -> 27` becomes `z = 7 -> 11 -> 17 -> 26`. Each map is
symmetric under reflection about its own centre label, `1/2` and `-1/2`, and negation swaps the two maps and the two
centres. The unit translation does not commute with full halving, so it does not conjugate the accelerated maps
`(3x +- 1)/2^v` (kind-pasteur check T1b); negation does, which is THM-4474 (E).

**(c) Classification of the alternation.** Write `E^+-_k(m) = {(x, (3x +- 1) 2^-k) : x in Z/m}`.

1. Table shift. `E^-_(k+t) = {(-x, y) : (x,y) in E^+_k}` for all `k` holds iff `2^-t = -1 (mod m)`; equality of the
   tables forces it by comparing two sources (kind-pasteur A.3–A.4).
2. Prime orders. By CRT, `2^t = -1 (mod m)` iff `2^t = -1 (mod p^a)` for every `p^a || m`. Since `(Z/p^a)^*` is
   cyclic of even order with unique involution `-1`, this holds iff `o_p = ord_(p^a)(2)` is even and
   `t = o_p/2 (mod o_p)`. Here `o_p = ord_p(2) p^j` has the parity and 2-adic valuation of `ord_p(2)`. The congruences
   `t = o_p/2 (mod o_p)` are compatible iff `o_p/2 = o_q/2` modulo `gcd(o_p, o_q)` for all `p, q`.
   - If `v_2(o_p) = v_2(o_q) = e`, both halves are `2^(e-1) (mod 2^e)` and `0` modulo the odd part of the gcd, so
     they agree.
   - If `v_2(o_p) < v_2(o_q)`, then modulo `2^(v_2(o_p))` the first half is `2^(v_2(o_p)-1)` and the second is `0`,
     so they disagree.
   Hence `-1 in <2> (mod m)` iff `ord_p(2)` is even for every `p | m`, with one common value of `v_2(ord_p(2))`.
3. Centre coset. `-1/4 in (1/4)<2>` iff `-1 in <2>`, since `1/4` is a unit.
4. Residue classes. If `p = 3 (mod 8)`, then 2 is a non-residue and `p - 1 = 2 * odd`, so `v_2(ord_p 2) = 1`. If
   `p = 5 (mod 8)`, then 2 is a non-residue and `p - 1 = 4 * odd`, so the valuation is 2. If `p = 7 (mod 8)`, then 2 is
   a residue and `(p-1)/2` is odd, so the order is odd. Primes `p = 1 (mod 8)` vary: 17 has 3, 41 has 2, 73 has 0,
   97 has 4. The alternating odd moduli are exactly the products of primes from one class `v >= 1`. Up to 101 they
   are `3, 5, 9, 11, 13, 17, 19, 25, 27, 29, 33, 37, 41, 43, 53, 57, 59, 61, 65, 67, 81, 83, 97, 99, 101`. `QED`

**The owner's split `{2,3,5,7} = {3,5} | {2,7}`, typed.** 3 and 5 alternate separately (`t = 1, 2`) but lie in
different classes, so `15` does not alternate; it is the least odd modulus of that kind. 7 is the least prime that
never alternates. 2 is the halving prime itself, and there are no tables mod an even modulus. Two further exact
facts:
- `sigma(14) = (1+2)(1+7) = 24 = (1+3)(1+5) = sigma(15)`, which makes `{2,7} | {3,5}` the unique split of the four
  primes into pairs with equal `prod (1+p)`;
- `(14, 15)` is the closest complementary divisor pair of 210 (session66 divisor bridge).

The untouchable residue profile is nearly multiplicative over `210 = 14 * 15` (S57 U-P2). Beyond these exact
statements no structure is claimed (MISTAKE-529).
