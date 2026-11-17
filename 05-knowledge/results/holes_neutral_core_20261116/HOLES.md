# 1, 3, 7, 21, the even prime between them, and 18 — in the grammar

Session klein-2026-11-16-core (2026-11-16). Owner: *"With the remaining infinite problem reduced to odd primes p and
values 7p, … study the behavior of the first two odd numbers, 1 and 3, and 7 times them, 7 and 21, and how they relate to
the only even prime sandwiched between them, and the number 18."* Concurrent (other machine, same day):
`ready_core_growth_and_mod18_20261117` reads 18 as the continuation coefficient `2E = 18` of the 15/17 rooted pair and
proves the phase `H ≡ 5 (mod 18)`, `p ≡ 11 (mod 18)` for the `g = 5` matched-ready `7p` lane; `one_port_neutral_continuations_20261117`
defines neutral addresses on an appended tape with at most two runs. Nothing below repeats those; the readings here are
the grammar's (THM-5266, THM-5267).

## 1. What the six numbers are in the count monoid

- **1 and 3** are the counts of the two rigid fibres: `F_1` = transitive tournaments, `F_3` = padded 3-cycles. Their
  response sets are the pure dyadic shifts `𝔅(1) = {2^j}`, `𝔅(3) = {3·2^j}` (THM-5266 (C)), and neither fibre admits a
  count-preserving move other than the trivial flips of adjacent ranks ([NEUTRAL.md](NEUTRAL.md)).
- **7 = 1 + 2·3 and 21 = 3 + 2·9.** The holes are the `k = 1` near-joins of the two rigid bases with a response borrowed
  from the next fibre: 3 is the response of every vertex of the 3-cycle (`𝔅(3) ∋ 3`, `𝔅(1) ∌ 3`), 9 the response of the
  vertices in the second block of `C_3 ⇒ C_3` (`𝔅(9) ∋ 9`, `𝔅(3) ∌ 9`). The pencils with the borrowed response miss
  only their first term: `1 + 2^k·3 = 7, 13, 25, 49, 97, 193` and `3 + 2^k·9 = 21, 39, 75, 147, 291, 579` are all counts
  for `k ≥ 2` ([receipts/holes_arith.out](receipts/holes_arith.out)); 13 and 25 are reached by other bases (13 only by
  sterile strong tournaments, THM-5263).
- **2**, the only even prime, is the `k = 1` doubling that both holes need (`2·3 = 6`, `2·9 = 18`) and the only count the
  near-join can never produce (counts are odd, Rédei). It is the centre of the first seam `(1, 3)` (THM-5267), and
  `7·(1, 3) = (7, 21)` is the seam with multiplier 7 whose centre 14 has both ends missing.
- **18 = 2·9 = 21 − 3**: the shift from the count 3 to the hole 21 is twice the borrowed response, as `6 = 7 − 1 = 2·3`
  is for the other hole. In the `7p` lane, `7p = p + 6p`: the holes are the two members whose shift `6p ∈ {6, 18}` is
  `2·(3p)` with `3p` the response of the wrong fibre; from `p = 5` on, `7p` is a count (35, 49, 63, 77, 91, …).

## 2. What is proved, what is read

- THM-1370 excludes 7 and 21 by the strong floors and the finite core (the concurrent packet says the same: "not by
  congruence"); THM-5266 (C) gives the additive sub-grammar's intrinsic reason, which is exactly the borrowed-response
  reading above: `7 − 1`, `21 − 3` are `2^k B` only for responses the rigid fibres do not have. PROVED.
- The lane reduction to primes and `7p` is THM-4094; `7p` must come from a strong component of count `7p`
  (its Lemma 5.2). KNOWN.
- "18 as `2E` of the 15/17 pair" and the `mod 18` phase are the concurrent packet's; the `2·9` reading here is
  arithmetic on the same number — two readings, no identification claimed. TYPED: arithmetic.
