---
id: THM-4641
title: "The Wilson / {2,4,8} / Fermat-multiplier dictionary, with every row typed: {2,4,8} are the powers of 2 adjacent to a power of 3 (Gersonides); they index the two-position relations of the repair algebra at -1 together with its padding relation (2 = 3 - 1, unit carry 4 = 3 + 1, Catalan carry 9 - 8 = 1), which repair descents but do not generate the collisions, and the integer cycles of unit denominator 0, -1, +1, {-5,-7}; Wilson's -1 is the product of a group's involutions, and for the Collatz multiplier's orbit <2> mod p the product is -1 iff ord_p(2) is even; the Wilson quotient is the sum of Fermat quotients (Lerch) and q_p(2^j) = j q_p(2), so the Fermat quotients of 2, 4, 8 are 1:2:3; the tournament count is T(p) = chi_8(p)(q_p(2)/2 - w_p) - A_p (mod p), so the Wilson and Fermat quotients of 2 enter through the same Fermat multiplier chi_8(p) = (2/p). The remaining numerical links are coincidences with stated baselines: 5 = 2^3 - 3 and 13 = 2^4 - 3 are the only Wilson primes among the primes 2^a - 3 with a <= 44; 2 is a primitive root mod 5, 13, 563 (baseline 0.052); p - 2 = 3, 11, 561 is absolutely Fermat (baseline about 0.09 with size-dependent rates); 563 is not a numerator at the -1 anchor"
status: >
  PROVED: rows 3, 5, the unit-denominator part of row 2b, and the reductions in rows 2 and 7 (elementary; row 5 upgrades THM-307's verified relation between
  the Euler and Fermat quotients to a proof and composes it with THM-307). KNOWN: rows 1, 3 (Wilson for finite abelian
  groups), 4 (Lerch 1905; logarithmic property of Fermat quotients) and the Wilson-prime search (Costa, Gerbicz and Harvey
  2014). FINITE-EXACT: the cycle census of row 2b and rows 6-10 (exact arithmetic; wilson_248_fermat_dictionary.py, 1,835 checks; THM-307's tournament
  decomposition re-derived from exact Burnside counts for every prime p < 48; audit A rechecked it for p <= 61 with A_p summed
  over partitions). HEURISTIC: row 11. Independently audited 2026-10-14 (audit A,
  04-computation/experiments/repair_algebra_20261014/audit_A/REPORT.md: CORRECT WITH FIXES; corrections applied: the carries are
  the two-position relations, not generators; the anchors -1, 1, 1/5 attached to {2,4,8} removed as numerology; the row-8
  baseline replaced by size-dependent rates; the count of non-Wilson primes 2^a - 3; Thakur's theorem stated for every q
  and not as support for the log log x count; Miller 1903 cited; MISTAKE-606).
session: mac-mini-2026-10-14-wilson
source: 05-knowledge/results/repair_algebra_minus_one_and_wilson_dictionary_20261014.md
scripts:
  - 04-computation/experiments/repair_algebra_20261014/wilson_248_fermat_dictionary.py (+ .out)
related:
  - THM-4640 (the repair algebra at -1), THM-4600 (anchors), THM-307 (Burnside-QR decomposition of T(p) mod p)
  - 05-knowledge/results/wilson_defect_trees_20261014c.md (Wilson lift trees; Lerch via the Teichmueller branch)
  - Gersonides (1343); Lerch (1905); Costa, Gerbicz and Harvey, Math. Comp. 83 (2014); Thakur, Wilson primes for F_q[t]
---

# THM-4641 — the Wilson / {2,4,8} / Fermat-multiplier dictionary

Notation:
* `W_p = ((p−1)! + 1)/p` and `w_p = W_p mod p` (the Wilson quotient).
* `q_p(a) = (a^(p−1) − 1)/p mod p` (the Fermat quotient).
* `f_p = (2^((p−1)/2) − (2/p))/p mod p` (the Euler quotient of 2).
* `χ_8(p) = (2/p)`.
* `T(p)` is the number of tournaments on `p` vertices up to isomorphism, and `A_p` is THM-307's remaining Burnside sum.

## The dictionary

| # | row | type |
|---|---|---|
| 1 | `|2^a − 3^b| = 1` (`a ≥ 1`) exactly at `(1,0), (1,1), (2,1), (3,2)`. The powers of 2 adjacent to a power of 3 are `{2, 4, 8}`. | KNOWN (Gersonides) |
| 2 | The three adjacencies `2 = 3 − 1`, `4 = 3 + 1`, `8 = 9 − 1` are exactly the padding relation `f_1(−1) = −1` (`(1)·s ~ s`) and the two two-position relations of the repair algebra at −1, the unit carry and the Catalan carry (THM-4640 (5)). The carries repair descending sequences but do not generate the primitive collisions (THM-4640 5(d)). The first-generation numerators at −1 are the reciprocal anchors `K((a)) = 2^a − 3 = 1/c_a`. | PROVED |
| 2b | A Collatz cycle with `k` odd steps and `A` halvings has cycle points `b/(2^A − 3^k)`, where `b` depends on the word. Gersonides' solutions are exactly the unit denominators `|2^A − 3^k| = 1`: `(A, k) = (1, 0)` gives `0`, `(1, 1)` gives `−1`, `(2, 1)` gives `+1`, `(3, 2)` gives the 2-cycle `{−5, −7}`. Every word with these `(A, k)` therefore closes on an integer cycle. So `{2, 4, 8}` indexes the integer cycles that exist for denominator reasons alone. The only other known integer cycle, through `−17` (7 odd steps, 11 halvings), needs `139 = |2^11 − 3^7|` to divide its numerator. Direct iteration over all odd `|n| ≤ 100001` finds exactly these cycles. | PROVED + FINITE-EXACT |
| 3 | In a finite abelian group the product of all elements is the product of its involutions (Miller 1903). In `F_p^*` this is Wilson's `−1`. For the orbit of the Collatz multiplier, `∏_(j<d) 2^j = 2^(d(d−1)/2) ≡ −1 (mod p)` iff `d = ord_p(2)` is even, and `≡ 1` otherwise. | KNOWN / PROVED |
| 4 | Lerch: `w_p ≡ Σ_(a=1)^(p−1) q_p(a) (mod p)`. Since `q_p` is a logarithm, `q_p(ab) ≡ q_p(a) + q_p(b)`, so `(q_p(2), q_p(4), q_p(8)) ≡ (1, 2, 3)·q_p(2)`. | KNOWN |
| 5 | `f_p ≡ χ_8(p)·q_p(2)/2 (mod p)`. Composed with THM-307, `T(p) ≡ χ_8(p)·(q_p(2)/2 − w_p) − A_p (mod p)`. The Wilson quotient and the Fermat quotient of 2 enter the tournament count through the same Fermat multiplier `χ_8(p)`. At Wilson primes `T(p) ≡ χ_8(p) q_p(2)/2 − A_p`; at Wieferich primes `T(p) ≡ −χ_8(p) w_p − A_p`. | PROVED (re-derived for `p < 48` from exact Burnside counts) |
| 6 | The anchor row: the primes `2^a − 3` with `a ≤ 24` are `5, 13, 29, 61, 509, 1021, 4093, 16381, 1048573, 4194301, 16777213`. Only `5 = 2^3 − 3` and `13 = 2^4 − 3` are Wilson. By the cited search there is no other Wilson prime below `2·10^13`, so none of the form `2^a − 3` with `a ≤ 44`. `563` is not of this form. | FINITE-EXACT + CITED; COINCIDENCE |
| 7 | 2 is a primitive root mod 5, 13 and 563. Every `2^a − 3` with `a ≥ 3` is `≡ 5 (mod 8)`, so 2 is a non-residue there and `4 | ord_p(2)` (PROVED). This partly explains 5 and 13. The frequency of "2 primitive" among primes `< 10^6` is 0.3738 (Artin's constant 0.3740), so all three by chance has probability about 0.052. | COINCIDENCE |
| 8 | `p − 2 ∈ {3, 11, 561}` is absolutely Fermat (`a^n ≡ a (mod n)` for all `a`): 3 and 11 are prime, and `561 = 3·11·17` is the least Carmichael number. For primes `5 ≤ p < 10^6`, `p − 2` is prime or Carmichael with frequency 0.1041, mostly twin primes. That global average understates small primes: the Hardy–Littlewood rates `2C_2/ln p` at 5, 13, 563 are 0.82, 0.51, 0.21, so the pattern has chance about 0.09 (audit A). The Carmichael numbers `c < 10^6` with `c + 2` prime are 561, 2465 and 656601. | COINCIDENCE |
| 9 | At the −1 anchor, 5 and 13 are numerators only as anchors (`(a)` and its root twin `(1, a)`). `563, 561, 1093, 3511` have no increasing `{2,3}`-representation. The search is complete: by THM-4640 (1), every canonical representation of `K` has exponents `≤ log_2(K + 3)`. | FINITE-EXACT |
| 10 | `684106401 = 3·11·17·401·3041` is Carmichael with `λ = 15200` and `λ + 1 = 27·563`. This is a single coincidence, observed in `carmichael_successor_blocks_20261014b.md`. | FINITE-EXACT; COINCIDENCE |
| 11 | If `W_p mod p` is equidistributed, the expected number of Wilson primes in `[x, y]` is about `log log y − log log x`: infinitely many in total, and about 1.6 expected in `(563, 2·10^13]`, where none was found (probability `e^(−1.6) ≈ 0.2`). Thakur, and Sauerberg–Shu–Thakur–Todd, proved there are infinitely many Wilson primes in `F_q[t]` for every `q` (the criterion `P'' = 0` is for odd characteristic). That criterion is algebraic: over `F_p` it rules out every degree `d` with `2 ≤ d < p`. So the function-field theorem does not support the `log log x` count; it only shows that finiteness is not forced by structure. | HEURISTIC + CITED |

## Proofs

**Row 2.** This is THM-4640 (1) and (5). `f_a(c_a) = c_a` for `c_a = 1/(2^a − 3)`.

**Row 2b.**
* A word `w` with `k` odd steps and total `A` has the unique fixed point `x = b_w/(2^A − 3^k)`, with `b_w = Σ_j 3^(k−1−j) 2^(S_j)`. It is an integer whenever `|2^A − 3^k| = 1`.
* The four solutions are Gersonides'. Enumerating their words gives `−1`, `+1` and `{−5, −7}`; the empty word (`k = 0`, `A = 1`) gives `0`.
* The census iterates `U` from every odd `|n| ≤ 100001` until a repeat and records each cycle's `(k, A)`.

**Row 3.**
* Inverse pairs cancel, so the product is that of the elementary abelian 2-subgroup of involutions. It is the unique involution if that subgroup has order 2, and 1 otherwise.
* The cyclic group `⟨2⟩` has an involution, namely `−1`, iff its order `d` is even.

**Row 4.** These are classical. The logarithmic property is `(ab)^(p−1) = (1 + p q(a))(1 + p q(b)) ≡ 1 + p(q(a) + q(b)) (mod p^2)`.

**Row 5.**
* `2^(p−1) = (2^((p−1)/2))^2 = ((2/p) + p f_p)^2 ≡ 1 + 2(2/p) p f_p (mod p^2)`. Hence `q_p(2) ≡ 2(2/p) f_p`, i.e. `f_p ≡ (2/p) q_p(2)/2`, using `(2/p)^2 = 1`.
* Substituting into THM-307's `T(p) ≡ f_p − (2/p) w_p − A_p` gives the stated form.
* The script recomputes `T(n)` for `n ≤ 12` from Burnside's formula against A000568, and for every prime `p < 48` it recomputes `A_p` and checks the decomposition exactly.

**Row 7.** For `a ≥ 3`, `2^a ≡ 0 (mod 8)`, so `2^a − 3 ≡ 5 (mod 8)` and `(2/p) = −1`. The 2-part of `ord_p(2)` is then that of `p − 1`, namely 4.

**Rows 6, 8–10.** Exact computation (`wilson_248_fermat_dictionary.py`).
* Wilson quotients are computed by multiplying `1, …, p−1` modulo `p^2`.
* Korselt's criterion is checked by factoring, and also by the definition `a^n ≡ a (mod n)` for every residue `a`.
* Representations are found by an exhaustive backward search. Its pruning `P ≥ −3^(e−1)` holds for every numerator whose exponents are below `e`.

## Reading

* **What is structure.**
  * `{2, 4, 8}` is not numerology: it is the complete list of powers of 2 next to a power of 3. Inside Collatz it plays two exact roles:
    * it is the set of two-term repairs at −1 (THM-4640);
    * it indexes the integer cycles that exist for denominator reasons alone: `−1`, `+1`, `{−5, −7}` (and `0`).
  * The padding relation and the two carries are these cycles' defining identities: `2 = 3 − 1` for the padding and the anchor `−1`, `4 = 3 + 1` for the unit carry and the `+1` cycle, and `9 − 8 = 1` for the Catalan carry and `{−5, −7}`.
  * The Wilson–Fermat link is classical and exact (Lerch). The tournament count sees both quotients through the single character `χ_8(p)`, which is the clean meaning of "Fermat multiplier" here.
* **What is coincidence.**
  * That the first two Wilson primes are `2^3 − 3` and `2^4 − 3` comes from small numbers. The other nine primes of that form with `a ≤ 24` are not Wilson, and neither is `2^29 − 3`, the only further one with `a ≤ 44`.
  * 563 has no representation at −1 at all.
  * The primitive-root and `p − 2` patterns have chance probabilities of about 5% and about 9% (size-dependent rates).
* **On "{5, 13, 563} is the whole list".**
  * Nothing in this dictionary, or in any known structure, forces `W_p ≡ 0 (mod p)` to stop.
  * The standard model predicts infinitely many Wilson primes, growing like `log log x`. In `F_q[t]` infinitude is a theorem, though by an algebraic criterion that does not test the count.
  * A proof that the list is complete is therefore not a realistic target: the statement is expected to be false.
