---
id: THM-4855
title: "Circulant tournaments as S-sequencings. (1) A Hamiltonian path of Circ(Z_q, S), with 1 in S, is a word over S of length q-1 with distinct partial sums mod q. Its positional reversal is an anti-automorphism iff every segment sum of the word and of its reversal lie in S together (segment-sum criterion, valid for every q). (2) For prime q this is the same as the difference word being a palindrome. (3) For every odd q >= 7, prime or composite, explicit words A1, A2, B, D give a non-palindromic Hamiltonian path, each with a named failing segment. This is an all-order arc-level proof along THM-4730 (B)'s original explicit route (three-interval with b = 2, descending pairs) plus a new family B. (4) For q = 3, 5 every S-sequencing is a palindrome; T5's are 1111, 2112, 2222, its three tilings. So the strong points of symmetry K1, C3, T5 are exactly the circulants of order <= 5, i.e. THM-4849's tournament cycles of period <= 5."
status: >
  PRE-AUDIT (mac-mini-2026-11-01-fractranTM).
  PROVED (elementary): (1), (2), (3), (4).
  FINITE-EXACT:
  - the segment-sum criterion agrees with the literal arc test on all 1,741,634 Hamiltonian paths of the circulants
    with q <= 13;
  - the A/B/D witnesses are valid and non-palindromic for all 131,071 connection sets with 1 in S at every odd q from
    7 to 35.
  CITED:
  - THM-4730 (A) (rotation lemma) and the audited strong-cycle reduction: a strong pure-blue tournament of odd order
    is a circulant along every Hamiltonian cycle; no even order;
  - affine isomorphisms of prime-order circulants (Turner 1967; Babai 1977; as used in THM-4849);
  - THM-4849 (A) (connection set read as a parity word).
  NOT NEW: the alphabet theorem itself. It is already proved at all orders by the audited prefix-rerouting argument
  (`pure_blue_atom_classification_20261020/UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md`). This file gives a second,
  word-level proof and closes THM-4730's open item, which said no all-order repair of its explicit paths was claimed.
source: mac-mini-2026-11-01-fractranTM (owner directive of 2026-11-01)
depends_on:
  - 01-canon/theorems/THM-4730-strong-pure-blue-alphabet-complete-points-of-symmetry-are-k1-c3-t5.md
  - 05-knowledge/results/merged_metagraph_purity_and_strong_core_20261020/PURE_BLUE_STRONG_CYCLE_REDUCTION.md
related:
  - 01-canon/theorems/THM-4849-circulant-tournaments-are-rational-collatz-cycles-shift-decimation-reversal-two-base-points-and-the-locking-theorem-behind-35-and-175.md
  - 01-canon/theorems/THM-2454-pure-blue-alphabet-and-center-law.md
  - 05-knowledge/results/pure_blue_atom_classification_20261020/UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md
  - 05-knowledge/results/fractran_turing_circulant_goldbach_20261101.md
scripts:
  - 04-computation/experiments/circulant_fractran_goldbach_20261101/circulant_sequencings.py
  - 04-computation/experiments/circulant_fractran_goldbach_20261101/witness_families.py
outputs:
  - 04-computation/experiments/circulant_fractran_goldbach_20261101/circulant_sequencings.out
  - 04-computation/experiments/circulant_fractran_goldbach_20261101/witness_families.out
---

# THM-4855 - circulant tournaments as S-sequencings

## Setting

Let q be odd. Let S ⊂ Z_q \ {0} contain exactly one element of each pair {x, -x}, and let T = Circ(Z_q, S), with
u -> v iff v - u ∈ S. Since x -> -x is an isomorphism onto Circ(Z_q, -S), we may assume 1 ∈ S. Write k = (q-1)/2,
and let χ be the indicator of S on 1, ..., q-1, so that χ(q - v) = 1 - χ(v).

A Hamiltonian path (v_1, ..., v_q) of T starting at v_1 = 0 is determined by its difference word
w = (d_1, ..., d_{q-1}) with d_i = v_{i+1} - v_i ∈ S. Its partial sums are distinct mod q. Call such a word an
**S-sequencing**. Every Hamiltonian path is a translate of exactly one path from 0, so H(T) = q · #S-sequencings.

## Statement

**(1) Segment-sum criterion.** Write r = rev(w) = (d_{q-1}, ..., d_1). The positional reversal ρ(v_i) = v_{q+1-i}
is an anti-automorphism of T iff

```text
for every segment [i, j] of positions 1..q-1:   sum_{t=i}^{j} w_t ∈ S   <=>   sum_{t=i}^{j} r_t ∈ S.
```

That is, iff the path is palindromic (blue) in the sense of THM-643 and THM-4730. A palindromic word always
satisfies the criterion.

**(2) Prime order.** If q is prime, a Hamiltonian path is palindromic iff its difference word is a palindrome.

**(3) Explicit witnesses, every odd q >= 7.** With 1 ∈ S, call u ∈ [1, q-2] a *rise* if u ∉ S and u + 1 ∈ S. Then
T has the following non-palindromic Hamiltonian path.

| case | condition | S-sequencing w | failing segment: sum on w, sum on rev(w) |
|---|---|---|---|
| A1 | 2 ∈ S and some s ∈ S, 4 <= s <= k+1, has s-1 ∉ S | 1^(s-2) 2 1^(q-1-s) s | suffix [s, q-1]: -1 ∉ S; 1-s ∈ S |
| A2 | 2 ∈ S, no such s (then S = {1, ..., k}) | 1 2 1^(q-4) 3 | prefix [1, k-1]: k ∈ S; k+1 ∉ S |
| B | 2 ∉ S and a rise u ≠ k | 1^(u-1) (q-u) (u+1) 1^(q-2-u) | [u, u+1]: 1 ∈ S; 2 ∉ S |
| D | 2 ∉ S and the only rise is k (then S = {1} ∪ [k+1, q-2]) | 1 ((q-3) 1)^(k-1) (q-2) | [1, 2]: -2 ∈ S; -1 ∉ S |

So no circulant tournament of odd order q >= 7 is pure-blue. A1 and A2 are THM-4730 (B)'s three-interval path with
b = 2. D is its descending-pairs path. B is new.

**(4) Orders 3 and 5, and the points of symmetry.**
- For q = 3 the only S-sequencing is (1, 1).
- For q = 5, S = {1, 2}, the S-sequencings are exactly 1111, 2112 and 2222. All are palindromes, so T5 is pure-blue,
  and these three words are its three tilings (fibre H/|Aut| = 15/5 = 3).
- With THM-4730 (A) and the strong-cycle reduction (a strong pure-blue tournament has odd order and is a circulant),
  the strong points of symmetry are exactly the circulant tournaments of order <= 5: K1, C3, T5.
- Under THM-4849 (A) they are the tournament cycles of period q <= 5:
  - {0} for q = 1;
  - the 3x+5 cycle {1, 2, 4} for q = 3;
  - T5's two 3x+23 cycles {5, 19, 40, 20, 10} and {7, 11, 22, 28, 14} for q = 5 (already in the cycle101 note).
- An S-sequencing is an ordering of the cycle's points by jumps T^s, s ∈ S. For example, T5's word 2112 visits
  10, 19, 40, 20, 5.

## Proof

**(1).** For i < j:
- v_i -> v_j iff v_j - v_i = sum_{t=i}^{j-1} w_t ∈ S;
- ρ(v_j) -> ρ(v_i) iff v_{q+1-i} - v_{q+1-j} = sum_{t=q+1-j}^{q-i} w_t ∈ S;
- substituting t = q - t', the second sum is sum_{t'=i}^{j-1} r_{t'}.

ρ is an anti-automorphism iff u -> v <=> ρ(v) -> ρ(u) for all u ≠ v. Since T is a tournament, it suffices to take
i < j.

**(2).** For prime q, every isomorphism between circulants on Z_q is affine (Turner; Babai). So an
anti-automorphism is x -> ax + c with aS = -S. Then:
- ρ is an involution, so a^2 = 1, so a = -1 (a = 1 would give a translation, not an anti-automorphism);
- ρ = (x -> c - x) on the path means v_{q+1-i} = c - v_i for all i, that is, d_{q-i} = d_i;
- so the word is a palindrome. The converse is (1).

**(3).** Each word is an S-sequencing. Its partial sums cover Z_q:
- **A:** [0, s-2], then s, ..., q-1, then q-1+s ≡ s-1.
- **B:** [0, u-1], then -1, then [u, q-2].
- **D:** 0, 1, then the pairs (q-2, q-1), (q-4, q-3), ..., (3, 4), then 2.

Its letters lie in S:
- B uses q - u ∈ S because u ∉ S.
- D uses q - 3 ∈ [k+1, q-2] because q >= 7.

The failing segments.
- **A1.** On w, the suffix [s, q-1] is 1^(q-1-s) s, with sum q - 1 ≡ -1 ∉ S. On r = s 1^(q-1-s) 2 1^(s-2):
  - the first s - 1 letters are s and s - 2 ones, because s - 1 <= q - s;
  - so the same suffix has sum (s-1) - (2s-2) = 1 - s;
  - this lies in S since s - 1 ∉ S.
- **A2.**
  - If s - 1 ∈ S for every s ∈ S ∩ [3, k+1], then S ∩ [1, k+1] is an initial interval.
  - k+1 ∈ S would force k ∈ S, which contradicts antisymmetry, so S = {1, ..., k}.
  - For w = 1 2 1^(q-4) 3, the prefix sums are j + 1 (2 <= j <= q-2); for its reversal they are j + 2
    (1 <= j <= q-3).
  - At j = k-1 the prefix sums are k ∈ S and k+1 ∉ S.
- **B.**
  - On w, [u, u+1] carries (q-u) + (u+1) ≡ 1.
  - On r = 1^(q-2-u) (u+1) (q-u) 1^(u-1), positions u and u+1 lie inside a run of ones:
    - in the leading run if u <= k-1;
    - in the trailing run if u >= k+1.
  - So the sum is 2 ∉ S.
- **D.** On w, [1, 2] carries 1 + (q-3) ≡ -2 ∈ S. On r it carries (q-2) + 1 ≡ -1 ∉ S.

**Exhaustiveness.**
- If 2 ∈ S, we are in case A1 or A2.
- If 2 ∉ S, then χ(2) = 0 and χ(q-2) = 1, so a rise exists in [2, q-3].
- Rises pair up as u <-> q-1-u, and k is the only fixed point of that pairing.
- If k is the only rise, antisymmetry forces χ = 1 0^(k-1) 1^(k-1) 0, that is, S = {1} ∪ [k+1, q-2], which is case D.

**(4).** Direct enumeration for q = 3, 5 (`circulant_sequencings.out`). The rest is (3) together with the cited
reduction.

## Verification

- **`circulant_sequencings.out`.**
  - All S-sequencings for odd q <= 13, one connection set per multiplier class.
  - The criterion (1) against the literal arc test on 1,741,634 paths.
  - (2) for every prime q <= 13.
  - The composite hostile Circ(Z_9, {1,3,4,7}): 66 palindromic paths have non-palindromic words. This is why (1),
    not the word test, is needed at composite q.
  - Counts. H(R7) = 175 with 9 palindromic words of 25. H(QR7) = 189, also with 9 of 27. At q = 13, about 1,120
    palindromic words among about 285,000.
- **`witness_families.out`** (pre-registered).
  - Every one of the 131,071 connection sets with 1 ∈ S, at every odd q from 7 to 35, receives a valid,
    non-palindromic witness by the case rule.
  - Case D occurs exactly once per q.

## Mechanism: why the cut is at 5

- A circulant's Hamiltonian paths are its S-sequencings.
- At q <= 5 there are too few letters for an asymmetric sequencing:
  - q = 3: S = {1}.
  - q = 5: the only sequencings are 1111, 2112 and 2222, all palindromes.
- From q = 7 on, S has at least three letters. One more letter is enough for the two-jump words above, and those
  cannot be palindromes.
- The Collatz side of THM-4849 sees only the connection set (the parity word), not these words.
  - "Short cycles" therefore describes K1, C3, T5 but does not explain them. That is a dictionary, not a mechanism.
  - That 3 and 5 are Fermat primes is NUMEROLOGY: order 17 is equally rigid, and none of its circulants is pure-blue
    (MISTAKE-698).
