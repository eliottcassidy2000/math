---
id: HYP-9590
title: "Run-window gap law: (A) for every n >= 13, the nu-pairs realised by the run-window lines of admissible words (THM-5020 (6)) are exactly the pairs allowed by HYP-9540 (B), minus three families - black (k, k) at n = 3k, blue (k, k) at n = 0, 1 (mod 3), and blue non-tight pairs with both entries odd at even n; (B) every pair in these three families is nevertheless realised by some line of the tiling metagraph at every n >= 9 (so HYP-9540 (B) holds); (C) the palindromic words never realise (k, k) at n = 0, 1 (mod 3)."
status: >
  OPEN.
  (A): FINITE-EXACT for 13 <= n <= 40 (exact census of admissible words, pairs_check.out). The "minus" part is
  PROVED for two of the three families:
  - black (k, k) at n = 3k is impossible for any word (unique maximum packing);
  - the both-odd blue pairs at even n are impossible for palindromes (THM-5020 (5), central-pair lemma).
  The "exactly" part (every other allowed pair is realised) is the conjecture.
  (B): FINITE-EXACT where complete line censuses exist: black (3, 3) at n = 9; blue (3, 3) and (4, 4) at n = 12. OPEN
  beyond. This is the residue of HYP-9540 (B) after THM-5020.
  (C): FINITE-EXACT for n <= 40 (it is part of the census). OPEN in general.
session: opus-2026-11-09-S88
source: 05-knowledge/results/word_law_markov_20261109/README.md
related:
  - THM-5020
  - THM-5021
  - HYP-9540
  - THM-4951
---

# HYP-9590 - run-window gap law and the residual pairs

**Setting.** As in THM-5020. A word is admissible if it avoids 11011 or avoids 00100. An admissible word w gives a
line with nu-pair {P(A), P(B)}, blue iff w is a palindrome. k = floor(n/3).

## (A) Gap law

**Statement.** For n >= 13 the unordered pairs {P(A), P(B)} of admissible words cover the pairs allowed by
HYP-9540 (B) except exactly:
- black (k, k) at n = 3k;
- blue (k, k) at n = 0, 1 (mod 3);
- blue non-tight pairs (i + j > k) with i and j both odd, at even n.

**Proved parts.**
- At n = 3k, P(A) = k forces A to contain the unique maximum packing {0, 3, ..., 3k - 3} of the positions, so
  P(B) <= k - 1.
- At even n, a palindrome never has P(A) and P(B) both odd (THM-5020 (5)).

**Evidence.** The census for 13 <= n <= 40 (PREREGISTRATION.md, Q1-Q2 outcome) misses exactly these families and
realises every other allowed pair.

**Falsifier.** An allowed pair outside the three families that no admissible word of some length n - 2 realises, or
a pair in the families that one does.

## (B) The residual pairs exist

**Statement.** For every n >= 9 each pair in the three families is the nu-pair of some line, black or blue as
stated. Such lines are then not run-window lines of admissible words.

**Evidence.** Black (3, 3) at n = 9 (121,953,098 lines) and blue (3, 3), (4, 4) at n = 12 (22,586 and 534,399,389
lines).

**Why it should hold.** The excepted pairs sit near the top (nu close to k on both sides), where most lines live.
At n = 9 the pair (3, 3) holds 91% of all lines.

**Falsifier.** An excepted pair with no line at some n >= 13.

## (C) Palindromic (k, k)

**Statement.** At n = 0, 1 (mod 3) no palindrome of length n - 2 has P(A) = P(B) = k.

**Likely proof route.** The maximum k-packings of [0, L - 1] at these lengths have at most two phases, and they are
exchanged by the mirror. A symmetric colour class containing one contains the other, and then leaves the second
colour too few positions.

**Falsifier.** A palindrome with P(A) = P(B) = k at such n.
