# Pre-registration: allowed nu-pairs from run-window words (HYP-9540 (B)) at orders 13-40

opus-2026-11-09-S88. Written and pushed before any (P(A), P(B)) pair census of words at n >= 13 was computed, and
before `pairs_check.py` was written.

## Background (proved this session; THM-5020 draft)

- **Notation.** R(w) is the run-window tiling of a word w of length L = n - 2: the union of the staircases on the
  windows [p, q + 2] of the 1-runs [p, q].
- **(L2).** nu(R(w) + 1) = P(B) for every w.
- **(L1).** nu(R(w)) = P(A) whenever w has no isolated 0 between two 1-runs of length >= 2 (no 0-junction).
- **The admissible words.** So every word w with no 0-junction gives a line with nu-pair (P(A), P(B)). Every word
  with no 1-junction gives one through its complement. Call these words admissible. The line is blue iff w is a
  palindrome.
- **Reduction.** HYP-9540 (B) would follow, for a given n, from a census of admissible words:
  - black: every pair 1 <= i <= j <= k = floor(n/3) with i + j >= k occurs as {P(A), P(B)} of a non-palindromic
    admissible word;
  - blue: every such pair with i + j > k occurs for a palindromic admissible word.

  The tight blue pairs are settled by the blue word law. The "only if" direction of (B) is proved.
- **In-sample data.** The complete line censuses (all lines n <= 9; all blue lines n <= 12) agree with (B).

## Predictions

- **Q1 (black, n = 13..40).** For every n in [13, 40] and every pair 1 <= i <= j <= k with i + j >= k, some
  non-palindromic admissible word of length n - 2 has {P(A), P(B)} = {i, j}.
- **Q2 (blue, n = 13..40).** For every n in [13, 40] and every pair 1 <= i <= j <= k with i + j > k, some palindromic
  admissible word of length n - 2 has {P(A), P(B)} = {i, j}.
- **Q3 (no new blue exclusions).** Consequently, for 13 <= n <= 40 the blue nu-pairs are exactly:
  - the non-tight pairs with i >= 1;
  - the boundary pairs of palindromic tight colourings;
  - (0, k).

A failure of Q1 or Q2 would show only that the run-window family is too small. It would not refute (B).

Outcome: to be appended below without editing the predictions.
