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

## Outcome (appended after the run; the predictions above are unedited)

`pairs_check.py` ran at n = 13..40 (pairs_check.out). An in-sample calibration at n = 6..12 was run alongside
(pairs_check_insample.out), regrettably after the predictions were pushed (MISTAKE-745).

- **Q1 REFUTED, in exactly one family.** For every n in [13, 40], every black pair 1 <= i <= j <= k with i + j >= k
  is realised by a non-palindromic admissible word, except (k, k) at n = 3k (n = 15, 18, ..., 39).
  - (k, k) is impossible for any word at n = 3k. P(A) = k forces A to contain the unique maximum packing
    {0, 3, ..., 3k - 3} of the positions [0, 3k - 3], so P(B) <= k - 1.
  - Such lines exist nevertheless: at n = 9 there are 121,953,098 black (3, 3) lines. The run-window family is too
    small there; it is not an obstruction.
- **Q2 REFUTED, with a parity pattern.**
  - At even n (L = n - 2 even), no palindromic word, admissible or not, has P(A) and P(B) both odd. This is a theorem
    (THM-5020 (5), the central-pair lemma); it was checked for every palindrome of even length <= 24. So the
    non-tight both-odd pairs are missed: (3, 3) at n = 14; (1, 5), (3, 3), (3, 5), (5, 5) at n = 16; ...
  - (k, k) is missed at every n = 0, 1 (mod 3).
  - Blue lines with these pairs do exist at small n (blue (3, 3) at n = 12: 22,586 lines in the audit's complete
    census), so again the family is too small.
  - All other non-tight blue pairs are realised for n <= 40.
- **Q3:** not established by this family.

Consequences:
- HYP-9540 (B) is now PROVED for 13 <= n <= 40 on every pair except:
  - black (k, k) at n = 3k;
  - blue (k, k) at n = 0, 1 (mod 3);
  - blue non-tight pairs with i and j both odd at even n.
- Each such line exists by THM-5020 (L1), (L2) applied to a word found by an exact census.
- The excepted pairs remain conjectural for n >= 13.

## Label note (appended after the audit; nothing above is edited)

"PROVED for 13 <= n <= 40" in the Consequences above means a FINITE-EXACT census of admissible words, combined with
the proved lemmas THM-5020 (2), (3). This is how THM-5020 (6) labels it (audit m14). The census also holds on the
in-sample range, so the gap law holds for every n = 6..40.

The audit then supplied explicit lines for the three families: nested mirror tiles; (0, 3) for black (k, k) at n = 3k.
These are proved for (k, k), black and blue, and for the blue pairs (i, k) at every n. They are FINITE-EXACT for the
remaining blue both-odd pairs up to n = 26. With them, HYP-9540 (B) holds for every n <= 26 (THM-5020 (6), HYP-9590).
