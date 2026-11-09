---
id: HYP-9590
title: "Run-window gap law: (A) for every n >= 6, the nu-pairs realised by run-window lines of admissible words (THM-5020 (6)) are exactly the pairs allowed by HYP-9540 (B), minus three families - black (k, k) at n = 3k, blue (k, k) at n = 0, 1 (mod 3), and blue non-tight pairs with both entries odd at even n (the exclusion is PROVED; the realisation of all other pairs is the conjecture); (B) every pair of these families is realised by some line at every n (PROVED for (k, k) black and blue and for the blue pairs (i, k); open for the blue both-odd pairs (i, j), j < k, beyond n = 26); (C) PROVED: no palindrome realises (k, k) at n = 0, 1 (mod 3)."
status: >
  OPEN, in two places: (A), the realisation half; (B), for the both-odd blue pairs (i, j) with j < k at even
  n > 26.
  (A):
  - The exclusion of the three families is PROVED (THM-5020 (6)).
  - The realisation of every other allowed pair is FINITE-EXACT for 6 <= n <= 40 (exact census of admissible words,
    reproduced by the audit).
  (B):
  - PROVED, by explicit lines (THM-5020 (6), the audit's constructions):
    - family (i), black (k, k) at n = 3k >= 6;
    - family (ii), blue (k, k), at every n >= 6;
    - the members (i, k) of family (iii);
    - also blue (1, k) at every n >= 4.
  - The remaining members of (iii), both-odd (i, j) with j < k at even n, are FINITE-EXACT for n <= 26 (the audit's
    structured family). OPEN beyond.
  (C): PROVED (THM-5020 (6); the k-packings S_j and the mirror). It is also checked on every palindrome of length
  <= 32 (audit).
  Audited with THM-5020: PASS WITH CORRECTIONS (m5, m6, m15 applied).
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

**Statement.** For every n >= 6 the unordered pairs {P(A), P(B)} of admissible words cover the pairs allowed by
HYP-9540 (B) except exactly:
- (i) black (k, k) at n = 3k;
- (ii) blue (k, k) at n = 0, 1 (mod 3);
- (iii) blue non-tight pairs (i + j > k) with i and j both odd, at even n.

**Proved: no admissible word realises the three families.**
- (i): the unique maximum packing at n = 3k.
- (ii): the k-packings S_j at n = 3k + 1 and the mirror exchange of S_0 and S_k.
- (iii): the central-pair lemma, THM-5020 (5).

For n <= 28 these pairs are missing for all words, admissible or not (audit, pairsx_*.out).

**Conjectural: every other allowed pair is realised.**
- FINITE-EXACT for 6 <= n <= 40.
- **Likely proof route.** Prefix a tight block word X with (100)^r. P is additive across the seam 00|11. This
  realises (r + P_A(X), r + P_B(X)) and avoids 11011. The top pairs need modified endings, such as (100)^(k-1) 10
  at n = 3k + 1.
- Appending 111 preserves "avoids 00100", and appending 000 preserves "avoids 11011". This suggests an induction on
  n.

**Falsifier.** An allowed pair outside the three families that no admissible word of length n - 2 realises.

## (B) The residual pairs exist

**Statement.** For every n >= 6 each pair in the three families is the nu-pair of some line, black or blue as
stated. These lines are not run-window lines of admissible words.

**Proved.** These are the audit's constructions (THM-5020 (6)).
- **Nested mirror lines.** Reverse exactly the mirror tiles (a, n - 1 - a), a < m. This gives blue (m, k) for every
  0 <= m <= k at every n >= 6. That covers family (ii) and the members (i, k) of family (iii).
- **Black (k, k) at n = 3k >= 6.** The nested mirror line with m = k, plus the tile (0, 3). That covers family (i).

**Open.** The both-odd blue pairs (i, j) with j < k at even n. They are realised for every n <= 26 by R(w) of a
palindrome with at most three mirror tiles flipped.

**Falsifier.** A both-odd pair (i, j), j < k, at some even n > 26 with no blue line.

## (C) Palindromic (k, k): PROVED

At n = 3k every word is excluded (the unique maximum packing). At n = 3k + 1:
- the k-packings of [0, 3k - 2] are S_j = {0, 3, ..., 3(j - 1)} u {3j + 1, ..., 3k - 2};
- the only disjoint pair is {S_0, S_k};
- the mirror x -> 3k - 2 - x exchanges S_0 and S_k, so a palindrome would need w_0 = w_(3k-2) in both colours.

The checks are `palindrome_kk_check.py` (the packing list for k <= 7; every palindrome to n = 31) and the audit's
`central.out` (to length 32).
