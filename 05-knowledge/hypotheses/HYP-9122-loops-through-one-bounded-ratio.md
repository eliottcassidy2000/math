---
id: HYP-9122
title: "Loops through 1 of every length with ratio below 3 (robust 1-escape)"
status: >
  OPEN HYPOTHESIS with FINITE-EXACT support to s=40. For every s>=1 there is
  a forward E-cycle through 1 with s multiplications and K halvings,
  2^K<3^(s+1). Together with the PROVED 1-escape lemma this handles every m
  with v_3(m-1)>=2 in Q2. For s>=2 every loop has ratio >= 13/9, and the
  observed minimal ratios lie in [1.517, 2.96], thin at s=11 and s=23. The
  robust form allows other exit words and is what the E-SCC induction
  actually needs.
source: collatz-procgen-20260922 (mac-mini)
depends_on:
  - 05-knowledge/results/collatz_procgen_20260922_choice_ladder.md
---

# HYP-9122 -- bounded-ratio loops through 1

Refutation form: an s for which every loop has 2^K >= 3^(s+1), together
with no alternative exit of bounded ratio.
