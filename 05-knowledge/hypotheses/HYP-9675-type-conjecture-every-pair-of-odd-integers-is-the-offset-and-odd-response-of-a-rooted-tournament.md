---
id: HYP-9675
title: "Type conjecture for rooted tournaments: for every pair of odd integers (c, r) with c >= 1 there is a rooted tournament whose response is c 2^j for some j and whose offset is r; then c 2^e + r is a Hamiltonian-path count for every e >= j + 1. The column r = 1 is HYP-9658 (A). Quantitative form needed for Chan-Pak: for every odd N other than 7 and 21 a splitting N = c 2^e + r with a root of type (c, r) and exponent below e"
status: >
  OPEN (CONJECTURE). FINITE-EXACT evidence (THM-5145, complete census of the rooted tournaments of order at most
  10): 243,055 of the 261,888 pairs with odd c <= 255 and odd |r| < 2047 are types (92.8 percent); every odd
  c <= 31 with every odd |r| <= 275; for c = 1 every odd r in [-791, 647]. At order 9 the same count is 87,182.
  Known obstructions: none for the existence of a type; the exponent is constrained (offset <= response; types
  whose ray would contain 7 or 21 need a larger exponent, e.g. (1, 3) has least exponent 5).
  The quantitative form is equivalent to the Chan-Pak conjecture (every path count is H(L o 0) for some root).
  Pre-registered predictions on the order-10 census: three confirmed, one failed (the offsets 649 and 825 have no
  cone root of order at most 10).
  AUDIT: AUDIT-PENDING (pre-audit checkpoint)
session: mac-mini-2026-11-11-rays
source: 05-knowledge/results/dyadic_cones_two_block_20261111/README.md
related:
  - THM-5145
  - HYP-9658
  - THM-5126
  - THM-5109
---

# HYP-9675 - the type conjecture

Source: [the note](../../05-knowledge/results/dyadic_cones_two_block_20261111/README.md), section 6.

**Conjecture.** Every pair of odd integers `(c, r)`, `c >= 1`, is the type of a rooted tournament: there are a
tournament `B` and a vertex set `O` with `sum_{y in O} U_B(y) = c 2^j` and `H(B) - sum_{pairs of O} W_ef = r`.

**Falsifier.** A proof that some pair is not a type. **First tests.** The pairs missing at order 10 (for example
`(1, 649)`, `(1, 825)`, `(245, 1)`) at order 11.

**What it would give.** For each pair, all but finitely many points of the progression `c 2^e + r`. It does not
give Chan-Pak without a bound on the exponent.
