---
id: HYP-9134
title: "The smallest Mahler-bridge instance: no xi > 0 has all xi (2187/1024)^j in Z + [4726, 4727]/1163 (equivalently, no positive integer's 3x+1 parity vector is eventually made of the blocks 0111101110 and 1101100111)"
status: >
  OPEN HYPOTHESIS. By THM-4469 (Theorem M') the two forms are equivalent.
  Collatz, no-divergence (T1) and the Periodicity Conjecture each imply it.
  It is beyond every proved Mahler-type exclusion: those reach interval
  length 1/p (p = 2187), and this interval has length 1/1163 = 1.88/p.
  Heuristically no xi exists, since alpha |I| = 2187/(1024*1163) < 0.002.
  The block class is supercritical (density 7/10 > log_3 2), has positive
  entropy (0.1 bit per step) and bounded discrepancy. Evidence (FINITE-EXACT):
  of the starts in [1, 2^22], 16 begin with two such blocks and none with
  three.
source: collatz-procgen-20260922 session, implication-atlas lane, 2026-09-24
depends_on:
  - 01-canon/theorems/THM-4469-mahler-bridge-adjacent-block-pairs.md
related:
  - 05-knowledge/hypotheses/HYP-9127-cube-swap-cubic-theta.md
  - 05-knowledge/hypotheses/HYP-9123-supercritical-strips-periodicity.md
---

# HYP-9134 -- the smallest Mahler-bridge instance

**Statement.** Let `alpha = 3^7/2^10 = 2187/1024` and `I = [4726, 4727]/1163`,
where `1163 = 3^7 - 2^10`. No real `xi > 0` has `xi alpha^j in Z + I` for all
`j >= 0`.

**Equivalent form (THM-4469).** No positive integer has a `3x+1` parity vector
that is eventually a concatenation of `B = 0111101110` (`R_B = 4726`) and
`B' = 1101100111` (`R_B' = 4727`).

**Refutation form.** An explicit `xi`, or equivalently a positive integer whose
orbit follows `B`/`B'` blocks forever. Such an orbit diverges, growing by
`alpha = 2.136` per block.

**Why it matters.**
* It is the smallest instance where Collatz no-divergence and a generalized
  Mahler Z-number problem are the *same* statement.
* A Mahler-side proof would settle a HARD slice of T1: supercritical, positive
  entropy and bounded discrepancy. The Collatz-side methods (capacity,
  repetition, Theorem D, the 2-adic Tschakaloff–Hankel arguments of
  Theorems Y and H1) do not apply, because the class is uncountable. It is
  the positive-entropy twin of the zero-entropy cube-swap number HYP-9127.
* Proved Mahler-type exclusions (Flatto–Lagarias–Pollington, Dubickas,
  Bugeaud) stop at length `1/p`. Beating `1/p` by any constant factor for
  one ratio `3^a/2^L` with this interval would already be new.

**Suggested attacks.**
1. A two-digit Pisot-type or discrepancy argument for `{xi alpha^j}` confined
   to a single interval of length `c/p` with `1 < c < 2`, for `alpha = 2187/1024`.
2. The 2-adic side: the parity vectors in `{B, B'}^N` form a Cantor set of
   dimension `1/10` in `Z_2`. One could show that its intersection with
   `Z_(>0)`, reached through the conjugacy, is empty, for example by a
   transversality count at scale `2^(10 j)` against the growth `alpha^j`.
3. Larger instances (`(16, 11)` with factor 1.59, `(20, 14)` with 1.28) are
   closer to `1/p` and may be easier on the Mahler side.
