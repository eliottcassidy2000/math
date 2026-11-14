# Exact frame decoding and the limits of natural witness selection

**Chan–Pak remains open.** This checkpoint builds two exact inverses,
distinguishes two growth exponents, and pinpoints a failed induction
rather than extending a numerical coverage census.

- [Terrace target inverse](RESETS_AND_DECODER.md). Given an exact count,
  decode the least-order normalized gap word or reject the entire first
  triangle terrace. The only rank-one corridors have lengths 0 and 13:
  their determinant is a negative-index Tribonacci term. Each creates
  a genuine factorization, so a prime path count cannot use either.
  Away from those resets, fixing the interior frame leaves at most two
  possible sums of the outer exponents. The same kernels give an exact
  conditional independence theorem and a uniform memory contraction.
- [Two-hub frame inverse](FRAME_DECODER.md). Compile an arbitrary fixed
  homogeneous-block skeleton into nonnegative surjection coefficients.
  Recover active widths from the whole count and, optionally, a given
  marked three-run overlap. Returns a canonical isomorphism class from a
  chosen bank, with inert widths normalized. It does not supply a
  universally successful bank.
- [Canonical selection](CANONICAL_SELECTION.md) and
  [exact append residuals](RESIDUAL_SELECTION.md). Count 63 alone rules
  out one-at-a-time descent through minimum-order scalar representatives.
  The two-hub append quotient has linear Hankel rank 17, allowing a
  compatible shortest-prefix choice for every actual residual. A
  successful global induction still must prove a physically available
  return for every unresolved prime and seven-times-prime target.
- [Pressure, and a lower bound using 2 and 3](PRESSURE.md). The old
  0.98 was a sufficient convergence parameter, not an observed or exact
  exponent. The new bounds for distinct terrace outputs are
  Omega(X^(11/30)) and O(X^0.952). The threshold for ordered corridor
  words is in (0.950, 0.952); it need not be the threshold for distinct
  integers. The lower bound decodes a mod-3 guard and a 2-adic width
  from the same physical output.
- [Long reversals and color channels](WRAP_CHANNELS.md). One long wrap
  around an arbitrary ordinal sum of terraces still gives only
  O(X^0.93) distinct values. A greedy scan of disjoint triangle windows
  gives the strict estimate. In contrast, iterated wraps require an
  unbounded ladder of color-sum coordinates in any linear interface.
  This provides an escape from the two-coordinate proof, not coverage.

Reproduce with `python replay_decoder.py`, `python replay_frames.py`,
`python pressure_replay.py`, and `python wrap_replay.py`.
Try `python decoder.py 385 573 25281` for the terrace inverse, or
`python frame_decoder.py 385 '1 0* 3 0 2* 2 1' --anchors 0 2 5 6 --marked 36`
for one marked two-hub frame.
