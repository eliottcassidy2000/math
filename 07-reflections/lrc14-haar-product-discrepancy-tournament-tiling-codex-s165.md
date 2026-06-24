# LRC14 Haar Product Discrepancy And Tournament Tiling

The new synthesis is that the 2D Haar product rule is the local square behind
several recent proof languages.

On a dyadic rectangle, `h_I(x)h_J(y)` gives the checkerboard

```text
[[+,-],
 [-,+]].
```

That is also the elementary 2-by-2 fixed-margin switch.  It is the move that
takes an anti-diagonal packet to a diagonal packet without changing row or
column margins.  So margins are exactly the quotient that forgets the mixed
orientation sign.

This makes the tournament tiling analogy much less decorative.  A fixed
Hamiltonian path is an observer axis.  The tiling square formed by two observer
cuts has a mixed sign, and that sign is the Haar product coefficient.  If the
proof remembers only row/column shadows, diagonal and anti-diagonal packets
are indistinguishable even though their mixed discrepancy is opposite.

That explains why HYP-2594's component count `K` is too pessimistic.  It counts
micro-boundaries before the mixed-product cancellations.  HYP-2595's
color-compatible resonance condition is what survives after the product signs
are respected.

It also clarifies the newest LRC14 packet work:

- HYP-2986 tope/cocircuit walls are the cyclic endpoint version of the same
  switch.
- HYP-2986 handoff arrows must preserve the mixed square or name its loss.
- HYP-2985 smoothing policies are admissible only when they do not smooth out
  the product mode without a boundary label.
- HYP-2981 Fejer certificates are higher-order mixed-sign detectors.
- HYP-2978 quotient guardrails are already telling us that row/column
  summaries are not enough.

The theorem-facing hope is a Haar-product discrepancy lemma for the colored CRT
placement layer:

```text
V*Sigma - actual_count is controlled by the independent color-compatible
mixed Haar switches, not by the raw number of continuous components.
```

This would give a conceptual route to the HYP-2595 bound
`Delta <= C*(k+c_GP)`: `k+c_GP` should count switch families after packet
compression, while `K` counts all visible boundary fragments before quotient
repair.

The strongest warning from the toy model is almost comically small:

```text
[[1,0],[0,1]] and [[0,1],[1,0]]
```

have the same row and column margins, but Haar coefficients `+2` and `-2`.
Any LRC14 quotient that cannot see that difference is not yet theorem-safe.
