# Boolean-Cube Balance as a Checksum

Session: opus-2026-05-30-S1

THM-351 changes the emotional shape of THM-346. The bucket balance is no
longer primarily a theorem about tournaments; it is a conservation law of the
Boolean cube that tournament quotients inherit.

That matters because the hard-looking line count

```text
2*internal + escaping = bucket_size * masks
```

is really a checksum for reversible binary motion. If every chosen mask is
nonzero, then xor has no fixed points. Internal half-lines therefore come in
two-element pairs, while escaping half-lines remain oriented boundary flux.
The quotient map can be arbitrary. It may be merged tournament class, even
graph class, good-cut bucket, an engineering feature bucket, or something not
yet named. The balance does not care.

This is a useful boundary between structure and bookkeeping:

- Bookkeeping: the row balance is forced by Boolean-cube reversibility.
- Structure: which buckets exchange mass, and how that exchange decomposes
  into spine/ribs/sea or even-graph motion, remains the live mathematical
  problem.
- Engineering: any implementation that computes quotient transport matrices
  can test row integrity against this identity before trusting downstream
  features.

The lesson is not just that Lean closed a formal gap. It identified the level
where the fact naturally lives. Once the theorem is stated at the Boolean-cube
level, the tournament-specific layers become applications rather than burden.

References: THM-346, THM-350, THM-351, HYP-1775, HYP-1778,
`05-knowledge/variables/tiling-bucket-balance.md`.
