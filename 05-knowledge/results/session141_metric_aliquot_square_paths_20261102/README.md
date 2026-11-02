# Session141: metric capacity, divisor fibres, and square-path interfaces

**Checkpoint1 · 2026-11-02.** Elementary proofs, recovered source mechanisms,
finite exact controls, and an independent audit. Status is claim-specific;
no novelty priority, convergence, or new plane-coloring bound is asserted.

## Accepted first checkpoint

- **PROVED / classical mechanism, new scoped application:** the first n for
  which the Erdős unit-distance maximum u(n)>kn is exactly the smallest
  plane unit-distance graph needing more than k functional successor layers.
  The dated AMP baseline gives N1=4, N2=11 and 25<=N3<=28. A single functional
  graph loses exactly m-n+t edges of an arbitrary finite graph, where t is
  its number of tree components; subdivision preserves this loss.
- **PROVED / recovered and self-contained:** faithful plane realizations of
  K3 Cartesian powers have n=3^d, m=d*3^d and chromatic number3. They require
  exactly d functional layers. Density can therefore grow without increasing
  chromatic difficulty. Honeycomb routing, Erdős edge maxima and
  Hadwiger–Nelson coloring are different consumers of the geometry.
- **PROVED all-size extension:** the solutions of sigma(n)=2n+2 with tau(n)<=17
  are exactly20,104,464,650,1952. The inherited proof covered tau<=12; the
  extension excludes the remaining exponent patterns. The next-layer hostile
  is130304, with tau18. This is not a global excess-two classification.
- **RECOVERED / exact source-preserving links:** s(20)=22, while the sum
  excluding1 is21. Twenty and26 share total divisor sum42, not aliquot sum;
  sigma^-1(42)={20,26,41}. The marked1–2–4/F8 relation retains rank and marking,
  and does not provide a global Collatz quotient.

## Read and reproduce

- [Synthesis and honest connection boundaries](SYNTHESIS.md)
- [Checkpoint1 scope](CHECKPOINT1.md) · [concept board](BOARD.md)
- [Unit-distance proofs](unit_distance/PROOF.md)
- [Aliquot proofs and literal source lists](aliquot_twenty/PROOFS.md)
- [Independent first audit](independent/checkpoint1/AUDIT.md)
- [Verification](VERIFICATION.md) · [source status](SOURCE_STATUS.md)

```sh
python3 -B replay_checkpoint1.py
sha256sum -c SHA256SUMS
```

Curvature and square-path refinements are not promoted by this checkpoint.
Collatz convergence, the exact chromatic number of the plane, LRC(14), and
unrestricted aliquot questions remain open in this work.
