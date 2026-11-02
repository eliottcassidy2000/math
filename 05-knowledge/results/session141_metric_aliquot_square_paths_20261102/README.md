# Session141: metric capacity, divisor fibres, and square-path interfaces

**Checkpoint2 · 2026-11-02.** Elementary proofs, recovered source mechanisms,
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

## Accepted second checkpoint

- **PROVED:** triangle inequality gives an upper ambient radius, not a lower
  distance bound. In curvature -kappa^2, separated bounded-multiplicity binary
  branching requires kappa*L>=log2. A proper noncrossing quasi-isometric
  embedding with equal geodesic edge length exists for every kappa*L>log2.
  The threshold is a sharp infimum; equality is unresolved. This construction
  is for the abstract binary tree, not the full native Collatz graph.
- **PROVED / refinement:** every nonempty unit-distance extremizer is
  connected. For n>=3 its exact functional deletion defect is u(n)-n and
  cycle rank u(n)-n+1. The earlier general formula still applies to
  nonextremal disconnected point sets.
- **PROVED + FINITE-EXACT:** common directed blocks give exact joint path-union
  sizes; conflicts occur precisely at reversed joins between singleton blocks.
  The recovered n25 square-guard polynomial extends at n26 to
  `2^296(536870561+320z+29z^2+2z^3)`:24 directed paths,12 compatible pairs,
  two triples. The n25 formula and general subset expansion are inherited,
  not new. At n24 there is no square-sum Hamilton path; finite spectra through
  n30 are retained without claiming an all-order law.

- **PROVED cross-lane bridge:** a compatible path union with b common blocks
  and J joins has cycle rank |J|-b+1 and reorientable successor loss
  max(0,|J|-b). Fixed arrows need m-n+z deletions, retaining the sink count z:
  the starting18 pair/triple cost1/3, their reversed ending18 families cost2/4.
- **PROVED normalization:** coefficient content equals the gcd of binomial
  moments. Reversal and exact joint multiplicities give contents2^276 and
  2^296 at25/26; the exponent jump20 has a precise free-orientation account,
  not an aliquot identity.

## Read and reproduce

- [Synthesis and honest connection boundaries](SYNTHESIS.md)
- [Literal common-block / successor-defect bridge](CONNECTIONS.md)
- [Checkpoint1 scope](CHECKPOINT1.md) · [Checkpoint2 scope](CHECKPOINT2.md) · [concept board](BOARD.md)
- [Unit-distance proofs](unit_distance/PROOF.md)
- [Aliquot proofs and literal source lists](aliquot_twenty/PROOFS.md)
- [Independent first audit](independent/checkpoint1/AUDIT.md)
- [Curvature and metric proofs](curvature_triangle/PROOF.md)
- [Square-path common-block proofs](square_paths/PROOF.md)
- [Connected-extremizer refinement](unit_distance_refinement/PROOF.md)
- [Binary content, reversal and exact normalization](square_paths_content/PROOF.md)
- [Independent refinement audit](independent/refinements/AUDIT.md)
- [Second independent audit](independent/checkpoint2/AUDIT.md)
- [Verification](VERIFICATION.md) · [source status](SOURCE_STATUS.md)

```sh
python3 -B replay_checkpoint2.py
sha256sum -c SHA256SUMS
```

The first-checkpoint scientific evidence remains byte-frozen.
Collatz convergence, the exact chromatic number of the plane, LRC(14), and
unrestricted aliquot questions remain open in this work.
