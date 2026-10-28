# From arithmetic code constraints to exact corner mass

2026-10-28. Research packet, with elementary proofs, exact finite controls and a
historical recovery atlas. No Collatz, planar Jacobian, square-sum conjecture or
Hadwiger–Nelson closure is claimed. No external novelty claim is made.

**Result:** the seven-prime Hamming carrier has a continuous, source-labelled
extension whose low-order marginals agree across arithmetic syndromes but whose
extreme tails differ exactly. For its zero and a nonzero syndrome, every
three-coordinate marginal is uniform on the same cube; nevertheless

\[
\mathbb E S^4=259/10\quad\hbox{versus}\quad139/10,
\qquad
\Pr(S>13/2)=1/10321920\quad\hbox{versus}\quad0.
\]

The general **rank-to-corner law** computes the entire distribution from a
weighted code enumerator and a box spline. Conditioning on r independent
parity checks amplifies the surviving top corner by exactly 2^r relative to
independent signs. The information discarded by low-order projections becomes
an explicit tail contribution, not merely a warning about a quotient.

- [Proofs and normalization](THEORY.md): the new continuous construction,
  Fourier/moment formulas, exact tail identity, and separator implementation.
- [The exact 56-object correspondence](FIFTYSIX.md): oriented Fano triangles,
  Frobenius conjugates and labelled matching kernels.
- [Recovered 7/56 atlas](ATLAS.md): old proofs, corrected temptations, genuine
  Fano/affine/E7 maps, and the separate Borwein threshold mechanism.
- [Applications and open targets](APPLICATIONS.md): square-sum port decoding,
  fixed-source Collatz, Kummer/planar-Jacobian carriers and six-colour contexts.
- [Inheritance and audit](AUDIT.md): source pins, negative controls, reproducible
  finite universes and limits.
- [Exact checker](check.py), [receipts](receipts.json), [separate verification
  path](audit.py), and [audit receipt](audit.json).

Incoming Session85 supplies the finite nonlinear recovery certificate, the
15/16-prime tables, linear scope-loss and signed incidence. These are credited
as inherited, not new results here. Its parity-plus-coprimality example and
70/364 hostile are rechecked. Session42 already proved square-clock twin
reconstruction and its first ambiguity at36; our duplicate recovery is a
control and not a new theorem. The Borwein thresholds and corner formula are
published results, with the incoming THM-4837 audit fixing attribution.

The atlas inventory scans **15,990 tracked prose paths / 14,914 distinct Git
blobs** at `4e88a9401adb3d1960053263d9254e8582149053`. It locates 12,691 files
mentioning7 and1,684 mentioning56. This is lexical search coverage, not a claim
that every matching file was read or every occurrence is mathematically
meaningful. Exact paths/line numbers are in [occurrences.json](occurrences.json).
