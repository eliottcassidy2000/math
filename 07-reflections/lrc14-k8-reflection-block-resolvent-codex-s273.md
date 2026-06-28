# LRC14 k=8 Reflection-Block Resolvent

Source: codex-2026-06-27-S273  
Anchors: HYP-3139, HYP-3132, HYP-3122, HYP-3085, HYP-3133, T1204, LTI-265,
LTT-163, OPEN-Q-108

The useful surprise in this pass is that the HYP-3132 solvable quartic is not
just an analogy for the hard `k=8` row.  It is visible as a concrete matrix
subproblem.  The pairwise sector co-emptiness matrix for `consec_8`, restricted
to sectors `1..5`, is exactly symmetric under `s->6-s`; the reflected pairs
`(1,5)` and `(2,4)` are the inner 2x2 shell page, sector `3` is the center
coupling, and sector `6` is boundary leakage.

That changes the proof-search shape.  The older question "can the quartic be
handled?" becomes the smaller question "can the inner reflected 2x2 shell be
bounded, then corrected by center/boundary ceilings and the known phi4 sign?"
This is the kind of quotient the repo has been asking for: it says exactly
which coordinates survive and which coordinates are being deliberately kept as
sidecars.

Niche connections pulled in:

- HYP-3085 supplied the exact pairwise `S2` covariance route and the warning
  that the core is reflection-symmetric, not circulant.
- HYP-3122 supplied the negative `kappa4` / `phi4` sign at `k=8`, turning the
  higher-order correction into a directional stabilizer.
- HYP-3132 supplied the De Moivre biquadratic vocabulary and the discriminant
  `9` rational collapse.
- HYP-3133 supplied the controlled-forgetting warning: A000568 is a legal
  middle shadow only after endpoint/child/boundary payloads are named.
- The older Worpitzky/Eulerian path-moment lane remains useful as coefficient
  vocabulary, but this session demotes it from load-bearing proof engine to
  expansion language.

Assumption challenge: the tournament vertices are not runners, arcs, or raw
matrix entries.  The vertices are proof pages: exact core block, inner
biquadratic fold, SPEC certificate, A000568 middle shadow, Worpitzky vocabulary,
antisymmetric nonmax block, boundary leakage, and raw spectrum.  The quotient
preserves the `k=8` bounded-core proof obligation only when the center and
boundary pages are retained; it destroys the actual LRC predicate if only the
inner shell is kept.

Candidate labelled packet theorem:

```text
For the bounded-core k=8 row, the pinned-sector co-emptiness packet decomposes
under s->6-s into inner-shell, center-coupling, antisymmetric, and boundary
pages.  If the inner-shell page satisfies the HYP-3132 biquadratic ceiling, the
center and boundary pages satisfy their rational leakage ceilings, and the
HYP-3122 phi4 sign supplies the S3/S4 correction, then the k=8 dip bound follows.
```

The next proof move should be exact and finite: derive the inner 2x2 shell
bound as a rational inequality, then show the center and sector-6 leakage terms
cannot restore the missing mass.  If that can be written in the HYP-3107
coverage-extremality interface, the LRC14 covering route has a much smaller
remaining surface.
