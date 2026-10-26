# Session71: Dudeney, Gaussian norms, and native Penrose alignment guards

2026-10-26. **PROVED** elementary interfaces below, **CITED** external dissection
optimality, and **FINITE-EXACT** controls are kept separate. No new proof of
Dudeney's four-piece minimum or broad equivalence of open problems is claimed.

## Main answer

Fermat's two-square theorem relates to the problem through **norm and coordinate
compatibility**, not as the reason three pieces are impossible. The strongest
structural connections found here are exact maps that retain the information
subsequent operations can observe:

1. **Gaussian squaring ↔ a canonical Pythagorean/Farey graph.** The inherited
   map preserves determinant/Lorentz adjacency; norms alone do not.
2. **Whole-side matching diagrams ↔ two named involutions.** Both involutions,
   node types and endpoint order are needed. Their product alone loses even
   connectedness.
3. **Marked Penrose placements ↔ the Boolean conjunction test system.** Actual
   nonoverlapping kite/dart moves realize every hidden-bit conjunction, giving
   the same exact correlation obstruction as Session70's tournament family.
4. **Length gluing ↔ guarded partial-affine maps.** Sign, ordered offset and
   admissible interval form a closed exact carrier; cancellation alone loses
   legality.

[Full cross-domain synthesis and proofs map](session71-root/BRIDGES.md).

![The exact marked Penrose attachment control](session71-root/penrose-guards.svg)

## Arithmetic: a precise bridge and a decisive boundary

The two-square criterion separates5=1²+2² from44=4·11. Yet triangle/square
pairs scaled to either area have the same four-piece optimum in the cited
category, while square-to-itself uses one piece. Area-norm membership cannot
therefore determine the minimum.

For the side-two triangle, the square side isθ=3^(1/4). Recovering the earlier
PS-D2 rigid-ear theorem gives an orientation-free edge-length obstruction:
each square side needs an exposed piece-edge length outside Q(√3). That is
four exposed edge occurrences, **not necessarily four different pieces**.

A sharper target-coordinate theorem distinguishes two conventions:

- Real Cartesian x/y coordinates need a number field of degree at least4;
  Q(θ) attains it. No totally real field works.
- Complex vertex coordinates need degree at least8; Q(θ,i) attains it with
  Galois group D4. No abelian complex number field works.

The four roots±θ,±iθ carry the exact square-dihedral permutation action,
but field automorphisms are not rigid motions of the entire plane. These
statements concern target vertices, not the field of Dudeney's cut coordinates.

There is a golden counterpart: equal-area similar Penrose kite/dart copies
also require complex degree at least8. The explicit field Q(ζ5,√φ) attains
it but is **not Galois**, unlike the triangle/square example. Native Penrose
tiles are not obstructed: they have area ratioφ, rather than equal areas.

- [Norms, Pythagorean maps, prior44/5 recovery](session71-norms/REPORT.md)
- [Independent field and labelled-rank proofs](session71-norm-audit/REPORT.md)
- [Golden norm review](session71-norm-audit/PENROSE_P5_REVIEW.md)

## Geometry: what gluing must remember

The Penrose construction uses disjoint finite marked placements and the same
named translation for either hidden kite orientation. Subset moves have guard
value∏b_j. Singleton guards recover k bits; the signed-law response space has
dimension2^k. Proper-subset marginals of parity mixtures agree, but a global
move distinguishes them. This is not a theorem about infinite tiling extension.

The dissection length carrier is(ε,c,I), representing x↦εx+c on its exact
admissible interval. A balanced twelve-step chain with totalsθ and2 has
identity formal composition but empty domain. The analogous angle-charge
carriers are Z6 for triangle/square and Z10 for pentagon/square: similar
corner-count consequences, **not isomorphic charge groups**.

Finally, flat Euclidean seam transport develops consistent frames but does
not ensure nonoverlap. A five-square path returns its last square onto its
first. This is the geometric counterpart of keeping Petersen's local run
options while still needing a separate global connectivity test.

- [Matching, Penrose, and Euclidean carrier proofs](session71-penrose/REPORT.md)
- [Partial-affine guards and angle charges](session71-dudeney/REPORT.md)
- [Independent matching and geometric audit](session71-matching-audit/AUDIT.md)

## Forty-four, five, and source versions

The displayed paper inventory has44 square cut types and5 triangle cut types,
with46 listed pairs surviving its filters. This is distinct from !5=44,
the earlier pure-black tournament count44, and rooted-complement target counts
5 and44. No structure-preserving bijection between those entire sets is proved.

The requested [arXiv v3](https://arxiv.org/html/2412.03865v3) has an official
[v4 correction notice](https://arxiv.org/abs/2412.03865v4). The current primary
publication is [ITCS2026, DOI10.4230/LIPIcs.ITCS.2026.47](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.ITCS.2026.47).
Our ledger records exact served-version observations and figure-label cautions;
it does not claim a full independent certification of the paper's case proof.
Its theorem covers polygonal pieces without flips, not the flipped/curved
extensions discussed there.

[Incoming repository corrections](INCOMING_RECONCILIATION.md) preserve the
updated HYP-9416 scopes without rewriting frozen historical source bytes.

## Reproduction and remaining frontier

Run `python3 -B verify_packet.py --replay` in this packet. The envelope checks
curated files, source provenance, producer and independent normal/optimized
receipts, and explicit hostile controls. Dependencies and exact finite universes
are recorded in the lane reports. The drawing is an illustration of separately
verified exact coordinates, not a floating-point proof.

The remaining challenge is an economical geometric carrier that retains
realization, positivity, marks, permitted motions and global nonoverlap without
merely storing every raw assembly. The exact subcarriers here provide constraints
and reusable components; they are not a universal dissection classifier.
