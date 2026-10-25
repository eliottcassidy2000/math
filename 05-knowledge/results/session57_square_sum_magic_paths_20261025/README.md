# Session57: square-sum paths and square-magic incidence

25 October 2026. Third checkpoint of an ongoing research session, prepared for independent publication review. Earlier packets are retained. The new affine-frame and prime-plane counting work is reserved for a later checkpoint. No nondegenerate integer magic square of squares, global impossibility theorem, or external novelty claim is made.

## Main conclusion

Square-sum paths do give exact certificates for a magic square when the original nine square values become **edge labels on a numerically decorated additive carrier**. Arbitrary paths on the original nine entry vertices do not give that equivalence. The decoration and nine-value distinctness are essential, not technical details that can be discarded.

Read [SYNTHESIS.md](SYNTHESIS.md) for the hypotheses, proof connections, attribution boundaries, obstructions, and next questions.

## Exact carriers and path counts

- [Original fixed-total carrier](carriers/THEOREM.md), [inherited connections](inheritance/CONNECTION_BOARD.md), and [audit](audit-carriers/REPORT.md): eight synchronized additive paths preserve the eight magic lines. The 22/28 line-intersection density is combinatorial, not an arithmetic criterion.
- [Finite-field separation](character-graphs/THEOREM.md) and [audit](audit-character-graphs/AUDIT.md): F17 and F59 separate original-entry square-sum Hamiltonicity from magic arrangements. A fixed characteristic-three affine magic plane has 432 layouts, not 432 real layouts.
- [Decorated K4,4](biclique-bridge/THEOREM.md), [path counts](biclique-path-counts/README.md), [nine-edge linear frame](two-path-frame/THEOREM.md), and [audit](audit-biclique-paths/AUDIT.md): two paths cover nine distinct formal edge labels; 204 optimal pairs have nine union edges and five shared edges. This does not mean two paths on the original nine vertices. Three paths are needed to cover all sixteen identified K4,4 edges.
- [AP-K3,3 bridge](ap-biclique-bridge/REPORT.md) and [audit](audit-reduced-frames/AUDIT.md): two paths sharing one edge cover the nine identified cross edges, with the AP and diagonal conditions retained. There are 36 optimal pairs, organized into six partner hexagons.

## Integer arithmetic, published families, and geometry

- [Integer bicliques](integer-bicliques/REPORT.md), [known construction and normalization](integer-bicliques/PRIMARY_CONSTRUCTION.md), and [audit](audit-integer-family/AUDIT.md): exact first occurrences are N=274 for K2,4 and N=673 for K3,3. Bare positive integer K4,4 **exists**. The least N lies in **[10001,1112833]**; the upper endpoint is not a global optimum.
- [Lagrange-family obstruction](lagrange-obstruction/proof.md): no nonzero rational member of this particular published family has both required parallelogram balances. The [AP addendum](lagrange-obstruction/ap-addendum.md), with [independent audit](audit-integer-family/AP_AUDIT_ADDENDUM.md), excludes every distinct AP triple in its original A-part. This obstructs selection of original rows/columns, not arbitrary re-embedding or regrouping.
- [Square-root geometry](square-matrix-geometry/geometry.md), [audit](audit-geometry-variant/AUDIT.md), and [known-point addendum](audit-geometry-variant/KNOWN_POINT_ADDENDUM.md): the bare X44 is a six-dimensional normal complete intersection of degree 512 with dualizing sheaf O(2). Its specified decorated pullback has 128 sign components, each the familiar six-quadric magic surface in P8. Bare-point existence does not imply a nondegenerate decorated point.
- [Quartic reduction repair](norm-quartic-repair/repair.md), [direct polynomial verification](norm-quartic-repair/POLYNOMIAL_CHECK.md), and [audit](audit-quartic-elliptic/audit.md): Choudhry's intermediate printed cross equality was a diagonal norm equality. The construction is repaired, and all sixteen square identities of the printed polynomial formula are independently checked, including boundary specializations. This is not a refutation of Theorem 3.

## Other exact interpretations and bridges

- [Unequal square totals](square-total-variant/RESULT.md): nine distinct positive squares have six row/column totals 4001² and diagonal totals 5151² and 2920². This example is **not magic**. Equal line totals that are themselves squares are impossible for a nonzero integer-square magic grid, since T=3e².
- [Congruent-curve bridge](congruent-curve/PROOF.md) and [audit](audit-quartic-elliptic/audit.md): a guarded AP of three centers in one fixed-gap square set is equivalent to the nine-square condition, and to an AP in eligible doubled elliptic x-coordinates. Common scaling multiplies area and gaps by L², not by 1. The integer-gap exclusion d≤20000 is bounded, not global.
- [Exact SVG diagrams](diagrams/README.md): [edge-label repair](diagrams/edge-label-repair.svg) and [square totals are not magic](diagrams/square-totals-not-magic.svg).

## Verification and scope

The packets pair proof arguments with exact finite or symbolic certificates. Some producer scripts use Python assertions and must be validated in ordinary mode; running them with -O is not effective verification. Their independent exception-safe audits retain checks in both modes. The quartic/elliptic transformed scripts are robustness checks, supplemented by independent symbolic derivations and source comparison, not independent implementations of every producer calculation.

Source attribution and version/hash receipts are retained. Full downloaded articles, preview PNGs, compiled search executables, caches, and rendering logs are excluded. Historical manifests can mention such excluded audit inputs; that is provenance, not publication of those inputs. Earlier Session55 modular results retain their stated scope; finite fields are not integer lifts, and extension-field encodings are not arithmetic modulo q.
