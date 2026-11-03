# S144 correction lineage

The live THM-4727 patch changes only the tangent-sweep completeness claim, the false higher-scale existence question, the stronger tangent-weight gate, and invariant-hyperplane terminology. Its original source is pinned below. The general fiber-degree and W=1 ramification limitations are made explicit in the correction banner, without pretending the entire legacy theorem was reaudited.

- Original commit: `5e1cedf452c7cc1ec447205939a9d6cf4c0e8743`.
- Original path: `01-canon/theorems/THM-4727-mixed-weight-torus-reduction-of-keller-maps-ramification-exponent-twisted-planar-jacobian-equation-and-the-dimension-profile.md`.
- Original SHA256: `eebf30e950564595a72f48b62304184930f7a48b3b6503473daebef96c9288df`.
- First failed inference: affine orbit scale in the seed construction was treated as exhaustive after closure under composition.
- Decisive witness: `G^2` has determinant 4, geometric degree 9 and orbit-scale degree 7.
- Strongest survivor: the displayed tangent-sweep seed family has affine scale and the quotient determinant identity holds.
- Repair: retain composition data and the scale divisor; classify primitive solutions only after fixing the equivalence relation.

The complete original bytes remain in Git history. No reserved theorem number or empty claim was introduced.
