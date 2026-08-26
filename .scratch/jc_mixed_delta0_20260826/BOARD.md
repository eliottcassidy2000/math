# Planar Jacobian exact-`M=9` concept board

Status: session scratch; claims below are not canon until independently audited.

## Inheritance

- Closest proved mechanism: `THM-4147-generic-exact-weight-nine-planar-jacobian-monodromy-exclusion.md`, with exact source critical ideals, Newton packets, prime carrier, and response squeeze.
- Canonical hostile example: coefficient specialization can contract a Newton face or a source resultant; `THM-4180` and `THM-4183` show that one must specialize first and replace the lost face explicitly.
- Corrected near miss: `MISTAKE-519`; the repeated-top and `K=0` coefficient faces are already closed, so they are not current frontiers.
- Least-used sidecar: the replacement diagonal coordinate `W=s*p` on `Delta=0`, together with a direct source resultant before projection.

## Live concepts

1. **Anchor — mixed off-antidiagonal `Delta=0`.**  Universe
   `eta*zeta*(eta+zeta)!=0`, forced `K=2848/45`, with `Phi,Theta` arbitrary.
   Test whether the contracted polygon has complete packet
   `(8,8,4,2,2,2,1)` and whether the source critical length stays `25`.
2. **Coefficient-wall atlas.**  Locate every remaining exact-`M=9` wall after
   the proved `Delta=0`, `K=0`, pure, and repeated-top pieces; do not infer
   coverage from theorem titles.
3. **Replacement-face principle.**  Track source, target, map, preserved
   ramification index, and lost coefficient when `Q*p^5*(Delta+eta*s)`
   contracts to a monomial.  Candidate replacement is
   `Q*p^4*(epsilon+eta*s*p)` of index four.
4. **Source critical endpoint.**  Compute `Res_s(A,B)` after `Delta=0`
   specialization.  Units at both ends are required; a generic degree alone
   does not close all `Phi,Theta`.
5. **Entry / exact-`M=10` niche.**  Inventory the first coefficient that the
   exact-`M=9` source universe cannot absorb and seek a typed obstruction,
   rather than assuming the same packet persists.
6. **Wildcard — dual critical projection.**  Compare `(A,B)` with `(A,C0)`
   and the Hessian bridge.  Any disagreement is a chart-loss detector, not a
   new theorem.

## Connection contract for the anchor

- Source: exact-`M=9`, `b=d=0` normalized Keller source at `Delta=0`.
- Target: a complete boundary monodromy packet plus an affine critical-length
  bound sufficient for the planar response contradiction.
- Map: rational source coordinates `(X,T)->(s,p)=(XT,T+X^2T^2)` followed by
  Newton normalization of `(s^2-p)(1-QH)-Qs^2/2=0`.
- Preserved predicate: source smoothness/critical locus on `p(p-s^2)!=0`
  and labelled boundary ramification.
- Destroyed information: the two coordinate fibres and the vanished
  `Delta`-face coefficient.
- Required sidecars: explicit universal-fibre ledger, replacement diagonal
  face, residue-field degrees, and source-first resultant endpoints.
- Cheapest decisive test: exact symbolic resultant plus independent Newton
  hull/face reconstruction at hostile `Phi,Theta` controls.
