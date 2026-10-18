# Graphs, Pythagorean triples, and the isosceles boundary: research board

**Session:** codex-graph-pythagorean-20261018. **Status:** research routing;
individual claims retain their statuses in the linked proofs and receipts.
Started from fetched `ceb77fcd33a0` in an isolated worktree. The original
checkout's untracked files were preserved. No literature-priority claim.

## Inheritance before derivation

- Closest proved arithmetic mechanism: [THM-3333](../../01-canon/theorems/THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md),
  the Gaussian-square lift and retained parity/scale; [THM-3334](../../01-canon/theorems/THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md),
  the difference-one parabolic spine.
- Closest graph mechanism: [THM-261](../../01-canon/theorems/THM-261-petersen-root-orthogonality.md),
  Petersen as disjointness on the edges of K5; [THM-4116](../../01-canon/theorems/THM-4116-boundary-state-gluing-and-ap-odd-shell-tree-synchronizers.md),
  exact ordered boundary-colouring vectors and their dot-product gluing law.
- Canonical hostile: K4 is three-edge-colourable but has no nowhere-zero
  Z/3Z flow. All cubic graphs have edge/vertex ratio 3/2, including K3,3,
  the cube and Petersen; density cannot detect snarks.
- Corrected near miss: MISTAKE-507 separates pair-indexed roots from Lie
  algebra isomorphisms, Z/3Z flows from F2^2 flows, and arithmetic analogies
  from maps preserving an obstruction. MISTAKE-558 warns against crediting
  the snark property when a weaker obstruction already suffices.
- Least-used relevant sidecar: the actual common order of boundary edges,
  the colour gauge, and boundary-extension multiplicity. THM-4113 additionally
  retains cyclic boundary order and noncrossing Kempe pairings.

## Portfolio and live concepts

| Lane / object | Predicate and operation | Lost information / decisive probe |
|---|---|---|
| Anchor: K5, K3,3, Petersen | Find a structural triple, test Seidel spectra and switching | Planarity and edge-colourability can change under switching; test small controls |
| Anchor: primitive triples | Euclid lift from square-sized shores of a marked clique | Bare clique order forgets partition; root must be checked explicitly |
| Niche: cubic boundary vectors | Tensor swap and sum/difference preserve exact compatibility | Scalar norm loses which boundary states overlap; compare Petersen with prism |
| Niche: maximum cuts | Maximum-cut sides are matchings and isolates | Normalize without arbitrary choice of cut; compare all maximizing cuts |
| Wildcard: Pell and conference families | Compare positive set substitution and Seidel-conference extension | Cubic degree is not preserved; test the first descendant |
| Wildcard: asymptotic geometry | Does approaching the isosceles point predict snarks? | Test connected planar triangle expansions and flower-snark cut defect |

Cards used from META-PATTERNS: **Search the statement before the method**,
**Type every analogy and every implication**, **Controlled forgetting and
unlabeled quotients require a sidecar**, and **Find the hidden second
coordinate in a nearly true theorem**. No new universal method card is
promoted from this one session.

## Connections and decisions after the first probes

1. The difference-one spine and the isosceles Pell spine are different
   arithmetic limits. Keep both; neither classifies graphs on its own.
2. The Seidel radii of Petersen, K5, K3,3 really are 3,4,5. An explicit
   switch relates Petersen to an isolated vertex plus L(K3,3), while the
   inherited disjoint-edge operation relates K5 to Petersen. This is an
   operation-level bridge, with an explicit obstruction-loss boundary.
3. Marked clique partitions give a faithful arithmetic construction:
   (x-y, 2 sqrt(xy), x+y). The root partition is 4+1 in K5; its cut is
   K4,1, not K3,3. Complete bipartite K3,3 gives the degenerate 0,6,6.
4. The graph-side promising invariant is boundary compatibility. If C is
   the labelled Tait-colouring count and s=|E|-MaxCut, maximum-cut pole
   norms satisfy FG=6^|V|/3^s. The normalized tensor overlap is therefore
   q=C^2 3^s/6^|V|. The proof and exact tests belong in the synthesis/audits,
   not in a numerical analogy. Its zero detects uncolourability, but
   approaching zero needs hostile testing.
5. Petersen's two five-cycle poles give equal nonzero orthogonal vectors.
   Investigate the minority-colour pair as an edge of K5; this preserves
   actual boundary gluing rather than merely matching a scalar.

## Artifacts

- [Arithmetic audit](graph_pythagorean_20261018_arithmetic_audit.md)
- [Graph audit](graph_pythagorean_20261018_graph_audit.md)
- [Inheritance and invariant proof](graph_pythagorean_20261018_inheritance.md)
- [Final synthesis](graph_pythagorean_20261018_synthesis.md)

None of these developments is a new proof of Kuratowski's theorem or the
Petersen-minor-free cubic colouring theorem. Any such claim must transport
minor/subdivision or boundary-extension data through the actual operations.
