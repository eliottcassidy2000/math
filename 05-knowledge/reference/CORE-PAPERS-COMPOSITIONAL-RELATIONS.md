# Core-paper sidecar: knots and compositional relations

This sidecar preserves detailed records routed compactly from
`CORE-PAPERS.md`. Freshness dates are per entry.

## Brittenham--Hermiller -- *Unknotting number is not additive under connected sum*

- **Primary / freshness:** [arXiv:2506.24088v2](https://arxiv.org/abs/2506.24088),
  revised 2025-09-15. **PREPRINT v2; PDF and stated examples checked
  2026-07-24.**
- **Imported role:** gives the first counterexample to connected-sum
  additivity:
  `u(T(2,7))=u(mirror(T(2,7)))=3` while the mirror-paired sum has unknotting
  number at most five. Gordian adjacency propagates this shortcut to infinite
  families, and repeated constructions make the gap arbitrarily large.
- **Repo consumers:**
  [THM-2176, continuation profile and interaction cocycle](../../01-canon/theorems/THM-2176-gordian-continuation-profile-and-interaction-cocycle.md)
  and the earlier
  [common-intermediate LRC synthesis](../../07-reflections/puzzle-atlas-adaptive-prime-spectral-quantum-endpoint-cocycle-opus-20260724.md).
- **Does not prove:** the exact unknotting number of
  `T(2,7)#mirror(T(2,7))`; the paper leaves it between its known lower and
  upper bounds. It does not supply a finite complete knot invariant, a
  connected-sum congruence for `u`, or additivity after homogenization.
  THM-2176's cocycle, universal-profile theorem, and pure-bypass
  decomposition are repo-derived consequences.

## Brittenham--Hermiller -- *Unknotting number and connected sums: The knots 4_1 and 5_1*

- **Primary / freshness:** [arXiv:2601.18757v1](https://arxiv.org/abs/2601.18757),
  submitted 2026-01-26. **PREPRINT v1; checked 2026-07-24.**
- **Imported role:** introduces the term **symbiont** for a pair with strict
  unknotting defect and proves the certified examples
  `u(4_1#9_10)<=3<1+3` and `u(5_1#8_2)<=3<2+2`, with the stated mirroring
  conventions. It also propagates examples through minimal unknotting
  sequences.
- **Repo consumer:** [THM-2176](../../01-canon/theorems/THM-2176-gordian-continuation-profile-and-interaction-cocycle.md),
  which identifies symbiosis as the support of a weighted 2-coboundary rather
  than a tournament edge.
- **Does not prove:** that the candidate trefoil partner `10_6` has
  unknotting number three; the paper explicitly retains `u(10_6) in {2,3}`.
  It does not classify symbionts or determine all composite unknotting
  numbers.

## Zakharov -- *An isoperimetric inequality for word overlap*

- **Primary / freshness:** [arXiv:2602.20143v2](https://arxiv.org/abs/2602.20143),
  revised 2026-03-03. **PREPRINT v2; checked 2026-07-24.**
- **Imported role:** for two length-`n` word families with no directed
  suffix--prefix overlap, proves the sharp upper bound
  `mu(A)mu(B) <= (1/n)(n/(n+1))^(n+1)`. The reusable proof move is a disjoint
  first-hit decomposition, a shift-density convolution, and optimization at
  a positive generating-function root.
- **Repo consumer:** the
  [Gordian/min-plus relation synthesis](../../07-reflections/the-gordian-object-is-a-min-plus-profile-not-a-tournament-codex-20260724.md).
- **Does not prove:** an LRC word model, a cyclic-overlap inequality, or a
  tournament theorem. A cyclic LRC application must first retain cut,
  rotation, reflection, and endpoint-owner data and verify the convolution
  on the actual packet.

## Rybin--Zhang--Luo -- *XX^t Can Be Faster*

- **Primary / freshness:** [arXiv:2505.09814v2](https://arxiv.org/abs/2505.09814),
  revised 2025-05-16. **PREPRINT v2; checked 2026-07-24.**
- **Imported role:** exploits transpose structure to give a `4 x 4` block
  scheme with eight recursive symmetric calls and 26 general products. Its
  discovery pipeline generates rank-one candidates, enumerates exact linear
  relations, and solves a second minimum-cover problem over the targets.
- **Repo consumer:** the
  [Gordian/min-plus relation synthesis](../../07-reflections/the-gordian-object-is-a-min-plus-profile-not-a-tournament-codex-20260724.md),
  as a proposed exact span-and-cover search architecture.
- **Does not prove:** any knot, tournament, or LRC statement. Exact algebraic
  span does not retain LRC phase owner, scale, or analytic uniformity.

## Krenn--Gu--Soltesz -- inherited vertex coloring of perfect matchings

- **Primary / freshness:** [*Questions on the structure of perfect matchings
  inspired by quantum physics*, arXiv:1902.06023v2](https://arxiv.org/abs/1902.06023),
  and the [maintained problem page](https://mariokrenn.wordpress.com/graph-theory-question/).
  **ORIGINAL PROBLEM / dynamic status page; checked 2026-07-24.**
- **Imported role:** a perfect matching inherits one endpoint color at each
  vertex, and the weight of a coloring is the sum of matching-weight products
  over the entire fiber inducing it. This is the precise sum-product model
  behind the repo's whole-fiber cancellation comparison.
- **Repo consumer:** the
  [Gordian/min-plus relation synthesis](../../07-reflections/the-gordian-object-is-a-min-plus-profile-not-a-tournament-codex-20260724.md).
- **Does not prove:** that the matching object is a tournament, a min-plus
  Gordian model, or an LRC certificate. Pairwise orientation loses edge
  absence, multiedges, endpoint colors, matching multiplicity, and complex
  amplitudes. Current claims on the maintained page require their own
  paper/formal-proof audit before entering the proved dependency graph.
