# Independent review of the realized refinement obstruction

23 October 2026, port-rigidity agent. **PASS.** Read the full README and
checker, then reproduced the default consumer: 297 checks, receipt LF
SHA256 `5a3413d4159279dbeeee833f3e76bc57fa12905aa48d153322976db5c418a1e1`.

The explicit difference sets define actual strong self-converse circulants:
opposite pairs are disjoint, step 2 gives a Hamilton cycle, translations
give transitivity, and negation reverses arcs. Literal cyclic-triple and
strong-four enumeration establishes the common seed (6,21,112).

The all-depth obstruction is inductive, not an extrapolation from a bounded
number of rounds: a uniform seed in a common regular degree produces the
same next seed. In the common transitive-source extension, each prefix
position is a matching singleton and the core is one equitable cell.
Retaining the old color ensures that the checker reaches genuine partition
stability. The argument for every prefix length is independent of the
finite controls k=0,...,8.

The two Hamilton counts come from independent forward all-start subset DP
and backward fixed-start recursion, with direct permutation controls.
Adjacency-consistent backtracking exhausts root-fixing automorphisms;
the trivial stabilizer together with translation transitivity gives
Aut=13. A Hamilton path in a source-prefix extension must traverse the
prefix in its unique order and then the core, so the unequal H values
persist and prove nonisomorphism at every extended order.

The third edge-incidence moment 234 versus 221 repairs these two examples,
and the document correctly declines a universal claim. The stronger DRT
obstruction in the concurrent minimal-signature packet has every such
moment equal; its transitive-five repair is a different sidecar. These are
compatible boundaries, not contradictory completeness claims.

Reviewed checker LF SHA256:
`8fb9a7bd0034c9a6b574381f21b463ccaa391055a4a19ac8929f36fe99e379a3`.
No implementation change was requested.
