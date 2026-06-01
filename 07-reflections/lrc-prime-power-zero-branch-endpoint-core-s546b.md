---
source: codex-2026-06-01-S546b; codex-2026-06-01-S548
status: theorem-level formalization (THM-391) + supporting computation
tags: [lonely-runner, p-adic, prime-power, cover-core, trienerment, tournament-analysis, gabor, zero-branch, endpoint-core]
---

# Prime-power zero branches kill unit witnesses but do not hold the cover core

**Prompt (user):** this long session, consider prime power zero branch cover cores.

This is the narrow cover-core complement to the p-adic tree picture in HYP-2035.
The object is not a runner, a sector, or a time chamber.  It is a branch of the
p-adic proof tree:

```text
Z_q(V) = { v in V : q divides v }, where q = p^d <= n.
```

If `Z_q` is empty, THM-369 gives the unit witnesses `u/q`, `(u,q)=1`: no speed
lands at the observer, and every nonzero residue has distance at least `1/q >=
1/n`.  If `Z_q` is nonempty, those unit witnesses are killed.  The question was
whether the killed branch itself can become the hard cover core.

The answer from S546 was no in the exact audited rows, and S548 turns the local
reason into THM-391.  A q-divisible speed `s` contributes local danger
intervals around each unit point `u/q` with radius `1/(n s)`.  For fixed `u`,
these intervals are concentric.  Distinct unit points are separated by at least
`1/q`, while every half-width is at most `1/(n q)`, so the unit-point clusters do
not interact in the LRC range.  The outer interval has an unprotected endpoint,
peels away, and the same argument repeats.  A prime-power zero branch therefore
kills a rational witness but cannot store a nonpeeling endpoint-protection core
by itself.

The formal surprise is that primality never enters the proof.  THM-391 works
for every `2 <= q <= n` and every nonzero q-grid center set.  Prime powers are
structural labels in the p-adic tree, not special local interval geometries.

## Computation

Artifacts:

- `04-computation/lrc_prime_power_zero_branch_core_s546b.py`
- `05-knowledge/results/lrc_prime_power_zero_branch_core_s546b.out`
- `01-canon/theorems/THM-391-lrc-zero-branch-star-core-peeling.md`
- `04-computation/lrc_zero_branch_star_theorem_s548.py`
- `05-knowledge/results/lrc_zero_branch_star_theorem_s548.out`

The script audits selected `n=6,8,14,18` rows: AP walls, the `n=6` `Z_3` bridge
row, the `n=18` no-`9` row, the `18`-gate and `36`-double-gate rows, and a deep
`3`-zero-branch row.

For every audited branch:

```text
Audited local branch cores: 86, nonempty=0
Audited full cover cores: 9, nonempty=0
```

The THM-391 verifier then checks the stronger q-grid star theorem and the
explicit peel-layer formula in `3255` bounded exact cases:

```text
all bounded cores empty and all peel layers match the formula
```

The n=18 signal is the useful warning.  The p-adic zero-flow side is huge:

```text
n18_AP_wall:                 250199979305648
n18_18_gate_skip8:           250199979240064
n18_36_double_gate_skip8:    250199979240064
n18_deep_3_zero_branch:      187649984495616
```

But every cover core still peels to empty.  So HYP-2025/HYP-2017's zero-flow
abundance is real arithmetic structure, not a local cover obstruction.

Representative branch facts:

```text
n18_no_9_zero_branch: q=9 zero=0 status=open-witness
n18_AP_wall:          q=9 zero=1 radius=1/162 core=0
n18_18_gate_skip8:    q=18 zero=1 radius=1/324 core=0
n18_36_double_gate:   q=18 zero=1 radius=1/648 core=0
```

The gate ladder is visible: speed `18` kills the `1/18` unit targets, speed `36`
kills them more narrowly, but neither creates a local endpoint core.  The debt is
exported to descendant event layers.

## Tournament Analysis

Tournament vertices are prime-power branches `q=p^d` plus the optional non-prime
`q=n` gate diagnostic.

Pairwise observable:

```text
(covered?, local_core_size, zero_count, local_radius)
```

Switch/gauge: the branch with more zero debt, larger local radius, or a nonempty
core beats the other.  Exact equality is a t(r)ienerment tie.  Tie Hamiltonian
path: increasing `q`.

In every audited row, the branch trienerment collapses to a transitive debt
ladder:

```text
ties=0, strict_c3=0, SCCs=(1,1,...,1), Hamiltonian paths=1
```

That collapse is information, not failure.  It says prime-power zero branches are
ranked invoices, not cyclic repair mechanisms.  The tournament becomes
meaningful only when the vertices are refined to the exported obligations:
`(q, endpoint descendant)`, `(q, event owner)`, or product-tree pairs
`(p^a, r^b)` in rank >= 2.

## Assumption Challenge

The naive tournament has runners as vertices.  This branch view abandons that.

Preserved predicate:

- Does a THM-369 rational unit witness `u/q` survive?

Destroyed data:

- runner identity away from divisibility by `q`;
- circular order inside the chamber;
- which endpoint becomes the next proof obligation;
- cross-prime coupling;
- the Gabor phase of the killed witness.

This destruction is the point.  It isolates a single proof channel and asks
whether that channel alone can be hard.  It cannot: the zero branch is only a
gate.  The hard object, if it exists, must live in the exported labels.

## Gabor Connection

The previous observer-anchored Gabor result made the good face into a marked
zero-column target: loneliness deletes near-observer runner ties, while the
Gabor lift creates an all-zero observer column.

The p-adic zero branch is the arithmetic analogue of that move.  A q-zero branch
creates a zero at the unit column `u/q`; the local time-frequency packet around
that zero is a nested star, so it carries no cover core.  This suggests a more
interesting hybrid vertex:

```text
(prime-power branch q, unit harmonic u, exported endpoint owner)
```

That vertex remembers both the arithmetic zero column and the descendant
boundary label.  It is a plausible next place where Gabor trienerments stop
being transitive.

## External Touchpoints

The web/literature scan reinforced the split between local branches and global
proof currencies:

- Tao's remarks reduce attention to integer speeds of bounded size and improve
  the trivial lower bound, which fits the repo's finite invoice perspective:
  https://arxiv.org/abs/1701.02048
- Giri and Kravitz frame LRC spectra as distances of one-dimensional subtori
  from the central point, matching the idea that our branch is a projection of a
  richer torus object:
  https://arxiv.org/abs/2304.01462
- Fan and Sun's amended spectrum work treats common factors as significant in
  small cases, aligning with the branch occupancy signal but not identifying it
  with a cover core:
  https://arxiv.org/abs/2306.10417
- The shifted LRC/zonotope covering-radius line shows that cover structure is a
  serious global formulation; S546 says a single prime-power branch is too local
  to be that global cover:
  https://arxiv.org/abs/2506.13379
- Bedert's Riesz-product improvement is another warning that the current
  frontier uses global harmonic/probabilistic mass, not one branch in isolation:
  https://arxiv.org/abs/2511.16636

## Next Experiments

1. Replace vertices `q` by `(q, endpoint descendant)` and orient by which branch
   protects, exports, or abandons the endpoint.
2. Add event-owner labels to the branch trienerment and compare SCC leakage
   against HYP-2029's symbolic bad-word cycles.
3. Run the same audit on rank-two product-tree cases such as `n=12,20,30`, where
   cross-prime branches can couple.
4. Build artificial circular-arc systems with centered branch stars plus one
   off-center descendant layer to find the smallest nonempty core.
5. Add Gabor labels `(q,u,m)` and test whether zero-column branches become
   cyclic only after harmonic/endpoint memory is included.

## Synthesis

Prime-power zero branches are not counterexample cores.  They are kill switches
for simple rational witnesses.  The moment a branch is covered, the proof
obligation moves downward: to endpoints, event owners, product-tree coupling, or
Gabor-labeled descendant columns.  The tournament of bare branches is therefore
the wrong final object but the right diagnostic: it tells us where not to look,
and it names the labels that must be added before cyclic trienerment structure
can appear.
