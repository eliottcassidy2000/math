# Independent review of full frozen-observer rigidity

**Verdict: PASS** for `triples/PROOF.md`, Theorem9 and Corollaries10–11, with the
stated labelled, fixed-incidence, integer/rational coefficient assumptions.
Reviewer: Session292 manuscript lane, independent of the construction/proof
producer. No full re-audit of external realization-fibre theorems is claimed.

## Hand proof review

1. `Q_st` really contains H(T[S]) as its coefficient at every nonempty labelled
   interior subset S. Thus equality implies identical cyclic-triple predicates,
   since the two triple counts are1 and3 over the stated coefficient ring.
   This would not hold for equality only modulo2.
2. Let D be the flipped-pair graph. A vertex outside a D-component has unchanged
   incidences with every vertex in that component. If it separates the ends
   of a flipped edge, exactly one of the two corresponding triples is cyclic.
   Equality of the triple predicates forbids this. Connectivity propagates
   uniform adjacency, so every D-component is a module in both interiors.
   All intercomponent arcs agree. Singleton components cause no exception.
3. The degree-two coefficient at u,v in the port-start inventory is
   `r_u A_uv+r_v A_vu`. Its change is `(r_u-r_v) Delta A_uv`. On a flipped
   pair, Delta A_uv is±1, so equality forces r_u=r_v. Every retained port
   therefore sees each D-component uniformly. Only degree-two port data
   are needed; demanding full port-start equality is stronger than necessary.
4. For a fixed labelled subset of one module and a prescribed number k of
   visits, the count of ordered nonempty path runs is the sum over ordered
   partitions of the product of each block's H. This uses no unproved
   endpoint invariance: uniform outside incidences make internal endpoints
   irrelevant to connection with the unchanged external skeleton.
5. Every boundary-to-boundary path decomposes uniquely into maximal module
   runs. Its *admissible* quotient sequence is determined by common quotient
   arcs, and its labelled resources are partitioned across the runs. Summing
   proves coefficientwise equality of each eliminated boundary entry, including
   closed excursions with the same retained start/end. Direct boundary edges
   are fixed separately. Arbitrarily many retained ports are allowed: the
   argument is entrywise and requires the degree-two guard for each port.

Thus no seam can alter even one boundary coefficient while keeping these
all-subset and port guards fixed. The result is a rigidity theorem, not a
physical construction or a theorem that scalar H identifies tournaments.

## Cubic compression dependency

Read the corrected§3 of
`01-canon/theorems/THM-4803-intransitivity-is-the-relation-modulo-its-modular-tree-and-a000568-from-the-prime-cycle-index-census.md`.
It explicitly retains full ordered path-cover/H invariance under equality of
the labelled cyclic-triple hypergraph. Its original automorphism-order and
ordinary component-isotype invariance claims are explicitly retracted. The
corollary uses neither retracted assertion.

Equality of cubic Q_st data persists on restriction to every labelled subset.
Applying the retained H invariant there gives the full subset-count hypothesis
of Theorem9. This proves Corollary10 **relative to the cited structural input**.
The source's linked correction file was not available in this checkout during
this review; the canonical corrected statement was read. We therefore do not
claim to have reconstructed or audited that external/module classification.

Corollary11's iff with augmented cyclic-triple equality also checks. Full
boundary equality gives the interior triples and quadratic port coefficients.
For a fixed port incidence and interior pair, flipping that pair changes the
port triple exactly when its incidence bits differ, exactly the same condition
as a quadratic coefficient change. Triples with two or three retained ports
are fixed already. The reverse implication is Corollary10. This remains a
relative, fixed-incidence statement, not an unlabelled reconstruction claim.

The degree bound is3 for the interior observer and2 for each port observer;
the number of labelled coefficients grows with the interior. There is no
constant-memory or fast-enumeration consequence without further work.

## Independent hostile and finite controls

`rigidity_checks.py` builds all directed simple paths by **literal permutations**,
not the producer's dynamic program or inverse algebra. It groups tournaments
by their labelled cyclic-triple counts, then tests all degree-two-compatible
port incidence masks. It checks common modules, port uniformity, full start
and end inventories, and every ordered pair of allowed port cuts.

- Interior order3: all8 tournaments, all16 same-C3 pairs,136 ordered port pairs.
- Interior order4: all64 tournaments, all336 same-C3 pairs,5328 ordered port pairs.
- Interior order5: all1024 tournaments are generated;595 of9912 same-C3 pairs
  are tested by deterministic1/17 sampling, covering17428 ordered port pairs.
  This is **not** an exhaustive order5 pair test.

All tests pass. The directed-triangle versus converse with a port beating only
one vertex is the sharp lower-degree hostile: degree-one responses match,
while a flipped pair changes its degree-two start coefficient. The theorem
correctly excludes it. Removing the fixed-incidence hypothesis, dropping
labels, or retaining only a final scalar is outside the reviewed result.

Run from repository root:

```sh
python3 05-knowledge/results/session292_triple_bridges_20261113/audit/rigidity_checks.py
python3 -O 05-knowledge/results/session292_triple_bridges_20261113/audit/rigidity_checks.py
```

Normal and optimized outputs agree with `rigidity_checks.json`. Final source
hashes belong in the coordinator's final manifest after the producer finishes
integrating this review. No discrepancy was found in the new proof.
