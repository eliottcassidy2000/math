# Tournament composition: the 81/27 boundary

2026-10-26. The user's 12-vertex example is the first possible failure of
automorphism-order invariance under equality of labelled cyclic triples.
The investigation gives exact structural and arithmetic classification
results, with an operation-aware starting point for the metagraph atlas.

| Result | Scope and source |
|---|---|
| Decorated canonical module trees classify ordinary tournaments; automorphisms retain the child-colour stabilizer | **PROVED** from **CITED** modular decomposition; [framework](CLASSIFICATION.md) |
| The p-part of automorphism order is constant in a C3-fibre below 4p, sharply for every odd prime p | **PROVED**; [fibre F2](fibre/RESULT.md) |
| At order 12 there is exactly one bad C3 hypergraph up to isomorphism; 128 labelled lifts, four ordinary classes, Aut orders 81/27/27/81 | **PROVED + FINITE-EXACT**; [fibre F3](fibre/RESULT.md) |
| All ordered path-cover counts survive the C3 map; H=462915 throughout this fibre, but H/Aut changes | **PROVED**, displayed value **FINITE-EXACT**; [fibre F4](fibre/RESULT.md) |
| All tournaments attaining the maximum p-part are Sylow-adapted towers of circulant p-tournaments on base-p digit blocks | **PROVED**; [arithmetic classification](arithmetic/RESULT.md) |
| Maximum p-valuation is v_p(n!); for strong tournaments it is v_p(n!) minus one exactly when the base-p digit sum is two | **PROVED**, n>=3 for strong case; [arithmetic bounds](arithmetic/RESULT.md) |
| Order-join is a relocation Cartesian product; free-arc complement transports it to an associative reversed join with marked hinges | **PROVED**; [metagraph M1–M3](metagraph/RESULT.md) |
| One blue loop already has two composition-distinguishable half-edges; the first two-factor example is order 1+4, with resulting H=11/13 | **FINITE-EXACT minimality**; rooted response convolution **PROVED**; [metagraph M4–M5](metagraph/RESULT.md) |

For the emblematic pair set `A=C3→1`, `B=1→C3`. Then `C3[A,A,A]` admits
the rotation of its three equal children, while `C3[A,A,B]` does not. The
internal groups are `C3³` in both. Thus the groups are `C3³⋊C3` and `C3³`.
The observer preserves the cyclic triples but forgets the child colours.
For a fixed observer hypergraph C, the valid aggregate identity is
`sum_[T] 1/|Aut(T)| = |R(C)|/|Aut(C)|`; uniform orbit division is unwarranted.

At n=12, base three gives 110. The unrestricted maximum 3-part is 243,
attained by a 9-block joined to a 3-block; the strong maximum 3-part is 81.
These are prime-specific extrema, not a theorem identifying the maximum
total automorphism order by combining separate prime maxima.

The [random recovery](RANDOM_CONNECTIONS.md) freezes an actual 24-file sample
and follows five of its signals. The useful connections are digit carries,
action cycle indices, and the distinction between independent module choices
and transport constraints. A suggested Kummer identification has no supplied
map and is explicitly stopped. The Sylvester reflection's erroneous recurrence
was repaired where found.

## Correction lineage and inheritance

This packet repairs the automorphism and strong-component-isotype clauses of
[THM-4803](../../../01-canon/theorems/THM-4803-intransitivity-is-the-relation-modulo-its-modular-tree-and-a000568-from-the-prime-cycle-index-census.md),
its [source](../a000568_intransitivity_families_20261026/INTRANSITIVITY.md),
and the mistaken identification of the `A_n` orbit split with prime converse
choices. It preserves the published C3-realization theorem and the earlier
finite census. The original working-tree base is commit `fc3768ac5fcd`;
old statements remain recoverable there. The mistake mechanism is recorded
in [MISTAKES.md](../../../01-canon/MISTAKES.md).

The closest mechanism is canonical modular decomposition; the canonical
hostile is loss of a repeated-child stabilizer; the corrected near miss is
finite invariance below the first possible failure; the least-used sidecars
are rooted half-edges and permutation cycle indices. Prior enumeration and
the order-join metric law are inherited, not claimed as new.

## Audit and reproduction

[AUDIT.md](AUDIT.md) records independent proof review and the finite universes.
Every check uses exact integers and explicit runtime gates, including under
Python optimization. No external package is needed. From the repository root:

```text
python3 -B 05-knowledge/results/tournament_composition_classification_20261026/verify.py --replay
python3 agents/check_docs.py
```

The replay verifies the manifest before computing, uses temporary output
files, compares normal and optimized results to frozen receipts, and reports
failure on any mismatch. It does not rewrite the committed expectations.

**OPEN:** an effective atlas of all prime tournaments, a small closed family
of rooted responses for general substitution and coloured metagraphs, and
LRC(14). The proved recursive classifier does not remove those obligations.
