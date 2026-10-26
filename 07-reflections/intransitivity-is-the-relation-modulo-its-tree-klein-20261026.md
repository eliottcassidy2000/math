# Intransitivity is the relation modulo its tree

> **HISTORICAL CORRECTION (2026-10-26):** the Aut-invariance assertion below is
> false, first at twelve vertices (orders 81/27). H-invariance and the cited
> fibre theorem survive. Current proof: [sharp 4p repair](../05-knowledge/results/tournament_composition_classification_20261026/fibre/RESULT.md).

klein, 2026-10-26, sixth directive of the session ("fully enumerate A000568 from the tournament
metagraph families, merge in work from many months prior regarding tournaments as information,
binary intransitivity in relation itself"). Reflection; truth lives in
`05-knowledge/results/a000568_intransitivity_families_20261026/` and THM-4803.

## What the directive turned into

Three phrases, three objects. "Families" became the modular decomposition tree: every class of
A000568 is a tree of linear nodes (strong components in order) and prime nodes (a prime quotient
with blocks), and the species identity session56 wrote down last week computes A000568 once the
prime quotients are supplied with their cycle indices. Supplying them (nauty, `n <= 10`) is the
enumeration part. "Intransitivity in the relation itself" became the hypergraph of cyclic triples,
`C3(T)`, and the exact statement of what it knows: the relation up to the symmetry group of the
tree, `∏ S_{k_N} × (Z/2)^p`. "Tournaments as information" became the Hartley fibre of that
observable, which MISTAKE-231 had already told us is the only honest information measure here.

## What was already known and had to be found first

The fibre theorem is published twice (Boussaïri–Ille–Lopez–Thomassé 2004 for prime tournaments;
Boussaïri–Chergui–Ille–Zaidi 2018, Corollary 57, for the count). I rediscovered the count by brute
force before finding the citation, and the right record is: cited, verified on every class to
`n = 7`, and extended by the counts. The lesson repeats the session's earlier ones: the repo's
own files (THM-1960's seed census, session56's species, MISTAKE-231's frame, THM-4739's
"second-order information is the cyclic-triangle 2-form") already pointed at the object; the
literature had the theorem; what was missing was the census and the connection.

## What is new, by my reckoning

- Lemma C with an elementary proof: the strong components are the connected components of the
  intransitivity hypergraph. (Equivalent to the root case of the shared-strong-modules result.)
- `H` and `|Aut|` are intransitivity invariants; scores are not. This explains THM-931's wild
  pairs at once.
- The counts: intransitivity structures `1, 1, 2, 3, 9, 32, 247, 3328, 93612, 4774417` and the
  labelled `1, 1, 2, 11, 208, 8563, 680426, 101747483, 28904318404, 15824793243139`, both absent
  from the OEIS, both reproduced by a species equation with `Aut±` (automorphisms and
  anti-automorphisms) — the dual-symmetry group is the right group for intransitivity counting.
- The information reading: one bit for almost every relation, `log₂ n!` for the transitive one,
  `1.22` bits on average at `n = 10`; level-resolved, the loss falls monotonically with `c3`.
- The prime census to `n = 10` (A280772 extended; A003505, the 1975 table, is wrong at 7 and 8).

## What I would do next

The `n = 10` cycle-index types (not collected) would give `n = 11` predictions for both
sequences without enumerating the 903,753,248 classes. The realizability characterization of
3-uniform hypergraphs (BCIZ 2018) suggests counting all realizable hypergraphs directly as a
hypergraph species; `c(n)` would then be a check. And HYP-9351's port question is the place where
"the 2-form determines the relation up to the tree" has to be pushed one level up, to the
metagraph lines; I have not done that.
