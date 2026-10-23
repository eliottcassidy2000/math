# Independent addendum: full-support walks recover the C3-structure

Verdict: PASS for the proof communicated by the context producer and independently checked below. This is a literature-linked consequence, not an independent novelty claim.

## 1. Full-support detector

Fix a labelled n-vertex tournament, n≥3, and a triple S. Give each vertex of S multiplicity n−1 and every exterior vertex multiplicity1. The corresponding commutative directed-walk coefficient is positive exactly when S is cyclic.

If S is cyclic, take any directed Hamiltonian path and insert n−2 complete directed rounds of the triangle at an occurrence of one of its vertices. Each round returns to that vertex and adds one occurrence of every vertex in S; outside multiplicities remain1.

If S is transitive, the n−3 exterior occurrences divide a walk into at most n−2 runs lying in S. Each run is acyclic and visits each S vertex at most once. It cannot supply the required n−1 occurrences of any S vertex. This establishes both directions without asymptotic or genericity assumptions.

## 2. Literature bridge verified

Independently opened Boussaïri, Ille, Lopez and Thomassé, “The C3-structure of the tournaments,” author-hosted primary PDF:
https://perso.ens-lyon.fr/stephan.thomasse/liste/bilt.pdf

Theorem2 states that two tournaments on the same vertex set have equal cyclic-triple families exactly when one can be obtained from the other by successive inversions of homogeneous intervals. The paper explicitly defines each inversion using an interval of the then-current tournament. Corollary1 gives the same-or-dual conclusion for indecomposable tournaments.

The imported theorem is not being reproved or claimed as new here.

## 3. Walk-series consequence

For a uniform substitution, splitting a walk into maximal block runs proves ordinary formal substitution of the internal walk series into the quotient walk series. Commutativity of vertex-content variables makes reversal preserve an internal tournament's walk series. Therefore each homogeneous-interval inversion preserves the full commutative walk series.

Combining this with the full-support detector gives equivalence of:

- full-support walk-series equality;
- full walk-series equality;
- equality of the labelled cyclic-triple family;
- homogeneous-interval inversion equivalence;
- equality of every vertexwise nonempty tournament-substitution Hamiltonian response.

The final equivalence also uses the audited Stirling inversion from transitive probes. The order1 and2 cases are immediate and should be separated from the triple detector. This proves the earlier support restriction does not create a stricter equality relation on tournaments, although the substitution formula directly consumes only full-support coefficients.

## 4. Independent exact checks

`independent_c3_walks.py` / `c3_walk_receipt.json` checks every labelled tournament through n5:

- 10,504 exact full-support triple detectors;
- 3,830 homogeneous-interval inversions;
- every full walk coefficient with individual exponents0..2 is constant within each cyclic-triple class.

The numbers of labelled cyclic-triple structures for n1..5 are1,1,2,11,208. Additional deterministic controls perform560 detector tests on112 tournaments of orders6 through8. Normal and optimized outputs are identical.

These tests corroborate the elementary detector and conventions. Equality of the entire infinite walk series uses the interval-inversion theorem and maximal-run proof, not the finite truncation.
