# Independent audit: all-order pure-blue atom classification

**Outcome: PASS.** The elementary prefix-rerouting proof establishes that
no odd circulant tournament of order at least seven is pure-blue.
Together with the previously audited odd-circulant reduction for nontrivial
strong pure-blue tournaments, the nontrivial strong pure-blue atoms are
exactly **C3 and the regular C5**, up to isomorphism. This is an all-order
proof, not extrapolation from finite enumeration.

## Mathematical audit

The audit checked every displayed path increment, permutation cycle
decomposition, sign implication, and modular reindexing formula.
In the prefix lemma the central vertex m must be fixed: the two supports
then lie in the disjoint sets {0,...,m−1} and {m+1,...,2m}.
Consequently the product of the two path-reversal maps has the same
even order as the prefix permutation. The independent direct violated-arc
table also passes and makes the automorphism-order lemma unnecessary for
the particular classification.

The three residual sign patterns are exhausted. Reindexing uses ±2,
which are units for every odd order, including composite orders.
The smallest boundary m=4 is valid: residual pattern B first maps to A
and requires a second reindexing. The four order-seven path certificates
and both order-five normalized positive controls were checked literally.
The positive C5 control does not assume that regularity or circulancy
implies purity.

The strong-to-odd-circulant reduction imports Camion Hamiltonicity and
the cited Lu crossing theorem. Those source dependencies belong to the
previously audited reduction; the new prefix proof does not use the
Hamiltonian-cycle counting bound, the order-131 cutoff, the finite
circulant sweep, or the new path-rigidity classification.

## Independent executable evidence

`independent_prefix_obstruction.py` contains no producer imports or data
dependencies. It independently constructs and verifies:

- Every normalized skew connection set, with 1 in S, at odd orders
  9 through 31: **32,760 presentations**, zero failures.
- The actual directed Hamiltonian path and a failed anti-automorphism
  arc in every presentation; also the prescribed failed arc from the
  direct certificate table.
- The prefix permutation's even order, disjoint supports, and the exact
  order of the product of reversal maps.
- All four order-seven certificates.
- Every Hamiltonian path for the order-three control and both normalized
  order-five controls: respectively 3, 15, and 15 paths, all pure-blue.
- All three residual sign patterns at every m from 4 through 200:
  **591 structured cases**, zero failures.

Reproduce in a directory containing the standalone script:

```sh
python3 independent_prefix_obstruction.py --max-n 31
```

Matching results are recorded in `independent_prefix_obstruction.out`
and `independent_prefix_obstruction.json`. These finite checks validate
the implementation and boundaries; the symbolic argument supplies the
all-order quantifier.

## Consequence and scope

Combining the result with the separately audited SCC center law and
the all-order rigid alphabet {singleton,C3} gives the exact pure-blue
node series, including the empty object,

    P(z) = (1 + z + z^3 + z^5) / (1 − z^2 − z^6).

The series consequence uses those additional structural dependencies;
it is not inferred from the prefix computation alone. The classification
does not assert that counts or moment constraints suffice to realize a
merged metagraph arrangement, and it does not equate pure-blue with
path rigidity: C5 is pure-blue and nonrigid.

## Separate geometry corollary: the central C5 has three distinct branches

This additional argument also passes. For T=A→C5→Aop with |A|=r≥1,
the original five center scores are r+2, while the left scores are at
least r+5 and the right scores are at most r−1. Free-flipping the
nonpath arcs sends an interior score s to n−1−s, with only the global
first vertex receiving +1 and the global last vertex receiving −1.
The resulting center is therefore exactly the five vertices of score
r+2; every outer vertex has score at most r or at least r+4.

The induced tournament on this intrinsically identifiable center is one
of the three free-flip partners of C5. Independent enumeration of all
5! permutations gives (Hamiltonian-path count, automorphism order)
(11,1), (15,3), and (13,1), respectively, for the center paths
01234, 02341, and 02413. The different path counts distinguish them
even modulo converse. Thus the three global neighbors are distinct.
The r=0 case is the same direct order-five calculation.

The standalone verifier includes this calculation under `center_five`.
Combined with the audited fiber multiplicativity and pure-blue
classification, center-C5 nodes have three distinct blue neighbors;
the other pure-blue nodes have degree one. For n≥4, all these
neighbors are mixed, as follows from the separately audited
triangle-rigid independent-set law or from the strong/non-strong
reduction and the small-order controls. This line-model statement
does not include the degenerate n=2 cube or the n=3 exception.

## Exact count, line-budget, and size-changing-bijection corollaries

The full note `PURE_BLUE_COUNTS_AND_BRANCH_GEOMETRY.md` also passes
symbolic audit. If a_0=1, a_j=0 for j<0, and
a_j=a_(j−1)+a_(j−3), then

    p_(2m)=a_m,
    p_(2m+1)=a_m+a_(m−1)+a_(m−2),
    L_(2m)=a_m,                 J_(2m)=0,
    L_(2m+1)=a_m+a_(m−1),      J_(2m+1)=a_(m−2).

Here L and J count mass-one and mass-three pure-blue nodes.
The recurrence p_n=p_(n−2)+p_(n−6) holds starting at n=6.
The displayed counts at orders 3 through 16 are correct.

For n≥4, exactly L_n+3J_n blue lines meet a pure-blue node, and
all terminate at mixed nodes. Subtracting this number from
2^(floor((n−1)^2/4)−1) gives the exact number of mixed-to-mixed
blue lines, including loops and multiplicities. At n=3 the unique
blue line has two pure-blue endpoints, so its two half-edges cannot
be counted as two lines.

Finally, the identity p_(2m+1)=a_(m+2)=p_(2m+4), for m≥0, has a
valid explicit bijection. Replace the singleton, C3, or C5 center
by appending the suffix 11, 3, or 31 to the left rigid word,
respectively. Every rigid word of weight at least two ends in
exactly one of these three disjoint suffix classes. Removal recovers
both the center and the original left word. The map changes order
by three; no graph isomorphism or preservation of complement
neighbors across different orders is asserted.
