# Independent finite audit: prime cardinalities do not preserve square-sum graphs

Status: **PASS.** The audit enumerates exactly the 105 unordered pairs from
the supplied 15-prime support, tests membership of each pair sum in the exact
integer square set, and reconstructs components independently with union-find.
It does not re-prove which primes divide the Monster or Baby Monster orders.

The complete graph has 14 edges:

    2–7, 2–23, 2–47, 3–13, 5–11, 5–31, 5–59, 7–29,
    13–23, 17–19, 17–47, 23–41, 29–71, 41–59.

It is connected on all 15 vertices, hence a tree. Its leaves are exactly
3,11,19,31,71, so it cannot have a Hamiltonian path. A spanning path in a tree
would use every one of its 14 edges and would leave exactly two leaves.

The three induced blocks contain respectively:

- {17,19,23,31,47}: edges 17–19 and 17–47, with 23 and 31 isolated.
- {29,41,59,71}: edges 29–71 and 41–59.
- {2,3,5,7,11,13}: edges 2–7, 3–13 and 5–11.

These are not isomorphic to the connected G25 blocks of the same respective
sizes 5,4,6. Cardinality equality is therefore not an adjacency-preserving
identification for the proposed square-sum relation.

Removing {29,41,59,71} gives exactly the two components

    {2,3,7,13,17,19,23,47},
    {5,11,31}.

Restoring the removed primes restores the connecting path 23–41–59–5 and
the pendant path 7–29–71. All claimed edges and component statements in the
corrected root note agree with the independent census. The retained history
of the rejected initial hand-list is appropriately labelled; the audit
relies on the corrected complete edge list, not that informal prediction.

The source note, checker and receipt are pinned in `prime-certificate.json`.
`audit_prime_support.py` has no dependencies beyond the standard library and
uses explicit checks. Normal and optimized outputs agree byte-for-byte.
