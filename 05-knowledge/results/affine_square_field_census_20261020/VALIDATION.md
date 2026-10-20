# Validation and publication receipt

2026-10-20. **VERIFIED.** Companion to AFFINE_SQUARE_FIELD_CENSUS.md.
The all-q family theorem is a written elementary proof; the full q=25
automorphism census is FINITE-EXACT.

## Independent audit findings

1. The finite universe is exactly four odd unordered-pair orbits. The verifier
   constructs the complete group H of order75 and computes its pair orbits,
   obtaining75,75,75,75. It checks each of16 graph and16 tournament objects
   against every element of H, verifies graph/tournament type conditions,
   and rejects duplicate objects.
2. Rooting at0 is lossless because translations act transitively on every
   object. Any isomorphism can be post-composed with a target translation to
   fix0; every automorphism has a unique such factorization. Therefore rooted
   automorphism counts multiply by25, and rooted nonisomorphism implies full
   nonisomorphism.
3. The equitable-partition algorithm exhausts every possible image of each
   individualized vertex and rejects only unequal invariant count classes.
   The independent MRV constraint search uses adjacency in both directions
   and common-neighbor counts. Both check the entire adjacency relation at
   every leaf; their returned sets of permutations agree exactly, not just
   their cardinalities. Neither uses a node cap or pruning heuristic capable
   of discarding an isomorphism.
4. All nontrivial rooted automorphisms equal the rooted normalizer stabilizer
   lists. Empty and complete graphs are handled by the direct S25 proof, not
   by an infeasible25! enumeration. The only equal-degree cross-orbit test,
   masks3 versus5, is independently rejected by both searches and by their
   different common-neighbor histograms.
5. Full sign-evenness follows from rooted signs plus translations of odd
   order. The sign computed is edge-orientation reversal, not the sign of an
   unordered-edge permutation or parity of vertex degree.
6. The full normalizer assertion is proved for q=25 using the characteristic
   translation group and the irreducible C3 centralizer. It is deliberately
   not extrapolated to the infinite family.
7. In the infinite family, the semilinear action on directed classes modulo8
   is sufficient: one source orbit proves the source count, and explicit
   odd-sign stabilizers eliminate every nonempty graph. The induced action
   has order16; the actual field Frobenius need not have order2. This avoids
   an otherwise tempting but false identification of the full Frobenius
   order with its quotient action.
8. The sign-witness table was checked for both p≡3 and p≡5 mod8. The q=25
   verifier additionally compares the compressed sign formula with actual
   labelled-edge signs in288 stabilizing cases. The family proof retains
   the necessary odd pair-orbit size qr.

The two isomorphism searches are algorithmically independent but share the
same explicitly defined field arithmetic and adjacency input. The separate
signed-cube proof provides an independent route to the global one-versus-one
count; no claim of completely independent input encodings is made.

## Completed checks

Status: PASS.

Universe: All 16 H-fixed tournaments and 16 H-fixed graphs on F25; all rooted isomorphisms for nontrivial normalizer representatives; exact signed four-cube actions for p mod 8 in {3,5}.

Main script SHA-256: `56aaece4963399158fc77d98efdbc4a5e782d78d13e719d02a7152c4f14f4881`.

JSON output SHA-256:
`51a0b5a465d84ca63cacf728267ae13d8025f6a140d262541355f172cecaf7f8`.

Normal and optimized JSON files are byte-identical. Checks use explicit
exceptions, not Python assert statements.

Reproduce from this directory:

    python3 verify_affine25.py --output verification.json
    python3 -O verify_affine25.py --output verification.optimized.json
    cmp verification.json verification.optimized.json

`probe.py`, `probe.json`, `independent_probe.py`, and `audit.*` are inherited
or intermediate research artifacts and are not in the publication inventory.
The consolidated verifier contains both audited searches.

## Scope boundaries retained

- Ordinary isomorphism counts admitting the stated embedded H, not all odd
  symmetry subgroups of a given order.
- Finite fields, not cyclic rings of prime-power order.
- A proved infinite family with v₂(q−1)=3, not every prime power.
- Full automorphism orders proved computationally only at q=25.
- No conclusion that H is pronormal in S_q.
- No general claim that scalar source fusion and graph sign escape cancel.

The closest next test is the eight-coordinate fixed cube at v₂(q−1)=4,
retaining the actual Frobenius residue in its signed action.
