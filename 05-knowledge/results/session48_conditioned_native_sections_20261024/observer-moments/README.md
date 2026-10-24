# Observer moments: closed proof and exact checks

Read `THEOREM.md` for the all-size proofs and exact scope.

Main results:

- Uniform marked tournament/Hamiltonian-observer records satisfy
  `E Q=n^2-3n+4` and `E c3=(binom(n,3)+n-2)/4` for `n>=2`.
- Fixed-observer class fibres have size `H/Aut`; literal converse merging
  gives `H/Aut` or `2H/Aut`. Every original order-eight atlas fibre passes
  an independent Hamiltonian-count/automorphism audit. Its triangle and
  score-square checksums are `32,505,856` and `92,274,688`.
- Arbitrary coherent free-host completion has
  `E Q=2f+average_O ||b_O||^2`, with no measure assigned to an empty domain.
- Four explicit two-owner records give the nonreal-root polynomial
  `(1+z)^21(4+10z+10z^2)`. The converse-symmetric disjoint union, weighting
  ALL forward source arcs, gives `2z^2(1+z)^22(5+2z+5z^2)`.
- These are marked counting laws, not isomorphism-enumeration formulas,
  capacity transport, routing, or self-adjoint spectral statements.

Reproduce from this directory with Python 3.10+ and the standard library:

```sh
python3 check.py --output fresh-normal.json
python3 -O check.py --output fresh-optimized.json
cmp fresh-normal.json normal.json
cmp fresh-optimized.json optimized.json
```

For a fresh-directory, manifest-only reproduction:

```sh
python3 verify_packet.py
```

The latter verifies hashes, copies only listed packet files into a new
temporary directory, and executes both normal and optimized checks there.
No Git operations, external packages, network, or parent workspace files are
used. Snapshots of inherited prose, the canonical decoder, explicit source
matrices, and the original atlas have `.txt` extensions and are never imported
as code. `sources/SOURCE_PINS.json` records provenance; `MANIFEST.json` pins
every deliverable byte except the manifest itself.

The original atlas's `canon_a` encoding is MSB-first binary on lexicographic
owner pairs. This differs from the checker’s internal LSB-first census code;
the dedicated atlas decoder is explicitly separate. It is checked against
four frozen explicit source matrices before reading atlas statistics.
