# 66, 70, five and fifteen: lattice content, class inversion, and guarded bridges

2026-10-28. **PROVED elementary deductions with explicitly CITED classical
inputs; FINITE-EXACT controls; OPEN transfers.** Same-agent mathematical
review and independent computational formulations; no independent-agent
proof audit or literature-priority claim. No universal Collatz, colouring,
Jacobian, or LRC theorem is asserted.

The strongest extension in this packet is an all-size exclusion, not another
list of coincidences:

> For every integer k>=1, **72k²−6 is not idoneal**. The primitive strictly
> reduced forms `(6k−1, ±4, 12k+2)` are distinct inverse classes of
> discriminant `−4(72k²−6)`. When q=12k²−1 is prime, this is exactly the
> source 6q on the square-sigma Berggren lift, whose aliquot edge is
> `(72k²−6) -> (72k²+6)` with total `(12k)²`.

This joins the previously separate Berggren/aliquot and idoneal lanes with
an explicit map and a uniform obstruction. It requires no completeness
theorem for the idoneal list and no conjecture about the supply of primes q.
See [proofs, P4](PROOFS.md#p4-an-explicit-non-idoneal-form-on-the-entire-lifted-ray).
The proof extends to the parent family8t²−6: for positive integer t it is
idoneal exactly at t=1. The lifted Berggren ray is the subsequence t=3k.

A second bridge gives a geometric interpretation of the common total144:

| Source | σ(N) | index `[SL₂(Z):Γ₀(N)]` | genus X₀(N) | class group of discriminant −4N |
|---:|---:|---:|---:|---|
|66=2·3·11|144|144|9|C4 × C2|
|70=2·5·7|144|144|9|C2 × C2|

Here σ(N) counts index-N sublattices of Z², and the modular index counts
those with cyclic quotient. They agree for squarefree N. For arbitrary N
the exact missing information is lattice content:

    σ(N) = sum_{d² | N} ψ(N/d²),    ψ(N)=N product_{p|N}(1+1/p).

At a balanced mixed source p²qr, the content layers are exactly1 and p.
This identifies precisely what a cyclic-quotient model omits when a prime
exponent is promoted. It also explains why genus alone cannot recover
class inversion:66 and70 agree on the two displayed modular statistics,
but the genus-character channel of the former loses one bit.

The incoming **Session84** packet already proves the complete
square-preserving prime-step law, the balanced square-sigma parametrization,
the eight-vertex promotion cube at66, and the full guarded Berggren maps.
Those results are inherited, not claimed as new here. In particular:

    σ(N) square and p prime:
    σ(pN) square iff p=3 and 3 does not divide N.

Its all-exponent necessity uses the classical Ljunggren square-repunit
theorem. This packet independently checks those conclusions and extends
the comparison to the promotion cube at70, lattice content, modular genus,
and the uniform non-idoneal form. [Incoming integration](INHERITANCE.md).

The requested history is organized in [the research atlas](ATLAS.md).
The pinned occurrence inventory searches every tracked `.md`, `.tex`, and
`.org` file outside agent/dot directories and this packet, recording exact
paths, blob identities and matching line numbers. It includes historical
and copied evidence but distinguishes these from reviewed truth sources.
It is **search coverage, not a claim that every occurrence was read or is
mathematically relevant**. Deleted historical versions and binary PDFs are
outside that inventory. [Search scope and counts](inventory_summary.json).

The most useful inherited maps are retained with their limits:

- Fifteen marked K6 matchings ↔ A5 involutions ↔ icosahedral half-turn axes;
  the marked point is essential to the first displayed construction.
- Fifteen matchings ↔ genus-two Lagrangian kernels; beyond six labels the
  missing Hamming-code structures matter.
- K5's two-subsets with disjointness form the Petersen graph; a bipartite
  sign gauge connects sums and differences on trees and K3,3, but cannot
  extend over Petersen's odd cycle.
- The Monster/Baby exponent pattern has the literal height cycle
  `4→2→5→3→7→4`, with13 feeding into2. Dividing the normalization by2
  changes it to `2↔4`.
- The fifteen supersingular primes use a Frobenius/level predicate. The
  connection with Monster/Baby Monster is stronger than cardinality, while
  the five known Fermat primes are not a proved complete list.

Reproduce from the repository root:

```sh
python3 -B 05-knowledge/results/aliquot_66_five_fifteen_20261028/check.py
python3 -B -O 05-knowledge/results/aliquot_66_five_fifteen_20261028/check.py
python3 -B 05-knowledge/results/aliquot_66_five_fifteen_20261028/inventory.py --revision f19c0ededc61205810b94384b51e84743c674ad1
```

The two checker modes must agree with [receipts.json](receipts.json).
[Audit and scope](AUDIT.md) distinguishes finite controls, proof dependencies,
the rejected shortcuts, and the remaining research questions.
