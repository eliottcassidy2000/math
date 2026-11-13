# Witt-inspired graded carriers for tournament seams

Source-informed elementary bridge, independently reviewed (see ../audit/REVIEW.md).

- [Proof and scope](PROOF.md): marked necklace carry and parity obstruction;
  square-zero determinant cycle gas; endpoint resolvent and seam elimination.
- [Sources and prior repository mechanisms](SOURCES.md).
- [Verifier](verify.py) and [receipt](receipt.json).

Reproduce from repository root:

```sh
python 05-knowledge/results/session291_cluster_ports_20261113/witt/verify.py
python -O 05-knowledge/results/session291_cluster_ports_20261113/witt/verify.py
```

Both modes gave identical receipts. Tests compare a permutation-determinant
square-zero quotient against independently counted subset Hamiltonian paths
for all 1,099 labelled tournaments through order five and 80 seeded examples
each at orders six and seven. Explicit full-path enumeration checks endpoint
marginals. Five small marked necklace orbit censuses check the prime-power
formula, plus scalar odd-prime parity checks. These support implementation,
not a claim of Chan–Pak coverage. The Schur-complement identity is proved
algebraically, not separately exercised by the current script.

Main actionable distinction: norm/necklace arithmetic alone has no physical
realization theorem; the graded resolvent does retain exact physical endpoint
and no-repeated-vertex information, but its exponentially growing carrier has
not been uniformly compressed. This is a precise target for stable-port
quotients, rather than evidence that such a quotient already exists.
