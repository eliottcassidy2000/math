# Session70 graph and labeling bridge

[REPORT.md](REPORT.md) contains the proofs, recovery/correction ledger,
connection contracts, finite-test universes, and primary-source scope.

- **PROVED:** exact perfect-matching minimum for unrestricted named
  linear-forest attachments; `(2r−1)!!` states.
- **PROVED:** inherited factorial graceful fibers have exactly `k!`
  Boolean overwrite-legality classes, not `k!2^k`; actual signs are invisible
  for that carefully specified language.
- **PROVED:** Petersen-to-K5 port transitions; all six Petersen2-factors
  are two5-cycles. The K3,3 subdivision requires its suppressed-edge marks.
- **PROVED:** K3,3 is graceful; K5 is not. Neither is a tree, and both
  are nonplanar and Hamiltonian.
- **RECOVERED / REPLAYED:** actual square-sum graph Q_n, endpoint-aware
  path-cover operations, and the m12 fixed-color graceful counterexample.

No unrestricted tournament minimum, graceful-tree existence theorem, or
minimality for the restricted square-sum attachment language is asserted.

Reproduction requires only Python3:

```sh
python3 verify_graph_alignment.py > /tmp/graph70.json
python3 -O verify_graph_alignment.py > /tmp/graph70-opt.json
cmp receipt.json /tmp/graph70.json
cmp /tmp/graph70.json /tmp/graph70-opt.json
python3 verify_graph_alignment.py --hostile  # intentional nonzero exit
```

`VALIDATION.json` records isolated temporary-directory normal/optimized
replays and intentional hostile rejection. `SOURCE_PROVENANCE.json` pins
the read-only inherited source snapshots to base commit
`fb872272125510c317b8114dd3f697b433c595f9`. Snapshots preserve provenance,
not precedence over subsequent corrections.
