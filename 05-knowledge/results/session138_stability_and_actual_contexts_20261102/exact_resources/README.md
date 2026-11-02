# Session138 exact-spectrum resource lane

## Results

* `PROOF.md`: one literal exterior arc reversal yields universal actual-context
  responses `6E(U)-H(U)` and `6E(U)+H(U)`, with full polynomial operators.
  This supplies unbounded strong fixed-gap pairs, a source-certified inverse,
  quantitative input-height sparsity, and the minimal four-exterior-vertex
  isolator classification. No unrestricted successor or spectrum theorem.
* `RESOURCE_BUDGET.md`: a bounded actual `(h,L)` source bank fails even
  profinite universality, using its finite affine-dyadic normal form and the
  modulus `2^r-1`; quantitative height/interior lower resource bounds follow.

## Reproduction

```
python3 -B verify.py > normal.json
python3 -B -O verify.py > optimized.json
cmp normal.json optimized.json
python3 -B budget.py > budget.normal.json
python3 -B -O budget.py > budget.optimized.json
cmp budget.normal.json budget.optimized.json
```

Both programs are standalone Python standard library. They do not import any
repository Python module, graph atlas, discovery probe, or remote service.
The all-order arguments are in the manuscripts; the receipts identify their
complete finite universes and controls. Actual adjacency certificates are
retained in `normal.json`. The inherited mathematical dependencies and their
read-only source hashes are in `SOURCE_PINS.json`.

`probe_flip.py`, `search_contexts.py`, and `profile_probe.py` are exploratory
scratch, excluded from `PUBLIC_FILES.txt` and from the theorem's reproduction
contract. No Git mutation was performed by this lane.
