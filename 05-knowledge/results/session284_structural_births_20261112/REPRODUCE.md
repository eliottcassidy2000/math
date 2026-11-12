# Reproduce Session284

From this directory, using standard Python3 only:

```bash
python3 architecture/verify.py
python3 births/verify.py
python3 extension/verify.py
python3 wildcard/verify.py
python3 bridge_verify.py
python3 audit/independent.py
```

Repeat each command with `python3 -O`; matching receipts are supplied or recorded by the independent review. Implementations raise explicit failures and do not depend on assertions removed by optimization. The wildcard literal9! permutation check is intentionally slower than the symbolic grammar.

The tests are hostile/positive controls, not proofs by extrapolation. Scope:

- Architecture: layer recurrence through5; all parent orbit controls through4; W5 prime children.
- Births: all1023 binary words through length9 for modules; all255 through7 for literal induced prime inventories; formula controls at arbitrarily parameterized larger indices.
- Extension: exact marked-template and Hamiltonian semantics controls; labelled batch probability laws for fixed parents, not a pooled isomorphism measure; explicit Gray-code interface witnesses.
- Wildcard: literal group and endpoint controls, exact host grammar and embedding receipts; inherited path-profile formulas are dependencies, not freshly proved here.
- Bridge: word coefficient controls, with membership equivalence established by the two family theorems.
- Independent: separate adjacency-matrix implementations of extension moments, root-forgetting identity and pointwise stabilizers; review specifies precise scope.

`audit/SHA256SUMS` binds the reviewed source snapshot. Verify from this directory with `sha256sum -c audit/SHA256SUMS` when available. Release receipt records any repository-wide documentation warnings separately from these mathematics checks.
