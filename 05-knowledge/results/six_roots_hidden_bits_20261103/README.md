# Six roots, hidden bits, and the price of retaining a port

2026-11-03. **PROVED** exact rooted-ear formula and port-memory lower bound
for E_n; elementary deleted-tree boundary and golden-unit coordinates;
**FINITE-EXACT** explicit source and finite-field controls; **OPEN** universal
odd Hamilton-path spectrum completeness and marked ternary reentry.

The numerical survey led to a concrete improvement and an all-size limit:

- The root65 can use the strong source E8. Its first ear can give195 or197,
  whereas the earlier seven-vertex witness's least first-ear count was203.
  This changes the chosen root witness, not the vertices of a supplied source.
- On E_n, n>=4, all n marked vertices have distinct responses to the same
  ear attachment. Preserving that operation needs at least ceil(log2 n)
  bits for the mark even when the base tournament and its order are known.
- The65 actual Hamilton paths of E8 admit a projective F64 chart. It carries
  exact63- and65-cycles of path addresses. These are not vertex symmetries:
  E8 is rigid, and its marked responses remain distinct.
- The six roots have an exact origin as the boundary of two adjacent
  deleted nodes in a ternary tree. Their resemblance to square-sum
  exceptional orders is an arithmetic comparison, not a transferred
  Hamiltonicity proof.

[PORT_MEMORY](PORT_MEMORY.md) contains the proofs and the actual source
comparison. [BOUNDARIES](BOUNDARIES.md) treats exceptional islands, the
Mersenne7→21→63 identity, and the golden-unit sign/direction carrier.
[NUMBER_ATLAS](NUMBER_ATLAS.md) connects older work across arithmetic,
geometric encodings, tournaments, Collatz, and Zeckendorf.
[HANDOFF](HANDOFF.md) gives the research board, hostile tests and next targets.

The discovery index scans18,234 tracked prose files at a pinned revision;
it is not a claim to have mathematically reviewed every literal occurrence.
The same-agent audit uses a second counting algorithm and graph search,
and all failure checks remain active under Python optimization.

## Reproduce

From the repository root:

```text
python3 05-knowledge/results/six_roots_hidden_bits_20261103/verify.py
python3 05-knowledge/results/six_roots_hidden_bits_20261103/audit.py
python3 05-knowledge/results/six_roots_hidden_bits_20261103/projective.py
python3 05-knowledge/results/six_roots_hidden_bits_20261103/probe67.py --distance 2 --max-order 10
python3 05-knowledge/results/six_roots_hidden_bits_20261103/probe67.py --distance 3 --max-order 10
python3 05-knowledge/results/six_roots_hidden_bits_20261103/survey.py
```

The last command is the larger lexical census, using the revision fixed in
the script. It excludes this packet. `occurrences.json.gz` retains the raw
discovery inventory, including historical and noisy hits. `source_pins.json`
pins the selected mathematical routes separately.
