# Reproduce the new constructions

From the repository root, using Python3 and its standard library:

    python3 -B 05-knowledge/results/session263_spin_flip_20261111/transfer/controls.py
    python3 -B 05-knowledge/results/session263_spin_flip_20261111/transfer/one_query_controls.py
    python3 -B 05-knowledge/results/session263_spin_flip_20261111/selectors/selector_probe.py
    python3 -B 05-knowledge/results/session263_spin_flip_20261111/converse/construct.py
    python3 -B 05-knowledge/results/session263_spin_flip_20261111/fourier/filters.py

Each writes its adjacent JSON receipt. The controls concern the NEW operations:

- Transfer: exact cone bounds, inverse decoding, state distinctness and
  tomography for all2,047 binary words through depth10; nine literal composites
  and four hostile rejects. Inherited transfer matrices are inputs.
- One-query: exact radix decoding and all-count distinctness for those2,047
  states; four actual small graph encodings and a negative-remainder hostile.
- Selectors: all63 binary apex words of length0–5; all225 two-apex return
  joins with word lengths0–3; three nested nine-vertex gates. Direct cycle
  lists and independent subset path counts are compared for the new graphs.
- Converse: exact six-prime valuation certificate, finite literal examples,
  caps, averaging obstruction controls, and component-level repeated bilateral
  filtering. The large filtered graphs are NOT exhaustively recounted.
- Fourier: independent subset-DP and direct permutation enumeration of every
  induced subset for C3, F8 and the partial-filter hostile; exact length product,
  parity-jet, primitive-only failure, finite closure heads and10010 architecture.

The proofs state the all-size results. The finite ranges are not their proofs
unless the assertion itself concerns an explicit finite block/certificate.
No float computation, random sampling result, inherited audit or external
snippet is a load-bearing all-size premise.

`fourier/discovery.py` optionally repeats only the exploratory10010-filter
search. It is not needed to verify or use the explicit constructor.
`FINAL_CHECKS.json` records the final replay and documentation checks.
`MANIFEST.json` records SHA256 hashes of this packet, excluding itself and
publication receipts whose contents depend on the final commit.
