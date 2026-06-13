---
id: HYP-1969
status: OPEN
source: codex-2026-06-01-S505
related:
  - HYP-1812
  - HYP-1942
  - HYP-1950
  - HYP-1960
  - HYP-1961
---

# HYP-1969: LRC proof routes split into finite-sieve, kernel, zonotope, and endpoint-pressure currencies

## Statement

Current Lonely Runner proof progress is organized by four proof currencies:

1. **Finite product-sieve currency.**  Reduce any minimal counterexample to a
   bounded product, then force enough prime divisors of that product to exceed
   the bound.
2. **Spectral/kernel currency.**  Use nonnegative Fourier, Fejer, or
   Riesz-product tests on safe/unsafe indicator functions to improve global
   loneliness bounds or certify mixed threshold intersections.
3. **Zonotope currency.**  Reinterpret LRC variants as covering-radius
   questions and build finite certificates or finite-checking bounds.
4. **Endpoint-pressure currency.**  Treat a counterexample as an all-protected
   endpoint incidence system and search for private rows, pressure leaves, or
   labelled SCCs.

The 2025-2026 fixed-case proofs show that the first currency is currently the
winning frontier technology.  The repo's best chance of turning this into a
general or semi-general proof program is not to imitate the product sieve only
as a black-box computation, but to translate the sieve's improper tuples into
endpoint row-cover cores and pressure-DAG peel certificates.

## Evidence

`04-computation/lrc_proof_landscape_s505.py` encodes the literature/proof routes
as vertices and compares them with two Tournament Analysis gauges.

The declared Tournament Analysis object is:

```text
pairwise observable = feature-wise comparison between two proof routes
actionability gauge = majority(repo_signal, boundary_overlap, proof_engine,
                               reproducible, frontier_signal)
credence gauge      = majority(trust, direct_lrc, scope_strength,
                               reproducible, generality)
tie path            = chronological route order
```

Stored output:

```text
05-knowledge/results/lrc_proof_landscape_s505.out
```

Key fingerprints:

```text
Actionability tournament:
  score histogram {0:1, 1:1, ..., 12:1}
  SCC sizes [1,1,1,1,1,1,1,1,1,1,1,1,1]
  directed 3-cycles 0
  Hamiltonian-path count 1
  top route ST13, then Repo26, T9T10, MSS25, R8, R9

Credence tournament:
  SCC sizes [8,1,1,1,1,1]
  directed 3-cycles 9
  Hamiltonian-path count 99
  top route T9T10, then ST13, MSS25, BS07, R9, R8

Gauge flips:
  29 of 78 route-pair edges flip between credence and actionability.
```

This split is meaningful.  Mature older proofs and peer-reviewed reductions win
credence edges, while the newest product-sieve and endpoint-pressure routes win
actionability edges.

External proof anchors used in the atlas:

- Barajas-Serra prove the seven-runner case via distance-graph/chromatic-number
  tools.
- Tao proves finite-checking and improves the asymptotic lower bound over the
  trivial union bound.
- Malikiosis-Santos-Schymura improve finite checking to a linearly-exponential
  bound using zonotopes.
- Rosenfeld, Trakulthongchai, and Sungkawichai-Trakulthongchai build the
  current fixed-case product-sieve line through thirteen total runners.
- Bedert's Riesz-product paper and Jensen's mixed-threshold paper give the
  strongest nearby spectral language for endpoint-kernel pressure.
- Shifted-LRC zonotope work is useful but must be quarantined: the 2026
  counterexamples show shifted analogies are not automatically valid for the
  original LRC.

## Predictions

1. A useful repo implementation of the product-sieve proof should expose the
   set `J(k,p)` or `I(k,p,l)` as an endpoint incidence object, not only as a
   residue tuple count.
2. Product-sieve eliminations should have a pressure interpretation: a tuple
   removed by a low-level ansatz should correspond to a pressure leaf or private
   endpoint row after quotienting.
3. Jensen's unequal-threshold two-function formula should detect the first
   nonuniform endpoint pressure before a full high-order kernel is needed.
4. Bedert-style Riesz products are the best candidate for HYP-1812's
   nonnegative endpoint-kernel certificate.
5. Zonotope certificates remain valuable as finite-checking and geometry
   witnesses, but shifted-LRC counterexamples warn against importing shifted
   conclusions without preserving the exact unshifted endpoint object.

## Next Tests

- Implement a small `J(k,p)` reproducer for `k=7,8` and record, for each
  improper tuple, the endpoint rows it tries to cover.
- Add a mixed-threshold safe-pair integral check to the endpoint-pressure
  scripts and compare its sign with pressure-DAG source/sink layers.
- Build Fejer/Riesz kernels on the exact closed safe arrays from
  `lonely_runner_endpoint_protection_s360.py`.
- Compare product-sieve surviving tuples against HYP-1961 pressure-DAG layer
  vectors: surviving tuples should be the first place where pressure leaves are
  hard to find.

## See Also

- `04-computation/lrc_proof_landscape_s505.py`
- `05-knowledge/results/lrc_proof_landscape_s505.out`
- `07-reflections/lrc-proof-route-landscape-s505.md`
- HYP-1812, HYP-1942, HYP-1950, HYP-1960, HYP-1961
