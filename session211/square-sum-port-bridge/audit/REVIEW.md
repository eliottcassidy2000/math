# Independent square-sum bridge review

**Verdict: PASS.** Reviewed proof SHA-256: `c9d7c4586d2ce61c3813aa74ee4b182d94ff072ed0e31a98133973f74c583199`.

## Approved mathematical scope

- Source-pinned distinction among square-sum vertex orders, square path counts, tournament Hamiltonian-path counts, and normalized tournament indices.
- Nonexistence of a universal tournament assignment satisfying `H(F(Q_n))=2p(Q_n)+1`, witnessed independently by `p(Q23)=3` and `p(Q25)=10` and the inherited exclusions `H=7,21`.
- The exact degree-two insertion/two-port spanning-path-forest bijection, with the no-predecessor-Hamilton-path hypothesis explicit.
- The scalar recurrence obstruction for input `(old path count,new degree)=(0,2)`, without excluding order-dependent or richer boundary-state recurrences.
- Exact finite path census `Q1,...,Q25`, with singleton convention stated; all14 onset paths and nine designated-port controls.
- Finite replay of the inherited Pythagorean construction: all11 primitive triples `1<=s<t<=60`, at each of orders `t^2-1,t^2`, giving22 explicit spanning-path checks. These are construction checks, not counts of all paths at those larger orders.

The packet does not prove the cited all-order square threshold, tournament-spectrum completeness, or literature novelty.

## Reproduction results

The two producer scripts were copied into an independent replay directory and run under normal Python and optimized `-O` mode. Both census outputs are byte-identical to frozen `RESULT.json`; both macro outputs are byte-identical to frozen `MACRO_RESULT.json`. All checks survive optimized mode.

An additional independently written script, `independent_square.py`, enumerates every simple path from every starting vertex, without memoization, leaf-root selection, or feasibility pruning, for every `Q_n` with `2<=n<=25`. Its only terminal reduction selects the path orientation with first endpoint smaller than the last. All counts and all14 actual onset path edge orders match the producer after reversal canonicalization. Its normal and optimized JSON outputs are identical.

All five historical source copies match both their pinned Git blob IDs and their SHA-256 values. Every17 payload entry in the producer's `SHA256SUMS` passes. No shared-checkout edit was used to establish these results.

## Final producer artifact hashes

- `PROOF.md`: `c9d7c4586d2ce61c3813aa74ee4b182d94ff072ed0e31a98133973f74c583199`
- `verify_square_bridges.py`: `570ba60e9c3fa5dd69a0ccc49d85d4aeccc73dc693a5f0793b797d40f7634097`
- `RESULT.json`: `8b711b85f47efaa483c2e678805fa3f00afa08ab3a98c45c4843836849f7295f`
- `verify_pythagorean_scaling.py`: `7f3094c0e82cf76243cb308892beff21bebd7475356ee449d1813b50abb53b44`
- `MACRO_RESULT.json`: `490c6033525766a2964716d0f78a28f34f5f17cb4e23d35597111acf0c863c19`
- `SOURCES.json`: `e801b9bee0f40409731b249d702dc99ea63e86b8fc5f2701c522df93c6177d6e`

The independent script and its normal/optimized outputs are approved for inclusion beside this report. Their numerical findings are controls supporting the scoped mathematical proof, not replacements for it.
