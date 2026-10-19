# Canonical mixed (2,4,1) maximum-span exclusion

**Transport setup:** The oversized gzip certificate is stored losslessly as raw-byte parts. Before running the audit or replay commands below, run `python3 reconstruct_transport.py` from this directory. See [the transport instructions](TRANSPORT_README.md) for original hashes, part verification, and safe reconstruction.

This theorem excludes morphology (2,4,1), up to duality, only for a root selected by maximum majority-order span and then maximum common-incomparable count. It is not an any-root exclusion.

All 91,604 degree-valid cores are covered: 21,501 analytic exclusions and 70,103 residual cores. The canonical screen removes 274,819 orders and wholly removes 15,711 previously nonempty residual cores. Exact positive certificates handle 415,884 remaining orders. The known six-thin theorem removes 12,716 short-span orders.

Run `python -O audit.py` for provenance/coverage or `python -O replay_complete.py` for a fresh complete standard-library mathematical replay. Search uses `run.py` with NumPy/SciPy. Every record was independently replayed before commit; no cone obstruction occurred. No support-minimality, full seven-thin or novelty claim is made.
