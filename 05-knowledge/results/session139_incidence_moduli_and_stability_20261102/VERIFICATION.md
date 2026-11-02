# Session139 checkpoint-one verification

**PASS in the stated scopes.** This release inherits the accepted Session138
inputs; it does not claim to rerun every earlier theorem or to close 10..747.
The source proof's drafting-status sentence is retained verbatim; its final
independent acceptance is recorded here and in [frozen scope](CHECKPOINT1.md).

## Sharp coefficient for every n>=748

- [Proof](stability/PROOF_N748.md): retained binary baseline, square completion,
  and one simultaneous heavy/light wedge budget.
- [Independent audit](independent_threshold/AUDIT.md): checks quantifiers,
  nonempty retained matching, both parities, strict inequalities, nonnegative
  boundaries and the inherited equality classification.
- Producer certificate: 20 explicit exact integer/Fraction checks.
- Independently authored certificate: 22 named checks, importing no producer
  implementation. Both normal and optimized Python outputs are byte-identical.
- Threshold 747 is deliberately rejected by the endpoint certificate; this
  is a control of the method, not a counterexample at order 747. Setting the
  producer's h to zero is also rejected in both execution modes.

Uniformity comes from the proof's monotonicity arguments, not finite sampling.
The reviewed producer proof and executable match the independent audit's
SHA-256 pins. No invalid termwise Pfaffian unit-coefficient lemma is used.

## Weighted spectral loss and historical recovery

The [independent bridge audit](independent/AUDIT.md) reconstructs the full
bivariate rational conjugacy, determinant polynomial and Hamiltonian sum with
standard-library polynomial arithmetic, not SymPy or producer imports.
It verifies the two-parameter `a^2*b^2/54` H difference and the exact rational
witness, including the `M23` repair. Conjugacy fixing the all-one vector proves
all rooted scalar walk responses coincide, not just a sampled list.

The same frozen audit contains auxiliary exhaustive binary F2 controls on all
labeled tournaments at n=4,5,6 and sharp ordinal witnesses through 20. This
retained evidence is **not promotion of the separate pending all-moduli or
root-cocycle continuation**. False H coefficient 1/55 and false Paley7 count
35 are rejected with and without optimization.

The release's [historical checks](riemann_five_eleven/verify_history.py) use
SymPy only for exact polynomial algebra and an independently written
Held–Karp dynamic program for the rational witness's Hamiltonian sum. They
also check the corrected transfer/tetranacci/pentanacci discriminants,
Bernoulli denominator repetition, and Paley11 character partial sums with
natural/optimal ranking counts 33/35. These are bounded arithmetic checks,
not a literature-priority certificate or an RH/GRH result. Deliberately false
tetranacci, Paley optimum and incidence values are rejected under both modes.
The original exploratory `probe.py` and its printed output are retained as
non-load-bearing producer evidence; they do not themselves assert checks.

[Source pins](SOURCE_PINS.json) verify 20 historical Git blob references and
six inherited proof inputs against base commit
`291b4d9757842c02aa451beb65acfade2848609b`. Historical prose is inherited only
with the corrections and type boundaries in [recovery](riemann_five_eleven/HISTORICAL_RECOVERY.md).
No older canon or central mistakes-ledger file is edited by this checkpoint.

## Fresh relocated replay

Run from this packet directory:

```bash
sha256sum -c SHA256SUMS
python3 replay.py --output /tmp/session139-checkpoint1-replay.json
```

`replay.py` copies only the four named kernels into a fresh temporary directory,
then performs **22 executions: eight successful scientific runs and fourteen
expected hostile rejections**, including normal and `python3 -O` modes.
The four positive output pairs are byte-identical; frozen producer and audit
outputs must match exactly. The report records source hashes and stdout hashes.
Python 3 and SymPy (version recorded in the history receipt) are required;
the two stability certificates and symbolic bridge audit use only the standard
library. [Fresh assembled receipt](REPLAY_CHECKPOINT1.json).

The older `independent/SHA256SUMS` is an unchanged archival source receipt with
original absolute scratch paths. The root `SHA256SUMS` is the portable release
manifest; it is the one used by the command above.

## Documentation and publication boundary

The native documentation checker cannot pass directly in this sparse checkout
because its baseline documents are absent. The unmodified same-version checker
passes with read-only navigation overlays on a clean, materialized Session137
ancestor; all 245 intervening Session138 paths are explicitly accounted for.
No baseline files are changed or hydrated. The exact result is
[recorded here](DOC_OVERLAY_VERIFICATION.json): startup 179966/180000 bytes,
frontier 450/450 lines, maintained warning sentinels and router smoke PASS.
New packet links, explicit source pins, frozen-copy integrity, conflict markers,
and the owned-path staging boundary are checked separately.

The local release commit and remote transport/publication receipt are distinct.
A prepared local commit is not evidence of a remote update; the publishing
agent must verify the actual remote commit and every owned blob before reporting
publication. Later continuations must retain this checkpoint's original scope.
