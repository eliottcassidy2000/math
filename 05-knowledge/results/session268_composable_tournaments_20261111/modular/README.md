# Count tournament families through their compositional responses

**PROVED:** exact Burnside pushforward of ordinary tournament-class inventories by full ordered path-cover response, with a reciprocal-automorphism mass sidecar. Together with the spectrum lane's tomography, this response is exactly the observational equivalence for Hamiltonian counts under uniform substitution contexts.

**FINITE-EXACT:** all tournaments whose canonical prime quotients are C3, with arbitrary ordered strong-component composition, through order 16. This is an infinite subclass, not all tournaments.

| Order | Ordinary classes | Full response types | Distinct H | Labelled tournaments |
|---:|---:|---:|---:|---:|
| 8 | 182 | 63 | 50 | 4,399,360 |
| 10 | 1,558 | 350 | 293 | 3,049,625,600 |
| 12 | 14,225 | 2,076 | 1,757 | 3,306,847,667,200 |
| 14 | 135,960 | 12,904 | 11,248 | 5,167,669,361,254,400 |
| 16 | 1,343,039 | 83,073 | 74,149 | 10,997,137,501,622,272,000 |

At order 16 this is a 16.17-fold reduction from ordinary classes to the exact composition states; the largest single response bucket holds252 classes. The full inventory through 16 ran in 9.728 seconds on this environment, including independent scalar checks and small permutation controls, before JSON serialization. This is a workload measurement, not a comparison against nauty or an optimal implementation. The run used69,883 distinct C3 response contractions and103,793 ordinal-sum response contractions.

The important arithmetic is the diagonal Burnside term: a response bucket containing m actual child classes contributes (m^3+2m)/3 ordinary C3 substitutions, while its mass u contributes u^3/3. Treating a response as if it were an isomorphism type is wrong. The inherited two-child 81/27 example gives four ordinary classes in one output response.

Full q, not scalar H, is necessary. At order 6 the generated family already contains two profiles with H=9 that a C3[child,1,1] context sends to177 and153. Pure substitution preserves full q equivalence; a nonuniform cross-module seam need not, as the interacting-automata lane demonstrates.

## Files

- `PROOF.md`: identities, aggregation theorem, exact family grammar, contextual equivalence, complexity and losses.
- `census.py`: standard-library response-inventory compiler. No graph-isomorphism work in its production census.
- `receipt.json`: readable count summaries and hash; `full_receipt.json.gz.b64`: lossless full H histograms through 16 and full response buckets through 12.
- `independent_controls.py`, `controls.json`: actual labelled orientation universes through 5, canonical-form class grouping, and an explicit scalar-H loss witness.

## Reproduce

    python3 -B 05-knowledge/results/session268_composable_tournaments_20261111/modular/census.py --max-n 16 --out /tmp/modular-receipt.json
    python3 -B 05-knowledge/results/session268_composable_tournaments_20261111/modular/independent_controls.py

Normal and Python -O runs agreed on every output field except elapsed time. All output fields except elapsed time are deterministic. `--bucket-max-n` defaults to12, so extending the run bound does not save every larger full profile. Existing canonical modular decomposition, the general species equation and historical prime cycle-index census are credited in the proof rather than claimed anew. Chan–Pak and a complete general prime classification remain OPEN.
