# Session26: randomized bridges and exact observation

Three independently audited packets distinguish arithmetic reachability, continuation observability, and the information retained by a distribution. These are scoped elementary results with exact computational controls, not literature-priority claims or solutions of LRC(14).

- [One-apex continuation observer](session26-tournament/PROOF.md): exact minimal continuation-state counts modulo every integer, an integral observer lattice, and the shared last-bit CRT constraint. [Independent audit](session26-tournament-review/REVIEW.md).
- [Runner clearance height curvature](session26-runners/RESULT.md): signed height atoms recover the full clearance distribution and maximum; the even two-runner family has unbounded normalized-moment recurrence order. [Independent audit](session26-runners-review/REVIEW.md).
- [Constrained arithmetic readers](session26-wildcard/RESULT.md): exact Moore-machine size from reachable and unobservable subgroups, with quadratic digit ideals and Smith invariants. [Independent audit](session26-wildcard-review/REVIEW.md).

## Reproduction

Keep the sibling directories intact. From each producer directory, run its Python checker: tournament `python3 verify.py` and `python3 verify_crt.py`; runners `python3 probe.py` and `python3 family_check.py`; wildcard `python3 probe.py`. Python optimized-mode controls are also included where supplied.

Independent checks: from tournament-review run `python3 independent_verify.py`; from runners-review run `node independent.js` and `python3 compare.py`; from wildcard-review run `node verify.js`. Saved outputs and source hashes preserve the precise audited inputs. Finite runs corroborate, but do not replace, the all-parameter proofs.

## Provenance and packaging

The random manifests and exploration annotations record the sampled repository paths and original commit. Copied archival sample files and tournament random-sources.txt are intentionally omitted; recover those research inputs from the recorded repository commits. Their omission does not affect the self-contained proof or verification scripts. Bare repository paths in research prose refer to repository-root source citations, not paths relative to this packet.

The producer observer text's original 'independent audit pending' label is historical: the included source-bound audit supplies its completed verdict. Original audit paths identify the research workspace; the corresponding sibling packet files carry the published bytes. Producer and reviewer checksum files preserve original audit provenance rather than a promise that every sampled research input is distributed here.

The observer independent checker's final proof-hash lookup is packaged using a portable sibling path; its mathematical computation is unchanged. No theorem IDs are reserved or canon statuses changed by this packet. Cut-partition asymptotics and colored wiring are separate work, not results claimed by this first checkpoint.
