# Session26: randomized bridges and exact observation

This independently audited research session connects continuation observers, attachment counts, feedback strata, exact metagraph wiring, and the information retained by a distribution. These are scoped elementary results with exact computational controls, not literature-priority claims or solutions of LRC(14).

Start with the [final synthesis](SYNTHESIS.md) for proved results, scope limits and next questions; the [exploration record](EXPLORATION.md) tracks the hypothesis/test cycles. [Independent final audit](session26-synthesis-review/REVIEW.md). Sampling provenance: [root manifest](ROOT_SAMPLE.json), [103-path summary](SAMPLING_SUMMARY.json), and the lane manifests retained below.

Run `python3 session26-synthesis-review/check_portable.py` from this packet directory (or pass the packet directory as its argument). This checks final document hashes, all 24 synthesis/exploration links, sample-manifest counts, exact scalar identities and the boundary table. It explicitly does not reverify the omitted archived source copies. The original check.py and original outputs in that review directory preserve the earlier workspace audit and require those omitted source copies; use check_portable.py for distributed reproduction.

- [One-apex continuation observer](session26-tournament/PROOF.md): exact minimal continuation-state counts modulo every integer, an integral observer lattice, and the shared last-bit CRT constraint. [Independent audit](session26-tournament-review/REVIEW.md).
- [Runner clearance height curvature](session26-runners/RESULT.md): signed height atoms recover the full clearance distribution and maximum; the even two-runner family has unbounded normalized-moment recurrence order. [Independent audit](session26-runners-review/REVIEW.md).
- [Constrained arithmetic readers](session26-wildcard/RESULT.md): exact Moore-machine size from reachable and unobservable subgroups, with quadratic digit ideals and Smith invariants. [Independent audit](session26-wildcard-review/REVIEW.md).

## Reproduction

Keep the sibling directories intact. From each producer directory, run its Python checker: tournament `python3 verify.py` and `python3 verify_crt.py`; runners `python3 probe.py` and `python3 family_check.py`; wildcard `python3 probe.py`. Python optimized-mode controls are also included where supplied.

Independent checks: from tournament-review run `python3 independent_verify.py`; from runners-review run `node independent.js` and `python3 compare.py`; from wildcard-review run `node verify.js`. Saved outputs and source hashes preserve the precise audited inputs. Finite runs corroborate, but do not replace, the all-parameter proofs.

## Provenance and packaging

The random manifests and exploration annotations record the sampled repository paths and original commit. Copied archival sample files and tournament random-sources.txt are intentionally omitted; recover those research inputs from the recorded repository commits. Their omission does not affect the self-contained proof or verification scripts. Bare repository paths in research prose refer to repository-root source citations, not paths relative to this packet.

The producer observer text's original 'independent audit pending' label is historical: the included source-bound audit supplies its completed verdict. Original audit paths identify the research workspace; the corresponding sibling packet files carry the published bytes. Producer and reviewer checksum files preserve original audit provenance rather than a promise that every sampled research input is distributed here.

The observer independent checker's final proof-hash lookup is packaged using a portable sibling path; its mathematical computation is unchanged. No theorem IDs are reserved or canon statuses changed by this packet. The initial checkpoint covered the three observer packets above; the independently reviewed continuation below supplies the cut-partition and wiring results.

## One-apex continuation: cuts, spectral laws, and wiring

- [Cut and Ferrers formulas](session26-tournament/CUT_ADDENDUM.md): exact weighted inversions, extremal boundaries, and exclusion of Hamiltonian count 13 throughout the one-apex stratum. [Concrete H=13 witness outside that stratum](session26-tournament/H13_WITNESS.md).
- [Apex classification and spectral law](session26-tournament/APEX_CLASSIFICATION_AND_SPECTRAL_LAW.md): exact isomorphism counts and the strong-word law W(X)=(log 2)X+O(sqrt(X)), giving a positive lower bound on missing odd values. The ordinary/converse-merged asymptotics count strong cores, not all padded tournaments. [Source-bound independent audit](session26-cut-review/REVIEW.md).
- [Blue incidences and forced wiring](session26-blue/RESULT.md): the feedback-number inequality for every metagraph edge, exact one-apex masses and blue incidences, and the resulting stratum constraints. [Independent audit and enumeration](session26-blue-review/AUDIT.md).

Reproduce the producer finite controls with `node verify_cut.js`, `node verify_asymptotic.js`, and `node verify_wiring.js` in session26-tournament. Independent blue controls are `python3 check.py` and `python3 check_n8.py` in session26-blue-review. The cut audit records separate independent computations, but no independent cut-verifier source is distributed; do not confuse the producer scripts with those independent tests. All-order conclusions rest on the reviewed proofs.

## Attachment transforms, fixed feedback strata, and the last internal edges

- [Attachment transforms](session26-attachments/proof.md): exact attachment means recover ordered path-cover profiles, with precise conditioning and moment obstructions. [Independent audit](session26-attachments-review/review.md).
- [Fixed feedback strata](session26-fixed-apex/proof.md): asymptotic uniqueness, rigidity, class counts and forced wiring for each fixed positive feedback bound. [Independent audit](session26-fixed-apex-review/AUDIT.md).
- [Last one-apex internal wiring](session26-last-wiring/RESULT.md): the sharp order-eight boundary, exact internal edges and endpoint classes. [Independent exhaustive audit](session26-last-wiring-review/AUDIT.md).

Producer reproduction: `python3 audit.py` in session26-attachments, `python3 check.py` in session26-fixed-apex, and `python3 catalogue.py` in session26-last-wiring. Independent attachment and fixed-apex checks run with `python3 independent.py` in their review directories. For last-wiring review, a C++ compiler is required: from its review directory run `g++ -O3 enumerate.cpp -o /tmp/last_apex_enumerate`, `/tmp/last_apex_enumerate > selected.txt`, then `python3 audit.py`. The enumerator exhausts all fixed-path masks through order eight, not merely the producer's one-apex generation.

Two packaged reviewer scripts replace original workspace-absolute paths by sibling-relative paths; computations are unchanged. Original review checksums preserve audit provenance and therefore identify the pre-portability script bytes. The attachment review explicitly rebinds its final status and Markdown-only proof updates. Earlier independent output intentionally records the original proof hash; portable-replay.out and portable-replay-optimized.out record the final proof hash and unchanged successful calculations. The last-wiring audit's historical absolute-path note is superseded by this portable packaging.

## Final counting instrument and repaired historical claims

- [Fixed-core attachment-word counter](session26-counter/proof.md): a supplied finite-core presentation has an exact nonnegative weighted-automaton count, retaining empty slots and entrance/exit ports. [Independent source-bound audit](session26-counter-review/REVIEW.md). Run `python3 audit.py` in session26-counter and `python3 check.py` in session26-counter-review; the latter's source-directory path is packaged portably. Original review checksums retain the original absolute-path script provenance. The packaged portable check.py has SHA256 `fa6363d31d75384731226f920f6a1263b479077059aec661ee2083e28cb6404c`, not the original hash in that historical checksum file.
- [Narrow correction packet](session26-corrections/README.md): reciprocal truncation is not reciprocal oddness; numerical Coxeter coincidences are not bridge proofs; the bounded-denominator LRC census generalization remains retracted. The lcm family has exact maximum 1/8 while defeating every prescribed finite rational-denominator census. This is not a general LRC proof. Run `python3 check_math.py` in session26-corrections. The four reviewed source repairs are applied at their original repository paths and logged as MISTAKE-656; manifest and patch retain the preimages. Copied original/proposed source trees and the auxiliary banner-only patch are omitted; FREEZE-SHA256.txt preserves their original audit provenance.

