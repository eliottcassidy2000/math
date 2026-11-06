# Independent balanced n14 audit

PASS, 2026-11-06. Reconstructed code292187853813044375413825338 independently, verified strong connectivity by reachability from each of14 vertices, and counted directed Hamilton paths with a recursive suffix recurrence. The count is11390625=15^6. All14 start counts agree with the main verifier, whose independent forward/end and start subset recurrences and arc-incidence normalization are correct.

Replayed both checkers from outside this directory; outputs byte-match verify.out and independent.out. The discovery source is byte-identical to the frozen audited n13 oracle source and its SHA256 matches DEPENDENCIES.json. The saved target line records restart1525,step42 and the verified code. No discovery replay is needed to certify this explicit witness.

The dominance-join inference m(15^k)<=14 floor(k/6)+5(k mod6) is valid, hence the inherited subadditivity limit has upper bound7/3. The proof clearly avoids claiming m(15^6)=14, global optimality, exhaustive order14 classification, search convergence or a positive linear-rate lower bound. It preserves the distinction between single-arc search and native path-complement adjacency. The prior n13 witness remains valid; this is an improved construction, not a correction. No main logarithmic preprint theorem is needed.

This review covers the complete written proof, main verifier and standalone bundled independent checker. Ready for frozen packaging.
