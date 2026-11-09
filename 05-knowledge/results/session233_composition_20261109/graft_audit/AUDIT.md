# Session233 coprime moving-core composition audit

Verdict: PASS for the stated universal consecutive-window theorem; no Chan–Pak completeness result. Audited source: `/workspace/scratch/07d5532f2fce/session233-graft`, against repository checkpoint `2ecb04e29` and the inherited Session231 full-interface results.

## Mathematical review

The bounded coin-box proof is valid for coprime integers a,b>=2, K>=b-1, L>=a-1. Fixing the unique residue of j modulo a, the quotient intervals meet because K>=b-1. Their first and last excluded members lie outside [C,M-C], C=(a-1)(b-1). Thus the inclusion holds uniformly, not merely for sampled rectangles. The explicit selector has legal bounds. Reflection proves the corresponding upper obstruction as well. For a=45,b=37 the conductor is 1584, and 1583 is unattainable.

The two independently checked finite flanked identities and the matched endpoint bits imply the all-order family statements by an adjacent-move induction. Composition of two pure-T changes multiplies them by the opposite block's h. Therefore the T difference is -90i+74j; after j=L-u the closure counts are exactly B+2(45i+37u). All other coordinates are unchanged. Different position pairs can collide; no unjustified isomorphism-distinctness claim is made.

The inherited coordinates agree with the available THM-4960 definition of full interface, in which T is the varying cut coordinate and h,q,A,G,S,E are fixed. Saturation persists since both q factors are divisible by their h factors. Common singleton contexts have h=1, so preserve pure-T difference and the scalar step; this does not confer strongness for arbitrary unprotected contexts.

The protected construction is strong: the first deletion singleton is selected and dominates every later deletion vertex, r reaches it, every deletion vertex reaches v then r, and v reaches r. Counting orders gives (1188K+10)+(12L+11)+2=1188K+12L+23. The width M-3167, thresholds K>=36,L>=44, the K=36 specialization, and both numerical examples are correct. Oddness follows from actual tournament path-count parity, and is independently checked at every tested rectangle baseline.

The universal conclusion follows from exact identities and the coin-box proof; finite tests do not prove universal coverage on their own. The interval baselines remain exponential, and no placement/overlap mechanism for arbitrary targets is supplied. The stated limitation is accurate. The 37-fold improvement concerns this explicit asymptotic width/order family, not an optimum or literature-priority claim.

## Independent calculations

`audit.py` does not import or execute the producer verifier. It independently enumerates all directed ordered vertex tuples for each core and combines subset split counts; checks every tournament orientation; recomputes the finite flanked word identities; checks literal Hamilton counts of both actual core star closures; and verifies the ordinal law on 32 actual cut-labelled triangle/pair compositions.

A separately implemented residue selector and direct reachable-set enumeration check 714 coprime bounded rectangles with a,b<=21, including sharp omitted lower/reflected upper boundaries. For each of the five producer parameter pairs, word powers are computed by repeated ordinary composition rather than producer exponentiation. Every physical position pair in each rectangle is checked for scalar count and common interface, and every target in each advertised interval is selected and physically recomposed. All five lower-endpoint hashes agree with the producer receipt.

Normal and optimized receipts agree byte-for-byte. `python3 -O audit.py --sentinel` intentionally fails with RuntimeError and exit code 1, recorded in SENTINEL.txt. All checks use explicit exceptions, never optimizable assertions.

These are direct graph controls only for small cores and small ordinal examples. Large witnesses use exact transfer counts; no literal subset DP on 43,000 vertices is claimed.

## Provenance clarification

The parent resolved the historical shorthand citation Session229 `dynamic/PROOF.md` to `05-knowledge/results/session229_combined_constructions_20261109/dynamic/PROOF.md`. The parent supplied this exact path after the initial audit; it was not materialized in the auditor’s checkout and was not independently inspected. The same seven-coordinate law and definitions are also displayed in canonical THM-4960. The minor provenance question is resolved; the stable theorem verdict remains PASS. No required fix identified.

Reproduce in an assembled packet with sibling directories `graft_audit` and `graft` by running `python3 graft_audit/audit.py` and `python3 -O graft_audit/audit.py`. For another layout use `--source PATH_TO_GRAFT`. The sentinel is `python3 -O graft_audit/audit.py --sentinel` and must exit nonzero before reading sources.

Artifacts: audit.py, CHECK.json, CHECK.optimized.json, SENTINEL.txt.
