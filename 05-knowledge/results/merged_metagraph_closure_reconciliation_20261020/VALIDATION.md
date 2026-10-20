# Validation and provenance

Date: 20 October 2026. Verdict: **PASS within the stated scopes**, with a separately audited published theorem as an explicit input. No finite check is extrapolated to prove an unbounded statement.

## Mathematical proof audit

The [independent receipt](INDEPENDENT_AUDIT.md) accepts:

- the full closed-class classification using the already source-audited Benhocine–Wojda theorem, Redei parity, and exact small cases;
- closure inventory and colour-resolved all-shore inequalities;
- positive strong-core outdegree and reachability of a directed cycle at n>=6, allowing loops and without claiming connectivity;
- both explicit asymptotic bounds and the leading coefficient 2;
- the actual anti-automorphism hostile, with the distinction between invalid criterion and unrefuted structured witness family retained.

The published bypass theorem is a cited input, not newly reproved literature. Its primary citation and precise n>=7, p=2 scope appear in the theorem note. The independent audit does not recertify every historical proof in THM-4730 or every growth-packet result.

## Portable exact replays

From this packet directory:

    python3 verify.py
    python3 independent_audit.py

The outputs match `verification.json` and `independent_audit.json`, respectively. The first checker uses explicit permutation enumeration and direct arc comparisons. The second independently uses adjacency bitsets, DFS Hamiltonian counts, separately constructed path cubes, and direct tests on directed arcs. Neither imports the other. Failure checks are explicit exceptions, not removable assertions.

Validated universes:

- All 8, 64 and 1,024 labelled tournaments at n=3,4,5; all Hamiltonian words; full-permutation isomorphism canonicalization of the fully closed survivors.
- Positive controls C3 and nonregular Q5; negative control regular cyclic C5 with five open Hamiltonian paths.
- All 36 directed arcs of the order-nine criterion hostile.
- Every normalized skew connection set at odd n=7..23, all eligible unordered three-interval parameter pairs, and each exceptional descending-pair path: 2,044 connection sets, 41,992 genuine violated-arc witnesses.
- Every fixed-path tiling at n=3..6. All node shores through n=5 for each of total, blue and black incidence; singleton, strong, nonstrong, empty and full shores at n=6.
- Exact integer recurrence and inequalities through n=100 in the independent checker, with no floating-point proof step.

Normal and Python -O outputs agree. The independent checker is copied byte-for-byte from its final portable audit version, with a fixed reviewed-source digest as provenance metadata; its output is byte-for-byte identical to the independent retained output. No workspace or repository path is required to run either checker.

## Provenance and proof boundaries

Initial incoming snapshot: `70ca1d8210dfd7db3a32a2687f74a331dfccf910`; targeted final refresh: `eb6c3e9818e00fb2e768d6a4e7eb41ec731ee60f`.

- The closure coordinate comes from concurrent THM-4729; first-cut/nonstrong results overlap independent earlier proofs.
- Benhocine–Wojda was already used and source-audited in the independent rigid-strong packet. The new connection here is its full closed-class consequence.
- The odd-circulant classification used by the corrections is the separate explicit-arc prefix proof, not a repaired version of the incoming sum test.
- THM-4730's rotation/even-order route receives separate credit; the incoming odd proof gap is not concealed by the true final classification.
- The corrected leading asymptotic follows from the existing exact recurrence and the explicit endpoint/middle bounds given here.
- The chain-root growth packet remains independently credited; newly populated THM-4731/HYP-9346 are preserved as source-reported finite work, not imported as new proved dependencies.

The add-only result packet has no scarce theorem-ID reservation. Proposed canon corrections separately preserve all historical body text and original status, inserting dated precedence notices and updating current frontmatter. They require exact-prior fingerprint checks; no incoming file is silently overwritten.

The final refresh found both canon targets changed. Their overlays were
rebuilt from the new exact priors, retaining the incoming proof-gap warnings
and self-line correction verbatim. The duplicate ledger correction was
removed. No original approved publication manifest was changed. The new
independent order-nine hostile in the remote audit was checked directly;
its broader census, the new self-line counts, and the weighted-core proofs
were not independently rerun here.
