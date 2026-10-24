# Independent audit: counting lenses and native metagraph boundaries

**Verdict: PASS for the stated inherited counting formulas, finite computations, operator-separation claims, and proposed historical warning.** This is an independent audit of `session52-counting-lenses`, not a new enumeration theorem, new scalar-count algorithm, native wiring reconstruction, or an order-eleven native prediction. No shared repository files were edited and nothing was published by this reviewer.

## What was independently checked

The producer's complete README, source code, historical-warning patch, correction note, provenance and receipts were read. The pair-cycle polynomial and anti-automorphism proof are valid as stated, including both necessity and sufficiency. The exponent

    sum_i L_i/2 + sum_(i<j) gcd(L_i,L_j) + f*k

uses non-singleton cycles L_i congruent to 2 modulo 4 and f in {0,1}. The two missing restrictions in the historical draft are genuinely independent: (1,1) and (4) both have zero anti-fixed tournaments. Positive controls (2), (2,2), (6), (1,6) have respectively 2, 16, 8, 16 anti-fixed tournaments.

`independent_audit.py` imports no producer functions and uses a different computational route: it writes literal affine orientation equations over GF(2) and counts their solutions by Gaussian elimination, without decomposing the pair action into cycles. Every vertex-cycle type at every order 0 through 20 is tested against the claimed criterion and exponent, then averaged with its conjugacy-class denominator. Independent literal labelled-tournament orbit generation through order five reproduces ordinary and self-converse class counts. Exact integer characteristic polynomials of the directly built arc-flip quotients agree with the producer's signed-polynomial multiplicities at these orders. These are stronger controls than matching only the two polynomial endpoint values.

The independent audit has **5,598 explicit checks**, unchanged under Python -O; normal and optimized receipts are byte-identical. The producer has **399,751 explicit checks**, also with byte-identical normal and optimized receipts.

Both scalar vectors through order 20 match separately retrieved public sources: the [OEIS A000568 b-file](https://oeis.org/A000568/b000568.txt) and the [OEIS A002785 entry](https://oeis.org/A002785). The latter explicitly identifies self-converse tournaments and starts at n=1; the empty order-zero convention is supplied separately. The A002785 b-file fetch failed and was not used. `inputs/external-oeis-prefixes.json` records the retrieved prefix, URLs and this limitation. Separate retrieval is source independence, not a claim that OEIS enumeration is theoretically unrelated to Burnside.

## Operator and quotient typing

The n=3 operators are independently constructed in transitive/cyclic order:

    native fixed-HP full-free-bit complement: [[0,1],[1,0]], eigenvalues 1,-1;
    labelled arc-cube S_n quotient:           [[2,1],[3,0]], eigenvalues 3,-1.

This is already a decisive counterexample to transporting the arc-quotient spectrum to native line incidence. Both classes are self-converse, so converse merging does not remove the discrepancy. The matrix weighting and operation differ, not merely the names of vertices.

Literal fixed-path enumeration also confirms the first non-descent order is four. Tilings 0 and 1 have the same tournament class, but their free-bit complements 7 and 6 have different converse-merged classes. Hence a native line target is not a function of the bare tournament class; a Hamiltonian-path presentation or aggregated line incidence remains necessary. Counts of classes are nevertheless well defined and the Burnside node-count identities remain valid.

The signed polynomial proof correctly averages character traces, and its t=-1 specialization is correctly restricted to a graded antipodal trace and a fixed-class count. No differential, topological chain complex or Betti-number theorem is silently supplied. The producer does not repeat the inherited THM-588 sentence that incorrectly calls *uncentered* cyclicity an eigenfunction. When that statement is used elsewhere, the exact eigenfunction is c(T)-binom(n,3)/4. At n=3, four times that vector is (-1,3), whereas c=(0,1) is not an eigenvector. Likewise, the isolated THM-587 proof phrase averaging dim(W^sigma) should be read as averaging traces, not imported literally. These inherited wording defects do not affect this packet's self-contained valid arguments and were not expanded into unrelated canonical edits.

## Correction lineage and scope

The warning patch is approved. It prepends a dated HISTORICAL / SUPERSEDED boundary, redirects to the already authoritative THM-283 and its October 21 endpoint-monodromy repair, and leaves the entire old body intact. It does not recast an inherited correction as a new theorem. The proposed MISTAKES entry distinguishes warning the stale draft from correcting THM-283's already repaired explanation. Its final packet link must be replaced by the actual publication destination, and any maintained MISTAKES front banner must stay at byte zero.

The parent synthesis should keep A000568(11)=903753248 typed as a known classical baseline. If it uses the separately audited correction E_11=152690 with U_11=11*A000568(11)+E_11, the resulting **9,941,438,418 is a backsolved consistency target**, not an independently evaluated native upper trace or a novel prediction. The present lane neither computes E_11 nor evaluates the native U_11 sum. No such computation is certified by its partition formulas.

## Portable reproduction and source closure

All eight producer manifest entries and its manifest are frozen under `inputs/producer/`. All twelve provenance source files were hash-verified against their original declared locations and frozen under `inputs/source-snapshots/`, including the corrected theorem, historical false draft, earlier proof/review, and historical growth packet. The snapshots exist for audit provenance, not to reinstate superseded claims.

From this directory:

    python3 -B independent_audit.py --receipt /tmp/audit-normal.json
    python3 -B -O independent_audit.py --receipt /tmp/audit-optimized.json
    python3 -B closure_replay.py

The closure replay copies **only manifest-listed files** into a fresh temporary tree, verifies every hash, and invokes producer and independent auditors with isolated Python, in both ordinary and optimized modes, from an unrelated working directory. Each generated receipt must exactly match its frozen reference. It also validates all source-snapshot hashes, applies the warning to the temporary historical source and verifies byte-for-byte body preservation and redirect targets, then mutates an expected external count and requires rejection in both Python modes.

`closure-replay.json` records the actual replay results. The manifest excludes itself and that replay receipt to avoid self-referential hashes; neither is a runtime input of the mathematical checkers. No undocumented sibling directory, untracked helper, optional package, live network request, or shared repository file is needed by the four relocated mathematical executions.
