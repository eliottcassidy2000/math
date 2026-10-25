# Independent P11 stress audit

Verdict: PASS. Frozen producer MANIFEST.json SHA-256: 0855158ced20c231100eeb5d8faf9094337e2a28dd26ad74befb437997814f26. No producer or repository files edited.

## Mathematical audit

For a doubly regular tournament n=2r+1=4lambda+3, a deficiency-two Hall certificate deletes r|S|-e(S,T) >= 2r arcs. Equality is exactly N^-(T) subset S. Cauchy applied to the incidence degrees gives the stated expansion bound. Direct expansion confirms D(t)=r(r-3)+lambda(t-1)(2r-2-t); it is strictly positive throughout 1<=t<=2r-2 for r>=5. The two endpoint cases are precisely two sinks and two sources. A matching-deficient deletion admits a certificate with deficiency at least two; enlarging T to deficiency two is valid and gives a required deletion subset of size at least n-1. Thus equality exhausts all size-(n-1) Hall blockers. Distinct source pairs and sink pairs give n(n-1) masks: other vertices retain positive degrees, so the corresponding two sources or sinks are recoverable from the residual digraph, and the two classes cannot coincide.

At n=7 the interior equality endpoints t=1,4 are genuinely attainable. Independent enumeration finds 7 certificates/masks in each extra branch, plus 21 in each source/sink branch. This verifies the sharp exceptional equality, not a general HP-classification theorem.

The order-five multiplier fixes no arc and has eleven size-five arc orbits. Hence invariant deletion sizes are multiples of five. Empty deletion has HP (already follows from any verified five-deletion survivor). All eleven invariant five-deletions survive; among 55 invariant ten-deletions only star0 blocks. Since stars give beta<=10, the minimum-blocker family's fixed-point count is zero for beta<10 or one for beta=10. The full ten-blocker family has exactly one fixed point. The specified F1024 nonzero Frobenius action fixes gcd(3,1023)=3 elements, contradicting either equivariant identification. No bare cardinality 1023 contradiction follows.

The extension to nontrivial Frobenius powers is sound: gcd(4^j-1,1023)=3 for j=1,...,4. For a genuine order-five semilinear multiplier a*z^(4^j), the norm-one condition puts a in the image of z -> z^(1-4^j), yielding precisely three nonzero fixed points. This does not exclude arbitrary actions on a 1023-element carrier.

## Reconstructed computation

The independent audit reconstructs arcs from quadratic residues and verifies regularity and all common-in-neighbor counts. It enumerates all 497,420 P11 and 2,002 P7 Hall certificates by scalar rectangle costs, comparing complete histograms and masks. It independently uses forward subset endpoint propagation for HP existence, rather than importing the producer's backward endpoint recurrence.

- P11: 110 Hall masks, 11 non-Hall stars; 121 HP blockers; all 1,210 single-arc restorations have HP.
- P7: 56 Hall masks, 7 stars; 63 HP blockers; all 378 single-arc restorations have HP.
- Every supplied star-complement cycle is checked directly.
- Order-five partition reconstructed independently; all 11 invariant five-sets and 55 invariant ten-sets checked. All 65 surviving path witnesses are verified directly, and the unique failure is star0.
- Independent normal and -O results agree byte-for-byte.
- Producer replay passes from arbitrary cwd with normal/-O equivalence and built-in mutation hostiles.
- Additional missing-file, nested-extra-file, and modified-content closure hostiles are rejected in normal and -O modes.

## Scope boundary

This certifies the general Hall theorem and the specified symmetry refutation. It does not determine beta(P11), enumerate all ten-arc HP blockers, refute the bare number 1023, or prove general inclusion-minimality of stars/source/sink pairs. Inclusion-minimality here is finite: 121 P11 examples and 63 P7 examples. No blocker, denial, or uncertain action remains.

Reproduce: python -B audit.py and python -B closure.py from this directory; also run each with -O. Results are independent.json, independent-O.json, closure.json, closure-O.json, and replay.log.

## Portable release packaging

The complete necessary producer packet is frozen under producer/ and pinned to the manifest hash above. Audit and closure scripts resolve that local copy, with no sibling-directory dependency. MANIFEST.json lists explicit relative paths and SHA-256 digests. replay.py validates the recursive exact file set and all digests before computation, then reruns the independent audit, closure checks, and original producer replay in normal and optimized modes. Its --self-test copies only manifest-listed files into an isolated directory, actually runs that copy in both modes, and proves corrupt, missing, and nested-extra files are rejected in both modes. Run python -B replay.py --self-test; optionally repeat with -O. Producer files and mathematical claims are unchanged.
