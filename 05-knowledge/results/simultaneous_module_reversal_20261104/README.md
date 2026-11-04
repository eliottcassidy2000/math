# Simultaneous disjoint module reversals: exact query isolation

**PROVED S1-S3.** In any fixed named tournament, simultaneously reversing
any nonempty selection of pairwise disjoint nontrivial modules changes its
full named two-port response, already modulo four. More strongly, for r
such modules there are r fixed actual queries recovering every one of the
2^r reversal vectors. At query i the exact integer value is

    F_e(q_i)=A_i+e_i*d_i,          d_i=2 (mod4).

Give the two ports arbitrary internal incidences on block i, and the root's
incidences on the exterior. Every other block remains a module, and its
complete ordered path-cover profile is converse-invariant. Its reversal
therefore leaves this query exactly unchanged. The inherited E3 separator
and odd-factor T1 transfer make d_i nonzero. This is an exact isolation
argument, not an assumption that full response differences add.

**All cubic and quartic terms are retained.** The complete tables contain40
responses for two-level reversal families of old orders4,6,8. A genuine
mixed quartic interaction is -2; it vanishes on the decoder queries.
The incoming cover normal form independently checks all16,432 quartic
coordinates, including zeros, and20 whole-converse pairs. The earlier
two prime twelve-core complete named cubics remain linked and pinned.

**FINITE-EXACT:**37,210 compiled actual-query counts; eight complete
Hamilton-word polynomial checks;40 independent forward-DP counts; and
all75 labelled quotients through order4 inflated by pairs, covering1098
reversal states and4306 decoder queries. All normal/optimized receipts
agree. Full response vectors are compared literally, not just by hashes.

**OPEN boundary:** arbitrary nested or overlapping reversal families and
unrelated-pair reconstruction. The cyclic three-pair nested control defeats
uniform modulo-four quotient probes, but its full responses separate251/249.
No exact full-response collision is found. This packet does not prove
completeness of the Hamilton spectrum outside7 and21.

- [PROOF.md](PROOF.md): statements, mechanisms, quantifiers and exact hostile.
- [REVIEW.md](REVIEW.md): validity gate and reproducibility scope.
- [BOARD.md](BOARD.md): inheritance, concept board and stopping boundary.
- [INCOMING.md](INCOMING.md): integrated source mechanisms and scope.
- `full_response.json`: full named coefficients and actual interpolation counts.
- `certificate.json`, `audit.json`: exact decoder and independent cover controls.

From the repository root:

    python 05-knowledge/results/simultaneous_module_reversal_20261104/verify.py --replay

The nonmutating replay runs both producers normally and with python -O in
a temporary directory and compares every output byte. It needs the existing
Python and g++ installations. The all-order theorem rests on the proof,
not on the bounded controls; the audit is same-session, not a separate agent.
