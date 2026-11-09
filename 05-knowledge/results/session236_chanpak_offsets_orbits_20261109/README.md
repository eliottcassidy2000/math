# Tournament offsets and the fifteen-matrix orbit

**Session236, 2026-11-09. Independently audited. Chan-Pak completeness remains OPEN.**

[THM-5080](../../../01-canon/theorems/THM-5080-port-orbit-polyhedral-cones-sharp-growth-and-all-word-tower-ray-rigidity.md)
proves three new all-order results:

- The fifteen-generator ternary orbit has an exact four-ray conical hull;
  its nonempty-word orbit has an exact five-ray hull. These are minimal cones,
  not descriptions of every reachable integer state.
- Every move after the first expands h by at least 292/193, sharply. Every
  represented h has a logarithmic-length word and an O(log h)-vertex witness.
- In the complete 57-matrix ternary alphabet, only tower powers preserve
  (3,1,2), at every word length. This replaces the previous length-at-most-four
  check with a proof and blocks one specific fixed-ray digit strategy.

**Keep three meanings separate:** undirected graph offsets count dominating
sets; binary tournament H=I(Omega,2) uses the 56-matrix alphabet; the fifteen
matrices evaluate I(Omega,3). Neither graph equinumerosity nor the new cone
geometry supplies a statistic-preserving transfer or a completeness proof.

Actual tournament sink-padding followed by one final seam gives C+A*2^m,
where C=H(T)H(B) and A retains rooted-response data. Near-join transport is
inherited from THM-4865/HYP-9453; mixed-periodic-word drift is inherited from
THM-4882/Session198. The offset audit corrected initial novelty wording; all
formulas survive. The positivity restriction concerns one final bridge only.

Concurrent S19's audited multi-bridge tails are a different grammar: its
order-seven core yields 2^K-159 for K>=10. That does not contradict this packet's
single-final-bridge barrier. S19's finite coverage and its withdrawn finite-bank
cofiniteness proposal do not prove our fifteen-matrix orbit conjecture either.
See [S19 audit](../kps_s19_fixed_port_tails_20261109/audit/AUDIT.md).

## Proofs and verification

- [Orbit cones, sharp growth and fifteen-generator proof](orbit/PROOF.md)
- [Complete 57-generator all-word proof](last_letter/PROOF.md)
- [Recovered offset mechanisms and exact limits](offsets/PROOF.md)
- [Independent geometry audit](orbit_audit/AUDIT.md)
- [Independent offset/provenance audit](offsets_audit/AUDIT.md)
- [Reproduction](REPRODUCE.md)

The geometry audit rebuilds all 1,105 marked blocks by actual gluing and direct
cycle/path enumeration, checks all 615 cone inequalities, and independently
recovers 53,732 states through h=3000. This reproduces old finite data, not a
coverage extension. The offset audit checks 116 padding identities, 1,682
bridges, 511 apex counts and all 131,070 binary period words through length16.

**Session end:** completed and independently audited. The remaining obligation
is uniform integer reachability for every sufficiently large admissible first
coordinate. Neither a completeness proof nor a counterexample is asserted.
