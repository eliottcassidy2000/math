# Independent review of manuscript lane

**Verdict: PASS, with scope qualifications retained.** Reviewer: S292 triples lane. This review does not audit the external RH manuscript or its Lean development.

- Reran `manuscript/verify.py` normally and with `python3 -O`; both passed. Independently reconstructed the seven-zero witness by rational Gaussian elimination on the seven moment equations (fixing the final coefficient to7), rather than using Lagrange denominators. The displayed W and first difference3991680 agree.
- Checked the finite-observer theorem directly. With the convention M_w is the ordered product, left multiplication prepends a letter. Equality V_j=V_(j+1) therefore gives invariance under every letter. Since nonzero V_0 has dimension at least1, the bound d−1 is correct; zero output columns are the explicitly handled degenerate case. The block direct-sum bound d1+d2−1 is sufficient (not asserted optimal). This is a bound on full modeled continuation outputs, not arbitrary unmodeled physical guards or mere scalar observations without a carrier.
- Checked the eight-mode Vandermonde determinant and characteristic recurrence. Seven values do not suffice in the full rational span, even for a positive/negative split into nonnegative coefficient profiles. This does not claim either profile is a tournament, and does not contradict S291's single-length separation for its much smaller4096-profile set. The shifted witness is scoped to m=2..8 and separates at9.
- Independently differentiated F_S(t): F'_S=(A−E log t−2E)/(2sqrt(t)); the unique stationary point is the global real maximum. The ray restriction uses max(log(rad n),t_*). The discrete exponents need not attain that maximum. Prime support alone does not fix t or the Nicolas observable. The factors15/4 and7/2 are exact.
- Quartet-polynomial expansion, roots, strip coordinates and compact error bound are correct. The example disproves only inference from the stated symmetry/compact-convergence hypotheses; it is not a zeta-function counterexample.
- Source use is bounded and appropriately typed: external RH/independence claims are not imported, no Lean build is asserted, Nicolas is cited rather than derived. I checked this boundary in PROOF/SOURCES, not the independent validity of those external papers.

- Reviewed the added arbitrary-eight-distinct-real-length upgrade: the Rolle induction on nonzero exponential terms is correct. In particular the eight wheel-surviving sample lengths are valid and the exact rational elimination determinant is nonzero. Replayed the updated verifier both normally and optimized; transcripts match. This identifies the known eight-mode profile, not an arbitrary physically guarded family.

No corrections required by this review. Final reviewed hashes:

```json
{
  "PROOF.md": "74037ef702a3d616d47488a1be22e20023b249172b7a3ac026687e065fe95af7",
  "verify.py": "f1ce3451ae7553946b3906e12bb07b9689fa2f0da00d2f9030ddbdd94c3aafc3",
  "SOURCES.md": "4401872670ea861aa916fd62b69e8b02ec89a58f7bf29326defef485e03dd171"
}
```
