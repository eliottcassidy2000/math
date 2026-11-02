# Final Session137 weighted-stability manuscript receipt

**Verdict: ACCEPT.** Read and checked the final publication manuscript against
the previously accepted independent stability proofs. No new mathematical
claim or changed constant was found.

## Reviewed manuscript

- Repository-relative path:
  `05-knowledge/results/session137_overnight_proof_frontiers_20261102/discriminant/WEIGHTED_STABILITY.md`
- Size: **4,798 bytes**.
- SHA256:
  `6bd43447a41e092721d0705069721e46f164eed2428c1655333209eb5fd5ce7c`.

## Scope checked

1. The augmented-skew principal-minor argument gives the two signed
   cofactor bounds for nonempty remainders. Only the unsigned u bound is
   extended to the empty remainder; the false signed empty bound is excluded.
2. The weak gap satisfies `delta>=c_n V` for n>=2, with the exact n=2 base.
   Its recurrence is applied inductively for n>=3, not as a definition of
   an undefined c_1. The closed form and finite-sum normalization agree.
3. The strong gap satisfies `Delta>=beta_n V` for n>=3. The root averaging,
   numerator increment, lower bound `6/[n(n-1)]`, and asymptotic `24/n^2`
   all match the accepted proof. The n=2 counterboundary is retained.
4. The identity `Delta=V` at n=3 and the two symmetric-pair n=4 sharpness
   example have the correct H,u,v,V values and constants. No sharpness claim
   above n=4 is made in this Session137 manuscript.
5. The rounding bound uses unordered-edge distance and applies in the
   strong-bound domain n>=3. The factor two for full matrix L1 distance is
   correct; transitivity of the rounded matrix is not inferred.
6. The verification counts and the inherited-versus-new provenance agree.
   The informal THM3729 reference resolves to
   `01-canon/theorems/THM-3729-rooted-pfaffian-response-and-sign-root-deletion-average.md`;
   its rooted principal-minor identity is consistent with this argument.

## Existing evidence retained unchanged

The matching independent proof is
`independent/WEIGHTED_POLYNOMIAL_AUDIT.md`, SHA256
`8f0188bbc655e2e3632c8ff4b30c823948869995e5e3423daa0879274f5462d3`.
The earlier cofactor proof is `independent/WEIGHTED_STABILITY_AUDIT.md`,
SHA256 `b1dafa9de34c491f9ab33574937ac8febc883fe090415409b03c4baf9758bb0c`.

The existing normal/optimized polynomial receipts both have SHA256
`5a953fb038f72f19b256e0962706907efdc48f15f4003b243b267c896293fa7b`:
856 rational matrices, 3,627 weak and 3,621 strong root inequalities,
27 ternary order-three identities, and coefficient checks through n=200.
Their bytes were checked; no large replay was needed for this unchanged
manuscript reconciliation.

No sharpness claim beyond order four or additional triangle refinement is included.
The frozen 89-file independent payload was not edited. Only transient review-location metadata is omitted from this publication copy.
