# Independent positive-approximation promotion readback

**PASS.** This approval covers only the constructive actual-source approximation packet. Exact target occupancy remains **OPEN**. Union witness recognition and bounded-error counting have separate proofs, audits, and promotion approvals.

## Exact approved changes

The independently regenerated promotion diff matches exactly. The proof/README changes are audit status, the correct digital-graph sibling link, and packet-local reproduction commands. No source convention, formula, run estimate, discrepancy constant, constructor, residual, finite control, or occupancy limitation changed.

- Promoted `PROOF.md` SHA256: `d87d754a536d550f2d1b02949db667cf6c697dd524bb49fd4a35057eead8c546`.
- Exact `PROMOTION.diff` SHA256: `7ae3958521e7ad6f8a34162189752038d3644cd2b54239b6c9b101956a02e188`.
- Unchanged `verify.py` SHA256: `49b32ffefd693c41ba89137126a738e73b388343642f0c74389719d2914a074d`.

`RELEASE.md` was separately read. Its summary correctly describes approximate actual indices rather than exact prescribed counts. The smaller grid-aligned bound is explicitly qualified by k>=1; the k=0 singleton and k=1 real-cutoff hostile remain correctly treated in the proof. Final summary SHA256: `ee816d4d3533485ff84854a0347761fc86cc375131e7f6850df2523c0fe27394`.

## Frozen artifacts

Every source file matches its frozen input hash, either unchanged in this release or archived under `originals/`. The verifier, source ledger, and both validation receipts are unchanged. All nine previous manifest-listed independent audit artifacts and their original audit manifest were copied byte-for-byte. The corrected real-rectangle bound and its correction lineage remain intact.

## Relocated replay

Normal and optimized author runs from the release directory both reproduce the frozen receipt exactly, SHA256 `e96482f1b12006f85d45fa4ea85c1257419b055b18e042fb7c630909ec297fd2`.

Normal and optimized independent runs from release `audit/` both reproduce the frozen receipt exactly, SHA256 `fd083354de87b2c984de8b155aaa1f87ae6ecfa9a573275de1269696100593d6`.

These replays retain the exact actual-index, run, endpoint-CDF, rectangle, and J-residual controls. No new mathematics or code is introduced by the promotion.

The audit manifest is refreshed to include this report. The publisher must refresh and verify the aggregate release manifest afterward. No GitHub write or remote publication is certified by this local approval.
