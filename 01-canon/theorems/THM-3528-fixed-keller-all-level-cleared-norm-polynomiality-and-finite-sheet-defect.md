---
id: THM-3528
title: "Fixed Keller all-level cleared-norm polynomiality and finite-sheet defect identity"
status: >
  RESERVED / PROVISIONAL PROOF UNDER AUDIT.  Candidate statement: the raw
  cleared norms of the fixed sporadic Keller map are polynomial and carry the
  THM-3522 packet at every level.  The old-L multiplicity of the next raw
  rung equals the nonnegative valuation on the regular finite inverse sheet.
  This would not prove L-coprimality at later levels, newest-image status,
  irreducibility, separability, distinct nonproperness components, an
  arbitrary-map law, or any general Jacobian-conjecture claim.
source: codex/all-level-cleared-norm/2026-08-16
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-3498-level-four-old-boundary-cancellation-and-degree81-discriminant-gate
  - THM-3522-fixed-keller-five-face-renewal-propagation
related:
  - THM-3523-fixed-R6-finite-sheet-unit-and-next-old-L-clearing
  - THM-3527-fixed-R7-finite-sheet-unit-and-next-old-L-clearing
---

# THM-3528 -- reserved all-level cleared-norm gate

**RESERVED / PROVISIONAL PROOF UNDER AUDIT.**

The proposed proof combines three already-proved mechanisms:

1. over `U=Spec(Q[a,b,c,L^-1])`, the inverse cover is finite etale, so the
   norm of every source polynomial is regular on `U`;
2. a complete packet `A(e,m)` gives valuation `-e/2` on each of the two
   divergent old-`L` sheets, while the regular finite sheet has valuation
   `s>=0`; hence

   ```text
   v_L(N(P))=-e+s,
   v_L(L^e N(P))=s>=0;
   ```

3. polynomiality then activates THM-3522 and propagates
   `(e,m)->(7e-2m,3e-2m)`.

The proof, cone audit, normalization conventions, and exact scope are being
written and independently hostile-audited.  This reserved file has no proved
downstream dependencies until its status is promoted.
