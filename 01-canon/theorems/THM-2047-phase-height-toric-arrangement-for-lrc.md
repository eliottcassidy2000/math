---
id: THM-2047
title: The phase-height toric arrangement is an exact carrier for the lonely-runner max-min
status: RESERVED / IN PROGRESS. The phase-height dictionary and vertex mechanism are proved below; the proposed deletion-localization consequence for Wall A is not proved.
source: codex-2026-07-21-LRC-arrangement-audit
depends_on:
  - THM-1002
  - THM-1017
  - THM-2043
related:
  - HYP-8830
  - HYP-7310
---

# THM-2047 -- the phase-height toric arrangement

For a finite set of positive integral speeds `S`, put

`E_S={(t,delta) in (R/Z)x[0,1/2] : delta<=||v t|| for every v in S}`.

The boundary walls are the translated toric characters

`v t-delta in Z` and `v t+delta in Z` (`v in S`).

The following claims are reserved here and proved in the completed version:

1. the height-`delta` fiber of `E_S` is exactly the weak-safe phase set
   `G_delta(S)`, hence `M(S)=max{delta:E_S intersects height delta}`;
2. a top vertex is an intersection of oppositely signed owner walls (allowing
   the self-cusp), which recovers the exact pair-sum denominator law of
   THM-1002;
3. the owner, sign, phase, and height labels are therefore a lossless LRC
   carrier, while the unlabelled relation lattice, a standard toric-complement
   invariant, or any fixed local residue packet need not be;
4. deleting one speed deletes its two signed walls. Whether the resulting
   deletion/restriction cell structure can force the AP core required by
   THM-1017 remains open and is the intended Wall-A experiment.

This reservation explicitly does **not** claim that Shi-arrangement region
counts, Orlik--Solomon Betti numbers, or a standard toric-complement volume
equal the lonely-runner safe measure.
