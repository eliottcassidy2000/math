---
id: HYP-8846
title: "LRC14 pointed transport on the rank-eleven relation planes"
status: >
  OPEN / namespace reserved. THM-2052 reduces every hypothetical
  counterexample to a finite atlas of rational planes. HYP-4346 gives an exact
  rank-drop theorem for two independent directions with one common scale-only
  template, but its conclusion is existential in the plane and does not prove
  the original direction lonely. The missing theorem is a pointed
  owner-labelled transport from that escaped direction back to the target, or
  an independent relation raising the code to rank twelve.
source: codex-2026-07-21-DC2-LRC14-termination
related:
  - THM-2050
  - THM-2052
  - HYP-2108
  - HYP-4336
  - HYP-4346
  - HYP-8841
---

# HYP-8846 -- pointed transport, not an unpointed plane escape

This file reserves the exact next namespace after the THM-2052 reduction.
The known input and the quantifier gap are:

```text
known:  target v lies in a rational plane ker(W);
known:  a genuinely rank-two plane cannot send every independent direction
        to one of finitely many common proportional templates;
missing: transfer a safe direction in the plane back to the specified v,
         preserving a strict phase or an owner-labelled Euler endpoint.
```

The intended terminal is a pointed trichotomy for the target row `v`:

1. some peel has the HYP-2108 endpoint functional `P_w>=0`;
2. active owner data supplies a bounded relation outside `W`; or
3. a finite wall-crossing chain carries a strict/Euler certificate from a
   rank-rigidity escape direction to `v`.

Clause 3 is open. An unlabelled statement that *some* direction of the plane
is loose has the wrong quantifier and is not an LRC14 proof.
