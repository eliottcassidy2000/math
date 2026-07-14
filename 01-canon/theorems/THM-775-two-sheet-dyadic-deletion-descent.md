---
id: THM-775
title: Two-sheet dyadic deletion descent and the exact Z/4 seam
status: CLAIMED (number reserved; elementary proof written and independent audit in progress)
source: codex-2026-07-14-S9
depends_on:
  - THM-765
  - THM-769
  - THM-772
related:
  - THM-774
  - HYP-6820
---

# THM-775 — Two-sheet dyadic deletion descent

This stub reserves the theorem number for the next recursive reduction of the
two-exception equality packet.  Let

```text
A=2U union {x,y},       |U|=10,       x,y odd
```

be primitive and tight, with THM-769's exact folded cover over `G_U`.
THM-772 already proves that `U` is primitive and divisor-complete for every
modulus `2,...,12`.

The intended theorem is:

1. either every nine-speed deletion `U\{u}` is primitive, or there is a unique
   odd `u in U` such that

   ```text
   U=2V union {u},       |V|=9,       gcd(V)=1;
   ```

2. in the second case, for every `sigma in G_V` the four lifts are covered by
   an exact disjoint `2+1+1` ownership tiling: `2u` owns one parity class of
   two sheets and `x,y` own the two sheets of the other parity;
3. the quotient `V` again contains a multiple of every `m=2,...,12`.

The proof carrier is a saturated sheet hypergraph.  For an imprimitive
deletion of gcd `d`, a closed `1/13` tooth of the deleted core speed and a
closed `2/13` eligibility tooth of an odd exception cover the `d` sheets.
Each has capacity at most one half.  Equality forces `d=2`, producing the
dyadic seam.  The seam lifts to `Z/4`; nearest-integer colours give the exact
`2+1+1` partition.  Unit fractions then transfer all divisor pins to `V`.

A proposed corollary iterates this mechanism: every failure of hereditary
primitivity in a quotient is another dyadic seam, so the only gcd pathology is
a tower of `Z/2` sheet extensions ending at a hereditarily primitive,
divisor-complete quotient.  Until the capacity, strictness, and iteration have
passed the independent audit, cite this file only as a claimed theorem.
