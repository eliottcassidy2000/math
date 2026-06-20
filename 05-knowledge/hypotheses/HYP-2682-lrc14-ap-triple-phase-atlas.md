---
id: HYP-2682
title: LRC(14) rank-one AP-triple phase atlas
status: OPEN; claimed for exact scout
source: codex-2026-06-20-S53
depends_on:
  - HYP-2681
  - HYP-2680
  - THM-548
  - HYP-2679
  - HYP-2677
  - HYP-2676
related:
  - HYP-2639
  - HYP-2637
  - HYP-2606
  - OPEN-Q-108
---

# HYP-2682 - AP-Triple Phase Atlas

## Claim Being Tested

The rank-one far relation

```text
u - 2v + w = 0
```

is too coarse to close the signed three-far proof obligation from HYP-2680.
For AP triples `F_m=(m,m+1,m+2)`, the exact relation lattice rank is constant,
but S51 and S52 already show the signed order words and cube-root packets change
with the offset `m`.

This hypothesis claims the missing finite datum is a phase/support address:

```text
(bounded core B, offset m, residue address, seven packets A..G,
 order sums R1,R2,R3, pair-tax shadow, Eisenstein imbalance).
```

The target is not a new scalar invariant.  The target is a finite atlas
coordinate that tells the simultaneous-peel proof when the rank-one branch is
a safe finite-resonant model and when it must pay the height-weighted signed
relation-lattice bound.

## Challenged Assumption

Do not assume exact relation rank determines the sign or size of
`Delta_3(B;F_m)-Phi_3(B)`.  The quotient to rank preserves the Freiman
branching predicate but destroys phase.  The quotient to the cube-root packet
preserves cyclic order information while destroying individual wall ownership.
The atlas should retain enough wall/support address to test whether the lost
ownership is harmless.

## Startup Evidence

A quick exact probe over `m=15..80` for three bounded cores found multiple
packet sign words even though every far triple has the same relation
`(1,-2,1)`.

For `B=(0,4,6,8,10,12,14)`, the sign word over

```text
(R1, R2, R3, total residual, pair-tax shadow)
```

has six observed types:

```text
+++++ 21
+-+++  2
++++- 9
++-+-  1
-+++- 19
-++-- 14
```

The top direct and total-residual rows occur at small offsets (`m=15,16,28,29`
in the quick probe), while later offsets appear safer for the tested cores.
This suggests the rank-one AP-triple branch may split into a finite small-offset
atlas plus a tail phase estimate, but no proof is claimed.

## Next Computation

Create `04-computation/lrc14_ap_triple_phase_atlas_codex_s53.py` and store
`05-knowledge/results/lrc14_ap_triple_phase_atlas_codex_s53.out`.

The scout should:

1. scan named bounded cores across AP triples `(m,m+1,m+2)`;
2. record exact seven-packet signs, order sums, pair-tax shadow, and
   Eisenstein imbalance;
3. group phase words by `m mod 7`, `m mod 14`, and offset range;
4. run selected all-core banks for small `m` where the quick probe sees the
   highest direct `p0`;
5. include Tournament Analysis with vertices chosen from phase words,
   residues, proof lenses, and core families rather than only runners.

No LRC(14) proof is claimed.
