---
id: HYP-2138
status: OPEN classification route; S592 verifies n=6,8,10 and isolates n=14 prime-3 shell residual
source: opus-2026-06-03-S592
related:
  - THM-401
  - THM-402
  - HYP-2135
  - HYP-2153
  - HYP-2137
  - HYP-2136
  - HYP-2132
  - HYP-2088
  - HYP-2083
---

# HYP-2138: non-transversal sporadic tight configs come from composite C=2n-1 shell strata

## Claim

Non-transversal sporadic floor-tight LRC configurations occur precisely when
the THM-401 shell modulus

```text
C = 2n - 1
```

has nonunit antipodal shell strata.  If `C` is prime, every missed shell is
unit-visible and the support ledger should force a clean transversal or an
explicit edge witness.  If `C` is composite, nonunit shells can hide the
unit-inverse witness and create sporadic floor rows.

## Evidence From S592

`04-computation/lrc_worryset_sporadics_s592.py` checks the small tight census:

```text
n=6,  C=11 prime: sporadic is still a transversal flip.
n=8,  C=15 composite: two non-transversal sporadic tight rows.
n=10, C=19 prime: AP only.
```

For `n=14`, `C=27=3^3`; the known `V*` row misses/doubles only gcd-3
shells, so the classification residual localizes to the prime-3 shell branch.

## S602 Additive-Chain Refinement

`lrc_p0_collapse_additive_chains_s602.py` checks the operational
`p0`-collapse predicate in targeted primitive boxes and recovers the two known
`n=8`, `C=15` non-transversal floor rows:

```text
(1,2,3,4,5,7,12)       misses nonunit shell 6, with 12=5+7;
(1,4,5,6,7,11,13)      misses nonunit shell 3, with 13=6+7.
```

Both still collapse to the unit-boundary skeleton, and both are two-seed
addition chains.  So the composite-shell branch is not only a shell-support
phenomenon; it carries an additive-chain normal form that should be part of
the n=14 prime-3 ledger.  HYP-2153 records this as the `p0`-collapse
classification subproblem.

## Relation To HYP-2135

HYP-2141 supplies the labelled support calculus.  HYP-2138 is the proposed
classification law for one of its residual fibers:

```text
unit hole      -> inverse 2/C witness;
prime C        -> all holes are unit holes;
composite C    -> nonunit holes require gcd/CRT/descent labels;
sporadic floor -> nonunit shell fiber survives to the floor.
```

This is complementary to the prime-2/doubling/C-prime residual.  The n=14
proof route is therefore not one prime-local story:

```text
prime 2: C-prime / multiple-of-14 / doubling seam;
prime 7: solved cyclotomic block Q(zeta_14)=Q(zeta_7);
prime 3: composite C=27 shell-sporadic branch.
```

## See

`07-reflections/lrc-sporadics-iff-2n-1-composite-prime-3-at-n14-s592.md`,
`04-computation/lrc_worryset_sporadics_s592.py`,
`05-knowledge/results/lrc_worryset_sporadics_s592.out`,
`04-computation/lrc_p0_collapse_additive_chains_s602.py`,
`05-knowledge/results/lrc_p0_collapse_additive_chains_s602.out`,
HYP-2153, HYP-2141, THM-401, THM-402, HYP-2088, HYP-2083.
