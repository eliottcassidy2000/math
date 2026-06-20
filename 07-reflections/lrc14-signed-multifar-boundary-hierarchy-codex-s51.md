# LRC14 Signed Multi-Far Boundary Hierarchy

codex-2026-06-20-S51

## Summary

The signed two-far curvature from THM-548/HYP-2679 is the second member of a
Newton hierarchy.  The decorrelated `s`-far mixed term is

```text
Phi_s(B) = 7^-s * sum_{t=1}^s (-1)^(s-t) t! S(s,t) p_t(B),
```

where `p_t(B)` is the missed-sector profile of the bounded core.  The first new
target is therefore

```text
Phi_3(B) = (p1 - 6*p2 + 6*p3)/343.
```

## Correction

The exact Newton identity is over all subsets of the far set:

```text
p0(B union F) = sum_{S subset F} Delta_S(B).
```

It should not be stated as truncating at `|S| <= 6`.  The six-sector structure
enters through the fact that the coefficient formula only needs `p_0..p_6`,
not through literal disappearance of higher Newton differences.  Redundant
far-sector hits still create higher signed differences.

## Exact Evidence

The script verifies the boundary-value identity

```text
P_r(B)=p0(B)+sum_{s=1}^r binom(r,s)Phi_s(B)
```

against THM-548's direct profile formula for three cores and `r<=6`.

The three-far deviations are relation-lattice effects.  For the HYP-2675
dilated core:

```text
F=(15,16,17): dev=31753/11428760, relation -15+2*16-17=0
F=(15,18,21): dev=18443/3025260, relation -15+2*18-21=0
```

For the separated `consec8` triple `(17,23,31)`, the deviation drops to
`421/71900346`.

The all-core bank for `(15,16,17)` checked `3003` primitive cores.  All have a
triple relation and no pair relation; deviations split `1999` positive and
`1004` negative.  The largest absolute deviation was
`40633081/445721640`, but that row still had direct cap margin
`2476301/7147140`.

## Interpretation

Pairwise resonance is not the whole story.  Consecutive triples have no exact
pair relation at small height, but they have the exact three-body form
`u-2v+w=0`.  That is the first place where the signed multi-far proof must use
the full relation lattice rather than a pair atlas.

The geometry is encouraging: the large three-far discrepancies found here are
not direct cap threats.  They are proof-obligation threats: without a signed
relation-lattice bound they prevent a clean analytic closure of the true-wide
branch.

## Tail-Rank Update

The second S51 scout shifts the target from isolated residuals to signed order
sums.  For a whole far block `F`,

```text
p0(B union F)-P_|F|(B)
  = sum_s R_s(B;F),
R_s(B;F)=sum_{|S|=s}(Delta_S(B)-Phi_s(B)).
```

The consecutive far blocks show the pattern the user pointed at: opposite
bounded signs across adjacent Newton orders.  For the dilated core
`(0,4,6,8,10,12,14)`, far `(15,16,17)` has order signs `+-+`, far
`(15,16,17,18)` has `+--+`, and far `(15..20)` has `++-+-+`.  In the exact
four-far bank for `(15,16,17,18)`, `R2/R3` have opposite signs in `1644/3003`
primitive cores and `R3/R4` oppose in `2053/3003`.

Scaling the AP-relation block keeps the relation rank but changes the signs.
For the same dilated core, triples `(m,m+1,m+2)` at
`m=15,22,31,43,61,89` give sign words `+-+,-++,-++,+++,-++,-++`.
Four-blocks at `m=15,22,31,43,61` give `+--+,-++-,-++-,+++-,-++-`.
So the finite rank-one atlas must retain phase/support labels.

This matters because an absolute-value proof would spend the same margin twice.
The right analytic lemma should keep the signed order ledger until the last
step, with exact low-height relation rank routing the finite atlas branch and
high-rank blocks paying Abel/Koksma cancellation.

This is not a rank-only lemma.  The repo's older summand/multiplicand notes
warn that equal additive energy or equal relation rank can carry different
observer-visible signs.  The residual packets need typed labels: summand shell,
multiplicand clearance, support, visibility, and sign.

It is also not a fixed-box denominator lemma.  The adversarial THM-548 audit
found two-far examples whose shortest relation sits just outside the tested
coefficient box, so `resdist` at bounded height overstates independence.  The
multi-far denominator should be a height-weighted relation-lattice packet.

## Next Move

Use the incoming THM-548 simultaneous-peel architecture.  For `r=3`, decompose
the full row as

```text
p0(B union {u,v,w}) =
  P_3(B)
  + three one-far residuals
  + three two-far residuals
  + one three-far residual.
```

The one- and two-far pieces are now routed by THM-547/548.  The new analytic
piece is a signed Abel packet bound for

```text
Delta_3(B;u,v,w) - Phi_3(B)
```

stratified by the minimum low-height form `|m*u+n*v+l*w|`, with exact
low-height forms routed to finite Freiman/scale atlases.  Then generalize from
single triples to signed order sums `R_s(B;F)` using the Stirling coefficient
table and the tail-rank scout.
