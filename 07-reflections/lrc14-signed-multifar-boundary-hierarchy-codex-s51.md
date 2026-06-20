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

## Next Move

Prove a signed Abel packet bound for

```text
Delta_3(B;u,v,w) - Phi_3(B)
```

stratified by the minimum low-height form `|m*u+n*v+l*w|`, with exact
low-height forms routed to finite Freiman/scale atlases.  Then generalize the
same mechanism to `s>=4` using the Stirling coefficient table.
