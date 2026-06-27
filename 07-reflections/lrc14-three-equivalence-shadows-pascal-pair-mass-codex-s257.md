# LRC14 Three Equivalence Shadows And Pascal Pair Mass

The useful correction to the current LRC14 language is to keep three equality
notions separate.

This note is a Pascal/pair-mass companion to mainline HYP-3091's
three-sameness fiber on the lonely set, HYP-3092's verified pair-normalized
cap mass, and HYP-3093's equivalence-triad audit.  HYP-3091 gives the fiber
`Phi(S)=(covering | D,1/lmax,arc-spectrum | meas)`; HYP-3092 verifies the
cap-side pair mass; HYP-3093 names the forgetting-cost protocol.  This note
asks how the row-14 counts, cap defects, and Farey pair-mass coordinates enter
that fiber.

Equidistribution is the analytic shadow: a large speed sweeps through a seed's
safe set and removes the expected fraction.  This is the Node-3 story and the
Weyl/single-far story.

Equinumerosity is the cardinal shadow: all-base counts, Pascal slices,
Burnside-style counts, and row sizes such as `3003,2002,1001`.  It tells us how
many objects are present, not whether they carry the same proof obstruction.

Equidecomposability is the scissors shadow: after two rows have the same count
or the same cap, do they still have the same sector-pair avoidance matrix,
endpoint owners, component structure, level-7 lift status, and Farey
`(p,q),p+q,p*q` address?  HYP-2187 made this distinction for tournaments using
`H`, `beta1`, and packet data.  LRC14 needs the same discipline.

The prompt numbers expose a compact invariant:

```text
C(14,2) = 91
1001 = 11*91 = C(14,4)
2002 = 22*91 = C(14,5)
3003 = 33*91 = C(14,6)
4004 = 44*91 = 2*C(14,5) = C(14,4)+C(14,6)
```

The first three are real row-14 Pascal entries.  The fourth is not a row-14
entry; it is the affine completion of the `11,22,33,44` pair-mass scale.  That
is exactly why it is interesting rather than suspicious.  It appears in
`cap_9 = 45/91 - 1/4004`, while HYP-3090/HYP-3092 say
`cap_k=C(k+1,2)/91` is exact for `k>=10`.  So `4004` is the one-pair-mass
denominator where the clean triangular cap law first fails.

This suggests a new invariant family:

```text
triangular cap shadow  = C(k+1,2)/C(14,2)
Pascal pair mass       = row-14 subset count / C(14,2)
sector scissors packet = pairwise sector co-emptiness + endpoint owners
cap defect             = triangular shadow - actual cap or p0 bound
```

The fundamental object is not the cap scalar.  It is the decomposition that
explains when the cap scalar can be obtained by equidistribution, when it is
only a Pascal/equinumerosity count, and when a Dehn-like scissors defect
survives.

The most concrete test is small: add `pascal_pair_mass_unit`,
`triangular_cap_shadow`, `cap_defect`, `sector_pair_scissors_signature`, and
`farey_additive_lane_mod_91` to a packet audit.  Then compare rows with the
same base count or cap shadow.  If the route outcome splits by sector scissors
but not by count, this HYP is doing real work.
