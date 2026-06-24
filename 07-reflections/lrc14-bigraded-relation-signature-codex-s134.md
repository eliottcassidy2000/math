# LRC14 Bigraded Relation Signature

S130 through S133 taught us not to replace the denominator.  S134 adds the
older lesson from the summand and multiplicand graphs: do not replace the
relations either.

The summand graph is additive and dense.  It gives the relation shells:

```text
a+b=C,      a+b=c,      a+b=c+d.
```

The multiplicand graph is sparse.  It tests whether those shells are visible to
the LRC observer:

```text
C | w,      a in (Z/CZ)^*,      C=27 unit/nonunit shell.
```

The mistake would be to say "high additive energy" and stop.  AP and shifted AP
show why.  They can have the same raw AP-like sumset profile, but AP's folds
are observer-visible while the shifted AP pushes the same additive shape into
hidden balanced shells.  Same scalar additive structure, different LRC
meaning.

So the useful object is a bigraded relation signature:

```text
(Farey branch, C=27 shell type, visible folds, hidden collisions,
 multiplicand blockers, K_{p,q} incidence)
```

In that language, S133's split becomes more precise.  The `2/27` branch is not
merely "before K33"; it is the old `C=27` antipodal shell with unit/nonunit
gcd strata, plus a two-block multiplicand tag.  The `3/41` branch is not merely
"bigger product"; it is the first three-owner incidence packet.

This suggests a useful proof rhythm:

```text
q/Farey branch first,
C27 shell or K33 incidence second,
visible/hidden relation sign third,
then scalar estimates.
```

That order matters.  If we scalarize first, AP, GW, shifted AP, the `2/27`
petals, and the `3/41` near-miss all throw similar-looking shadows.  With the
typed relation signature, they sit on different ledgers.

The proof hope is now a little sharper: the remaining non-AP/GW atom should
either reveal a C27 shell defect that petal/lift rigidity can kill, or reveal a
three-owner K33 packet with enough sign-visible relation mass to construct the
HYP-2908 state lift.  Still no proof, but the map has fewer foggy rooms.
