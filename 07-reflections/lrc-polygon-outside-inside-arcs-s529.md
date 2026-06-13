# The polygon outside remembers the simplex inside (S529)

The user's distinction is the right one:

```text
simplex view:  tournament = one binary relation on each edge
polygon view:  outside cyclic gaps -> hidden inside chord orientations
```

The simplex picture makes `C(m,2)` arcs look free.  The regular-polygon picture
shows why LRC arcs are not free.  Once the points lie on a circle, every inside
chord is a consecutive sum of outside gaps tested against the half-turn.  Rotate
or reflect the outside necklace and the inside tournament rotates, reflects, or
opposes with it.  The dihedral group is the visible symmetry; the hidden arcs
are its chord-length channels.

For LRC the observer is a clasp vertex on the outside polygon.  Source is not
primarily an interior tournament condition.  It says the two outside gaps around
the clasp are both at least `1/n`.  The interior runner tournament is important,
but it is downstream from that outside clasp geometry.

This is clearest at the initial AP witness.  At `n=14`, speeds `{1,...,13}` and
time `t=1/14` put the runners on every non-observer vertex of the regular
14-gon.  The outside is perfectly simple: every gap is `1/14`.  But deleting the
observer leaves a large hidden interior:

```text
outside gaps: 14
inside runner chords: 78
D_14 channels after deleting clasp: 12,12,12,12,12,12,6
hidden diameter ties: 6
```

Those six diameter ties are the compactified boundary.  They are why the regular
polygon is not just a pretty extremal picture; it is the place where the outside
necklace stops deciding a strict inside tournament without a tie path.  The
boundary-source classes live exactly there.

So the LRC proof shape should be:

```text
outside: force the clasp gap to open
inside: account for hidden chord-channel debt when the clasp only closes on a wall
```

This ties HYP-1998, HYP-2001, and HYP-2003 together.  The round/A000016 body is
the open polygon body.  The n=14 regular-polygon AP row is a closed wall whose
hidden inside diameter debts are perfectly balanced.  Proving n=14 by this lens
means showing that no non-AP sweep can keep the outside clasp closed while
circulating those hidden endpoint/diameter debts forever.
