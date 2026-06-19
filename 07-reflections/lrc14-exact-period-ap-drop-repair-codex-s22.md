# LRC14 Exact-Period AP-Drop Repair

The useful new reframe from this pass is that the AP cap face is not only a
dihedral-mouth object and not only a squarefree-modulus object.  It is an
exact-period packet object.

HYP-2628 said the `Q=210` radical grid misses exactly AP drops `6` and `12`,
while the raw `Q=1260` carrier catches all one-drop cores.  The tempting slogan
was "1260 is bigger."  The computation says something sharper:

```text
drop 6 at Q=1260:  denominators 63, 420, 630
drop 12 at Q=1260: denominators 12, 315, 630, 1260
```

None of these reduced denominators divides `210`.  So the squarefree quotient
does not merely coarsen the two cusp mouths; it deletes their exact-period
packets.

This is a good instance of the user's "large absolute mass but tiny signed
mass" warning.  The safe mouths have positive rational length, but after the
wrong quotient their visible grid mass is zero.  Restoring the raw carrier does
not make the mouths large; it restores the few exact-period atoms that land
inside them.

The most clarifying detail is the asymmetry:

- Drop `6` is still the dihedral `D6` shadow from HYP-2569, but inside the
  `1260` carrier its first visible packets use repeated `3` and repeated `2`
  denominators.
- Drop `12` is the Goddyn-Wong neighbor.  Its dyadic packet `12=2^2*3` is
  completely invisible to the squarefree `210` grid, and all four denominator
  `12` units appear as safe residues in the raw carrier.

The caveat matters: drop `6` has an earlier `q=98` witness outside the Hill
carrier, so this is not a global minimal-denominator theorem.  It is instead a
statement about the order of proof quotients:

```text
raw product / exact-period packets first,
squarefree or coimage projection second.
```

That order now feels less like bookkeeping and more like a proof constraint.
The radical face is allowed only after the AP mouth packets have been accounted
for.  Otherwise the two most important cusp repairs vanish.

For the LRC14 route, this suggests a finite state machine:

```text
endpoint mouth -> reduced denominator -> Euler-copy mask -> coimage class
```

instead of either runner vertices or squarefree masks alone.  Runner vertices
forget the mouth; squarefree masks forget the packet.  Reduced denominators are
the small recursive detail that survives the handoff between geometry and
number theory.
