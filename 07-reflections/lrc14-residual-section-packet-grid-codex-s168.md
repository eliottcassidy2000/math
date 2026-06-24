# LRC14 Residual Section / Packet Grid Reflection

HYP-2996 makes the "residual section" language operational.  The important
shift is that F7 is no longer the complement of a checklist.  It is now a
section predicate: a packet must be hard, non-AP/GW, zero-open, and not
discharged by an owner-strip, cross-handoff, or nested-refinement exit in the
Haar-product grid.

On the default HYP-2963 bank, that predicate is empty.  There are `7235` hard
non-AP/GW packets, and every one has a named grid exit: positive open-Haar
owner strip, C27/unit-petal owner strip, covering nested refinement, or
K33/THM-572 cross handoff.  AP/GW are the only same-tile boundary cohomology
atoms.

This is still bounded evidence, not closure.  The next proof step is not to
claim "no F7" globally; it is to replace the bounded rows with family templates:

1. q-witness packets as exact `0`-cochain exits.
2. AP/GW as same-tile boundary cohomology.
3. C27/unit petals as owner-strip coboundaries.
4. K33 as cross-handoff state-lift debt.
5. Covering rows as nested boundary-moment descent.

Only after those templates exist should Fejer/Ramanujan certificates be grouped
by grid section.  That grouping is the likely way to turn the packet-grid audit
into a family compression engine rather than another finite table.
