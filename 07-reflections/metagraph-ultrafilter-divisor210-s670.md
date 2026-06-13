# Metagraph Ultrafilter and Divisor 210, S670

The productive correction is simple:

```text
not "the metagraph is an ultrafilter";
yes "the tiling cube has ultrafilters whose descent defects are the metagraph problem."
```

The divisor lattice of `210=2*3*5*7` is the right toy because it is exactly a
Boolean lattice.  Choosing one prime gives a principal ultrafilter: the upper
set is "divisible by p", the lower set is "not divisible by p", and every
complement pair `d <-> 210/d` is decided.  The Hasse diagram has 12 lower
internal edges, 8 crossing edges, and 12 upper internal edges for every prime
choice.

That is the clean version of the black/blue intuition.

The tournament tiling cube has the same raw behavior.  In a fixed Hamiltonian
path frame, the non-path tiles form a Boolean cube `Q_m`, and a principal tile
coordinate decides every complement-tiling pair.  This is also where the
MISTAKE-033 correction matters: blue/black means flip all tile bits inside
`Q_m`, not all arcs of the tournament.

But quotienting is brutal.  At `n=4`, every principal coordinate already has a
mixed tournament isomorphism class.  At `n=5`, the best coordinate still has 5
mixed classes and leaks 30 of 32 complement pairs.  At `n=6`, the best
coordinate has 35 mixed classes and all 512 complement pairs leak at quotient
level.  So the side choice is exact before quotienting and mostly invisible
after quotienting.

This makes the ultrafilter idea better, not worse.  It tells us precisely what
the proof obligation is: descend the side choice.  If it does not descend, add
the missing address.  That is exactly HYP-2240/HYP-2241 in another dialect.

This also reabsorbs S669.  Equinumerosity is too weak: fixed-path tilings and
labelled degree-even graphs have the same size, but degree-even isomorphism
classes are not tournament classes.  Royle-even graphs are the deeper count
match, but still need a predicate-preserving fiber functor.  The ultrafilter
language says the same thing in order-theoretic form: cardinal equality is a
shadow; upper/lower membership is proof data.

For LRC14, the most promising object is not a runner tournament but an
owner/carry/proof-obligation filter.  Vertices can be residues, cover arcs,
endpoint protectors, D/U/N obligations, fixed circle sections, wall events, or
carry cocycles.  A good quotient is one where a side choice remains pure:
floor rows on one side, strict rows on the other, with AP/Vstar/2AP handled as
named boundary atoms.

The wild next move is to search for quotient-pure filters.  Principal tile
filters leak on ordinary tournament isomorphism classes, but maybe owner-paired
filters, residue-shell filters, or isomorphism-invariant tile-family filters
descend.  That would turn "black line / blue line" into a certificate rather
than a picture.

The phrase I would hand to the next agent:

```text
Find the ultrafilter that descends.
If none descends, attach the smallest address that makes it descend.
```
