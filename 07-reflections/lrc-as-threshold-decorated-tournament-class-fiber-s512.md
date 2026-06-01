# LRC as a threshold-decorated tournament-class fiber

The user's hypothesis is right in the direction that matters: LRC should be a
statement about which tournament classes a runner system can realize.  S512
shows the necessary refinement.  It is not a statement about raw A000568
classes alone.  It is a statement about a decorated fiber over A000568.

The base layer is still A000568.  A tournament class is a possible complete
binary relation shape.  The half-turn runner clock maps positions on the circle
to a walk in that class space.  In open cells, the image is very restricted:
for total `n=5`, the open half-turn clock sees only four of the twelve
A000568 classes.  That is the old S24 circular-menu insight in sharper form.

The first twist is that LRC is often tight on equality walls.  Once the fixed
tie Hamiltonian path is included, the boundary compactification can add many
classes: for total `n=5`, S512 sees eleven of twelve classes.  So the object is
not just the open circular menu.  It is the circular menu plus boundary strata.

The second twist is more important.  Raw iso class is too coarse for loneliness.
For total `n=3` and `n=4`, the unmarked phase class, observer-marked phase
class, gap-rank class, and pair-deficit-rank class all have mixed fibers:
the same class can occur at good and bad LRC states.  That kills the naive
proof route.

The route that survives is threshold decoration.  Put colors on the fiber:

```text
ordinary tournament class
+ observer mark
+ which observer-adjacent gaps are >= 1/n
+ which pair-cells have zero danger deficit
```

With those colors, the small cases become class-avoidance statements.  For
total `n=3`, `gap_threshold_fiber` has two good-only classes and certifies all
`79/79` primitive speed sets sampled through max speed `16`.  For total `n=4`,
it has five good-only classes and certifies all `109/109` primitive speed sets
sampled through max speed `10`.  The pair-cell threshold fiber also certifies
all sampled sets.

This gives a proof template:

```text
runner speeds
 -> finite wall/cell walk
 -> decorated tournament-class fiber over A000568
 -> cannot avoid good-only classes
 -> LRC witness
```

For `n=14`, enumerating A000568(14) is not the plan.  The plan is compression.
Use the restricted circular/pair-cell/endpoint-pressure fibers already exposed
by S506-S509.  The hard quotient-ladder rows should live in a tiny subregion of
the decorated fiber, with dyadic/odd/product-sum labels telling us which
operation-grid branch we are on.  A proof would show that every allowed walk in
that subregion hits a good-only threshold fiber before the clock closes.

The cleanest future test is to compute, for the known `n=14` hard rows, the
sequence

```text
phase class
gap-threshold class
pair-deficit-threshold class
endpoint-pressure class
```

at every wall and midpoint of the positive corridor.  If all bad states lie in
a few mixed fibers and the corridor is forced to exit through a good-only fiber,
then the A000568 analogy has become a usable LRC proof mechanism.
