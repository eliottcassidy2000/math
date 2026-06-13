# H as a loneliness meter, after S506

The useful correction is that `H` is not one loneliness meter.  It is two
meters sitting on the same runner-clock movie.

The low end is unanchored.  In the half-turn circular tournament, `H=1` means
transitive, and THM-374 says this is exactly the open-semicircle state: all
points fit inside a half circle, so the complement has a gap greater than
`1/2`.  This is a real loneliness signal, but it is not the Lonely Runner
Conjecture's marked signal.  It says the pack is bunched somewhere on the
circle and leaves a huge empty arc somewhere else.

The LRC signal is marked and local.  A stationary runner is lonely at threshold
`1/n` when every moving runner is at distance at least `1/n` from it.  More
generally, a vertex in the circular order is locally lonely exactly when the
two adjacent gaps around it are both at least `1/n`.  That pushes toward
regular spacing, so it lives at the high-H end of the half-turn clock, not the
`H=1` end.

S506 made this visible by comparing exact clock cells with marked witness
times.  For small clocks, `H` is strongly anti-correlated with maximum circular
gap: initial `n=5` has Pearson `(H,max_gap)=-0.946`, primes `n=5` has `-0.902`,
initial `n=6` has `-0.912`, and spread `n=7` still has `-0.716`.  This supports
the slogan "low H means a big empty arc."  But the initial LRC witnesses
`t=1/n` for total `n=5..9` have all vertices locally lonely and high `H`
(`15, 41, 175, 629, 3267`).  The same scalar points in opposite directions
depending on which notion of loneliness is being asked.

The `n=14` and `n=18` rows are the clean warning sign.  At `t=1/14`, the
initial row has `H=24104937` and every vertex lonely; at `t=1/28`, it has
`H=1`, maximum gap `15/28`, and no lonely vertices.  At `t=1/18`, the initial
row has `H=115642276825` and every vertex lonely; at `t=1/36`, it has `H=1`
and no lonely vertices.  The ladder times stay high-H but can lose the marked
stationary witness, which means the scalar has to be paired with the safe-gap
mask.

This is where the LRC becomes a tournament-structure problem.  The scalar `H`
is a shadow of a richer object:

- the half-turn tournament from pairwise circular phases;
- the circular order and exact gap vector;
- the marked stationary vertex;
- the safe-gap mask and locally lonely vertex set;
- pressure/deletion graphs recording who blocks whom under perturbation.

So the right methodology is not "maximize H" or "look for H=1."  It is
Tournament Analysis: choose the pairwise observable, choose the binary switch,
record the tournament, and keep the marked data that the projection forgets.
For LRC, the forgotten data is exactly where the proof lives.

The pressure-DAG work now fits into this picture.  A pressure search returning
a DAG is not a failed search; it is a marked tournament certificate with a peel
order.  A pressure SCC would be the counterexample-shaped object.  High `H`
finds the regularized circular environment, while the pressure DAG and
safe-gap mask decide whether the marked stationary demand is actually met.
