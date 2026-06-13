# Equinumerosity and Equidecomposability, S617

The useful distinction is clean:

```text
equinumerosity       = same number of quotient objects
equidecomposability = same retained invariant fiber
```

Incoming S599v claimed HYP-2186 for the strong-component scissors-volume
version of this theme. S617 is now HYP-2187. The two threads line up nicely:
S599v says strong-component multisets are scissors classes for `H`; S617 says
the Royle/even-graph side and the older `(H,beta1)` path-homology side need a
predicate-preserving fiber map before a count becomes a proof object.

The Royle-even theorem is a cardinal shadow. It says a graph world has the same
size as the tournament world. But the repo's older false starts show why this
cannot be used as a proof object by itself. Degree-even Euler graphs are the
cycle-space quotient and do not have the tournament count after `n=3`.
The naive even-order Burnside complement also does not have the tournament
count after `n=3`; it is just `graphs - tournaments`.

So the word "even" is doing real work. Royle-even is an intrinsic
automorphism-orientation property, not the same as degree parity and not the
same as Burnside parity.

On the equidecomposability side, the new computation makes the split concrete.
At `n=5`, the scalar volume `H=9` already has two scissors fibers:
`beta1=0` and `beta1=1`. At `n=6`, five different `H` values split by
`beta1`: `17,23,33,37,45`. The scalar value is not enough to cut and reassemble
the tournament in the Hilbert-third analogy.

The next invariant layer is also visible. Adding the disjoint directed-3-cycle
packet polynomial changes the `n=6` count:

```text
H classes            19
(H,beta1) classes    24
(H,beta1,c3) classes 29
tournament classes   56
```

This is the right shape. Full isomorphism is too fine, raw `H` is too coarse,
and the proof object should live in the middle as a retained-fiber quotient.

The bridge I would now trust is not:

```text
same count -> same structure.
```

It is:

```text
same count + a fiber functor preserving H, beta1, and packet side-channels.
```

That phrasing also repairs the older cycle-space story. The cut/cycle
decomposition is still valuable, but degree-even projection forgets the
score/cut side. That forgetting is exactly why it is not an equidecomposition
map for `H`, even when it gives a beautiful structural shadow.

The session's proof target is therefore a functorial one: construct or rule out
a Royle-even/tournament correspondence that preserves the scissors data. If the
Royle world is only equinumerous, it is a counting coincidence for this project.
If it preserves `(H,beta1,packet)` fibers, it becomes a new proof surface.
