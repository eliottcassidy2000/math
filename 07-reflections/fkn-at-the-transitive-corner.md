# FKN at the transitive corner

The first clean lesson is negative: `H` is not globally dictator-like on the
tiling cube. The degree-1 Walsh share drops from `0.75` at `n=4` to `0.21` at
`n=7`. So if I ask for a literal Friedgut-Kalai-Naor theorem with the full
staircase cube as the ambient Boolean space, the model pushes back immediately.

But the second lesson is exactly the one the user wanted. The transitive tiling
is still a free reference state, and it has a unique strongest escape direction.
Among all single tile flips, the extreme tile `(n,1)` dominates every time. In
the shell-1 table, strip `s` gives `H = 2^s + 1`, so the largest strip
`s=n-2` is the unique best move. That is not a global dictator theorem; it is a
corner dictator-shadow theorem.

This is the right way to phrase the analogy. The tournament is not free. The
first layer only tells me which perturbation wants to grow fastest. The actual
interaction enters through the deletion-lattice Möbius defect

```text
mu_H(T) = Σ_{U subset V} (-1)^(n-|U|) H(T[U]).
```

That quantity is exactly the recursive version of "subtract off all lower-order
subtournament shadows and see what remains." It turns the user's
`A+B+C-D-E-F+G` sentence into a precise operator.

The parity split in `mu_H` is the best new surprise:

```text
near-transitive mu_H:
n=4 -> 0
n=5 -> 2
n=6 -> 0
n=7 -> 2
```

So the most FKN-like direction is not itself free of interaction. It carries a
small, parity-sensitive irreducible defect. That feels exactly right: the board
predicts the dominant escape coordinate, but the real tournament remembers that
the perturbation is creating a cyclic portal, not just toggling an independent
bit.

The conceptual reframe I want to keep is:

- the full cube is the wrong scale for the stability theorem;
- the bounded-radius ball around the transitive corner is the right scale;
- the first layer chooses the candidate junta coordinate `(n,1)`;
- the Möbius defect measures when the interacting tournament departs from the
  free assembly of its subtournaments.

That makes the next proof target much more believable:

> among radius-1 perturbations of the transitive corner, near-transitive is the
> unique maximizer of Hamiltonian growth; among larger low-radius perturbations,
> the obstruction to a pure junta picture is exactly the Möbius interaction
> defect.
