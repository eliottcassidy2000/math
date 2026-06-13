# Infinite Go, n+2 Recursion, and Partition Functions, S616

The MathOverflow construction is more than a cute infinite-game gadget. It is
a clean model of how a local finite delay becomes an ordinal boundary state.

Hamkins builds a doomed Go group with infinitely many hopeless cuts. Black can
choose a far-away knob, forcing an arbitrarily long ladder before white finally
kills the group. The group dies in finitely many moves, but black can make the
finite number as large as desired, so the life/death value is `omega`.

The update is the part that connects to the repo: give the main group `n+2`
liberties and black gets `n` serial announcements, hence value `omega*n`.
The two extra liberties are not fuel. They are the terminal boundary that lets
the group be killed at the end. The fuel is exactly:

```text
active fuel = liberties - 2.
```

That is the same shape as THM-291. In the tournament fixed-path model,
`n->n+2` adds two endpoint vertices. Setting the new boundary variables to
zero makes those vertices a source and a sink, so the Hamiltonian-path
partition function collapses back to the inner `n`-tournament. Turning on the
boundary variables adds the correction.

The exact count is the useful new bridge:

```text
new fixed-path boundary variables
  = binom(n+1,2) - binom(n-1,2)
  = 2n - 1.
```

That is precisely the LRC floor modulus `C=2n-1`. So `C=27` at LRC `n=14`
is also the boundary shell of the tournament `14->16` recursion. This makes
the carry/lift problem feel less accidental: a floor resonance is a boundary
state of a partition function.

The right abstraction is:

```text
partition function + retained boundary state.
```

For a finite cutoff `K`, the Go ladder fuel packet is

```text
G_K(q)=1+q+...+q^K.
```

With `r` serial fuel liberties,

```text
Z_{r,K}(q)=G_K(q)^r.
```

Taking `K` large, the scalar `Z(1)` is not the invariant we want. The important
datum is the number of serial geometric factors: the pole order at the critical
boundary, or the tropical/open-game shadow, becomes `omega*r`.

This improves the ongoing partition-function picture:

- Tournament OCF: `H(T)=I(Omega(T),2)`; forbidden `H=7,21` are unavailable
  partition-function evaluations.
- LRC: `Z_delta(y)=sum p_k y^k`; the target is the ground cell `p_0`, and
  `2n-1` shells are boundary observers for all-order cancellation.
- Infinite Go: `Z_{r,K}` records finite ladder choices, while the open-game
  limit reads serial fuel depth.
- Unit distance: arithmetic carriers supply the incidence/energy partition
  function before projection.
- Collatz/two-block: the residue automaton is a partition function over live
  determinant states.

So yes: partition functions are everywhere, but the lesson is not "count
everything." It is "keep the boundary state that tells the count which limit
it is in." Scalarizing too early is exactly how the hard part disappears.

Source anchor: MathOverflow answer 299504 by Joel David Hamkins,
`https://mathoverflow.net/questions/299491/is-there-a-position-in-infinite-go-for-which-the-life-of-a-particular-stone-has/299504#299504`.
