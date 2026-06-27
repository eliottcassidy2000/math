# LRC14 Danger-Count Moment Duals

This pass deliberately avoided the last lift-packet angle.  The new object is
the integer danger count

```text
N_S(t)=#{s : ||s t|| < 1/14}.
```

A counterexample has `N_S(t)>=1` everywhere.  That turns LRC14 into a finite
dual problem on the count alphabet `0..13`: any polynomial `P` with `P(0)=1`
and `P(n)<=0` for `n>=1` gives `safe_mu >= E[P(N)]`.

The computation is exact.  It sweeps all rational danger endpoints, records the
Haar distribution of `N`, and searches exact factorial-moment duals.  This is
not a low-degree miracle: AP/GW and the named hard rows are not separated by
degrees up through `6`.  But degree `8` or `9` starts to see the positive safe
mass in the hard rows, while AP and GW remain nonpositive until the full
degree-13 interpolant gives exactly `0`.

The useful result is the shape of a different theorem:

```text
positive degree <=9 count-dual expectation
  or AP/GW equality
  or C27/K33/HYP-2908 state lift.
```

This route forgets where the safe interval is.  That is both its strength and
its weakness.  It cannot replace endpoint-owner or lift-packet proofs, but it
can act as an upstream gate: if the count distribution already certifies
`safe_mu>0`, no boundary-gap or aperture work is needed.

Assumption challenge:

```text
Considered vertices:
  runners, danger arcs, exact endpoints, safe intervals, lift packets,
  danger-count states, factorial moments, polynomial duals, Fourier modes,
  boundary owners, proof obligations.

Chosen vertices:
  danger-count states and moment-dual certificates.

Preserved:
  cover predicate N>=1, exact Haar distribution of N, AP/GW equality
  visibility, and low-degree dual proof obligations.

Destroyed:
  safe interval location, endpoint ownership, Q/R split, and lift geometry.
```

So the missing theorem is no longer "find the interval."  It is "prove the
right finite family of dual polynomials exists on every labelled packet
outside the equality atoms."
