# Convolution lifts and hidden factor grids

The most useful correction to the coefficient-tiling picture is this:

```text
The coefficient row is not the object.  It is the shadow.
```

If

```text
f(x)=g(x)h(x),
```

then every visible coefficient is a diagonal sum from a hidden grid:

```text
a_k = sum_{i+j=k} b_i c_j.
```

That is a much sharper version of the tiling model.  The triangle of
coefficient signs is the boundary readout.  Reducibility asks whether a
nontrivial rectangle can live behind it.

This also cleans up the role of the fixed Hamiltonian path.  The path is not
where irreducibility lives.  It is the arbitrary deterministic order in which
we inspect split obligations:

```text
1 x (d-1), 2 x (d-2), ...
```

or in which we rank local primes:

```text
2, 3, 5, 7, ...
```

That is exactly the kind of tiebreak the repo keeps finding: necessary for a
finite algorithm, but not the invariant itself.

The computation was encouraging.  In the same degree-4 coefficient scout where
signs were completely mixed, adding the least small mod-p convolution blocker
certified `3058/3096` irreducibles with no false positives.  So the first
strong carrier is simply:

```text
which split rectangles survive after local reduction?
```

Then the caution appears immediately.  Eisenstein-style rows can look
completely split mod p, while the p-adic Newton hull proves irreducibility.
The residue face sees the grid after flattening coefficients into `F_p`; the
valuation face sees the height profile of the same grid before flattening.
They are not rivals.  They are orthogonal projections of the hidden lift.

That is why the next state should not be "a polynomial has sign pattern S" or
"a denominator q is blocked."  It should be:

```text
local gate -> surviving split obligations.
```

For polynomial irreducibility, the gates are residue primes, Newton polygons,
Cohn digit addresses, Singh value depths, and recombination traces.

For LRC14, the gates should be denominator/resource fibers, carry owners,
Pisano classes, and Bprime/private-owner deletion.  A counterexample-like row
is not merely a scalar blocked row.  It is a row whose survivor set refuses to
empty across all local projections.

The sentence worth keeping:

```text
Signs show the visible tiling; split survivors show the hidden factor grid.
```
