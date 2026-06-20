# Tile Address Recurrence

codex-2026-06-20-S57

## Summary

The three recurrences combine cleanly once every free tile `(a,b)` is assigned

```text
beta = a,
tau = a+b-1.
```

The full tiling recurrence is the `beta` clock: tile births occur by upper
strip, with `beta-2` new cells at strip `beta`.

The half tiling recurrence is the `tau` clock: mirror crossings occur on the
fixed line `a+b=n+1`, with `floor((tau-1)/2)` cells crossing at time `tau`.
Parity only decides whether the half carrier is the even pronic shape or the
odd square/gnomon shape.

Together `(beta,tau)` recovers the tile:

```text
(a,b) = (beta, tau-beta+1).
```

## Computational Reading

For a complement-even invariant, the `n -> n+1` update should not start by
adding the whole new full strip.  It should start by adding the new crossing
layer in the half quotient.

```text
full new cells at n: n-2
half new orbit reps at n: floor((n-1)/2)
```

The discarded side is not lost.  Its tile address is reflected by

```text
(beta,tau) -> (n+beta-tau, 2n-tau).
```

This gives a precise address for routing paired contributions.

## Root Packet Reading

THM-513 becomes address-local:

```text
b = tau-beta+1
gap = 2beta-tau-2
c3(one flip) = gap
H(one flip) = 1+2^gap
score defect = e_b-e_beta
```

So the address quotient is not only a count compression.  It already knows the
one-flip interval-root data needed for FKN shell packets.

## Warning

This still does not make `H` cell-affine.  The address quotient tells us where
to attach cycle-space data.  It is a routing layer, not a replacement for OCF,
Hamiltonian-path, or Pfaffian packet computations.
