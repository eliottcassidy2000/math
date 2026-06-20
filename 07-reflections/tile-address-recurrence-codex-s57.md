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
new half-address crossing coordinates at n: floor((n-1)/2)
```

The discarded side is not lost.  Its tile address is reflected by

```text
(beta,tau) -> (n+beta-tau, 2n-tau).
```

This gives a precise address for routing paired contributions.

Precision: the half-address coordinates are literal one-bit coordinates only
on the grid-symmetric/self-converse subcube.  For the full complement quotient,
nonfixed mirror pairs still carry unordered pair states.  The address quotient
is therefore a routing compression first; it becomes a bit compression only
after imposing grid symmetry.

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

## First DP Check

The follow-up scaffold `half_tiling_address_dp_codex_s57.py` enumerates the
grid-symmetric half cube directly.  It computes `c3`, score multisets, and `H`
for `n=3..8`.

The useful displayed point is `n=8`:

```text
half coordinates = 12
assignments = 4096
Hmax = 661
c3 at displayed leaders = 20
scores = (3,3,3,3,4,4,4,4)
```

This is a self-converse subcube computation, not a full complement-quotient DP.
