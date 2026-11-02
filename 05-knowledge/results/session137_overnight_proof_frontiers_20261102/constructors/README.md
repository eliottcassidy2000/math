# Actual-source geometry and exact constructor barriers

**The unrestricted odd Hamiltonian-path spectrum remains OPEN.**

[THEORY.md](THEORY.md) proves:

- Actual marked-source states generate the full integer lattice and full
  nonnegative cone; their closed convex hull is exactly `v>=0, g+t>=1`.
  Consequently every universal linear four-state operation has a nonnegative
  integer matrix. This closes S134's general fresh-block positivity caveat.
- [BOUNDED_INTERIOR.md](BOUNDED_INTERIOR.md) strengthens the face result:
  fixing **both h and L** gives a finite actual source bank followed by a zero
  tail, with exact output `N=2^k q0+h-2L`. Thus bounded height and interior have
  only `O(log X)` outputs, even with variable radices and phases. Exact dyadic
  face bases at h=1,3,5 and the first positive-interior fibre are supplied.
- An entire right-local migration history retains its original marked-port
  cost. If `q>=N`, any output at most N must be an attained multiple of the
  original F. H-neutral padding can therefore destroy every nonmultiple lift.
- A source-cost overflow quotient decides arbitrary finite right-local recipe
  banks at all depths. It independently recovers the inherited cap173 language
  using 2,592 states instead of 6,234 positive-saturation states.

All claims retain actual graph recipes. Positivity, parity and the convex hull
are expressly not realization criteria. The radix-80 actual family survives
by growing its parent height; no universal integer acceptance is claimed.

Run `python3 -B check.py`, `python3 -B face_bases.py`, and
`python3 -B bounded_interior.py`, also under `-O`.
See [VERIFICATION.md](VERIFICATION.md) for exact universes and controls.
