# LRC14 Moser/Fibbinary Automatic Gap Carrier

The useful conclusion is negative first: Moser-de Bruijn is not a free
`2`-adic LRC14 carrier.  It is a base-4 / even-bit-position carrier.  The moment
we apply the operation that the `2`-adic Littlewood paper makes salient,
`n -> 2n`, every nonzero Moser value leaves the language.  The missing datum is
small but essential: the parity of the bit position.

Fibbinary behaves differently.  Its no-adjacent-1 DFA is closed under
`n -> 2n`, and its failure mode under `n -> n+1` is exactly a carry-boundary
failure.  That makes it a better match for the Zeckendorf/path-normal-form
thread: it is not a scalar density estimate, but a local finite-state memory
for when carries are allowed to pass.

This changes how I would use the prompt's four outside ideas:

```text
2-adic Littlewood / Hurwitz by 2:
  ask whether a proposed carrier survives multiplication by 2.

Fermat-Catalan:
  treat as a p-adic valuation gate, not as a license to import arbitrary
  sparse automatic sets.

Ostrowski-Hadamard:
  use as a gap warning for primitive atoms, but do not confuse atom lacunarity
  with the sorted value language.

Moser/fibbinary finite automata:
  retain DFA state as packet data; never use raw membership counts alone.
```

The ACM DOI in the prompt is also a useful discipline check.  It is a SODA
paper on approximating large sticks and convex bodies inside polygons.  In LRC
terms that is a safe-component/inner-approximation analogy.  It says nothing
about the p-adic or automatic closure we need here, so I kept it as the
`stick_potato_safe_geometry` carrier below the arithmetic/automaton carriers in
the tournament.

## Assumption Challenge

I explicitly considered vertex sets based on runners, gaps, fixed circle
sections, section boundaries, wall crossings, residues, cover arcs, Fourier
modes, matroid circuits, proof obligations, automatic-language states, and
valuation gates.  The selected vertex set is proof-language carriers because
the local predicate is "does this quotient preserve the data needed by the LRC
packet route?", not "which runner wins?"

The quotient preserves only an automatic/valuation side condition:

```text
2-adic phase closure
carry-boundary state
eventual-periodic recurrence gate
```

It destroys:

```text
endpoint-owner geometry
exact safe interval length
Farey binding scale
C27/K33/state-lift labels
```

This is why the result is a guardrail rather than a proof.  A future HYP-2963
extension can attach an `automatic_gap_carrier` field, but the field must be
checked for route mixing.  If the same automatic fiber contains Fejer, endpoint,
C27, K33, and F7 routes, it is only a weather report.

## Proof Route Consequence

The next usable theorem shape is:

```text
Any sequence-shadow quotient used inside the LRC14 proof must either
  (a) retain product automaton state,
  (b) be closed under the relevant 2-adic shift,
  (c) expose its carry-boundary debt, or
  (d) route to an explicit valuation/SML residual sector.
```

Moser alone satisfies none of `(a)-(c)` for `n -> 2n` unless the even/odd phase
state is kept.  Fibbinary satisfies `(b)` and exposes `(c)`.  A genuine
Fermat-Catalan/SML import would have to satisfy `(d)` by naming an ultimately
periodic index coordinate, not by pointing to a sparse automatic value set.
