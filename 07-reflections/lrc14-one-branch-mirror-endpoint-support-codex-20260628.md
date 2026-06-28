# LRC14 One-Branch Mirror / Endpoint-Support Reflection

Date: 2026-06-28

The new angle is that HYP-3425's two-branch obstruction is hiding a simple
mirror symmetry.  Since the map `u -> 1-u` preserves the even-safe set and
swaps the branch-0 near-integer condition with the branch-1 near-half
condition, branch 1 is not an independent proof burden once branch 0 survives.

So the proof target can be compressed to a one-color interval cover:

```text
E_safe not subset B0_odd.
```

This is still hard, but it is cleaner.  It changes the shape from a two-color
bad-core cover into a standard one-dimensional interval-piercing problem.

The endpoint-owner audit is the real new leverage.  On `162` rows every
branch-0 survivor has labelled rational endpoints, and no survivor needs more
than three speed owners to explain its two boundary walls.  The tight canonical
row is even sharper: every survivor is just one odd owner plus the even owner
`42`.  That suggests a proof by local endpoint-owner triples, or by showing
that a full one-color cover would force an impossible alternation word on an
`E_safe` component.

Next pass: classify which owner triples can bound a survivor and which triples
would be required to destroy all survivors.  If the triple list closes, the
Helly lemma becomes a finite endpoint-owner theorem instead of a global mass
estimate.
