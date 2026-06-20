# Half-Tiling Parity Address Quotient

codex-2026-06-20-S56

## Summary

The user's half-tiling sequence is exactly

```text
h_n = floor((n-1)^2/4).
```

It is the fixed-line-plus-one-side fundamental domain for THM-280's grid
reflection, which is tournament complement up to relabeling.  The fixed line has
`floor((n-1)/2)` cells, so

```text
h_n = (binom(n-1,2) + floor((n-1)/2))/2.
```

## Recurrences

Even tournament sizes have the two-corner law:

```text
h_n = A+B-C = 2h_{n-1}-h_{n-2}.
```

Odd tournament sizes have the folded seven-slot law:

```text
h_n = A+B-C+D-E-F+G.
```

The `C` and `D` terms have equal size but different positions.  That is the
important geometric fact: cardinality cancels them, but a proof should not.

## Repo Connections

THM-442's full recursion `A+B+C-D-E-F+G` is the third finite difference of the
full staircase.  THM-550 now says the half model is the mirror-folded parity
refinement of that recursion, complementing the concurrent THM-549
complement-quotient theorem.

This should not be applied blindly to Hamiltonian path counts.  THM-442 already
warns that `H` is multiplicative/cycle-space, not cell-affine.  The right next
move is to combine the half-domain with THM-513's interval-root packets: mirror
orbits of roots, fixed-line roots, and their OCF packet effects.

The LRC analogy is also address-level only.  The same seven-symbol syntax shows
up in HYP-2681's three-far packets, but there it is a pair-tax shadow.  The
half-tiling model explains why signed packets recur; it does not turn them into
a direct inequality.

## Useful Warning

Even and odd tournament sizes have different half-carrier geometry.  Any
recursive tournament proof that treats the folded carrier as parity-blind is
throwing away corner data before the cycle-space calculation even starts.
