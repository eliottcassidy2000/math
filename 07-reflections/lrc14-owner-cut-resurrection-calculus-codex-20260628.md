# LRC14 Owner-Cut Resurrection Calculus Reflection

Date: 2026-06-28

This pass turns the strongest HYP-3409/HYP-3410 obligation into a finite
certificate language.  The useful graph is the cross-exit clause hypergraph:
rows are vertices, different-exit row pairs are obstruction clauses, and
endpoint-owner labels are the cut coordinates that hit those clauses.

The key correction is that the visible data does not support a singleton-owner
theorem.  The height leak and `13->26` owner leak have singleton cuts
`5:g1` and `1:g1`, but the newer `10->20` frontier has minimum cut size `3`,
five minimum cuts, and empty core.  A proof route should therefore target a
bounded owner-transversal theorem with an exit-pure cut code, not a universal
label.

Rebase integration: incoming HYP-3411 reframes the residue/equioscillation
half as `C6`-invariant and locates the hard core in the magnitude layer where
the covering flex breaks symmetry.  Incoming HYP-3412 tests compressed
special-function sidecars on the same expanded-bank leaks, and HYP-3413 ties
the Goddyn-Wong doubling switch to `q == 1 mod 3`.  This owner-cut calculus is
the finite certificate language for that broken layer: once the residue/Galois
quotient stops being exit-pure, the endpoint-owner clauses say which owner
coordinates must be retained.

The Menger/Schwarz/BDH/Bring reframes become proof-facing in this order:

```text
Menger      -> minimum owner transversal
Farkas      -> finite cross-exit clause current
SC          -> owner/accessory labels missing from turn words
BDH         -> variance prefilter, never final certificate
Bring       -> theorem-exit branch alphabet made single-valued on quotients
charal      -> child-deck stability of transversal number and terminal exits
```

Next concrete tests:

1. Extend HYP-3406 beyond `(72,20)` and stop at the first
   `residue+owner_support` mixed fiber, if any.
2. Run the owner-clause calculus on the first failure and record whether the
   minimum transversal remains bounded, increases, or needs a new sidecar.
3. Build endpoint deletion, mirror-swap, and `+14` child decks for the three
   known fibers and test whether transversal number and cut-code purity are
   stable.
4. Add terminal-exit routing after every cut certificate.  Owner cuts separate
   exits; they do not by themselves prove positive open mass or boundary
   equality.

Assumption challenge: tournament vertices were considered as runners, gaps,
fixed sections, section boundaries, wall-crossing events, residues, cover arcs,
Fourier modes, matroid circuits, owner labels, and proof obligations.  The
selected vertices are proof obligations because they preserve the predicate at
stake: theorem-exit purity on quotient fibers.  This destroys raw runner
identity, row order, irrelevant owner labels, and scalar/special-function
analogies, which is acceptable only because the cut certificate records exactly
which owner labels survive.
