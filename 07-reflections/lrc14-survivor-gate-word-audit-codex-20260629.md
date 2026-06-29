# LRC14 Survivor-Gate Word Audit Reflection

HYP-3438 changed the shape of the next LRC14 proof step.  HYP-3436 had already
shown that the two-color bad core

```text
E_safe cap B0_odd cap B1_odd
```

does not swallow the audited rows.  The new audit asks a sharper local question:
what is the exact word around each survivor gap inside a mixed `E_safe`
component?

The answer is less corridor-like than expected.  Across the `135`-row bank the
dominant survivor proof carrier is an edge gate, not an interior `B-S-B`
corridor:

```text
survivor_gates=8702
edge_singleton_parent_gate=5950
edge_survivor_residual=1080
owner_current_small_delta=794
mixed_owner_residual=878
```

An `edge_singleton_parent_gate` is a boundary survivor adjacent to exactly one
`(1,1)` bad-core block and touching the parent even wall.  That is important:
the canonical family is not asking us to prove a generic interior corridor
lemma first.  It is asking us to prove a boundary-gate lemma.

The incoming branch-mask ledger sharpens this rather than competing with it:

```text
survivor_branch_mask_hist={both:1064, branch0:3819, branch1:3819}
gate_adjacency_hist={left_bad_edge:3515, right_bad_edge:3515, two_sided:1672}
```

The symmetry of branch0 and branch1 says the residual problem is not a one-sided
branch artifact.  The hard cases are cover-delta and endpoint-gluing objects:
they need a cut, owner-current, endpoint-spine, or exact-period explanation.

The real creative payoff is the mod-35 canonical split.  For

```text
S_m = {1,2,3,4,5,6,7,8,9,10,11,13,84m}
```

the one-period probe `m=1..35` found:

```text
outer 13/7 edge gates  <=>  7 does not divide m
inner 11/5 edge gates  <=>  5 does not divide m
35 divides m           =>  no mixed gates; survivor is clean E_safe
```

So the old “canonical singleton tail” story refines into a residue theorem.
The factors `5` and `7` are not just numerology: they name which bad-core owner
pair is available as a boundary certificate.  When both owner pairs disappear
at `35 | m`, the mixed gate problem itself vanishes and the row survives by
clean even-safe components.

This gives a plausible proof decomposition:

1. Prove the canonical mod-35 gate law by exact endpoint arithmetic around the
   parent `E:84m` wall.
2. Use that law as the base case for HYP-3431 corridor-fence and HYP-3433
   endpoint-tail arguments.
3. Feed noncanonical residual gates into HYP-3437 overlap-tax cuts, HYP-3439
   rescue-core compression, HYP-3440 endpoint-cut vocabulary, or the older
   owner-current/two-adic/state-lift exits.

The quotient lesson is also now cleaner.  We considered runners, gaps, fixed
sections, section boundaries, endpoint walls, branch-bad owners, cover-pair
deltas, mixed even components, survivor gaps, owner-current labels, and proof
obligations as possible tournament vertices.  The useful vertices are
survivor-gate proof carriers because they preserve the branch-relocation
predicate while retaining the endpoint and owner payload needed for gluing.
They deliberately destroy raw runner order and scalar survivor mass.

That is the guardrail: a quotient may forget global runner order only after it
keeps the local gate word, adjacent B0/B1 covers, and parent even wall, or after
it names the sidecar that reconstructs them.  This matches the broader repo
lesson from irreducibility, unitals, Faulhaber moments, Pollock defects, and
tournament tilings: the theorem usually fails at the coordinate a quotient
quietly forgot.

Next useful computation: take the `878` mixed-owner residual gates and build
their owner-delta graph.  The target is a finite discharging table:

```text
mixed_owner_residual
-> overlap-tax cut
-> owner-current small delta after one relabel
-> endpoint-spine wall
-> exact-period/state-lift/signed-SPEC exit
```

If that table closes, the LRC14 proof state becomes much smaller: canonical
rows are handled by the mod-35 boundary/clean law, and noncanonical rows are
forced into named residual debt rather than raw survivor-measure optimism.
