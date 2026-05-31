# Natural/LRC Recursive Modes S378

**Session:** codex-2026-05-31-S378
**Script:** `04-computation/natural_lrc_recursive_modes_s378.py`
**Stored output:** `05-knowledge/results/natural_lrc_recursive_modes_s378.out`

The prompt asked for the old incomplete-tournament idea on natural numbers to
be brought back into contact with product-sum equations and with the Lonely
Runner recursion.  The useful correction is now sharp:

```text
addition shadow       = complete transitive order
multiplication shadow = sparse divisibility skeleton
```

So the additive graph is not interesting as a simple graph after labels are
forgotten.  Its information lives in the two-input fibers.  Multiplication,
unlike addition, keeps structure after projection because divisibility survives
as a sparse DAG.

## Product-Sum As Toy Model

Product-sum equations are the finite laboratory for this projection story.
They ask when a multiplicative seed can land at the same endpoint as an
additive fold after adding enough units.

The S378 table aligns product-sum arity `k` with the LRC speed count `k`, so
`n=k+1` is the LRC threshold denominator.  The minimum product-sum endpoint is
not a smooth function of `n`; it changes with factor-packing type.  Around the
frontier:

```text
n=14, k=13: m(k)=18, seed 1^10 + 2*3*3
n=15, k=14: m(k)=20, seed 1^11 + 2*2*5
n=16, k=15: m(k)=24, seed 1^13 + 3*8
```

That is a small but helpful warning: recursive state needs to remember the
shape of the factor packing, not just the endpoint value.

## LRC Initial Segments

Initial segments are the Dirichlet equality spine.  They are boundary-only for
the exact interval problem, and their unprotected boundary witnesses are the
unit residues modulo `n`.

S378 measured, for `4 <= n <= 22`:

- `phi(n)` and `phi(n)/(n-1)`;
- divisor count `tau(n)`;
- endpoint count and interval count;
- boundary witness count and boundary modulus;
- endpoint peel depth and terminal core size;
- additive gates, multiplicative gates, and divisor edges inside the speed set.

The key feeling: `n -> n+1` is not a gentle deformation.  It jumps when a prime
or composite structure changes the unit skeleton and divisor channels.  Prime
`n` has full unit skeleton.  Composite `n` has a thinner unit layer but more
possible quotient/descent channels.

## Hard Families

The feature vectors give a compact contrast.

The initial `n=14` spine is boundary-only, quotient layer `unit_mod_14`, with
six boundary witnesses.  The seven-ladder near-disproof has positive gap ratio
`0.005411`, five speeds divisible by `14`, `84` boundary witnesses, and a
higher-denominator quotient layer.  It is close in max-gap space but far in
endpoint-state space.

Similarly, the `n=15` ladder and mixed-gate examples have positive gaps and
large unprotected endpoint sets even though they intentionally overload divisor
channels.  This matches S373-S375: speed-first disproof routes protect one
visible layer and open another.

## Proposed LRC-TDA Features

The tournament side learned not to trust one scalar.  `H`, SCC defect, deletion
residue, path homology, support excess, and projection defects all became useful
because they expose different shadows of the same object.

The LRC analogue should expose at least:

```text
max_gap_ratio
forbidden_length
boundary_witness_count
boundary_modulus
unit_skeleton_size = phi(k+1)
divisor_gate_count
nonunit_residue_count modulo k+1
endpoint_peel_depth
terminal_endpoint_core_size
operation closure inside speeds
scalar-ramp distance
missed-cell repair deficit
```

This does not prove the conjecture.  But it changes the object from "a set of
speeds with a gap" into a structured point in a recursive feature space, which
is exactly what helped the tournament project stop staring at raw counts.

## Next Work

1. Add a reusable `lrc_tda.py` feature extractor.
2. Rank S373-S375 near-disproof families by the new feature vector.
3. Combine endpoint-protection incidence features with the S375 exposed-cell
   repair deficit.
4. Search for finite endpoint-protection patterns first, then solve speed sets
   realizing those patterns.
