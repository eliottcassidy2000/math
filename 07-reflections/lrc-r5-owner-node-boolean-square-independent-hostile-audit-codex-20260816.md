# Independent hostile audit: a five-state Gray path leaves a genuine Boolean square

**Status: FINITE-EXACT INDEPENDENT AUDIT ACCEPTS A SCOPED RANK-GREATER-THAN-ONE
SOURCE-SUPPORT REFINER.  LRC(14) remains open.**  Neither the submitted
Boolean-square script nor its common-owner parent was imported.  The audit
rebuilt the THM-2471/THM-2594 source service and actual THM-3514 endpoint
indicators, formed a new common-coordinate event sweep, and lifted `Q(13y)`
directly.  Its exact arrays and all submitted digests agree.

## Corrected carrier before any Fourier transform

MISTAKE-417 supplies the first gate: a refiner implied by `OWNER` is only a
delta-cell lift.  The source support is genuinely finer.  On the complete
base its realized support states, in cyclic order, are

```text
{0} -> {0,6} -> {6} -> {6,12} -> {12}.
```

Their toggle word and exact measures are

```text
toggles  = (6,0,12,6),
measures = (1,12,2,12,1)/28.
```

The putative sixth cut `{0,12}` never occurs.  Thus the six nonempty proper
subsets of `{0,6,12}` are an abstract cut completion, not the realized
carrier.  Every one of the five realized states nevertheless has the exact
partial-tournament census `(missing,one-way,two-way)=(1,2,0)`.

The owner window intersects four of these state types and excludes only the
middle type `{6}`.  Each visible state has one connected owner component and
exact measure `1/28`:

```text
{0}, {0,6}, {12}, {6,12}.
```

Complement pairs the first with the last and the middle two with each other.
With bits `(owner component, support multiplicity)`, complement is exactly
XOR with `11`.  This proves the four-state Boolean square over the rational
interval geometry before finite-field evaluation.

## Independent common-base construction

The source grid and endpoint residual grid embed in

```text
C = 9684279613402457983920,
L = 13C = 125895634974231953790960.
```

The audit independently checks the complete Lucas certificate for

```text
p=755373809845391722745761,
g=23,
rho=g^6=148035889,
```

and obtains the order-thirteen root
`266737884585332483769981`.

Its integrator differs from the candidate.  It materializes the thirteen
inverse branches of `Q(13y)` as exact intervals on `C`, sweeps the endpoint
sheet mask together with every source-profile breakpoint, and evaluates the
phase with exponent `742586*x`.  Five fully guarded rows independently match
the delete-then-restore guard construction.  The same-root term remains zero
pointwise.

Summing the four states recovers the corrected rank-one parent gamma hashes

```text
b5246eb2...a073d,
20c83a78...3e258.
```

This is an exact marginal identity, not merely agreement of inverse values.

## Rank gate and genuine mixed support

The state-by-residue matrices have ranks

```text
coupled = 4,
source-weight erasure = 4,
doubly centered = 3.
```

Rank three is maximal after removing the constant state direction.  Therefore
this object cannot be another separable boundary delta.  Both binary
marginals independently have rank two, so neither state bit is determined by
the other.

Only after those geometric and rank checks does the Walsh/Fourier census
become meaningful.  It is

```text
(total,DC,V4-axis,F13-axis,mixed)
  coupled = (52,1,3,12,36),
  erased  = (52,1,3,12,36),
  ANOVA   = (36,0,0,0,36).
```

All `3*12=36` centered mixed modes are nonzero.  At relation class `(1,0,6)`,
all four state values and all four Walsh channels are nonzero.  The trivial
Walsh channel equals `317699132065964946247468`, exactly the surviving
one-dimensional parent value.

The source-erasure control is important: unit numerical source weights still
give rank four.  Hence the result proves that the typed source-support
partition is load-bearing, not that the detailed source weights cause the
nonseparability.

## Tournament and recurrence interpretation

The natural five-object structure is a Gray path in the subset cube, not a
five-vertex tournament.  Each vertex is itself an oriented cut on the
three-root spine.  The visible four-state quotient is a Boolean square because
the owner removes the central path vertex and complement closes on the four
remaining types.  No Fibonacci, ternary-tree, or tournament recurrence follows
from this finite path alone.

This distinction also gives a useful sidecar criterion.  A proposed ancestry
coordinate is nonredundant only if, after adjoining it to the Boolean square,
it increases a rank, splits a state fibre, or changes a conditional spectrum.
Merely recovering the Gray-path order or one of the two square bits is
redundant.

## Connection contract and boundary

| field | exact answer |
|---|---|
| source | THM-2471 source root-service support and actual THM-3514 endpoint indicators |
| target | four owner-visible support states times thirteen relation residues |
| map | pointwise support -> `(component bit,multiplicity bit)` before integration |
| preserved | common base, source weights/roots, endpoint sheets/guard, complement, residue phase |
| destroyed | inverse ancestry sheets, exact address orbit, source/arrival distinction, `U_clock` |
| hostile | MISTAKE-417 occupancy/rank gate, source-weight erasure, both bit marginals |
| survivor | rank `4`, centered rank `3`, all `36` conditional mixed modes |

The result is one canonical `r=5` host computation.  It is not an exact
`C(a;X,m)`, source/arrival transport, `U_clock` chronology, uniform-row
statement, physical current, row exclusion, or LRC(14).

## Reproduction

```text
python -B 04-computation/lrc_r5_ufull_owner_node_boolean_square_independent_audit_20260816.py
python -B -O 04-computation/lrc_r5_ufull_owner_node_boolean_square_independent_audit_20260816.py
```

Normal and optimized transcripts are byte-identical.  The semantic digest is
`b0996af3f1760b2118187490c93e0e01b322cb57fd6b25d3bf3688778b7e664c`.
