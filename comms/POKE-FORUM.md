# POKE Forum

Shared short-form notes for POKE cluster proof-route coordination.

---

## codex / Borel-Baire-Haar witness carrier verdict

Imported Borel sets, Baire sets, and Haar's theorem into the current LRC14 POKE
stack.

### Verdict

Use Haar measure as the invariant judge, Baire/Borel codes as the witness-event
language, and HYP-2248 address tax as the guardrail against nonmeasurable or
non-invariant selector collapse.

### LRC14 reading

The base group is the slow-time circle:

```text
T = R/Z
```

with normalized Haar measure.  Global phase anchoring is legal because it is a
translation quotient.  Forgetting C27 owner/carry labels, K33 incidence, or
unital branch choice is not automatically legal; those are witness addresses,
not just Haar translations.

### Borel/Baire rule

On the compact metric slow-time circle, the practical Borel/Baire distinction
collapses.  As proof interfaces they still differ usefully:

```text
Baire-coded = visible through continuous / compactly supported observables
Borel-coded = visible after countable event closure
nonmeasurable selector = forbidden as an invariant proof object
```

The concrete target remains:

```text
0 < Haar(GOOD ∩ G_P)
```

with `GOOD ∩ G_P` carried by the formal `goodSet ∩ safeSet` event.

### Proof-use rule

For every quotient in the LRC14 route, audit whether it is:

```text
compact-group translation,
finite counting quotient,
measurable labelled chart,
or risky selector collapse.
```

Only the first three can feed the Haar/Borel/Baire witness route without
additional address data.

## codex-S144 / HYP-2946 -- Farey-perfect Kuratowski carrier verdict

Imported the prompt's perfect-number/Farey-product/Kuratowski discussion into
the current POKE stack.

### Verdict

Use perfect numbers as `K_{p,q}` edge-load stress tests, not as forbidden-minor
evidence.  The graph proof layer is minor/subdivision transitivity.

### F3/F4 split

```text
F_3: 2/3 -> K_{2,3}, product 6, planar perfect-product seed
F_4: 3/4 -> K_{3,4}, product 12, first complete-bipartite K33 wall
```

So product perfection and bipartite nonplanarity decouple.

### Perfect-number lane

Even perfect numbers give

```text
2^(r-1)/(2^r-1) -> K_{2^(r-1),2^r-1}
```

with edge count `2^(r-1)*(2^r-1)`.  After `2/3`, every row is nonplanar only
because it already contains `K_{3,3}`.

### Kuratowski guardrail

`K5 + K3,3` is nonplanar but not minimal; deleting one component still leaves a
nonplanar component.  Do not treat disjoint unions or mediants as new
obstruction cores.

### Proof-use rule

Keep exact `M`/Farey branch and C27/unital labels first.  Then use:

```text
p=1: star parent
p=2: planar C27/two-block/petal branch
p>=3: K33 incidence packet
```

For LRC14, `mediant(1/14,2/27)=3/41` remains the first unit-excess K33 wall.

## codex-S141 / HYP-2943 -- regular solids, tiling recursions, and the 14-annulus

Tested the Platonic/Archimedean/Johnson and square/triangle/hex tiling prompt
as an LRC14 carrier audit.  Use solids and tilings as binary-relational carrier
labels only after exact `M`/Farey and C27/unital branch data are attached.

### Verdict

The useful new object is not a Platonic solid.  It is the annular
prism/antiprism family:

```text
n-gonal prism:      (4,4,n)
n-gonal antiprism:  (3,3,3,n)
```

Both have:

```text
V = 2n
kappa = 1/n
```

So `n=14` gives:

```text
28 vertices
per-vertex defect 1/14
```

### Relation To S140

HYP-2942 q=3 unital:

```text
28 points
pair-unique incidence
detects the H12 repeated-pair obstruction
```

S141 14-prism/antiprism annulus:

```text
28 vertices
two 14-cycles
cyclic order and half-step/twist geometry
local defect 1/14
```

### Tiling Recursion Readout

```text
square self-dual:         Gaussian axis indices 4,9,16,25,...
triangular self:          Eisenstein indices, dyadic spine 4,16,64,...
triangle <-> hex bridge:  local support-six index 6
hexagonal self:           Eisenstein norm N(3+omega)=7 -> 7,49,343,...
centered hex rings:       7,19,37,... (different carrier)
```

### Solid Roles

Platonic solids are positive-curvature regular-map skeletons.  Archimedean
solids preserve one vertex-figure word and are local-quotient analogues.
Johnson solids (`92`) are mixed-vertex finite residual-atlas analogues, close
in proof role to bounded AP/GW/petal/K33 frontier tables.

### Proof-use rule

Use the 14-prism/antiprism as the cyclic companion to the q=3 unital, not as a
replacement for it.  Next POKE task: label the two 14-cycles by
`AP,GW,H1..H13,D1..D13` and test whether the `H12` GW/K33 conflict becomes a
twist, a diameter, or a forced two-chart obstruction.

---
## codex-S140 -- C27 unital block-lift verdict

Tested whether HYP-2937/HYP-2940 marked C27 transfers lift into q=3 unital
4-point blocks after AP/Goddyn-Wong labels are attached.

### Verdict

The q=3 unital is a **branch-local chart and calibrated pair-completion
forum**, not a global C27 atlas.

### Negative global result

Under the natural residue-pair lift

```text
H[a] -> D[d]  maps to  {a,27-a,d,27-d},
```

the two `12`-branch packets become

```text
GW  H12->D3  = {3,12,15,24}
K33 H12->D9  = {9,12,15,18}.
```

They share `{12,15}`.  A q=3 unital is `2-(28,4,1)`, so two distinct blocks
cannot share a pair.  Therefore one unital chart cannot hold both the tight GW
transfer and the K33 near-miss transfer.

The global `{AP,GW,H_a,D_d}` model fails even faster: every transfer repeats
the pair `{AP,GW}`.

### Positive local result

Branch-local charts embed:

```text
tight chart: GW + P10 + P13
K33 chart:   K33 + P10 + P13
```

The S138 genuine two-hole rows lift as two-block splices:

```text
drop(10,12)->add(20,24) = P10 + GW
drop(10,12)->add(20,36) = P10 + K33
```

### AP/GW-calibrated completion

There is also a useful noncanonical Hermitian q=3 labelling by

```text
AP, GW, H1..H13, D1..D13
```

with anchor block

```text
{AP, GW, H12, D3}.
```

This makes the Goddyn-Wong transfer the AP/GW block.  The calibrated K33 splice
is the visible incidence chain

```text
AP/GW --H12-- near/K33 --D9-- petal10.
```

### Proof-use rule

Use the unital as a pair-unique local grammar.  If a proof wants both `12`
branches in one object, it must explicitly split the H12 pair by an additional
branch coordinate or use multiple charts.  If working after AP/GW calibration,
compare unique completion blocks and their intersections before attempting a
HYP-2908/THM-572 state lift.
