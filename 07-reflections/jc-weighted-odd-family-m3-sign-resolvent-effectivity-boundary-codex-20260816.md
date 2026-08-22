# The first odd weighted quintic: common sign class, incomplete boundary

**Status: proved synthesis of THM-3517; the m=3 atlas is independently
audited.**  This reflection concerns the explicit cyclic weighted family
only.  It does not reconstruct THM-1605's historically reported family,
classify Keller maps, or change `JC(2)`.

**All-grade closure update:** reducing the `z`-reconstruction numerator
modulo the inverse polynomial gives a nonzero `w^(n-1)` coefficient in every
cyclic grade.  Hence `z` is primitive for all `ell`, so the common
three-coordinate square class is no longer only an `m=3` phenomenon.

**Independent-audit update:** a clean-room python-flint reconstruction now
matches every `m=3` eliminant, discriminant, and index hash, including the
191-term `z` quintic and 268-term index core.  A split `F_31` fibre gives
nonzero `x,y,z` Vandermonde determinants `(28,14,1)`.

## Inheritance pass

- **Closest proved mechanism:** THM-3494's trace-form rule: primitive views of
  one separable algebra have discriminants differing by basis-change squares.
- **Canonical hostile:** the base-field view `X=P`, whose quintic resultant is
  `(X-P)^5` and whose discriminant vanishes.
- **Corrected near miss:** MISTAKE-413—an odd divisor determines valuations,
  not the constant square class; the fixed cubic has `[-L]`, not `[L]`.
- **Least-used relevant sidecar:** THM-3448's component-wise infinity inertia,
  especially the distinction between a discriminant factor and an effective
  nonproperness component.

The first computation beyond the cubic closes with a productive failure:

```text
three primitive coordinate quintics
        -> one common class [L5]
        -/-> the complete Jelonek set V(C) union V(L5).  (1)
```

The missing implication fails for a structural reason, not a large-polynomial
accident.  The `C=0` component has local inertia `C3`, and a 3-cycle is even.

## Live concept board

| concept | invariant retained | first information lost |
|---|---|---|
| primitive coordinate | the whole degree-five field | the chosen power-basis index |
| discriminant square class | the sign character of root permutations | every even inertia cycle |
| raw discriminant divisor | ramification order with multiplicity | whether the degeneration is affine-effective |
| Jelonek divisor | actual source escape | permutation sign and labelled sheet matching |
| global `S5` monodromy | transitivity and primitive block structure | which local boundary subgroup occurs |
| historical family name | intended ancestry | the object itself when its formula is absent |

Every new observation changes at least two lanes.  The third-coordinate
quintic closes the primitive-coordinate lane, but the `C3` component refutes
the idea that this closure also closes effectivity.  Conversely, the two
Jelonek components identify the exact extra sidecar needed by the common
square class: component-wise inertia and escape residues.

## What the quintic computation actually says

For `m=3`, the inverse primitive satisfies

```text
T(w)=w^5-w^4+Pw-Q,
D5=Disc_w(T).
```

All three source coordinates generate `Q(P,Q,C,w)`, and exact resultants give

```text
Disc(E_i)=D5 J_i^2,             i=x,y,z.              (2)
```

After `P=BC,Q=AC^2`,

```text
D5(BC,AC^2)=C^4L5.                                    (3)
```

The quotient by squares turns (3) into `[L5]`.  Yet THM-3448 proves

```text
S_F=V(C) union V(L5).                                  (4)
```

The component dictionary is therefore

```text
V(C):   three escaping sheets, C3 inertia, even sign, invisible in (2);
V(L5):  two escaping sheets, transposition, odd sign, visible in (2). (5)
```

This cleanly separates four notions that had been drifting together in the
older cubic discussion:

1. coordinate primitivity;
2. equality of discriminant square classes;
3. support of the raw discriminant; and
4. effective nonproperness.

Only the first implication, `primitive views -> common class`, is automatic.

## The m=2/m=3 hostile pair

The strongest small hostile is not a random specialization but the adjacent
family member.  In both cases the pulled-back inverse discriminant has an even
power of `C`:

```text
m=2: C^2 L3,       but S_F=V(L3);
m=3: C^4 L5,       and S_F=V(C) union V(L5).           (6)
```

Thus even the **raw multiplicity plus square class** does not decide whether
`C=0` is an effective escape component.  At `m=2`, the factor is chart/index
data; at `m=3`, it is a real component with even local permutation.  The
needed sidecar is the Newton reconstruction ledger, not another quotient of
the discriminant.

## The family-level deduction

Reindexing the cyclic THM-3448 family by `ell=2m-3` makes every `C`-cycle odd:

```text
cycle length=2m-3,
raw C-order=2m-2,
sign=(-1)^(2m-4)=+1.                                  (7)
```

For all `m>=3`, `V(C)` is effective and invisible to the sign resolvent.  This
is a useful negative theorem: no amount of coordinate switching among
primitive elements can restore that component, because every such switch
multiplies the discriminant by a square.  To see `V(C)`, one must enrich the
observable, not search for a better primitive coordinate.

The global group remains `S_(2m-1)`.  The odd `C`-cycle and the `L`
transposition are local generators/constraints inside that global action;
calling the cover cyclic would discard the two-root incidence that proves
full symmetry.

## Connection contract

| field | exact content |
|---|---|
| source | `E_m^cyc`, its inverse root `w`, and actual coordinates `x,y,z` |
| target | coordinate discriminant classes and effective Jelonek components |
| map | resultant, discriminant, then target pullback `P=BC,Q=AC^2` |
| preserved | degree, primitivity, sign character, odd divisor parity |
| destroyed | even inertia, escape residues, component effectivity, labelled roots |
| restoring sidecar | Newton polygon plus reconstruction valuations on each root cluster |
| cheapest decisive hostile | compare the even `C^2` and `C^4` factors at `m=2,3` |

There is no intrinsic tournament here.  If coordinate views are made vertices,
their unsquared basis-volume ratios form a complete-graph coboundary as in
THM-3494; orienting the edges is only gauge.  The new nontrivial object is
instead the map from boundary components to conjugacy classes in `S5`.

## New frontiers exposed

1. **Closed third-coordinate gate; index recurrence next.**  THM-3517 now
   proves `z notin K` in every cyclic grade from the explicit remainder
   coefficient, and `S_n` maximality makes it primitive.  The open refinement
   is a recurrence or Newton-polytope law for the rapidly growing `z`
   power-basis index, not another primitivity test.
2. **Effectivity functor.**  Package each boundary prime as
   `(escape count, inertia cycle type, discriminant order, primitive index,
   residue rank)`.  The m=2/m=3 pair proves no projection to only parity and
   multiplicity can be faithful.
3. **Historical recovery.**  Locate a primary source containing the old
   `E_m` formula before comparing degrees or claiming conjugacy.  The mismatch
   `7,17,27,...` versus `7,13,26,...` is a stop sign, not an invitation to
   rename one family as the other.
4. **Coordinate/Jelonek stratification.**  On `V(C)` compute which coordinate
   projections collide among the finite sheets while the three-sheet cluster
   escapes.  This may explain the very different index-core sizes
   `17,26,268` without confusing projection collisions with ramification.
5. **Composition boundary.**  The fixed-map iterate changes the odd square
   class by norm cancellation.  The weighted atoms instead change the local
   cycle length with the seed.  Comparing these requires a divisor-orbit
   ledger, not multiplication of sign classes.

## Incoming fixed-tower corroboration

During this integration, the independent exact artifact
`keller_level_three_three_coordinate_primitive_independent_audit_20260816`
arrived.  Over `F_101` it gives a lawful degree-27 fixed-third-iterate fibre
where the `x,y,z` power-basis matrices all have rank 27, while retaining a
separate target where the `y` rank drops to 26.  This is strong FINITE-EXACT
evidence for the same primitive-coordinate/trace-square mechanism in the
fixed composition tower.  Subsequent independent proof audit promoted
THM-3519, which now proves the common fixed-map class `[-2J]`.  It remains a
corroborating related theorem rather than a dependency here and supplies no
weighted-family Jelonek data: primitivity supports (2), not the invalid
implication from (2) to complete effectivity.

The central lesson is compact: the three-cubic `-4*square^2*L` pattern does
extend to a three-quintic common square class, but its apparent boundary
completeness was a cubic accident.  At degree five the lost even monodromy is
already a whole irreducible nonproperness component.
