# LRC14 polynomial-method witness ledger

*codex-2026-06-27-S255.  Continuation after the S31ag paper bridge, the S254
hyperoperation grid address, and the earlier q-Pochhammer / full-modular-cusp
finite-principal-part work.  This note merges the incoming paper summary with
the user's request to keep the hyperoperation hierarchy and space-filling grid
inside the LRC14 proof advance.*

## Main synthesis

The incoming S31ag reflection is the right pivot: for 13 speeds / 14 runners,
the Sungkawichai-Trakulthongchai polynomial method fails for the same reason
the project sees the apex wall:

```text
k+1 = 14 = 2 * 7
Z_14 is not a field
CRT splits the failed interpolation into c=7 and c=2 lift obligations
```

The project already owns the `c=7` side as THM-573, the level-7 sieve:
`>=7` multiples of `7` imply `M > 1/14`.  What survives is exactly the
residual where that lift does not fire:

```text
primitive covering rows
with <= 6 multiples of 7
plus the c=2 / dyadic / analytic witness debt
```

This makes the paper's named bottleneck, computing `I(k,p,1)`, much more
concrete in project language.  It is the finite denominator-grid version of:

```text
Does the lonely set L(S) hit (1/p)Z?
```

The project's witness route replaces a large table by the largest-arc theorem:

```text
mu(L(S)) >= m0
components(L(S)) <= A0
=> largest arc >= m0/A0
=> every denominator d >= ceil(A0/m0) has a witness
```

That is the useful version of "Conjecture 7.1(13) equals our witness route."
The honest direction is important: the largest-arc route implies Conjecture
7.1(13), which implies LRC14.  Treating scalar LRC14 as already equivalent to
Conjecture 7.1 needs the finite-packet compactness and direct lonely-arc
bound that remain open.

## Where the q-cusp idea fits

The full-modular-group cue says a modular function invariant under the modular
group and meromorphic at the cusp has a q-expansion with only finitely many
negative powers.  In this LRC proof ledger, that becomes a disciplined
finite-exception rule:

```text
bad denominator debt must have a finite principal part
```

So a proof route is legal only if it proves a uniform denominator threshold
`D`, discharges each smaller denominator by exact packet work, or names the
remaining debt.  An infinite tail of unaccounted bad denominators is the LRC
analogue of an illegal infinite polar tail.

This does not make q-Pochhammer a direct LRC proof engine.  It supplies the
controlled-forgetting test: once the witness route proves a largest arc, the
only possible bad denominator debt is finite and must be explicitly recorded.

## Where the hyperoperation grid fits

The hyperoperation hierarchy on `(numerator, denominator)` and the older grid
tiled by `x+2` and `x*2` should remain attached to the ledger, but only as an
address chart:

```text
root:      p/q
sum lane:  p+q       additive / x+2 / endpoint-owner pressure
product:   p*q       product / x*2 / factor and v7-depth pressure
powers:    q^p,p^q   stress lanes
```

After THM-573, the vertical/product lane must track level-7 status, not merely
whether a speed is a multiple of `14`.  The proof-safe operation cell has to
retain:

```text
(p,q), p+q, p*q, power word, danger deficit, endpoint owner,
count_7_divisible, c=7 lift status, c=2 lift status, destroyed coordinate,
finite address, terminal exit.
```

The space-filling curve remains a scheduler through these cells.  It is useful
only if adjacent moves change a retained coordinate and preserve the witness
predicate or name the sidecar that pays for the loss.

## The improved LRC14 target

The remaining proof now has a shorter statement:

```text
For every non-tight primitive covering 13-tuple S with <=6 multiples of 7,
prove a uniform largest lonely interval in L(S).
```

The direct proof stack should be:

```text
1. THM-573 closes count_7_divisible >= 7.
2. THM-530/565 supplies a measure floor m0 on the relevant lonely/witness set.
3. Prove the direct 1/14 lonely-set component bound A0
   or reduce L(S) to the THM-565 maxgap>1/7 object with controlled loss.
4. Pigeonhole gives largest arc >= m0/A0.
5. Largest arc gives denominator-net witnesses for all large d.
6. This proves Conjecture 7.1(13) for the residual, hence LRC14.
```

The open technical item is still item 3.  Everything else is already in the
project in some form, or is elementary once the component bound exists.

## Tournament Analysis and challenged assumption

Tournament vertices should be proof obligations:

```text
largest_arc_denominator_net
direct_lonely_component_bound
lonely_measure_floor
crt_c7_level7_lift_THM573
crt_c2_dyadic_lift
continuous_I_substitute
finite_principal_part_bad_denominator_budget
hyperoperation_grid_address
polynomial_prime_field_packet
raw_I_table_enumeration
```

The observable is whether one vertex preserves the LR predicate, reduces CRT
debt, retains the denominator clock and lonely-set topology, bounds or
replaces `I(k,p,1)`, and names destroyed coordinates.  This gives a transitive
proof-routing tournament with the largest-arc denominator net first and raw
table enumeration last.

The challenged assumption is that runners, residues, or raw operation-grid
cells are the natural vertices.  They are not load-bearing enough here.  The
proof needs the objects that survive quotienting: CRT lift status, component
count, measure floor, denominator threshold, and terminal exit.

## New artifact

Created `HYP-3089` as the formal ledger:

```text
05-knowledge/hypotheses/HYP-3089-lrc14-polynomial-method-witness-route-ledger.md
```

Navigation hooks: `T1170`, `LTI-235`, `LTT-133`.
