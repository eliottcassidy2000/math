# LRC14 Lebesgue Wall Certificate (S676)

Prompt: consider Lebesgue measure and work on the LRC `n=14` proof.

The useful distinction is now very crisp:

```text
p_0 > 0  = open interval of lonely times
p_0 = 0  = endpoint-only wall
```

Lebesgue measure is therefore not the final proof object.  It is the filter
that throws away all strict rows and names the residual wall.  That wall still
needs set-level endpoint witnesses because the danger arcs are open.

## Exact Audit

S676 adds `04-computation/lrc14_lebesgue_wall_s676.py` and stores
`05-knowledge/results/lrc14_lebesgue_wall_s676.out`.

The script computes the exact rational circular subdivision induced by all
danger-arc endpoints.  For each speed `v`, the danger set is

```text
D_v(delta) = {t in R/Z : ||v t|| < delta}.
```

The safe measure is the depth-zero mass:

```text
p_0(V,delta) = measure {t : no D_v(delta) contains t}.
```

At `delta=1/14`, the three known normalized floor atoms are exactly
measure-zero:

| Row | `p_0` | endpoint witnesses |
|---|---:|---:|
| AP | `0` | `6` |
| Vstar | `0` | `6` |
| `2AP` | `0` | `12` |

AP and Vstar share the endpoint witnesses
`1/14,3/14,5/14,9/14,11/14,13/14`.  `2AP` has their two doubling-preimages.

The surprise worth keeping is that all three floor atoms open with the same
one-sided slope when the collar is relaxed:

```text
p_0(1/14 - eps) / eps = 23324/6435.
```

This held exactly for `eps=1/700,1/1400,1/2800,1/5600`.  It says AP, Vstar,
and `2AP` are not only equally tight; they have the same first-order
Lebesgue wall normal.

## Carry Tax As Measure

S666 had already shown local `+27` carries through weight `3` are strict in
maximin.  S676 upgrades this to positive safe measure:

- local carry probes: `1134`;
- walls: `3`, exactly AP, Vstar, `2AP`;
- positive-measure rows: `1131`.

This matters because `M>1/14` implies an interval of safe times by continuity,
but the direct `p_0` audit tells us how much interval mass has opened.  The
smallest local positive masses are not huge, but they are exact rational
certificates.

The one-swap AP scan through new speed `60` reinforces the same picture:
`12->24` is the only zero-measure swap, and every other one-swap row has
positive measure.  The smallest positive one-swap mass found is `1/1260` at
`12->36`.

## Proof Shape

The proof shape should be:

1. Normalize by scalar symmetry and the `Res_27` quotient.
2. Use exact measure to discard every positive `p_0` row.
3. Show the only normalized `p_0=0` walls are AP, Vstar, and `2AP`.
4. Verify those walls by endpoint witnesses.

This avoids the common trap: trying to prove a positive measure lower bound
for a problem whose extremal examples have measure zero.

The new no-leak target is stronger than the S666 owner-derivative target:

> every non-scalar, non-floor lift has positive Lebesgue safe measure.

HYP-2241's owner-private deletion bit is still the likely side channel, but
the conclusion should now be phrased in measure terms, not only in maximin
terms.

## Tournament Analysis

Vertices are proof routes, not runners.  I explicitly considered runners,
danger arcs, breakpoints, safe components, endpoint witnesses, `Res_27` rows,
carry vectors, owner obligations, and proof obligations as possible vertices.

The selected route tournament preserves the LRC predicate "positive interval
versus endpoint-only wall" while forgetting raw runner labels and most phase
order.  The route tournament is transitive:

```text
endpoint_wall_certificate
> local_carry_measure_tax
> res27_floor_quotient
> global_scalar_floor_orbit
> one_swap_wall_sieve
> positive_measure_interval
> raw_first_moment
```

The first moment sits last on purpose.  It is true and elegant, but it is
measure-blind to the wall.  The endpoint wall certificate sits first because
it is the only route that actually touches the residual after measure has done
all it can do.

## Next Move

Run the same exact `p_0` wall audit on coherent carry subspaces rather than
only Hamming balls:

- scalar unit lifts, as the globally coherent floor control;
- two-block carry patterns;
- HYP-2165 owner-route lifts of the 64 fixed classes;
- HYP-2241 private-owner fibers.

The target is a normalized theorem:

```text
no new p_0=0 wall in the Res_27 carry/owner fiber.
```

If that theorem lands, LRC14 reduces to the already visible endpoint witnesses
on AP, Vstar, and `2AP`.
