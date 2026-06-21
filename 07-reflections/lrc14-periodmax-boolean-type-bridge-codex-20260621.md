# LRC14 Period-Max / Boolean-Type Bridge

**Session:** codex-2026-06-21, HYP-2790.

The bridge test was useful because it failed in the right way.  After
`THM-563`, the one-far error is not an analytic tail but a finite signed
endpoint-period maximum:

```text
periodmax(B) = max_w w * Delta_w(B).
```

It was tempting to hope that the low-depth Boolean/type cut from HYP-2791, or
the KPS containment-type cuts from HYP-2752, also ranked this period pressure.
The exact scout says no.  Across `135` high-plateau bounded bases k=8..12,
every checked row satisfied the needed inequality, but the explanatory
coordinate was not the Boolean/type slack:

```text
global worst checked ratio = 13.280470214 < 15
B = (0,2,4,6,8,10,12,14), k=9, periodmax=86/49.
```

The AP dilation has zero Boolean/type deficit and zero containment deficit, as
it should.  That quotient sees it as the same extremal object as the consecutive
base.  But the period numerator doubles under the dilation, so the quotient has
destroyed exactly the phase-order data needed for the single-far oscillation.

For non-AP rows, endpoint arc count was the stable pressure signal.  The
correlations with `periodmax/margin` by k were:

```text
k=8  arc_count +0.6556
k=9  arc_count +0.7331
k=10 arc_count +0.6924
k=11 arc_count +0.4971
k=12 arc_count +0.7114
```

The containment-cut deficit changed sign/strength across k.  Sorted cell leaks
were also not stable enough.  So the Boolean/type basis remains a good bounded
plateau/cap carrier, but it is not the right carrier for the oscillatory
single-far numerator.

The next proof object should be the endpoint-period sum itself, in the same
direction as HYP-2792:

```text
sum_j sum_{t endpoint of A_j} +/- S_j({w t})
```

This is a signed generalized Dedekind sum over the endpoint orbit.  The target
is now to prove a uniform reciprocity/arc-pairing bound for that sum, then
compare it against `15*(cap_k-Plat(B))`.  A Boolean/type finite ledger would
have been convenient, but it would have hidden the `w`-multiplication phase
that actually makes the signed cancellation work.

Assumption challenge: I tested vertices as bounded bases, missed-sector type
atoms, containment atoms, sorted cells, endpoint arcs, and explanatory-lens
tournament nodes.  The quotient that preserves the LRC predicate best for this
branch is not the runner set or the sector-type law; it is the endpoint-arc
orbit under multiplication by `w`.
