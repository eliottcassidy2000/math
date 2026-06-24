# LRC14 Farey Graph/PZ Carrier Synthesis

The useful outcome of this pass is a hierarchy, not a new scalar invariant.
HYP-2931 already made the arithmetic non-negotiable: `q` is the binding scale
because `M(S)-1/14=(14p-q)/(14q)`.  The new pass explains where the tempting
extra data belongs.

On the unit-excess chain `p/(14p-1)`, `q` moves by `+14` and `p+q` by `+15`.
That is the additive `n+2` lane in Farey form.  The product `p*q` has constant
second difference `28`; it is the `n*2` / area ledger, useful for coimage and
packet growth but too globally inverted to order risk.  The power payloads
are not local addresses; they are stress tests for magnitude leakage.

The graph carriers also split cleanly.  The octahedron `L(K4)` has cycle rank
`7`, so it is the support-six current/curl carrier.  The Clebsch folded cube
and halved 5-cube are the covariance/cut side: triangle-free folded residual
masks on one side, complement/halved-cube cut density on the other.  These are
not substitutes for the Farey address; they are the packet labels to carry
after the address is fixed.

Paley-Zygmund remains useful only at the front door.  It can certify nonzero
mass from second moments, but the toy six-sector model shows too much loss
for the tight LRC14 cap.  The actual cap lane still points to HYP-2823's
degree-4 factorial moment region.

The role tournament therefore has proof carriers as vertices, not runners:

```text
q_binding_scale
> p_plus_q_additive
> octahedron_LK4
> p_times_q_product
> Clebsch_halfcube
> PZ_second_moment
> power_payloads.
```

The next proof step should try to make a reduced `|M14|<=6` atom carry
`(e,q,p+q,p*q)` plus octahedral/Clebsch packet data into either the non-tight
Farey class from HYP-2930 or the forbidden `H=7` state lift from HYP-2908 /
THM-572.  That keeps the creative sequence analysis, but routes it through
the two endpoints the repository already knows how to use.
