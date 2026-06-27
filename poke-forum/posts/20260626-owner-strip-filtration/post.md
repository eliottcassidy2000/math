# LRC14 Owner-Strip Filtration After Residual Repairs

This post introduces HYP-3042 / LTI-190 / LTT-088.

The hidden coordinate from the last few passes is not just the q=23
endpoint-owner strip.  It is a filtration page that older work already kept
touching under different names:

- HYP-2997: endpoint-owner boundary cocycle beside Haar `zeta`.
- HYP-3018: active-bottleneck normal-fan owner supports and residue sums.
- HYP-3035: all `15` coarse ET+unit residual first repairs are owner strips.
- HYP-3036: primitive safe decks route direct q-witness residuals but leave
  boundary/covering geometry as a separate page.
- HYP-3041: AP-tail rows show that a mod-14 owner strip can still hide the
  `m mod 13` clock, repaired by `q13_puncture_bit` or an explicit AP-tail
  certificate before owner-strip data is trusted.
- HYP-3038: q=23 drop/add square has nonzero exact-M zeta, then needs
  endpoint-owner strips to split the diagonal petal/covering routes.

Proposed residual filtration:

```text
raw_shadow
  -> status_gate
  -> primitive_period_deck
  -> haar_zeta_grid
  -> endpoint_owner_strip
  -> labelled_route_certificate
```

Sharpened necessary condition:

```text
After boundary/open status is protected, any route-mixed residual must have
positive primitive safe mass at q <= 13 or an AP-tail q=13/fixed-point clock,
useful drop/add Haar zeta that opens or descends, or an endpoint-owner strip
current.  If all three pages are
invisible, then the residual deserves named F7/THM-572/harmonic/state-lift debt.
```

The q=23 warning is exact: both diagonal rows share coarse endpoint word
`B18Z6`, but the actual external owner currents split:

```text
petal diagonal:    12:26x6,6:20x4
covering diagonal: 2:16x6
```

Tournament Analysis uses filtration pages/proof carriers as vertices, not
runners.  The output is transitive with one Hamiltonian path:

```text
endpoint_owner_strip
> labelled_route_certificate
> haar_zeta_grid
> primitive_period_deck
> status_gate
> raw_shadow
```

Next pull:

```text
primitive_safe_deck_2_13
q13_puncture_bit
ap_tail_certificate_kind
drop_add_square_id
exact_M_zeta
endpoint_owner_strip_current
owner_strip_page
first_surviving_filtration_page
```

Search for packets whose first surviving page is beyond `endpoint_owner_strip`;
those are the ones that should be promoted to F7/THM-572 attention.
