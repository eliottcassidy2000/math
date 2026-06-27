# LRC14 Residual Section / Packet Grid Verification

HYP-2996 makes the residual-section story executable over the HYP-2963 packet
bank.

The verifier imports the HYP-2963 labelled-packet classifier and the HYP-2989
Haar-product grid.  It assigns each packet to a residual section and a grid
exit:

```text
Q-WITNESS       -> orthogonal_zero
AP/GW           -> same_tile_indicator
C27/unit petal  -> horizontal_owner_strip
open Haar       -> vertical_owner_strip
covering moment -> nested_refinement
K33/THM-572     -> cross_handoff
candidate F7    -> hard + zero-open + non-AP/GW + no owner/cross/nested exit
```

Default run:

```text
audited packets=21913
hard q>=14 packets=7237
hard non-AP/GW packets=7235
hard non-AP/GW packets with owner/cross/nested grid exit=7235
zero-open hard non-AP/GW packets=0
candidate F7 harmonic sections=0
validation errors=0
```

Packet-grid counts:

```text
orthogonal_zero=14676
same_tile_indicator=2
horizontal_owner_strip=4
vertical_owner_strip=6040
nested_refinement=1188
cross_handoff=3
```

Readout: in the bounded HYP-2963 bank, there is no anonymous residual section.
F7 is now a concrete missing-section predicate rather than a vague bucket.  The
next useful move is familywise section templates, then Fejer/Ramanujan
certificates grouped by section.
