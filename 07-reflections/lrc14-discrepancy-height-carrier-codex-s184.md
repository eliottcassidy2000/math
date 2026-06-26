# LRC14 Discrepancy-Height Carrier Reflection (codex S184)

The useful move in this pass was to treat the recent automaton failure as a
carrier-design lesson.  HYP-3016 and HYP-3017 did not say automata are useless;
they said automata are missing the coordinate that decides whether a fiber is
boundary or open.  The S184 trident keeps three different kinds of missing
information separate:

```text
residue discrepancy     which denominator clocks are noisy
height / Mahler scale   how far the magnitude cocycle moved
Hensel singular status  which local p-adic lifts need guarding
```

On the bounded named plus AP single-swap bank, the combined trident has zero
mixed boundary/open fibers.  That is encouraging, but the result should be read
with discipline: the signature is almost exact on this bank (`2167` fibers for
`2173` rows).  The next proof improvement is not to celebrate that fine a key;
it is to compress the key until it stops being row identity while preserving
boundary/open purity.

The strongest insight is that a good proof carrier can be product-shaped.  We
do not need one scalar that recognizes LRC14 packets.  We need a small tuple of
independent clocks whose failure modes are typed:

```text
height failure      -> magnitude cocycle / Farey side channel
residue failure     -> Ramanujan, Haar, Fejer, exact-period side channel
p-adic failure      -> Hensel / finite lift / THM-572 side channel
topology failure    -> HYP-3015 barcode / endpoint-owner side channel
```

This fits the repo's no-free-slider rule.  A quotient is admissible only when
it can say what it forgot and which other tooth restores or annihilates that
forgotten coordinate.

Practical next steps:

1. Add the S184 trident sidecar to the full HYP-2963 classifier:
   residue L1 at selected denominators, Erdos-Turan proxy bins, height/Farey
   band, Mahler proxy, and Hensel `(root,singular)` counts at `2,3,7`.
2. Start with the mixed word `MFCMMCCFFFCCC`, then try coarse signatures until
   the AP/GW boundary rows no longer share a fiber with open rows.
3. Compare the resulting fibers with HYP-3015 barcode fields and HYP-2981
   Fejer interval manifests.  The best compression is probably not pure
   trident; it is trident plus the first certificate-facing topology field.
4. Build the incidence hypergraph of packet rows versus sidecar tokens.  If
   row arity stays bounded, a Beck-Fiala style discrepancy argument may give a
   route to "some coordinate must fire" without enumerating every row.

Assumption challenge: the vertices are not runners.  They are not even
automaton states.  The vertices in the useful tournament are proof carriers,
and the directed edge asks which carrier is allowed to forget less of the LRC
predicate.
