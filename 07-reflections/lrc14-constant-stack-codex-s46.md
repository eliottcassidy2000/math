# LRC14 Constant Stack - Codex S46

The useful correction from the incoming KPS HYP-2653c work is that the phrase
"the one open constant" was overcompressed.  It mixed two different currencies.

For local shell-full rows, the meaningful denominator is `p1(E')`, the mass of
base intervals missing exactly one inner sector.  That is why HYP-2670 and
HYP-2671 naturally produce `2/5`, `1/3`, and possibly `1/4` packet taxes.

For global far peels, the meaningful denominator is `w`.  There the object is
`C(k)=sup w*|Delta_w|`, and HYP-2653c says the old `~1.95` value was a
non-resonant sample.  Resonant scale clusters add bounded chunks, so the right
target is `C(k)=O(number of scale clusters) <= O(k)`.

The key synthesis is that these are the same endpoint discrepancy viewed after
different quotients:

```text
boundary-local quotient:  charge endpoint imbalance to p1(E')
wide far-peel quotient:   charge endpoint imbalance to scale clusters / w
```

This makes the route less magical.  The shell-damaged gate asks for the
dyadic-1 tower `{1,2,4,8}` to survive.  The shell-full quotient asks for
phase-packet cancellation below fixed boundary taxes.  The wide quotient asks
for a per-scale-cluster bounded variation/Koksma-style decomposition.

The finite-span numbers are also clarifying.  At the tight k=9 far-peel margin
`cap_9-Q(8)=129643/980980`, the old non-resonant `C=39/20` asks for span 15,
KPS's scale-count example `C~=2.93` asks for span 23, and the HYP-2655
multiscale sample asks for span 30.  That is an enlargement of the finite base,
not a collapse of the route.

The next proof should therefore not try to "prove the constant" in isolation.
It should prove the constant stack:

1. shell-damage gate at `426/35035`;
2. finite shell-full packet ledger below `2p1/5`;
3. new-speed packet tax below `p1/3`, with the `m=4` dyadic block isolated;
4. far-tail decay, plausibly below `p1/4`;
5. global per-scale-cluster `C(k)` and a finite base span in the 23-30 range.

Tournament Analysis: the vertices are proof currencies, not runners.  The
Hamiltonian path is

```text
shell_damage_gate
> finite_packet_tax
> new_speed_packet_tax
> far_tail_packet_tax
> scale_cluster_Ck
> raw_runner_vertices.
```

The challenged assumption was that a scalar constant determines the proof.
The better invariant is a quotiented currency: boundary mass locally, scale
cluster count globally.
