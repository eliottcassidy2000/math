# LRC14 Constant Stack - Codex S46

The useful correction from the incoming KPS HYP-2653c/HYP-2653d work is that
the phrase "the one open constant" was overcompressed.  It mixed two different
currencies, and one transient normalization was wrong as a proof target.

For local shell-full rows, the meaningful denominator is `p1(E')`, the mass of
base intervals missing exactly one inner sector.  That is why HYP-2670 and
HYP-2671 naturally produce `2/5`, `1/3`, and possibly `1/4` packet taxes.

For global far peels, the corrected HYP-2653d object is not `w*Delta_w`.
That quantity is a useful resonance detector, but it grows with scale because
`Delta_w` has a small nonzero dyadic-family floor.  The proof object is the
uniform tail cap

```text
sup_{max(E')>B} Delta_w(E',w) <= cap_k-Q(k-1),
```

with `max(E')<=B` checked finitely.

The key synthesis is that these are the same endpoint discrepancy viewed after
different quotients:

```text
boundary-local quotient:  charge endpoint imbalance to p1(E')
wide far-peel quotient:   bound absolute Delta_w after finite cutoff
```

This makes the route less magical.  The shell-damaged gate asks for the
dyadic-1 tower `{1,2,4,8}` to survive.  The shell-full quotient asks for
phase-packet cancellation below fixed boundary taxes.  The wide quotient asks
for a uniform cancellation estimate on the signed plateau deviation.

The finite-span numbers are still clarifying as a diagnostic.  At the tight k=9 far-peel margin
`cap_9-Q(8)=129643/980980`, the old non-resonant `C=39/20` asks for span 15,
KPS's scale-count example `C~=2.93` asks for span 23, and the HYP-2655
multiscale sample asks for span 30.  But after HYP-2653d those spans should not
be read as theorem targets.  They explain why the old normalization looked
plausible; the actual closer is finite pocket plus uniform `Delta_w` tail.

The next proof should therefore not try to "prove the constant" in isolation.
It should prove the constant stack:

1. shell-damage gate at `426/35035`;
2. finite shell-full packet ledger below `2p1/5`;
3. new-speed packet tax below `p1/3`, with the `m=4` dyadic block isolated;
4. far-tail decay, plausibly below `p1/4`;
5. global uniform `Delta_w` tail below `cap_k-Q(k-1)`, likely after a finite
   cutoff near `B=20`.

Tournament Analysis: the vertices are proof currencies, not runners.  The
Hamiltonian path is

```text
shell_damage_gate
> finite_packet_tax
> new_speed_packet_tax
> far_tail_packet_tax
> uniform_delta_tail
> raw_runner_vertices.
```

The challenged assumption was that a scalar constant determines the proof.
The better invariant is a quotiented currency: boundary mass locally, and an
absolute signed plateau-deviation margin globally.
