# LRC14 Two-Branch Obstruction / Helly Reflection

Date: 2026-06-28

The useful correction to the user's resonant-transparency picture is now even
sharper than HYP-3422, and it sits downstream of HYP-3424's covering-floor
duality transfer.  It is not enough to say the optimizer moves off the `14`
grid and resonant speeds become transparent.  The proof has to explain why the
two-adic lift has a legal branch after the even half chooses `u=2t`.

The exact obstruction is two-colored.  Branch `0` fails when an odd speed is
too near an integer in `o*u/2`; branch `1` fails when an odd speed is too near
a half-integer.  Therefore a real failure must sit in

```text
E_safe cap B0_odd cap B1_odd.
```

That is a much smaller and more proof-shaped object than "resonance."  The
audited rows all have positive leftover after removing this bad core.  The
tight row `{1..11,13,84}` is especially instructive: `E_safe` has measure
`107/245`, the two-color bad core has measure `314/735`, and the survivor is
exactly `1/105`, spread over four finite-ruler windows.

This suggests the next proof attempt should be local on components of
`E_safe`.  Try to show that no component can be fully covered by the combined
near-integer and near-half odd obstructions, or at least that one component
survives after a two-adic descent step.  Owner-current labels like `2:g2` are
still valuable, but only as exception names; they should not replace the
interval lemma.
