# LRC14 Parity Carry Bridge S654

The Basel carrier was a good nudge: do not collapse the product too early.
For zeta, the product side is the elementary-packet family
`zeta({2}^m)=pi^(2m)/(2m+1)!`; the moment side is `zeta(2k)`.  Newton moves
between them only because the packet side was retained.

The LRC14 version is almost embarrassingly literal.  A lifted speed over the
`C=27` quotient is

```text
v = r + 27k.
```

The residue `r` is not enough.  The carry `k` is the packet side:

```text
v mod 2  = r+k,
v mod 14 = r-k.
```

So the same bit of side information toggles even/odd parity and decides
whether the row has the required multiple of `14`.  This is exactly the user's
even/odd-pair intuition in LRC clothing.  Even numbers use an unordered pair
shadow; odd numbers need the ordered/doubled bridge.  Here the multiple-of-14
apex is the doubled bridge, and it appears only when

```text
k == r mod 14.
```

Pair-sums inherit the ordered version:

```text
14 | v_i+v_j  iff  k_i+k_j == r_i+r_j mod 14.
```

That is the useful conceptual move.  The residual is not "some lift might be
bad"; it is a carry-cocycle problem with one-coordinate and two-coordinate
faces.

The S654 computation pushed directly on the forced branch.  S612 searched
`L1<=6` and the Boolean lattice, which was already strong but could not touch
minimal apex carries for residues `7..13`.  S654 sets the minimal apex carry at
each AP and `V*` coordinate and computes exact loneliness.  Every one is strict.
The closest AP case is residue `13`:

```text
M = 28/365 > 1/14.
```

The closest `V*` case is also residue `13`:

```text
M = 2/25 > 1/14.
```

Then I added all one- and two-coordinate `+27` toggles around every minimal
apex lift.  Among `2054` screened rows, no near-floor and no below-floor row
appeared; all exact group minima stayed strict.  The Boolean parity profile
also behaved: the only exact floor row in the near-floor screen was the
zero-carry base row.

So the first carry bridge is not the monster.  If a new floor lift exists, it
must be more coordinated: a larger carry cocycle, a scalar AP/`V*` route, or an
owner/Cprime certificate.  That is real progress because it gives the next
theorem a local shape:

```text
minimal apex carries pay positive tax;
floor-preserving carries must be global.
```

This also explains why parity kept returning in the number-theory sessions.
Parity is not just a color.  In the `27` tower, parity is the mod-2 shadow of
the same carry that the mod-14 clock reads as apex obstruction.  Odd/even,
ordered/unordered, Basel packets/moments, and LRC lift/CRT are all saying:

```text
the scalar is downstream of a retained packet.
```

The theorem target I would hand to the next session:

```text
Prove an apex carry tax lemma over AP/V*:
if a primitive lift over either floor shadow has a minimal apex carry and all
other carries are locally small/noncoherent, then M>1/14.

Then identify the cohomology-like normal forms for carry patterns that can
avoid this tax.  Expect the only floor-preserving normal forms to be scalar
AP/V* carries or rows discharged by owner/Cprime certificates.
```

Artifacts:

- `04-computation/lrc14_parity_carry_bridge_s654.py`
- `05-knowledge/results/lrc14_parity_carry_bridge_s654.out`
- `05-knowledge/hypotheses/HYP-2230-lrc14-parity-carry-bridge.md`
