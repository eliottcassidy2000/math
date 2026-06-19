# LRC14 modular recurrence address: mod 6 proves, mod 30 addresses, mod 210 couples

**Source:** codex-2026-06-19-S18.

The useful correction from this pass is that mod `30` and mod `6` are not the
same kind of object.

Mod `6` is the HYP-2598 universal-center skeleton.  The centers `1/2`, `1/3`,
and `2/3` have gaps large enough to be honest large-spread cluster witnesses,
and the survivor count is exactly

```text
C(7,s) + C(9,s) - C(5,s).
```

Adding `5` gives the next clean integer sequence, but it is no longer a
standalone proof tool.  The denominator-5 centers are exact modular addresses:
they sort small parts and recurring intervals, yet their center gap is below
the `2/7` cluster threshold.  This is why mod `30` keeps appearing without
closing the proof.

The missing address is then not another isolated modulus.  The support-six
Fourier tail is intrinsically mod `7`, so the natural combined bookkeeping is
the squarefree mask on `{2,3,5,7}` - equivalently mod `210` when one wants a
single modulus.  The exact survivor rows are:

```text
{2,3}:     1, 11, 47, 109, 156, 146, 91, 37, 9, 1, 0, 0, 0, 0
{2,3,5}:   1, 13, 72, 223, 437, 581, 545, 366, 174, 56, 11, 1, 0, 0
{2,3,5,7}: 1, 13, 78, 280, 658, 1066, 1231, 1030, 623, 266, 76, 13, 1, 0
```

The last row should not be read as "mod 210 proves the theorem."  It says the
small-part recurrence and the support-six coimage tail now live in compatible
coordinates.

The concurrent HYP-2624 height-2 coimage wall audit makes this more concrete:
the mod-210 address is the coordinate system in which finite wall accounting and
the remaining k=10 repeated-residue signed tail can talk to each other.

The finite transfer state is the divisor mask of each speed in `{1,...,13}`:

```text
0000: 3, 0001: 3, 0010: 2, 0011: 2, 0100: 1, 0101: 1, 1000: 1.
```

That is the hidden recurrence relation.  It is not a scalar recurrence in the
runner list; it is a mask-polynomial / inclusion-exclusion recurrence coupled
to signed mod-7 coimage fibers.  This also explains the KPS primorial warning:
fixed moduli give fixed constants, while growing constants come from the
divisor lattice acquiring more small-prime bits.

Assumption challenge: I considered runners, residues, rational centers,
denominator sets, prime masks, support-six coimage classes, and proof
obligations.  Runner residues preserve too much geometry and hide the
recurrence.  The prime-mask quotient preserves denominator survival but loses
cluster-width data, so it must be paired with the signed coimage tail rather
than used alone.

Net result: no proof of LRC(14), but the modular search target is sharper:

```text
mod-6 universal skeleton
-> mod-30 recurrence address
-> mod-210 divisor-profile/coimage coupling
-> finite wall deletion
-> signed reciprocal tail.
```
