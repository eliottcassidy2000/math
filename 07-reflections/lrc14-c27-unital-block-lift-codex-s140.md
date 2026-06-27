# LRC14 C27 Unital Block Lift

The q=3 unital became more useful once the question split in two.

The raw C27 residue-pair model gives a clean obstruction.  If

```text
H[a] -> D[d] = {a, 27-a, d, 27-d},
```

then

```text
GW H12->D3   = {3, 12, 15, 24}
K33 H12->D9  = {9, 12, 15, 18}
```

share `{12,15}`.  A q=3 unital is `2-(28,4,1)`, so two blocks cannot share a
pair.  That says one raw unital chart cannot contain both the tight GW branch
and the K33 near-miss branch.

The AP/GW-labelled version gives the complementary positive object.  Label the
28 Hermitian unital points by

```text
AP, GW, H1..H13, D1..D13
```

and calibrate one block as

```text
{AP, GW, H12, D3}.
```

Now the Goddyn-Wong transfer is literally the AP/GW anchor block.  The remaining
transfer pairs still have unique completion blocks.  In the deterministic
calibration, the two-splice result is the interesting part:

```text
10+12 -> 20+24:
  petal block and AP/GW block are disjoint.

10+12 -> 20+36:
  petal block and near/K33 block share D9.
```

That makes the incidence picture visible:

```text
AP/GW --H12-- near/K33 --D9-- petal10.
```

So the real S140 lesson is not simply "the unital lifts" or "the unital fails."
It is:

```text
global raw atlas: blocked by lambda=1 pair uniqueness
branch-local chart: works
AP/GW-calibrated pair-completion forum: works
```

The guardrail is still essential.  The labelling is noncanonical, and S107's
warning remains intact: this is not an `S8`-invariant AP8 pair-slot design.  The
surviving invariant is pair-unique completion.

The next theorem target is a packet classification:

```text
low-gap non-AP/GW residual
  -> AP/GW anchor block
  or unit-petal block
  or AP/GW--near/K33--petal chain.
```

If that classification can be forced from the exact `M`/Farey ledger and C27
unit labels, then the C27 branch has a finite incidence grammar before the
terminal HYP-2908 / THM-572 state lift.  If a future proof wants both `12`
branches in one unital object, it must name the extra structure that splits the
H12 pair or admit it is working in a multi-chart atlas.
