# LRC14 Core-Gap Survival Bridge - S36 Reflection

Date: 2026-06-19
Session: codex-2026-06-19-S36
Claimed: HYP-2658 / T903

## Why This Session

The prompt asked to integrate other agents' ideas and keep pushing toward a
proof.  The live overlap sharpened while this session was in flight:

- HYP-2651 gave the fixed-observer positive-core gap atlas for THM-523.
- KPS HYP-2653 supplied the exact decorrelation engine and the honest
  uniform-discrepancy open problem.
- KPS HYP-2655 later refuted the small uniform-constant dovetail and moved the
  far route to joint plateau/Delta recursion.
- HYP-2654 and THM-541/542/543/544 proved the first AP-tail layers.
- S15's dyadic/apex-prime residue work reframed dyadic richness as a
  tiebreaker/address, not the determinant.
- HYP-2648 and HYP-2652 kept saying: retain address before scalarization.

So the local computation was reframed as a survival bridge: how does the far
decorrelation route look for the fixed-observer gap predicate, and what address
does the THM-543 exceptional row carry?

## Main Finding

For the core-gap predicate the far-speed limit should be survival, not filling:

```text
meas(G_{B union {w}}) -> (6/7) meas(G_B).
```

This is the fixed-observer sibling of HYP-2644's plateau.  In HYP-2644, a far
speed can fill a missing sector.  Here it imposes another avoidance condition,
so an independent far point survives with probability `6/7`.

The tested multi-far ledger has comfortable margin over the collar.  Even the
closest one-far base in the bounded scan gives

```text
(6/7)*(313/9702) = 313/11319,
313/11319 - 7/858 = 5737/294294 ~= 0.01949.
```

So the far limit itself is not tight.  The danger is finite resonance before
decorrelation, exactly the zone now controlled through the one- and
two-replacement AP-tail layers by THM-543/544.

## The THM-543 Bubble

The local scan rediscovered the mainline theorem's exceptional row:

```text
C20 = (1,2,3,4,5,7,8,9,11,12,13,20),
meas(G_C20) = 3859/420420 = 7/858 + 1/980.
```

THM-543 proves this is the unique one-replacement AP-tail row below
`426/35035`, namely `(a,b,r)=(6,10,20)`, and THM-544 then proves no
two-replacement AP-tail row falls below that same threshold.  The useful extra
data here is the component anatomy.  `C20` is just the collar with `10`
replaced by `20`; the four old collar components remain, and two new bubbles
appear:

```text
[29/98, 83/280]   length 1/1960   owners 7 -> 20
[197/280, 69/98]  length 1/1960   owners 20 -> 7
```

This is the "anti-coset everywhere" clue in interval form: doubling the missing
collar speed creates an owner-addressed symmetric bubble pair rather than a
generic tail.

## What Changed

Before fetching the other agents' work, this looked like a fresh local
hypothesis/tangent candidate.  After fetch, that namespace was already
productively occupied:
KPS HYP-2653 is the far decorrelation engine, HYP-2654 is the near-collar
template, and THM-543 proves the one-replacement AP-tail layer.  A second fetch
then claimed KPS HYP-2655 for multiscale decorrelation growth and THM-544 for
the two-replacement AP-tail layer, so this packet moved once more to
HYP-2658/T903.

The corrected role of this packet is HYP-2658/T903:

```text
near collar: THM-541 -> THM-542 -> THM-543 -> THM-544
genuinely far: fixed-observer 6/7 survival + HYP-2655 joint plateau/Delta recursion
glue data: endpoint-owner component ledger + HYP-2648 state words
```

## What Failed

Raw sumset excess failed again.  The old `B<=19` second row has excess `1`,
while the THM-543 exceptional row has excess `9`.  The scalar order is
backwards.  The component address sees the truth immediately: the exceptional
row has the collar components plus two tiny owner bubbles.

Raw far speed also failed as a determinant.  The best appended row among the
top 11-core bases uses `w=20`, but other bases have resonant best appenders at
`w=46`, `50`, `37`, and so on.  The far-discrepancy object needs a finite
resonance ledger before the asymptotic `6/7` survival estimate can be used.

S15's dyadic-richness result fits this rather than changing it: QR/NQR and
doubling-orbit structure can explain tiebreaks near apex-prime cover walls, but
inside this fixed-observer bridge it lives in the address layer alongside
state words and endpoint owners.

## Proof Target After This

1. Treat THM-541/542/543/544 as the proved near-collar AP-tail base layer.
2. Use the `10->20` component bubbles as the endpoint-owner model for later
   AP-tail/state-word template ledgers.
3. Prove the fixed-observer far survival recursion after finite mouth/replacement
   templates are removed:

```text
|meas(G_{B union {w}}) - (6/7)meas(G_B)| <= C/w.
```

4. Use the `0.01949` tested survival margin as the first quantitative target,
   but do not rely on a false small uniform constant; KPS HYP-2655 says
   multiscale discrepancies accumulate and must be paired with shrinking
   plateaus.
5. Retain endpoint-owner/state-word addresses until scalarization is forced.

LRC(14) is not proved.  But the proof shape is less blurry: proved collar and
two AP-tail layers first; genuinely far rows by joint recursion second;
addressed components as the glue.

## Tournament Analysis

Tournament vertices were proof quotients:

```text
proved_one_replacement_tail
proved_two_replacement_tail
exact_core_gap_components
drop6_collar_ledger
owner_bubble_ledger
far_survival_multiplier
joint_plateau_delta_recursion
finite_resonance_discrepancy
sumset_excess
raw_far_speed
```

Pairwise observable: which quotient preserves the lower-bound predicate and
explains the sub-`426/35035` near-collar rows.

Switch/gauge: keep fixed-observer positivity; judge far speeds after applying
the `6/7` survival quotient and the KPS HYP-2655 joint plateau/Delta recursion.

Hamiltonian path:

```text
proved_one_replacement_tail
> proved_two_replacement_tail
> exact_core_gap_components
> drop6_collar_ledger
> owner_bubble_ledger
> far_survival_multiplier
> joint_plateau_delta_recursion
> finite_resonance_discrepancy
> sumset_excess
> raw_far_speed
```

The tournament is transitive; no directed 3-cycles.  The useful vertices are
not runners or arcs.  They are proof obligations with enough address retained
to avoid false scalar orderings.
