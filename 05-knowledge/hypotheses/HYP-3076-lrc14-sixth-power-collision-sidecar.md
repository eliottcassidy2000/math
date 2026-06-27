---
id: HYP-3076
title: LRC14 sixth-power collision sidecar
status: SYNTHESIS / bounded exact scout; not a proof of LRC14 or of equal sixth-power finiteness
source: codex-2026-06-27-S244
tangent: T1159
technique: LTI-224
tournament_technique: LTT-122
script: 04-computation/lrc14_sixth_power_collision_sidecar_s244.py
result: 05-knowledge/results/lrc14_sixth_power_collision_sidecar_s244.out
related:
  - HYP-3075
  - HYP-3074
  - HYP-3073
  - HYP-3072
  - HYP-3071
  - HYP-3066
  - HYP-3062
  - HYP-3060
  - HYP-3058
  - HYP-3009
  - HYP-2963
  - HYP-2887
  - HYP-2636
  - HYP-2632
  - HYP-2618
  - HYP-2617
  - HYP-2614
  - HYP-2608
  - THM-538
  - THM-572
  - OPEN-Q-108
---

# HYP-3076: LRC14 sixth-power collision sidecar

## Claim

The two prompt equations should be merged into the LRC14 proof stack as two
different relation-lattice sidecars:

```text
a^6+b^6+c^6=d^6+e^6+f^6  -> native support-six collision
a^6+b^6=d^6+e^6            -> rank-lowered square-cube shadow
```

The first equation is a genuine six-term signed relation with three positive
and three negative slots.  It belongs beside THM-538, HYP-2614, HYP-2617,
HYP-2632, HYP-2636, HYP-2887, and the HYP-3071/HYP-3074 cycle-image/route-state
closure interface.  The second equation is still useful, but differently:
since `x^6=(x^3)^2`, it is a two-square collision in cube coordinates.  If it
is lifted into a six-slot LRC support relation, it must carry an explicit
canceling pair and be marked as degenerate padded support-six data.

This is not a proof shortcut.  It is a controlled-forgetting rule:

```text
equal sixth-power scalar
  -> retain arity, owner gcd, residue phase, support-six native/padded status
  -> then attach exact M, endpoint owners, topology, cycle image, discharge route
```

## Relation to S242

This sidecar extends the S242 HYP-3060 Desargues/Beal collision split.  S242
already separated binary sixth-power equality as a Gaussian owner/factor gate
from ternary sixth-power equality as a diagonal-current carrier.  HYP-3076
pushes that split down into the relation-lattice ledger: ternary collisions are
native six-slot support data, while binary equalities remain rank-lowered
square-cube shadows unless a separate packet makes the padding legal.

It also neighbors the S243 HYP-3075 Hurwitz-Markov-Pell reservation.  Both
threads say the same quiet thing in different dialects: rare scalar
coincidences become theorem data only after their arithmetic address is kept.
For HYP-3075 that address is continued-fraction/Pell/carry structure; for
HYP-3076 it is support arity, owner gcd, residue phase, native-versus-padded
status, and discharge route.

## Bounded scout

The exact S244 scout computes sixth-power residue masks at the moduli already
active in the LRC14 proof stack:

```text
7, 9, 13, 14, 21, 27, 28, 41
```

It also runs bounded equality searches:

```text
2-vs-2 bound 220: nontrivial_hits=0
3-vs-3 bound 80: collision_pairs=5, primitive_pairs=3
```

The first primitive `3-vs-3` collision in this bounded scan is:

```text
3^6 + 19^6 + 22^6 = 10^6 + 15^6 + 23^6
sum = 160426514
gcd = 1
```

Its residue profile is exactly the kind of sidecar we need:

```text
mod 7:  all six terms are unit residues, so sixth powers all equal 1
mod 9:  both sides have mask 0^1 1^2
mod 13: both sides have mask 1^2 12^1
mod 27: both sides have mask 0^1 1^1 19^1
mod 41: equality survives, but masks differ: (21^2,32) vs (5,10,18)
```

Thus small moduli expose why the relation is locally legal, while modulus `41`
keeps a magnitude-sensitive phase.  That mirrors the LRC14 warning that
mod-14/mod-27 shell data must remain attached to exact Farey scale and
magnitude cocycle before any proof quotient is trusted.

## LRC14 transfer

Add the following sidecar fields to packet ledgers that touch relation
lattices, support-six tails, power-lift guards, or route-state closure:

```text
sixth_power_collision_type:
  native_3v3_support6 | rank_lowered_2v2_square_cube | padded_support6_cancel

sixth_power_owner_gcd
sixth_power_residue_mask_mod7
sixth_power_residue_mask_mod9
sixth_power_residue_mask_mod13
sixth_power_residue_mask_mod27
sixth_power_phase_mod41
native_support6_flag
degenerate_padding_pair
rank_lowered_square_cube_shadow
power_collision_discharge_route
```

The intended theorem-facing use is:

```text
If a residual LRC14 relation-lattice packet contains a native 3-vs-3
sixth-power collision, then route it through finite low-height wall ledger,
support-six coimage/cycle image, or named THM-572/F7 debt.

If it contains only a 2-vs-2 sixth-power equality, first mark the rank-lowered
square-cube shadow and the padded canceling pair; do not count it as a native
six-term obstruction.
```

This complements HYP-3060's Beal common-owner gate.  HYP-3060 says a primitive
three-channel collision must not forget owner/factor coordinates.  HYP-3076
adds that equal sixth-power collisions must also declare whether their six
terms are genuinely native or artificially padded.

## Tournament Analysis

Vertices are proof carriers and sidecar fields, not runners or raw bases.

Pairwise observable:

```text
retained LRC predicate
native support-six status
owner visibility
residue phase
low-height wall compatibility
route-state legality
degeneracy protection
```

The S244 tournament is transitive:

```text
labelled_packet_sheaf
> native_three_vs_three_support6_collision
> sixth_power_residue_phase_mask
> route_state_closure_sidecar
> low_height_wall_ledger
> owner_gcd_common_factor_gate
> padded_support6_canceling_pair
> rank_lowered_two_vs_two_square_cube_shadow
> raw_equal_sixth_power_scalar
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1,1]
hamiltonian_path_count=1
```

Assumption challenge: the tempting vertices were bases in the Diophantine
equations, runners, and raw equal sums.  Those are too lossy.  The usable
vertices are sidecar fields and proof obligations.  This quotient preserves
whether a relation can legally enter the support-six/cycle-image LRC layer and
destroys the exact base identity unless it is retained in the owner fields.

## Next pull

Run this sidecar on an actual HYP-2963/HYP-3074 packet sample.  For every row
with a relation-lattice or support-six field, emit:

```text
relation_lattice_covolume
low_height_wall_ledger
native_support6_flag
sixth_power_collision_type
sixth_power_residue_masks
owner_gcd_common_factor_gate
cycle_class_image_status
route_state_closure_status
power_collision_discharge_route
```

The proof obligation is then clean:

```text
native 3-vs-3 wall -> finite wall / cycle image / THM-572-F7 route
rank-lowered 2-vs-2 wall -> padded-degeneracy guard, not terminal debt
```
