---
id: HYP-3233
title: Signed address chart-change sheaf for the Mobius, Eisenstein, and Legendre recursions
status: SYNTHESIS / proof-carrier extension; not an LRC14 proof
source: codex-2026-06-28
tangent: T1331
technique: LTI-331
tournament_technique: LTT-231
reflection: 07-reflections/signed-address-chart-change-sheaf-for-the-three-recursion-modes-codex-20260628.md
related:
  - HYP-3232
  - HYP-3231
  - HYP-3230
  - HYP-3216
  - HYP-3004
  - HYP-2902
  - HYP-2901
  - HYP-2899
  - HYP-2704
  - HYP-2685
  - HYP-2681
  - THM-553
  - THM-550
  - THM-549
  - THM-442
  - OPEN-Q-108
---

# HYP-3233: Signed Address Chart-Change Sheaf

## Claim

The recurring signed recursions

```text
full:       A+B+C-D-E-F+G
even half:  A+B-C
odd half:   A+B-C+D-E-F+G
```

are not three loose analogies.  They are local charts of one signed address
sheaf.  Each chart is a Mobius or character transform on a finite carrier, and
each chart is legal only while its local slot labels are retained.

The main guardrail is:

```text
A..G are not global letters.
They are local address slots in the current chart.
```

The same scalar sign word can mean different proof data depending on the chart:

- **Full Mobius chart:** `A+B+C-D-E-F+G` is the full staircase `B3`
  inclusion-exclusion from THM-442.  It is a cell-affine / third-difference
  identity, not a Hamiltonian-path or cap-risk recurrence.
- **Even Eisenstein chart:** `A+B-C` is the complement-folded `B2` edge from
  THM-549/THM-550.  It lives on the even pronic half-domain and matches the
  pure modulus fold in HYP-3232 when the cap kernel is still scale-covariant.
- **Odd Legendre chart:** prompt order gives `A+B-C+D-E-F+G`, but the corrected
  geometric Venn chart is `A+B+D-C-E-F+G`: corners `A,D,B`, edges
  `A+B-C`, `A+D-E`, `B+D-F`, center `G`.  The slots `C` and `D` have the same
  size and cancel only after scalarization; they must remain distinct proof
  addresses.

Thus the hidden similarity is not the literal sign string.  It is the repeated
operation:

```text
choose chart
-> name local generators
-> subtract overlap / destroyed coordinate
-> retain fixed-line, apex, exact-period, or phase sidecar
-> recurse only after the chart-change debt is explicit.
```

## Relation To The Current LRC Route

This extends HYP-3231 and HYP-3232.  HYP-3231 says the route is scale-normal:
normalize, expose the coordinate scale does not kill, attach the sidecar, and
recurse.  HYP-3232 says the three modes concentrate at the apex:

```text
Mobius:     moment order
Eisenstein: modulus fold
Legendre:   speed ratio / three-gap kernel
```

HYP-3233 adds the chart-change law: before any one of those modes is used, the
packet must record which signed address chart it is in and which coordinates
are destroyed by moving to the next chart.

The proof packet should therefore carry:

```text
scale_orbit
signed_address_chart
local_slot_basis
slot_size_vector
character_word
chart_change_map
cancelled_same_size_slots
fixed_line_or_apex_coordinate
apex_break_defect
denominator_exact_period_packet
moment_depth_target
destroyed_coordinate_sidecar
terminal_discharge_or_debt
```

## Similarities Across The Recursions

### 1. Every chart is an inclusion-exclusion carrier

The full formula is ordinary `B3` Mobius.  The even half is a degenerate `B2`
edge.  The odd half is simultaneously a prompt-order `chi_7` word
`++-+--+` and a geometric `B3` Venn word after the corner basis is changed from
`A,B,C` to `A,D,B`.  HYP-2899 gives the companion denominator side:
`phi=mu*id` is the divisor Mobius chart for exact-period packets.

The common object is therefore a product of address lattices, not a scalar:

```text
Div(D) x scale_orbit x parity_chart x Boolean/character_slot x moment_depth.
```

### 2. Negative terms are overlap taxes, not proof negativity

In the full chart, `D,E,F` are pairwise overlaps.  In the even chart, `C` is the
one overlap of the two retained half-carriers.  In the odd chart, `C` is the
`A cap B` edge while `D` is a new corner.  In HYP-2681, the same seven-packet
shape becomes a pair-tax shadow `H(1)-2(D+E+F)`, while the actual residual is
`A+B+C+D+E+F+G`.

So a minus sign means "this coordinate was counted twice in this chart."  It
does not mean "this packet is bad" until the chart and sidecar are fixed.

### 3. Same-size cancellation is a danger signal

The odd chart has `C,D` both of size `N-2`, so scalar size projection sees
zero at `N-2`.  That zero is not absence of structure.  It is the statement
that a corner and an overlap have been confused by scalarization.

This is the same warning as HYP-3231's scale ledger: a quotient is legal only
when the forgotten coordinate is reconstructed, dual-annihilated, descended,
boundary-stopped, or named as residual debt.

### 4. The apex is the chart-change singularity

HYP-3232 verifies that the cap kernel obeys the clean Eisenstein fold
`K^(2n)=K^(n)/2` while `n >= 2*max(speed)`, and breaks at the apex `n/2`.
That is exactly the point where the even edge chart stops being enough and
the odd Legendre / three-gap / moment-depth chart becomes visible.

For LRC14:

```text
speeds <= 7:  pure scale, low-depth, even-fold chart
speeds 8..13: apex-break defect, Legendre chart, binding constants
```

The proposed defect sidecar is:

```text
apex_break_defect(a,b,n) = K^(n)(a,b) - (2/n)*h(a,b),
```

with the convention that it is zero on the modulus-covariant side and named
debt on the antipode side.

### 5. Survival and cube-root packets are quotient charts

HYP-2704 shows that survival middle mass is not literally one of the seven
sign words.  It is a tail-weighted quotient that keeps live depth cuts
`G1,G5,G6` and forgets sector labels.  HYP-2681 shows that the cube-root
refinement keeps the cyclic phase coordinates

```text
S_omega = A + omega B + omega^2 C,
P_omega = D + omega E + omega^2 F.
```

Both are chart changes.  They are useful only when the destroyed labels and
phase addresses are recorded.

## Proof-Route Subtargets

1. **Even chart discharge.**  Prove that the even/pronic `A+B-C` chart plus
   the scale-normal fold handles every packet whose cap-kernel
   `apex_break_defect` is zero.
2. **Odd chart finite atlas.**  On the nonzero apex-break side, keep the
   corrected Legendre Venn addresses `A,D,B; C,E,F; G` and route them through
   HYP-3230's three-gap/Stern-Brocot kernel plus HYP-3216's moment-depth
   ladder.
3. **Product-lattice denominator split.**  Keep `Div(D) x chart` until
   divisor-killed denominators, CRT defects, and exact-period unit capacity
   are separated.
4. **Survival quotient lift.**  When using the HYP-2704 death-chain
   boundary, lift near-equality rows back to the seven-slot chart before
   bounding signed resonant deviation.

## Tournament Analysis

Vertices are proof charts, not runners or arcs:

```text
scale_normal_packet_sheaf
cap_kernel_modulus_covariance
odd_legendre_venn_chart
even_eisenstein_edge_chart
full_mobius_B3_chart
product_divisor_mobius_chart
cube_root_phase_chart
survival_depth_quotient
raw_scalar_sign_formula
```

Pairwise observable:

```text
(preserves LRC predicate,
 names local slots,
 exposes chart-change debt,
 detects apex or fixed-line singularity,
 keeps exact-period/moment-depth sidecars,
 formal proof maturity).
```

Switch/gauge: lexicographic comparison of that tuple.  This gives a synthesis
tournament with score histogram

```text
{0:1, 1:1, 2:1, 3:1, 4:1, 5:1, 6:1, 7:1, 8:1}
```

no directed 3-cycles, nine SCCs of size one, and one Hamiltonian path:

```text
scale_normal_packet_sheaf
-> cap_kernel_modulus_covariance
-> odd_legendre_venn_chart
-> even_eisenstein_edge_chart
-> full_mobius_B3_chart
-> product_divisor_mobius_chart
-> cube_root_phase_chart
-> survival_depth_quotient
-> raw_scalar_sign_formula.
```

The ranking is a workflow ordering, not an importance theorem.

## Assumption Challenge

Do not assume the tournament vertices are runners, arcs, or even fixed
subtournaments.  Alternate vertices considered here are:

- signed address charts;
- local Venn slots;
- fixed-line half-tiling cells;
- apex break events;
- exact-period denominator packets;
- cube-root phase modes;
- survival depth cuts;
- proof obligations.

Preserved predicate: the scale-normal LRC proof obligation plus enough local
slot data to reconstruct or legally forget the signed recurrence.  Destroyed
information: global letter identity, scalar size projection, raw sign word,
and packet ownership under chart changes.  Challenged assumption: one can
reuse `A..G` across modes without a chart sidecar.

No proof of LRC14 is claimed.  The contribution is a sharper recursion
interface: signed recurrences become legal proof carriers only as charted,
sidecar-completed packet maps.
