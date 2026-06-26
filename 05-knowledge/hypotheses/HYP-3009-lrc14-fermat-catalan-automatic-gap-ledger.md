---
id: HYP-3009
title: LRC14 Fermat-Catalan automatic gap extension
status: SYNTHESIS / packet-schema guardrail and proof-interface carrier; not a proof
source: codex-2026-06-25-S178
script: 04-computation/lrc14_fermat_catalan_automatic_gap_s178.py
result: 05-knowledge/results/lrc14_fermat_catalan_automatic_gap_s178.out
related:
  - HYP-3007
  - HYP-3008
  - HYP-3006
  - HYP-3003
  - HYP-3002
  - HYP-3000
  - HYP-2999
  - HYP-2998
  - HYP-2997
  - HYP-2996
  - HYP-2995
  - HYP-2963
  - HYP-2950
  - HYP-2944
  - HYP-2937
  - HYP-2702
  - HYP-2698
  - HYP-1920
  - HYP-1902
  - THM-572
  - OPEN-Q-108
---

# HYP-3009: LRC14 Fermat-Catalan Automatic Gap Extension

## Claim

The Fermat-Catalan / Moser-de Bruijn / fibbinary / gap-theorem cluster is useful
for LRC14 only as a labelled packet ledger.  It is not a scalar shortcut.

The proposed addition to the HYP-2963 packet schema is:

```text
automatic_language_class
fibbinary_carry_status
moser_even_bit_status
ostrowski_digit_system
lacunary_gap_ratio
power_lift_guard
fermat_catalan_residual
hurwitz_doubling_cf_state
visibility_potato_approx_guard
```

The guiding rule:

```text
forgetting an automatic/gap/power coordinate is legal only if the LRC predicate
is fiber-constant, the coordinate is reconstructible, a dual certificate
annihilates it, or the packet is routed to AP/GW, C27, K33, Res_27, covering,
or named F7/THM-572 residual debt.
```

## Computation

Script:

```text
04-computation/lrc14_fermat_catalan_automatic_gap_s178.py
```

Output:

```text
05-knowledge/results/lrc14_fermat_catalan_automatic_gap_s178.out
```

The script audits named rows by three integer languages:

```text
M = Moser-de Bruijn / even binary positions / carry-free base-4
F = fibbinary but not Moser / no adjacent binary ones
C = adjacent carry present
```

Selected readout:

```text
AP13:          word=MFCMMCCFFFCCC, counts M=3,F=4,C=6, tail_gap=13/12
GW_12_to_24:  word=MFCMMCCFFFCCC, counts M=3,F=4,C=6, tail_gap=24/13
K33_12_to_36: word=MFCMMCCFFFCCF, counts M=3,F=5,C=5, tail_gap=36/13
petal_10_20:  word=MFCMMCCFFCCCM, counts M=4,F=3,C=6, tail_gap=20/13
cover_12_84:  word=MFCMMCCFFFCCM, counts M=4,F=4,C=5, tail_gap=84/13
Res_27 probe: word=CCMCCFCCCCCCM, counts M=2,F=1,C=10, tail_gap=20/7
```

Two things matter here.

First, AP and GW share the same automatic language word.  So the Moser/fibbinary
shadow alone cannot distinguish the two equality atoms, just as Haar/Baire
boundary owners alone cannot distinguish the Goddyn-Wong hidden transfer.

Second, the non-AP/GW rows start to split into distinct automatic packets: K33
adds fibbinary tail, C27 petals add Moser tail, covering rows add lacunary
Moser tail, and the Res_27 probe is carry-dense.  This is not yet a theorem,
but it supplies a concrete field set for familywise routing.

This extends the incoming checkpointed HYP-3008 automatic-gap carrier and the
HYP-2937/HYP-2944/HYP-2950 family: C27/unital transfer blocks,
Farey-product perfect-number gates, and Borel-Baire-Haar witness labels are
upstream packet routes.  HYP-3009 adds the Fermat-Catalan, lacunary, and
visibility side labels those routes must carry before a sequence or power
analogy is allowed to scalarize them.

## Power-Lift Guard

The LRC unit-excess payloads

```text
M = p/(14p-1), q=14p-1, p+q=15p-1, p*q=14p^2-p
```

have small perfect-power stress points:

```text
p=2: q=27=3^3
p=4: p=4=2^2
p=8: p=8=2^3
```

The live rule is negative:

```text
perfect-power or power-sum coincidences are no-lift guards, not proof
certificates.
```

A Fermat-Catalan style equation can help only after its cyclotomic, p-adic, and
packet labels are retained.  Otherwise it is just another scalar shadow that
can mix AP/GW boundary atoms with C27, K33, and covering residuals.

## Gap-Theorem Guard

The Ostrowski-Hadamard gap theorem says a Hadamard-lacunary power-series support
cannot be analytically continued across the unit circle.  The LRC translation is
not literal analytic continuation; it is a guardrail:

```text
if a residual frequency or endpoint-debt support becomes lacunary, Fejer/Haar
smoothing may not erase it unless the lacunary address is carried or certified
irrelevant on that packet fiber.
```

This is most relevant for large-tail covering rows and q-witness bypasses.  A
large tail can look harmless to a scalar smoother while still carrying an
address coordinate that blocks a familywise certificate.

## Tournament Analysis

Vertices are proof carriers, not runners:

```text
labelled_packet_bank
fibbinary_zeckendorf_fsm
moser_square_x_plus_2y
ostrowski_hadamard_gap
fermat_catalan_power_guard
hurwitz_doubling_cf_state
potato_visibility_guard
raw_scalar_speed_word
```

Pairwise observable:

```text
preserves_lrc_predicate
finite_state_checkable
keeps_packet_labels
guards_power_or_carry_lift
connects_to_existing_routes
resists_scalarization
```

The conservative gauge is transitive:

```text
score_histogram={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3cycles=0
Hamiltonian path:
  labelled_packet_bank
  > moser_square_x_plus_2y
  > fibbinary_zeckendorf_fsm
  > ostrowski_hadamard_gap
  > fermat_catalan_power_guard
  > hurwitz_doubling_cf_state
  > potato_visibility_guard
  > raw_scalar_speed_word
```

This deliberately challenges the runner-vertex assumption.  The better vertex
sets are packet languages, gap supports, power-lift guards, continued-fraction
states, visibility guards, and residual proof obligations.

## Proof Targets

1. **Automatic carrier lemma.**  A HYP-2963 quotient that uses Moser,
   fibbinary, Zeckendorf, or recurrence shadows must keep the finite-automaton
   state or prove the LRC predicate is constant on its fiber.

2. **Moser sublanguage lemma.**  The even-bit / carry-free square channel is a
   strict sublanguage of fibbinary.  On this subfiber, the unique `x+2y`
   decomposition should prevent hidden square-carry lifts unless the missing
   odd-bit coordinate is explicitly reattached.

3. **Fibbinary/Ostrowski lemma.**  Endpoint debt that has no adjacent carries
   should admit a finite automaton compatible with HYP-1902/HYP-1920's forced
   odd fan plus even-bridge picture.

4. **Lacunary boundary lemma.**  Hadamard-lacunary residual supports require a
   boundary/frontier label before smoothing can treat them as harmless tails.

5. **Fermat-Catalan power-lift guard.**  Perfect-power payloads and power-sum
   equalities must be certified by cyclotomic/p-adic labels or routed to named
   residual debt.

The closest next implementation task is to add these fields to the labelled
packet classifier, then compute whether automatic words are route-pure within
the HYP-2963 bank.
