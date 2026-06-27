---
id: HYP-3028
title: LRC14 residual status-gate switchboard after coarse ET/Henselian convergence
status: EVIDENCE / proof-interface switchboard; not a proof
source: codex-2026-06-26-S192
tangent: T1109
script: 04-computation/lrc14_residual_status_gate_switchboard_codex_s192.py
result: 05-knowledge/results/lrc14_residual_status_gate_switchboard_codex_s192.out
related:
  - HYP-3026
  - HYP-3025
  - HYP-3024
  - HYP-3023
  - HYP-3020
  - HYP-3017
  - HYP-3016
  - HYP-3015
  - HYP-2963
  - THM-572
  - LTI-176
  - LTI-173
  - LTI-171
  - LTT-074
  - LTT-071
  - LTT-070
  - OPEN-Q-108
---

# HYP-3028: Residual Status-Gate Switchboard

## Claim

After HYP-3024, the best next proof angle is status-first rather than
route-first:

```text
coarse ET + Henselian unit gate
  should prove boundary/open status convergence,
  while the remaining open-route collisions become certificate scheduling debt.
```

This is not a proof of LRC14.  It is a proof-interface target that narrows the
load-bearing theorem to boundary/open purity inside automatic/residue fibers.
Full theorem-route purity can be handled afterward by q-witness, safe-stick,
Fejer/Haar, petal, K33/F7/THM-572, covering, or magnitude-cocycle certificates.

## Evidence

The S192 script parses the stored full-bank HYP-3024/S188 and named-bank
HYP-3026/S189 outputs rather than recomputing the expensive HYP-2963 exact
packet bank.

Key imported S188 facts:

```text
full_bank_coarse_et_unit_gate fibers=21702 mixed_route=15 mixed_status=0 max_mixed=4
full_bank_magnitude_cocycle   fibers=21909 mixed_route=0  mixed_status=0 max_mixed=0
target_word_coarse_gate       fibers=585   mixed_route=9  largest_mixed=4
```

Key imported S189 fact:

```text
named automatic_plus_barcode_shape fibers=11 mixed_status=0 mixed_route=0
```

The largest displayed S188 residual mixed-route fiber is already harmless for
LRC14 status:

```text
single swap 13->143  COVERING-MOMENT  M=11/144  mu=4307/194040  q0=14
single swap 13->31   Q-WITNESS        M=1/13    mu=194039/6015240 q0=13
single swap 13->59   Q-WITNESS        M=1/13    mu=53509/1635480  q0=13
single swap 13->115  Q-WITNESS        M=1/13    mu=13469/405720   q0=13
```

The Q rows have direct denominator witnesses.  The covering row is strict-open
with positive safe mass.  The route labels differ, but the LRC predicate does
not.

## Theorem Target

Prove the following before trying to prove full route purity:

```text
In every automatic/residue fiber, the coarse ET+Henselian-unit data cannot
identify an AP/GW boundary equality atom with a strict-open packet.
```

Once that is known, any residual route-mixed fiber is a scheduling problem:

```text
q<=13 direct witness
AP/GW zero-bar equality
strict safe-stick / barcode / Fejer / Haar certificate
unit-petal discharge
K33/F7/THM-572 state-lift debt
covering-moment certificate
magnitude-cocycle family formula
```

## Tournament Analysis

Vertices are proof obligations after the residual gate, not runners.  The
observable is:

```text
predicate, compression, route_split, arithmetic, topology, proof_cost, anti_address
```

The S192 fingerprint is:

```text
score_hist={0: 1, 3: 4, 4: 1, 5: 1}
directed_3cycles=7
scc_sizes=[6, 1]
hamiltonian_path_count=33
score_order=status_first_coarse_et_unit_gate >
  barcode_packet_zipper_address >
  exact_et_address_clock >
  henselian_unit_lift_rule >
  magnitude_cocycle_fallback >
  safe_stick_barcode_certificate >
  raw_automatic_word
```

The nontrivial SCC is useful: it says the post-status certificate layer is not
linearly ordered.  Magnitude, barcode, ET address, and Hensel unit data should
remain a switchboard, not a scalar ranking.

## Next Hook

Add a cached packet-ledger mode to HYP-2963/S181 so the 15 S188 residual
coarse-gate fibers can be listed without recomputing exact `M`.  Then attach
the first successful certificate tooth to each residual fiber.
