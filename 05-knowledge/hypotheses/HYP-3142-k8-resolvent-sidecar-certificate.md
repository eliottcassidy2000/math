---
id: HYP-3142
title: The k=8 hard node may close by an exact 4-moment resolvent sidecar
status: EVIDENCE / exact bounded-bank scout through B=14; not a proof
source: codex-2026-06-27-S273
tangent: T1207
technique: LTI-268
tournament_technique: LTT-166
script: 04-computation/lrc14_k8_resolvent_sidecar_scout_codex_s273.py
results:
  - 05-knowledge/results/lrc14_k8_resolvent_sidecar_scout_codex_s273.out
  - 05-knowledge/results/lrc14_k8_resolvent_sidecar_scout_codex_s273_B14.out
related:
  - HYP-3141
  - HYP-3140
  - HYP-3139
  - HYP-3138
  - HYP-3137
  - HYP-3136
  - HYP-3135
  - HYP-3134
  - HYP-3133
  - HYP-3132
  - HYP-3131
  - HYP-3129
  - HYP-3122
  - HYP-3119
  - HYP-3118
  - HYP-3113
  - HYP-3111
  - HYP-3110
  - THM-577
  - OPEN-Q-108
---

# HYP-3142: k=8 Resolvent Sidecar Certificate

## Claim

The current LRC14 frontier says the covering proof has collapsed to one
bounded-core node: the k=8 dip.  HYP-3132 explains that node as the
reflection-even De Moivre resolvent

```text
g(t) = (t-1)(t-2)(t-4)(t-5),
u=t-3,
g(u+3) = (u^2-1)(u^2-4) = u^4 - 5u^2 + 4.
```

HYP-3137 asks which generating-function payload survives scalar evaluation;
HYP-3138 tests the HYP-3132 even reflection fold as a legal quotient only
with an odd-coordinate resurrection table.  HYP-3139 sharpens the fold at
matrix level by splitting the k=8 reflection block into inner shell, center,
antisymmetric, and boundary-leakage pages.  HYP-3136 factors the covering
branch into the `Q` floor, `R-safe` floor, and signed `Rprime` coupling.
Incoming HYP-3140 turns that `Rprime` factor into a finite 14-sheet fiber-PGF
conditional first-moment packet.  Incoming HYP-3135 says that this branch must
keep the middle resolvent/SPEC payload before quotienting to a raw scalar.
Incoming HYP-3141 then recasts directed tournament edges as tip/tail
proof-information packets with an `edge_bounded_core_floor_exit` field.
HYP-3142 tests the terminal k=8 subpacket that these cards now need:

```text
k=8 terminal packet =
  exact 4-moment U4 cap certificate
  + biquadratic resolvent fold
  + Bravais full-residue-flatness sidecar
  + Savitch/coordinate-resurrection repair depth
  + A000568 edge-global-consistency guard.
```

The conjectural theorem is now concrete:

```text
For every primitive k=8 bounded-core shape E,
U4(E) <= U4(consec_8) = 2633/7350 < cap_8 = 2243/5880.
```

If this moment-majorization theorem is proved, the single remaining k=8 node
closes by exact rational arithmetic, with slack

```text
cap_8 - U4(consec_8) = 683/29400.
```

In the integrated LRC14 route, this is the exact bounded-core exit used after
HYP-3141 has made tip/tail edges into proof-information packets, HYP-3140 has
converted the `Rprime` scalar into a finite fiber-PGF conditional moment
theorem, HYP-3137 has named the generating-function payload, HYP-3138 has
named the reflection-fold adjoint repair, HYP-3139 has split the reflection
block into proof pages, HYP-3136 has reduced the multi-far floor to a finite
constant chase, and incoming HYP-3135 has preserved the signed SPEC / De
Moivre middle layer.

## Evidence

The executable scout
`04-computation/lrc14_k8_resolvent_sidecar_scout_codex_s273.py` computes, for
primitive bounded shapes `E={0<e_2<...<e_8<=B}`:

- the exact miss-count distribution `q_i`;
- the exact 4-binomial-moment relaxation `U4(E)`;
- the nearest PGF root and real-root count;
- Bravais residue DFT peak, entropy, and mirror defect;
- the fourth cumulant `kappa4`;
- a Savitch-style missing-sidecar depth;
- Tournament Analysis over certificate sidecars.

Stored run at `B=13`:

```text
primitive rows scanned = 1716
U4_over_cap_count = 0
worst row = E=(0,1,2,3,4,5,6,7)
U4 = 2633/7350
cap_8 - U4 = 683/29400
corr(U4,q0) = +0.997660
corr(U4,nearest_root) = +0.905012
corr(U4,bravais_peak) = -0.547962
corr(U4,kappa4) = -0.566836
```

Stress run at `B=14`:

```text
primitive rows scanned = 3431
U4_over_cap_count = 0
same worst row: E=(0,1,2,3,4,5,6,7)
```

The worst rows are all full-residue-flat in the Bravais sense:
residue counts `(2,1,1,1,1,1,1)`, DFT peak `1/8`, entropy `1`, mirror defect
`0`.  This is the older HYP-3113 Bravais lesson in a proof-facing place:
high k=8 coverage is not a Bragg peak; it is reciprocal-flat.

## Interpretation

This merges several previously separate motifs into one possible terminal
certificate.

- **De Moivre / HYP-3132:** the hard quartic is not generic; after reflection
  it is degree 2 in `u^2`, with square discriminant `9`.
- **Edge/fiber/GF payload and fold/block repair /
  HYP-3141-HYP-3140-HYP-3137-HYP-3139:** HYP-3141 says the proof edge must
  carry tail payload, tip payload, commutator defect, coordinate repair, and
  terminal exit.  HYP-3140 keeps the 14-sheet `Rprime` conditional first
  moment, while this note keeps the k=8 `U4` coefficient/cumulant payload.
  The even fold and reflection block are legal only after odd, center,
  antisymmetric, and boundary-leakage pages have resurrection tables or named
  debt.
- **Resolvent packet / HYP-3135-HYP-3136:** the exact `U4` inequality is the
  bounded-core terminal gate for the multi-far floor packet, not a competing
  scalar shortcut to HYP-3140's fiber-PGF floor.
- **phi4 / HYP-3122:** the consecutive row has the negative fourth cumulant
  stabilizer, but many positive-`kappa4` rows still close by `U4`, so phi4 is
  the explanation of the extremal stalk, not the whole certificate.
- **Bravais / HYP-3113:** the extremal stalk is residue-flat and mirror-even.
  Non-flat rows gain slack in the bounded scans.
- **Savitch / HYP-3118:** rows that fail a side signal still have shallow
  repair depth after the exact `U4` terminal gate is tried.
- **A000568 / HYP-3133-HYP-3134:** the packet must retain the edge/global
  consistency sidecar before quotienting away local child data.

## Proof Target

Turn the bounded-bank evidence into a theorem:

1. Prove an exact moment-majorization inequality:
   `U4(E) <= U4(consec_8)` for every primitive k=8 bounded-core shape.
2. Use the HYP-3132 reflection fold to rewrite the remaining fourth-order
   obligation as the coefficient inequality of `u^4 - 5u^2 + 4`.
3. Splice this inequality into the
   HYP-3141-HYP-3140-HYP-3139-HYP-3138-HYP-3137-HYP-3136-HYP-3135 packet theorem as the
   `bounded_core_U4_exit` field of the signed SPEC / Lee-Yang / edge ledger.
4. Prove strict improvement for non-full-residue-flat Bravais spectra, or
   route exceptions to a finite Hensel/CRT `2x7` resonance sidecar.
5. Attach `k8_U4_cap_certificate`, `biquadratic_resolvent_fold_status`,
   `bravais_full_residue_flatness`, `savitch_repair_depth`,
   `a000568_global_consistency_class`, and `terminal_exit_or_named_debt` to
   the k=8 packet ledger.

## Guardrail

This is not yet a proof of LRC14.  The script verifies bounded banks through
`B=14`; it does not prove the global moment-majorization theorem, nor does it
prove that every unbounded k=8 shape reduces to a checked bounded stalk.  Its
value is that the remaining theorem has become a single exact inequality with
a visible extremizer and multiple independent sidecar explanations.
