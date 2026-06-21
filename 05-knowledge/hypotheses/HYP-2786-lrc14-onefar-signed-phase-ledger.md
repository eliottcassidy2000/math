---
id: HYP-2786
title: One-far binding error has a low-mode signed phase ledger, not a finite residue table
status: OPEN proof target; exact Abel value plus signed Fourier scout
source: codex-2026-06-21-S75
depends_on:
  - HYP-2779
  - HYP-2782
  - HYP-2784
  - HYP-2785
  - HYP-2720
  - THM-546
related:
  - HYP-2745
  - HYP-2772
  - HYP-2783
  - OPEN-Q-108
---

# HYP-2786: One-Far Signed Phase Ledger

## Claim

For the wide binding family

```text
E = {0,1,...,k-2} union {w},   k=8..12,   w>14,
```

the one-far error

```text
Delta_w = p0(E) - Phi({0,1,...,k-2})
```

is controlled by a small signed Fourier head plus a signed analytic tail, not
by an absolute Koksma/BV envelope.  This refines incoming HYP-2784: even the true
arc-complexity `V` leaves a large absolute-bound gap, so the missing theorem
must exploit signed phase cancellation.

This is compatible with incoming HYP-2785, which refutes the tempting
residue-only `w mod 7` table for single-far `Delta_w`.  The useful quotient is
not a finite table in the far runner.  It is a low-head phase ledger, with the
remaining `w`-dependence routed to Dedekind-reciprocity/equidistribution rather
than to an absolute tail bound:

```text
n = 1,2,3 dominate the dangerous positive head;
n mod 14 in {1,2} carries most L1 pressure;
7 | n modes vanish by the apex-prime sector coefficient.
```

Odd support is usually the larger L1 envelope, but not a signed cone.  The
k=11 binding row is even-led (`top_mod14=2`, odd-share about `0.395`), matching
HYP-2720's warning that odd support must be used before scalarizing, with a
finite even-led exception ledger.

## Exact Scout

Script:

```text
04-computation/lrc14_onefar_signed_phase_ledger_codex_s75.py
```

Stored output:

```text
05-knowledge/results/lrc14_onefar_signed_phase_ledger_codex_s75.out
```

The script uses the exact Abel endpoint identity from the Thread-B/HYP-2784
engine for `Delta_w`, then compares it with signed Fourier heads.  In the
scanned binding band `w in [15,100]`, the max positive rows are:

```text
k=8  w=21  Delta/margin=0.114075  top_mod14=1  odd_share=0.810518
k=9  w=68  Delta/margin=0.081192  top_mod14=1  odd_share=0.884733
k=10 w=22  Delta/margin=0.104288  top_mod14=1  odd_share=0.622080
k=11 w=16  Delta/margin=0.083918  top_mod14=2  odd_share=0.395081
k=12 w=71  Delta/margin=0.042799  top_mod14=1  odd_share=0.735831
```

The exact error is only about `5%..16%` of the absolute THM-546 BV bound on
these risk rows.  The first two Fourier blocks already capture the signed
value well: over `w in [15,80]`, the worst `n<=13` head error is below `0.008`
of the relevant cap margin for every k=8..12.

## Proof Obligation

The next theorem should split the one-far binding branch into:

```text
P1. finite signed head n<=13, grouped by n mod 14;
P2. apex-prime vanishing of n==0 mod 7;
P3. finite odd/even exception ledger, with k=11/even-led as the first model;
P4. signed Dedekind/equidistribution tail small enough after P1-P3, not before.
```

This is narrower than an absolute Koksma proof and also narrower than a
residue-table proof.  It keeps HYP-2720's odd support idea in the right order:
odd L1 support supplies an envelope for most rows, but signed even-led packets
must be retained until the mod-14 head ledger is evaluated.

## Tournament Analysis

Vertices are the five binding rows at their scanned max positive `Delta_w`.
Pairwise observable:

```text
larger Delta_w/margin, then larger odd L1 share, then smaller BV ratio.
```

Switch/gauge: scalarize the exact Abel endpoint value only after splitting the
Fourier head by mod-14 phase and odd/even support, and only after separating
the HYP-2785 Dedekind/equidistribution tail from the finite head.  The
tournament is transitive:

```text
k=8,w=21 > k=10,w=22 > k=11,w=16 > k=9,w=68 > k=12,w=71
```

Fingerprint: score histogram `{0:1,1:1,2:1,3:1,4:1}`, no directed 3-cycles, one
Hamiltonian pressure path.

## Challenged Assumption

Rejected assumption: after HYP-2784, a sharper absolute `V`, raw BV constant,
or finite residue-only table is the remaining wide-bound target.

Preserved predicate: the exact one-far `Delta_w` and the cap-margin comparison.

Destroyed information: none that affects `Delta_w`; the split only changes
basis from absolute arc counts to signed Fourier phase packets.

What remains: prove the finite head ledger symbolically and attach the signed
Dedekind/equidistribution tail bound indicated by HYP-2785.  This is a direct
OPEN-Q-108 subproblem.
