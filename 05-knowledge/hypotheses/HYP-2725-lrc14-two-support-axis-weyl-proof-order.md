---
id: HYP-2725
title: LRC14 Weyl error needs a two-support-axis proof order
status: OPEN; synthesis with exact S69 triage
source: codex-2026-06-21
depends_on:
  - HYP-2724
  - HYP-2723
  - HYP-2722
  - HYP-2721
  - HYP-2720
  - HYP-2719
  - HYP-2718
  - THM-561
  - HYP-2717
related:
  - THM-559
  - THM-560
  - HYP-2698
  - THM-558
  - OPEN-Q-108
---

# HYP-2725 - Two Support Axes For The LRC14 Weyl Error

## Claim

The phrase "Weyl error is odd-support dominated" is correct only after
separating two different support axes.

1. **Relation-support axis.**  In the Fourier relation lattice,

   ```text
   corr(E) = sum_{0 != n in Lambda(E)} K(n),
   ```

   support means `|supp(n)|`.  On this axis, odd/even support parity is not a
   signed cancellation rule.  Incoming S9 work records the exact obstruction:

   ```text
   K(-n)=conj(K(n)),
   ```

   so relation reverse-pairs reinforce as `2 Re K(n)`.  The useful structure is
   support size, not parity: support-2 is the cut/2-body layer, support-3 Schur
   packets are the first cycle-like layer, and higher supports are the tail.

2. **Factorial-support axis.**  In the HYP-2718 origin-atom basis,

   ```text
   Q_0 = ProductCover - ActualCover
       = sum_j (-1)^j W_j
       = sum_even W_j - sum_odd W_j,
   ```

   support means the missed-sector subset size `j`.  On this axis, HYP-2720's
   exact S69 scout still supports an odd `L1` envelope for the high-height tail:

   ```text
   OddAbs = |W_1|+|W_3|+|W_5|
   EvenAbs = |W_0|+|W_2|+|W_4|+|W_6|.
   ```

The proof order should therefore be:

```text
relation-support filter first
-> factorial odd-L1 envelope second
-> finite signed-even-led ledger
-> evaluate Q_0.
```

## Exact Triage

The scout

```text
04-computation/lrc14_two_support_axes_codex_20260621.py
05-knowledge/results/lrc14_two_support_axes_codex_20260621.out
```

reuses S69's exact rational row bank and classifies rows by their factorial
support behavior.

Across the ten-row bank:

```text
odd_L1_signed_phase: 2
odd_L1_envelope: 5
finite_even_L1_ledger: 3
```

The three rows where the factorial odd envelope fails are not proof killers;
they are exactly the rows that should be routed to the finite low-height ledger:

```text
5+3 split
3+3+2 split
seven singleton carriers
```

The largest cap-risk row in that finite-even-led class is the singleton carrier
row, with

```text
|Q_0|/(cap-product) = 910424533191000/29642787224726503
                    ~= 0.030713189.
```

The two signed-odd phases are positive-`Q_0` high-relation two-carrier rows.
They are useful diagnostics: they show where odd support actually flips the
origin atom, but their cap-risk ratios are only about `0.0136` and `0.0120`.

## Proof Route

1. **Relation-support filter.**  Split the Weyl relation lattice by support
   size.  Keep the S9 correction: no reverse-pair cancellation.  Support-3
   Schur packets are the first binding layer.

2. **Finite Schur/low-height ledger.**  Route support-3 Schur packets,
   same-difference support-4 packets, and small-denominator relation classes
   to the HYP-2717/HYP-2714 finite ledger.

3. **Factorial odd envelope on the tail.**  After the relation-support filter,
   prove

   ```text
   EvenAbs_tail <= OddAbs_tail + finite_ledger_tax.
   ```

   This is the precise surviving form of the user's odd-support intuition.

4. **Origin atom evaluation.**  Only after those filters evaluate

   ```text
   Q_0 = even_support - odd_support.
   ```

   Signed odd dominance is not needed, and in the current data it is false for
   most negative `Q_0` rows.

## Fit With HYP-2720 Through HYP-2724

HYP-2720 is the factorial odd-support envelope.  HYP-2721 is the `Q_0`
Vitali atom-cone boundary view.  HYP-2722 is the generated miss-zeta word
compatibility filter that blocks cheap abstract `Q_0`-hiding cone moves.
HYP-2723 splits depth-law convex order from generated-word compression, and
HYP-2724 shows low-support relation counts explain only part of `corr`, so the
tail is not disposable.  HYP-2725 is the proof-order wrapper: first use
HYP-2719's relation-support filter, then HYP-2720's factorial odd envelope,
then HYP-2721/HYP-2722/HYP-2723's boundary atom-cone, generated-word, and
depth-law lenses, with HYP-2724 warning against replacing the tail by a finite
low-support count.

## Assumption Challenge and Tournament Analysis

This synthesis considered relation vectors, relation supports, missed-sector
support layers, residual masks, Fourier modes, Schur packets, same-difference
packets, rows, and proof obligations as possible vertices.

Productive vertices are proof obligations:

```text
relation_support_filter
support3_schur_ledger
factorial_odd_L1_tail
finite_even_led_packets
origin_atom_Q0
raw_odd_relation_parity
reverse_pair_cancellation
```

Pairwise observable: a vertex beats another if it preserves the HYP-2718
predicate while avoiding a known false cancellation.  The pressure path is

```text
relation_support_filter
> support3_schur_ledger
> factorial_odd_L1_tail
> finite_even_led_packets
> origin_atom_Q0
> raw_odd_relation_parity
> reverse_pair_cancellation.
```

The quotient preserves the origin-atom proof target and destroys the single
meaning of "support."  That ambiguity must stay explicit.

Challenged assumption: odd support should mean one parity rule everywhere.  It
does not.  Relation-support parity is refuted as a signed rule; factorial
support oddness remains a usable envelope after relation filtering.
