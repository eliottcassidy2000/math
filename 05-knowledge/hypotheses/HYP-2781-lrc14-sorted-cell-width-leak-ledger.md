---
id: HYP-2781
title: Independent sorted cell-width leak ledger for the HYP-2780 joint coupling
status: OPEN proof target; exact bounded scout
source: codex-2026-06-21-S75
depends_on:
  - HYP-2749
  - HYP-2770
  - HYP-2773
  - HYP-2778
  - HYP-2779
  - HYP-2780
  - THM-534
related:
  - HYP-2760
  - HYP-2761
  - HYP-2764
  - HYP-2751
  - HYP-2744
  - HYP-2726
  - OPEN-Q-108
---

# HYP-2781: Sorted Cell-Width Leak Ledger

## Claim

For the bounded, full-residue LRC14 stratum, the right compression object is
not a fixed cell, a fixed cyclic cell window, or a per-cell Huffer-Shepp
inequality.  The useful quotient first forgets the labels of the six nonzero
`Z/7` cells and compares the sorted survival-width vector

```text
sort_desc(W_1(E),...,W_6(E)).
```

In exact bounded scans, the AP/consec row weakly majorizes this sorted width
vector for all k=8 full-residue shapes and almost all k=9,10 shapes.  The only
target-range leaks seen by this scout are one-hole AP-extension shapes, and
each leak is killed by the deficit in the remaining cells / total `measS7`.

Incoming HYP-2780 is the stronger joint-coupling statement: sorted `W` partial
sums are the right object and the obstruction is localized at the single largest
cell.  HYP-2781 is an independent codex replication and leak ledger for that
coupling, with explicit witnesses and compensation obligations.

This is a proof route for the bounded piece of HYP-2773/HYP-2778, not a claim
that universal consec-max is true.  HYP-2778 shows universal consec-max is
false by k=12 and unnecessary for LRC14; HYP-2779 says the wide piece must be a
direct joint `p0` bound.  HYP-2781 therefore targets the finite bounded
compression certificate and should be read as supporting HYP-2780, not
competing with it.

## Exact Scout

Script:

```text
04-computation/lrc14_cell_width_majorization_codex_s75.py
```

Stored output:

```text
05-knowledge/results/lrc14_cell_width_majorization_codex_s75.out
```

The script enumerates anchored primitive full-residue shapes and computes exact
Fractions for the six cell widths `W_a`.  It also compares fixed cyclic
cell-window sums against AP and compares the sorted top-`L` prefixes, which are
dilation-safe because multiplication by `F_7^*` merely permutes the cell labels.

Bounded banks:

```text
k=8  span<=15  full-residue shapes=528
k=9  span<=14  full-residue shapes=688
k=10 span<=14  full-residue shapes=832
k=11 span<=14  full-residue shapes=620
k=12 span<=14  full-residue shapes=292
```

Sorted prefix results:

```text
k=8:  top1..top5 have 0 violations.

k=9:  exactly one violating shape, E=(0,1,2,3,4,5,6,7,9).
      It leaks top1..top4 by
      5/392, 11/1176, 37/4410, 11/2205,
      but the total deficit is -17/1176.
      top5 has 0 violations.

k=10: exactly one violating shape, E=(0,1,2,3,4,5,6,7,9,10).
      It leaks only top1 by 53/8820.
      The top2 prefix is already below AP by 5/3528,
      and the total deficit is -151/5880.

k=11: top1..top5 have 0 violations.

k=12: guardrail, not part of the bounded-consec target after HYP-2778.
      The sorted quotient has violations, including one shape with positive
      total surplus under this cell-width convention.
```

The k=12 line is useful negative signal.  It says sorted width majorization is
not a universal theorem, exactly matching HYP-2778's correction that the
project only needs a bounded cap certificate at the binding rows.

## No-Go: Fixed Cell Windows

Fixed cyclic cell-window sums fail badly and are not dilation invariant.  The
same exact run reports:

```text
k=8:  fixed L=1..4 windows have violations.
k=9:  fixed L=1..4 windows have violations.
k=10: fixed L=1..5 windows have violations.
k=11: fixed L=1..4 windows have violations.
k=12: fixed L=1..5 windows have violations.
```

This refines the mac-mini Route-3 obstruction.  Per-cell reflection symmetry and
Huffer-Shepp-style symmetrisation are real, but AP does not maximize individual
cells or fixed cell windows.  The bounded certificate must be aggregate and
dilation-quotiented.

## Proof Obligation

The proposed bounded proof has five parts:

```text
P0. Use HYP-2749 / HYP-2770 to reduce to the full-residue bounded stratum.

P1. Replace labeled cell data by the dilation quotient
    sort_desc(W_1,...,W_6).  This preserves total measS7 and destroys only the
    cell labels permuted by F_7^*.

P2. Prove sorted prefix dominance for all bounded target rows except the
    explicit one-hole AP-extension leak family.

P3. Prove a one-sink compensation lemma for the leak family:
    any excess in the largest cell is repaid by the next cell(s), and certainly
    by the total cell sum.

P4. Splice the bounded certificate into the THM-534/Delsarte cap layer, while
    leaving the HYP-2779 wide side to a direct joint p0 bound rather than a
    separable Q(k-1)+error estimate.
```

The likely exact formulation is not ordinary Schur-convexity of `measS7`; that
already fails in nearby aggregate tests.  It is a restricted finite certificate:
sorted width prefixes plus an addressed exception ledger.

## Tournament Analysis

Vertices are proof lenses.  Pair observable:

```text
(affine_safe, pins_ap, finite_leaks, formal_ready, composes_with_wide)
```

The switch is lexicographic comparison, producing a transitive tournament:

```text
relation_marginal_bound      score 7
one_sink_compensation        score 6
sorted_cell_majorization     score 5
conductance_bottleneck       score 4
residue_affine_moments       score 3
win_disconnected_split       score 2
fixed_cyclic_windows         score 1
per_cell_reflection          score 0

score histogram = {0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed 3-cycles = 0
Hamiltonian path count = 1
```

The ranking is intentionally conservative.  The sorted-cell quotient is a good
bounded compression, but relation marginals and the explicit one-sink
compensation ledger are closer to a formal proof.

## Challenged Assumption

Rejected assumption: tournament vertices / proof atoms can be fixed LRC cells
or fixed cyclic windows.

Preserved predicate: total `measS7=sum_a W_a` and the dilation-invariant
multiset of cell widths.

Destroyed information: which residue cell carries which width.  This is exactly
the quotient that `W_a(cE)=W_{ca}(E)` permits.

What remains to prove: the finite leak family is the only way to beat AP in a
top sorted prefix at k=8..10/11, and the leak is always repaid before the final
sum reaches the cap.
