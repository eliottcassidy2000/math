# LRC14 S68: Moderate Multi-Block Exact Budget

Date: 2026-06-21
Source: codex-2026-06-21-S68

## Context

HYP-2714 has localized the remaining sector-route gap to the moderate-span
balanced / multi-block branch.  Incoming HYP-2715--HYP-2717 changed the right
form of the missing lemma:

```text
not:  Product(E) is an upper bound for p0(E),
but:  |top-character error| <= cap_k - Product(E).
```

The useful quotient is the HYP-2716 Krawtchouk top character, followed by the
HYP-2717 carrier relation filter.  Exact low-gap and low-relation exceptions
should be routed to a finite ledger.

## New Computation

Script:

```text
04-computation/lrc14_moderate_multiblock_budget_codex_s68.py
```

Stored output:

```text
05-knowledge/results/lrc14_moderate_multiblock_budget_codex_s68.out
```

The scout uses exact integer boundary-sweep arithmetic and raw rational
comparisons.  The default run extends the previous fastcheck windows one step:

```text
k=8,  span <= 21
k=9,  span <= 18
k=10, span <= 17
```

All three windows have zero over-cap rows and zero rows within `0.02` of cap.
The global exact leaders remain the consecutive blocks.  The newly added
moderate-balanced leaders are far below cap:

```text
k=8:  E=(0,2,4,6,8,10,12,21), margin=1/6
k=9:  E=(0,1,2,3,4,5,6,7,18), margin=17499/140140
k=10: E=(0,1,2,3,4,5,6,7,9,17), margin=30437/194922
```

The script can be deepened by setting:

```bash
LRC14_S68_WINDOWS="8:24,9:22,10:20"
```

The first attempted deep run was too slow with the midpoint common-refinement
engine, so the script was patched to use a boundary sweep and raw integer-pair
comparisons.  The one-step exact extension is the pushed stable default.

## Proof Interpretation

S68 does not prove LRC(14).  It gives a cleaner role for the finite ledger:
moderate-balanced rows just beyond the prior window are not close to the cap,
so low-gap or low-height carrier resonances can plausibly be handled by exact
enumeration and targeted finite reductions.

The analytic residue should be:

```text
1. split E = {0} union_i (M_i+B_i);
2. write the scalar cover error as the HYP-2716 top character M_6;
3. expand M_6 in carrier Fourier modes;
4. route low-height exact relations n.M=0 and small |n.M| nonrelations to the
   HYP-2714 finite ledger;
5. prove a BV/Fourier tail estimate for the remaining high-height relation
   modes and high-denominator nonrelations under cap_k - Product(E).
```

The exact relation modes are unavoidable for integer carriers.  Therefore the
right target is not full-torus equidistribution, but a relation-height tail
bound after the top-character quotient.

## Tournament Analysis

The vertices are proof obligations, not runners or raw arcs:

```text
bounded_span_exact_window
single_far_comb_closed
far_count_domination
moderate_balanced_multiblock
carrier_product_error_bound
raw_gap_monotonicity
```

Pairwise observable: proven/finite status first, then cap margin, then
preserved proof data.  Switch/gauge: pass from raw rows to proof obligations
after the HYP-2716/HYP-2717 quotient.  The S68 scout tournament is transitive
with score histogram `{0:1,1:1,2:1,3:1,4:1,5:1}` and zero directed 3-cycles.

The productive Hamiltonian path is:

```text
far_count_domination
> moderate_balanced_multiblock
> single_far_comb_closed
> bounded_span_exact_window
> carrier_product_error_bound
> raw_gap_monotonicity
```

This ranking is not a proof of the open lemma.  It says which obligations
currently have the most usable margin and which are just diagnostics.

## Assumption Challenge

This session considered vertices as runners, gaps, fixed circle sections,
section boundaries, wall-crossing events, residues, cover arcs, residual masks,
Fourier modes, carrier relation classes, and proof obligations.  The useful
finite scout quotient is proof obligations; the useful analytic quotient is
carrier relation classes after the Krawtchouk top-character projection.

What is preserved:

```text
cap margin,
the scalar p0 cover target,
low-gap / low-relation routing data,
and the top-character error object.
```

What is destroyed:

```text
sector ownership,
individual residual-coordinate signs,
raw gap monotonicity,
and literal full-torus carrier independence.
```

Challenged assumption: the multi-block branch can be closed by proving
`Product(E) >= p0(E)` or by a scalar gap-monotonicity theorem.  HYP-2715
already refutes the product-envelope form, and S68 shows raw gap monotonicity
is weaker than the finite-ledger plus top-character-tail split.

