---
id: THM-769
title: Binding-scale residue packets and the two-tightener sheet-capacity reduction
status: CLAIMED (proof being written; the exact local residue packet and sheet-capacity inequality are elementary; no claim of sporadic-branch emptiness)
source: codex-2026-07-14-S4 (n=12 residue-forcing audit)
depends_on:
  - THM-593   # unit-residue pinning, conditional on no multiple of 13
  - THM-668   # pair-sum maximizers
  - LRCUpTo13 # lower-dimensional margin for the on-sheet core
related:
  - THM-617   # shift-pigeonhole precursor at the fourteen-runner threshold
  - THM-765   # hereditary primitivity / tooth decks
  - HYP-6820  # n=12 sporadic-branch audit
---

# THM-769 — Binding-scale residue packets and sheet capacity

This number is reserved for the following elementary reduction, whose full
proof and strict/closed endpoint bookkeeping are in progress.

Let `A` be a primitive twelve-speed set with `M(A)=1/13`, and let a global
maximizer be `p/(13s)` in lowest terms.  The phase residues lie in the closed
band `[s,12s]`, both endpoints occur, and the endpoint owners are divisible
by `s`.  Writing `E={v in A:s|v}=sU` and `F=A\E`, the off-sheet runners are
strictly interior at this maximizer and `p/13` is a local `1/13` maximum of
`U`.

At one global maximizer of `U`, its `s` lifts must be covered by the closed
`1/13`-danger teeth of `F`.  For `w in F`, put
`g_w=gcd(w,s)` and `D_w=s/g_w`.  The resulting necessary capacity inequality
is

```text
sum_(w in F) (floor(2D_w/13)+1)/D_w >= 1.
```

Every summand is at most `1/2`, with equality only at `D_w=2`.  Consequently
`|F|>=2`; if `|F|=2`, both tighteners have `D_w=2`, and primitivity forces
`s=2`.  The remaining two-tightener packet is therefore exactly ten even
speeds and two odd speeds, with complementary ownership of the two sheets.

Still missing from this file: the complete proof, the exact folded two-color
form of the `s=2` residual, and explicit guardrail examples showing why a
single binding maximum does not force a full nonzero residue system modulo
13.
