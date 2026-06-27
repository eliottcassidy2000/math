---
id: HYP-2940
title: LRC14 C=27 two-swap frontier splice
status: PROOF-INTERFACE / exact local frontier extension; not a proof
source: codex-2026-06-23-S138
related:
  - HYP-2937
  - HYP-2936
  - HYP-2935
  - HYP-2934
  - HYP-2932
  - HYP-2930
  - HYP-2920
  - HYP-2908
  - THM-544
  - OPEN-Q-108
---

# HYP-2940: C=27 two-swap frontier splice

S138 extends HYP-2937's `C=27` shell-transfer quotient from the
single-replacement AP/Goddyn-Wong atlas to the exact AP two-swap bank:

```text
start with AP = {1,...,13}
delete two AP entries
add two values <= 40
keep primitive rows
```

The purpose is not to prove LRC14 outright.  It tests whether interacting AP
holes create new low-gap atoms before the existing AP/GW/petal/K33 proof
routes can engage.

## Computation

The script
`04-computation/lrc14_c27_two_swap_frontier_codex_s138.py` stores output at
`05-knowledge/results/lrc14_c27_two_swap_frontier_codex_s138.out`.

It audits:

```text
primitive AP two-swap rows with added values <= 40: 27730
q_threshold >= 14 survivors:                         8674
exact M(S) computed for all survivors
```

The q-threshold filter is rigorous for the second-gap frontier.  If
`q_threshold(S)=q<=13`, then `t=1/q` gives `M(S)>=1/q>=1/13>2/27`, so no row
with `M<=2/27` is lost by filtering to `q>=14`.

The exact low frontiers are:

```text
M <= 3/41:
  AP                            1/14
  GW 12->24                     1/14
  near-miss 12->36              3/41

M <= 2/27:
  AP                            1/14
  GW 12->24                     1/14
  near-miss 12->36              3/41
  swap 10->20                   2/27
  swap 13->26                   2/27
  drop(10,12)->add(20,24)       2/27
  drop(10,12)->add(20,36)       2/27
```

Thus the first K33 wall is stable under the two-swap test: no new row appears
at or below `3/41`.  The `2/27` layer is not stable, but its two new rows are
not wild.  They are exactly the splice of the unit petal `10->20` with the two
known `12`-branch packets:

```text
10->20 plus GW 12->24
10->20 plus near-miss 12->36
```

## C=27 Transfer Readout

The two genuine two-hole rows have marked shell transfers:

```text
drop(10,12)->add(20,24):
  H[10:g1, 12:g3] -> D[7:g1, 3:g3]

drop(10,12)->add(20,36):
  H[10:g1, 12:g3] -> D[7:g1, 9:g9]
```

So the new layer is a product-like splice:

```text
unit hole 10 -> unit double 7
times
12-branch nonunit transfer (g3->g3 for GW, g3->g9 for near-miss)
```

This is the new theorem target.  The second-gap two-hole rows should be killed
by showing that the unit-hole petal component forces a `2/27` witness unless
the `12` component collapses back to the tight GW transfer.  If the `12`
component instead moves to the gcd-9 shell, the row is still loose at `2/27`
and should feed the K33/HYP-2908 state-lift branch.

## Guardrail

The C=27 transfer quotient remains only a carrier.  In the `q>=14` two-swap
bank, several low transfer labels recur in safely looser rows, including:

```text
perfect transversal
H[12:g3] -> D[3:g3]
H[12:g3] -> D[9:g9]
unit-hole petal labels
```

Therefore exact `M(S)` / Farey branch must stay attached.  The useful ordering
is:

```text
exact M/Farey branch
> q-threshold prefilter
> marked C27 shell transfer
> two-hole interaction flag
> petal or K33 discharge
> raw runner-residue tournament
```

## Tournament Analysis

S138 uses proof channels as vertices, not runners.  The pair observable is a
role vector:

```text
theorem-scale retention,
shell predicate retention,
two-hole interaction visibility,
finite-certifiability,
state-lift fit,
anti-scalar guard.
```

The resulting tournament is transitive, with fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
c3=0
scc=(1,1,1,1,1,1,1,1)
hp=1
```

Hamiltonian path:

```text
exact M/Farey branch
> q>=14 threshold prefilter
> marked C27 shell transfer
> two-hole interaction flag
> unit-hole petal discharge
> gcd3-to-gcd9 K33 packet
> AP-tail theorem S124/S35
> raw runner-residue tournament
```

## Proof Target

Bounded evidence suggests the following sharpened local lemma:

```text
In the AP two-swap neighborhood after the q>=14 reduction, every row with
M<=2/27 is either inherited from the S136 single-replacement frontier or is
the splice {10,12}->{20,24/36}.
```

The global proof version would need to replace the ceiling `40` by a structural
argument:

```text
Any genuine two-hole low-gap residual must contain a unit-visible C27 hole.
If the unit hole is not the 10->20 petal, q<14 or a larger Farey gap appears.
If it is the 10->20 petal, the second hole is forced into the 12-branch,
where GW is tight and the gcd9 near-miss is loose/state-liftable.
```

This gives a concrete next move toward LRC14: prove a two-hole splice lemma,
then attach it to the existing AP-tail theorem and the HYP-2908 K33 endpoint.
