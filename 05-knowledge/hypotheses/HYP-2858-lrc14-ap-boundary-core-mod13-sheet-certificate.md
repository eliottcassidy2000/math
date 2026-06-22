---
id: HYP-2858
title: LRC14 AP boundary core has finite mod-13 and sheet-count certificates -- {t,2t,...,12t,V} does not need V/t -> infinity
status: PROVED for the AP boundary-core family; mod-13 certificate verified on 387840 rows; sheet constants exact
source: codex-2026-06-22-S90/S91
related:
  - HYP-2581
  - HYP-2581d
  - HYP-2581f
  - HYP-2652
  - HYP-2838
  - HYP-2853
  - HYP-2855
  - HYP-2856
  - HYP-2857
  - HYP-+2863
  - THM-523
  - THM-524
  - THM-526
  - THM-527
  - OPEN-Q-108
---

# HYP-2858 -- AP boundary-core mod-13 and sheet-count certificates

## Claim

The coordinated-growth AP boundary core

```text
S(t,V) = {t, 2t, ..., 12t, V}
```

is not a genuine finite-`V` Part-A obstruction.  For every positive `t,V`,
`S(t,V)` has a direct LRC14 witness `tau` with

```text
min_{s in S(t,V)} ||s tau|| >= 1/14.
```

Thus the suggested Node-1 condition `V/t -> infinity` is not needed for this
AP core.  The remaining Node-1 difficulty is the non-AP wide/coordinated-growth
case, not the pure AP family.

## Proof

First divide by `g = gcd(t,V)`.  Scale invariance lets a witness `tau_0` for
the reduced set `{t_0,2t_0,...,12t_0,V_0}` give the original witness
`tau=tau_0/g`.

Assume then `gcd(t,V)=1`.

### Case 1: denominator `13t` certificate

If there is an integer `a` such that

```text
13 ∤ a
and
||a V / (13t)|| >= 1/14,
```

then `tau=a/(13t)` works.  Indeed, for `j=1,...,12`,

```text
j t tau = j a / 13,
```

and since `13 ∤ a`, the residues `ja mod 13` are nonzero, so
`||j a/13|| >= 1/13 > 1/14`.

It remains to prove that such an `a` exists except for the tail family below.

If `13 ∤ V`, multiplication by `V` permutes `Z/(13t)`.  Let

```text
M = {r mod 13t : 14 * dist_mod(r,0) >= 13t}.
```

Then `|M| = 13t - 2*ceil(13t/14) + 1 > t`.  The forbidden preimages
`13 | a` form only `t` residue classes, so some `r in M` has preimage
`a` with `13 ∤ a`.

If `13 | V` and `t>1`, write `V=13V'`.  Since `gcd(t,V)=1`, `gcd(t,V')=1`
and `13∤t`.  Choose a middle residue modulo `t`, then pull it back by
`V'` and choose one of its `13` representatives modulo `13t` not divisible
by `13`.  This again gives the denominator-`13t` certificate.

### Case 2: reduced tail `{1,...,12,13m}`

The only reduced case not covered by the denominator-`13t` certificate is
`t=1`, `V=13m`.  Here the explicit witness

```text
tau = m / (13m+1)
```

gives the exact tail-law margin.  For `1 <= j <= 12`,

```text
||j tau|| = min(jm, 13m+1-jm)/(13m+1) >= m/(13m+1),
```

and

```text
13m * tau = 13m^2/(13m+1) == -m/(13m+1)  (mod 1),
```

so the `V` runner has the same distance `m/(13m+1)`.  This is at least
`1/14`, with equality only at `m=1`, the tight AP `{1,...,13}`.

Therefore every AP boundary core is certified directly.

## Sheet-Count Proof For The Covering-Forced Pure Dilation

The user's follow-up isolates the purest kps-S4 hard core:

```text
S(b,V) = {b,2b,...,12b,V},   gcd(b,V)=1,   V == 0 mod 14.
```

This is already covered by the mod-13 proof, but it gives a more transferable
finite-`V` quotient.  Write

```text
tau = (n+u)/b,   n=0,...,b-1,   u in [0,1).
```

The dilated block is lonely exactly when

```text
u in G_12 = {u : ||j u|| >= 1/14 for j=1,...,12}.
```

Exact interval decomposition gives

```text
meas(G_12) = 6617/194040,
arcCount(G_12) = 12,
widest arc = [1/14, 13/168].
```

### Comb Floor

The THM-523/THM-518 comb argument specializes to

```text
L(S) >= (6/7) meas(G_12) - b*arcCount(G_12)/(7V).
```

So this lower bound is positive whenever

```text
V/b > arcCount(G_12)/(6*meas(G_12))
    = 388080/6617
    = 58.648935...
```

This is a good large-`V/b` cash-out, but it deliberately loses the regime where
the hard core was supposed to be most dangerous.

### Sheet Quotient

For fixed `u in G_12`, the offsets

```text
frac(V n / b),   n=0,...,b-1,
```

are exactly the `b` equally spaced points on the circle because `gcd(V,b)=1`.
The extra runner `V` is unsafe on an interval of length `1/7`, so it can kill
at most `b/7 + 1` of the `b` sheets.  Hence at least

```text
6b/7 - 1
```

sheets survive for every fixed good `u`.  Integrating over `G_12`,

```text
L(S) >= meas(G_12) * (6/7 - 1/b).
```

In the covering-forced primitive subcase `b` is odd and nontrivial, so this is
positive in the hard-core rows.  More importantly, it proves that the apparent
finite-`V` obstruction is not a `V/b` problem: the sheet coordinate itself
keeps a positive fraction of witnesses.

### Transferable Lemma Shape

The reusable statement suggested by this proof is:

```text
For bE plus h parked runners each coprime to b,
L >= meas(G_E) * (1 - h/7 - h/b)
```

before any slow-fast `V -> infinity` limiting step, with small `b` left to a
finite check.  This points at non-AP coordinated-growth clusters: the useful
vertices are sheet labels and parked-runner damage intervals, not the raw
growth ratio.

## Computation

Script:

```text
04-computation/lrc14_ap_boundary_core_certificate_codex_s90.py
```

Transcript:

```text
05-knowledge/results/lrc14_ap_boundary_core_certificate_codex_s90.out
```

Sheet-count script:

```text
04-computation/lrc14_pure_dilation_sheet_count_codex_s91.py
```

Sheet-count transcript:

```text
05-knowledge/results/lrc14_pure_dilation_sheet_count_codex_s91.out
```

The script verifies the constructive certificate over the box
`t <= 80`, `V <= 120t`, skipping duplicate rows `V in {t,...,12t}`:

```text
rows = 387840, failures = 0
```

Case minimum margins:

| case | count | minimum |
|---|---:|---:|
| denominator `13t` | `387120` | `1/14` |
| tail `{1,...,12,13m}` | `720` | `1/14` |

Representative certificates include:

| row | witness | margin |
|---|---:|---:|
| `{1,...,12,13}` | `1/14` | `1/14` |
| `{1,...,12,26}` | `2/27` | `2/27` |
| `{1,...,12,182}` | `14/183` | `14/183` |
| `{20,40,...,240,261}` | `19/260` | `19/260` |
| `{20,40,...,240,2007}` | `1/260` | `1/13` |

## Consequences

This repairs one specific reading of the old OPEN-Q-108/HYP-2581d AP-family
warning.  The family `{t,2t,...,12t,V}` is infinite and can make some
slow-fast margins look asymptotically tight, but it does not require
equidistribution, three-distance machinery, or `V/t -> infinity` to prove
LRC14.  A finite modular witness already exists.

The result does not close the full Node-1 finite-`V` problem: non-AP
coordinated-growth clusters can still require the `GOOD ∩ G_P` floor plus the
arc-count/rhoK approximation.  It narrows the target by removing the canonical
AP hard core.

Post-rebase integration: incoming HYP-2855 identifies the q-uniform
three-distance floor, HYP-2856 gives an explicit `3/pi^2` Farey lower bound,
and HYP-2857 gives the Abel/Dirichlet signed-tail route for the
quasi-independence floor.  Incoming HYP-+2863 independently closes the same
boundary-core family through the `rho_K` discretization window
`V/t > 12`.  This certificate is complementary: it says the canonical AP
boundary family can be handled by exact modular residues, sheet damage counts,
or the discretized `rho_K` route.  That leaves the three-distance/Farey/signed
Fourier machinery for the genuinely non-AP coordinated-growth residue.

## Tournament Analysis And Assumption Challenge

Pairwise observable: exact certified margin or positive lonely measure for the
AP boundary core.

Switch/gauge: a proof quotient wins if it proves the AP boundary core
uniformly for all finite ratios.

Hamiltonian path:

```text
sheet_quotient_damage_count > mod_13t_residue_certificate > tail_13m_exact_witness > scale_invariance_gcd_reduction > comb_floor > slow_fast_arc_error_budget > raw_V_over_t_growth > runner_vertices
```

Challenged assumption: the useful vertex is not the raw runner family nor the
growth ratio `V/t`.  The quotient preserving the LRC predicate is the
denominator-`13t` residue class with the `13∤a` AP-safety constraint and the
middle-residue condition for `V`, or the sheet label `n` after `tau=(n+u)/b`.
The sheet quotient preserves actual finite witnesses and destroys only the
irrelevant raw growth scale.
