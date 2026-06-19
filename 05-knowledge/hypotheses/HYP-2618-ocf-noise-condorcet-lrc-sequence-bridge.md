---
id: HYP-2618
title: OCF noise/Condorcet spectra as an LRC(14) sequence bridge
status: PARTIALLY-CONFIRMED
source: codex-2026-06-19-S15
depends_on:
  - HYP-2617
  - HYP-2615
  - HYP-2614
  - THM-538
related:
  - HYP-2616
  - HYP-2608
  - OPEN-Q-108
---

# HYP-2618 - OCF Noise/Condorcet Spectra as an LRC(14) Sequence Bridge

## Claim

The OCF evaluation

```text
H(T) = I(Omega(T), 2)
```

is canonically a hard-core partition function at activity `2`, equivalently a
biased independent-set density at `p=2/3`:

```text
mu_p(independent) = (1-p)^m I(Omega, p/(1-p)),
H(T) = 3^m mu_{2/3}(independent),   m=|V(Omega)|.
```

It is not, by itself, a nontrivial two-copy noise-stability functional at
`rho<1`.  Nontrivial biased noise stability needs the ordered-pair spectrum of
independent sets `(I,J)`, not only the OCF evaluation `H`.

Thus the useful reframe is:

```text
raw paradox mass
-> compatible packet address
-> activity-2 signed/compatible partition evaluation.
```

For tournaments read as majority relations, the forbidden values `{7,21}`
become forbidden compatible Condorcet-cyclicity spectra: not merely "some
Condorcet cycle is impossible", but "certain small inventories of compatible
odd-cycle paradox packets cannot be realized by any majority tournament."

The LRC(14) transfer is direct.  HYP-2617 supplies the finite support-six packet
address table (`159` projective mod-7 coimage classes).  The proof should now
delete/account for low-height wall classes and bound the signed reciprocal tail
over the remaining non-null packet addresses, rather than improving the raw
absolute Minkowski volume.

## Evidence

Script:

- `04-computation/ocf_noise_condorcet_lrc_sequence_bridge_codex_s15.py`
- output: `05-knowledge/results/ocf_noise_condorcet_lrc_sequence_bridge_codex_s15.out`

Exact OCF scan over all labeled tournaments `n<=6`:

| `n` | labeled tournaments | distinct `H` values | max `H` | missing odd values up to max |
|---:|---:|---:|---:|---|
| 1 | 1 | 1 | 1 | `[]` |
| 2 | 2 | 1 | 1 | `[]` |
| 3 | 8 | 2 | 3 | `[]` |
| 4 | 64 | 3 | 5 | `[]` |
| 5 | 1024 | 7 | 15 | `[7]` |
| 6 | 32768 | 19 | 45 | `[7, 21, 35, 39]` |

Formal OCF alpha candidates:

```text
H=7:  (1,1,1), (1,3,0)                 none realized through n<=6
H=21: (1,0,5), (1,2,4), ..., (1,10,0)  none realized through n<=6
```

Exact 3-voter majority profile scan:

| alternatives | profiles | unique majority tournaments | all tournaments? | `H` values |
|---:|---:|---:|---|---|
| 3 | 216 | 8 | yes | `[1,3]` |
| 4 | 13824 | 64 | yes | `[1,3,5]` |
| 5 | 1728000 | 1024 | yes | `[1,3,5,9,11,13,15]` |

So, for `m<=5`, even 3 strict voters realize every tournament.  In general,
McGarvey supplies the eventual majority-realization theorem, so the tournament
forbidden spectra are electorate-forbidden spectra.

Noise guardrail:

```text
same OCF evaluation H=23 at n=6,
left alpha  = (1,11,0),
right alpha = (1,9,1),
different ordered-pair/noise spectra.
```

Therefore a nontrivial `rho<1` noise-stability invariant cannot be recovered
from `H` alone.  The canonical normalization is activity `2` / `p=2/3` density;
diagonal same-state normalizations can force activity `2`, but they are not the
full two-copy stability.

LRC sequence splice from HYP-2617:

```text
projective support-residue classes = 159
zero-residue histogram z=0..5 = [80,42,22,10,4,1]
null class counts d=6..13 = [142,113,80,43,34,3,3,3]
```

Named support readout:

```text
AP d=7          |S_d|=0.93653539, abs/signed=15.9731
resonant 21    |S_d|=0.11706692, abs/signed=125.46
wide 68 d=8    |S_d|=0.21741,    abs/signed=38.0849
k10 wall 22    |S_d|=4.1e-16,    abs/signed=1.13e16
```

The k10 wall is the model example: raw absolute mass is huge, but the retained
coimage packet is null in the relevant dimension.  That is the LRC analogue of
not trusting raw paradox mass in OCF.

## What This Changes

This answers the OCF/noise prompt with a useful correction:

- **Yes:** OCF is a partition/density functional at the specific activity
  `x=2`, equivalently product density `p=2/3` after scaling by `3^m`.
- **Only degenerately:** it can be viewed as `rho=1` biased stability, or as a
  diagonal same-state mass after tuning `(p,rho)`.
- **No:** it is not a genuine nontrivial two-copy noise-stability functional
  determined by `H`.

For LRC(14), the practical lesson is to stop optimizing the absolute support-six
mass.  The durable sequence is the packet address:

```text
OCF: independent odd-cycle packet spectrum alpha
LRC: projective mod-7 coimage class plus low-height wall status
```

The next proof computation should classify which of the `159` coimage classes
remain after height-1 and height-2 wall deletion in the binding rows `k=8,9,10`.

## Tournament Analysis

The tournament vertices in the S15 quotient are not runners, arcs, or voters.
They are quotient choices:

```text
coimage_sequence_tail
> condorcet_alpha_spectrum
> hard_core_activity_2
> biased_density_p23
> diagonal_noise_embedding
> full_noise_stability
> raw_absolute_mass.
```

This transitive proof-path preserves the relevant predicate: which packet
address survives to the signed compatible sum.  It destroys witness-time and
individual-voter geometry, which is acceptable for this support-six tail but not
for constructing an actual lonely time.

## Status

Partially confirmed by exact finite scans and identity algebra.  The LRC(14)
proof is still open; this packet sharpens the next computation and rejects the
wrong noise-stability overclaim.
