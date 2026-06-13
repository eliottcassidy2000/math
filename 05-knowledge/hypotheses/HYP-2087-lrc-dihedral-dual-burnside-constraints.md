---
id: HYP-2087
status: SUPPORTED - finite dihedral constraint audit; proof use open
source: codex-2026-06-03-S573
related:
  - HYP-2086
  - HYP-2085
  - HYP-2082
  - HYP-2081
  - HYP-2080
---

# HYP-2087: dual Burnside constraints sharpen from cyclic stabilizers to the dihedral cosine sector

HYP-2085 corrected the naive cyclic Burnside count for LRC time slots: the
lonely subset is not usually invariant under all translations, so the legal
cyclic group is the stabilizer of the whole lonely word.

The sharper constraint is that every LRC lonely word has an additional
symmetry:

```text
1_X(t) = 1_X(-t).
```

This is forced by distance to the nearest integer:

```text
||v t|| = ||-v t||.
```

So the legal finite Burnside stabilizer is not merely cyclic.  It is dihedral:
the cyclic stabilizer plus its reflection coset.  If the cyclic stabilizer has
size `k`, the dihedral stabilizer has size `2k`.

## Evidence

`lrc_dual_burnside_constraints_s573.py` repeats the HYP-2085 `n=14` grid audit
on

```text
N = lcm(14,12,15,16,23) = 38640.
```

All audited primitive binary lonely words have trivial cyclic stabilizer, but
nontrivial dihedral stabilizer:

```text
AP_wall:                 6 lonely slots -> 3 dihedral orbits
V_star_wall:             6 lonely slots -> 3 dihedral orbits
near_AP_apex:          470 lonely slots -> 235 dihedral orbits
S562_packet_n14:       276 lonely slots -> 138 dihedral orbits
S562_packet_n14_lift:  264 lonely slots -> 132 dihedral orbits
no_small_pinch_proxy: 4810 lonely slots -> 2405 dihedral orbits
random_low_resonance: 4966 lonely slots -> 2483 dihedral orbits
```

The all-odd probe has one lonely reflection-fixed half-turn, so the Burnside
fixed-point correction is visible:

```text
all_odd_probe: 5845 lonely slots, reflection fixed lonely total 1,
dihedral lonely orbits = (5845 + 1)/2 = 2923.
```

The Fourier dual also tightens.  Time reversal kills the odd/sine sector:

```text
odd-sector L2 = 0.000e+00 on every audited row
imaginary Fourier L2 = numerical roundoff only, about 1e-32
```

Thus any LRC counterexample or boundary obstruction must survive after
projection to the real cosine character sector.  Sine/odd character mass is not
just unhelpful; it is forbidden for binary lonely words.

## Interpretation

HYP-2086 localized the hard LRC regime to the Burnside fixed side: boundary,
self-converse, resonance-maximal objects rather than generic round/open
orbits.

HYP-2087 sharpens the time-word version of that containment:

```text
cyclic quotient:      legal only after stabilizer restriction
dihedral quotient:   always contains time reversal
dual quotient:       only cosine/even characters survive
```

For `n=14`, the binary word is still too coarse for a proof, but it gives a
hard constraint on any richer labelled event word.  The next audit should ask
which runner-owner, pair-sum, wall-endpoint, and endpoint-core labels remain
fixed by reflection, and which apparent symmetries exist only after forgetting
ownership.

## Tournament Analysis

Vertices were speed-set time words.

Pairwise observable:

```text
(density, cyclic fold, dihedral orbit count, odd-sector defect, top cosine mode).
```

Switch/gauge:

```text
harder = lower density, more cyclic folding, fewer dihedral quotient orbits,
zero odd-sector defect, stronger top cosine mode.
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3_cycles=0
sccs=[1,1,1,1,1,1,1,1]
hamiltonian_paths=1
```

## Assumption Challenge

Possible vertices considered:

```text
time slots, translations, reflections, Fourier modes, wall-crossing events,
runner-owner events, pair-sum events, and endpoint-core states.
```

Chosen quotient:

```text
whole binary lonely time words under their legal dihedral stabilizer.
```

Predicate preserved:

```text
which grid times are lonely, which translations stabilize the word, which
reflections stabilize the word, and which Fourier sectors are allowed.
```

Information destroyed:

```text
runner ownership, wall endpoint labels, pair-sum pinch provenance, continuous
interval geometry, and endpoint-core labels.
```

Challenged assumption:

```text
cyclic translations are the whole legal finite Burnside group.
```

They are not.  Time reversal is always legal, and the dual Burnside constraint
is the cosine/even projection.

## Files

- `04-computation/lrc_dual_burnside_constraints_s573.py`
- `05-knowledge/results/lrc_dual_burnside_constraints_s573.out`
- `07-reflections/lrc-dihedral-dual-burnside-constraints-s573.md`
