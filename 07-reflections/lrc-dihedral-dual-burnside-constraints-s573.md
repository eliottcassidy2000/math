---
source: codex-2026-06-03-S573
status: reflection + finite dihedral Burnside audit
tags: [LRC, Burnside, Fourier, dihedral, reflection, cosine, n14]
---

# Dihedral Dual Burnside Constraints

The previous Burnside pass corrected the cyclic action: the full `Z/N` action
does not act on the lonely subset, only on ambient time.  The legal cyclic
object is the stabilizer of the whole lonely word.

The sharper thing is that cyclic symmetry is still not the whole story.

For every speed `v`,

```text
||v t|| = ||v(-t)||.
```

So every binary lonely word satisfies

```text
1_X(t) = 1_X(-t).
```

That gives a free reflection.  The legal finite group is at least dihedral:
cyclic stabilizer shifts plus the reflection coset.

## Constraint

If `K_cyc` has size `k`, then the legal dihedral stabilizer has size `2k`.
Burnside is now

```text
|X / K_D| = (|X| + reflection-fixed lonely slots) / (2k),
```

with one reflection term per cyclic stabilizer shift.

For the usual audited `n=14` rows, `k=1` and the reflection has no lonely fixed
slot, so the quotient simply halves the lonely times into `t,-t` pairs:

```text
AP/V*: 6 -> 3
S562 packet: 276 -> 138
S562 dyadic lift: 264 -> 132
random low resonance: 4966 -> 2483
```

The all-odd probe shows the fixed-point correction is real, not cosmetic:

```text
5845 lonely slots with one lonely half-turn fixed point -> 2923 orbits.
```

## Dual Meaning

The dual representation statement is cleaner than the orbit formula:

```text
the odd/sine sector is forbidden.
```

The FFT still sees the clock families from HYP-2085, but now every allowed
frequency occurs as a cosine pair.  Imaginary Fourier mass and odd-sector L2
mass are zero up to numerical roundoff.

This sharpens HYP-2086.  If the hard LRC regime is the Burnside fixed boundary,
then the time-word obstruction must live inside the dihedral fixed/cosine
sector.  A putative counterexample cannot spend any obstruction mass in the
anti-fixed sector, because the binary LRC predicate annihilates that sector
before any proof begins.

## Handoff

The next useful object is not the binary word.  It is the labelled event word:

```text
runner owner, pair-sum denominator, G component, wall endpoint, endpoint core.
```

The test should be:

```text
which labels remain reflection-fixed, and which symmetries appear only after
forgetting ownership?
```

If a boundary core is self-converse only after labels are erased, that is proof
debt.  If it is label-level dihedral fixed, then the obstruction is confined to
an even smaller fixed-point family.

**Artifacts:** `04-computation/lrc_dual_burnside_constraints_s573.py`,
`05-knowledge/results/lrc_dual_burnside_constraints_s573.out`, HYP-2087.
