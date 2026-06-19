---
id: HYP-2613
title: LRC(14) relative support-6 permanent count
status: OPEN
source: codex-2026-06-19-S11
depends_on:
  - THM-538
  - HYP-2608
  - HYP-2612
  - HYP-2610
  - HYP-2611
related:
  - THM-410
  - MISTAKE-078
  - OPEN-Q-108
---

# HYP-2613 - LRC(14) Relative Support-6 Permanent Count

## Claim

The remaining LRC(14) wide-spread lemma should not be attacked by a bare
successive-minima bound on the full relation lattice.  The useful object is a
relative, signed count:

1. Split off the bounded core relation lattice.
2. For supports that touch a large offset, count lattice points on each fixed
   support hyperplane `sum n_j e_j = 0`.
3. Keep the exact THM-538 signed kernel `K(n)`, whose exact support-6 layer is a
   six-sector root-of-unity permanent.
4. Enumerate finite low-height subset-sum resonance walls separately.

If this hypothesis is right, HYP-2608(a) should close by a proof of the form

```text
bounded core finite certificate
+ finite low-height subset-sum walls
+ relative signed hyperplane tail
< cap_k - wide ceiling
```

for `k=8,9,10`.

## Evidence From S11

This refines the concurrent HYP-2612 anti-coset count by measuring the
support-by-support signed layer and the exact residue permanent constants.

Script:

- `04-computation/lrc14_support6_relative_count_codex_s11.py`
- output: `05-knowledge/results/lrc14_support6_relative_count_codex_s11.out`

The direct support-6 envelope computation confirms that the free product
majorant was the wrong object, but also shows that a merely coupled absolute
Minkowski count is still too blunt:

- AP `k=8`: exact-support-6 envelope at height `H=28` is `32305.4`, while the
  actual correction is `0.302731`.
- Resonant one-stranger `{0,1,2,3,4,5,6,21}`: type-II support-6 envelope at
  `H=28` is `17190.3`, while the wide shape has `meas(S7)=0.217687`.
- High one-stranger `{0,1,2,3,4,5,6,211}`: type-II support-6 envelope is still
  `413.173` at `H=28`, despite exact `meas(S7)=0.195351`.

The exact signed support-6 layer is much smaller:

- AP `k=8`: signed exact-support-6 layer through `H=12` is `0.0471542`.
- `{0..6,21}`: type-II signed support-6 layer through `H=12` is `-0.00596474`.
- `{0..6,211}`: type-II signed support-6 layer through `H=12` is about
  `-7.63e-8`.
- `k=9` wide sample `{0..7,68}`: type-II signed support-6 layer through `H=10`
  is `0.000987`.
- `k=10` wide sample `{0..8,22}`: type-II signed support-6 layer through `H=8`
  is `-0.002595`.

The exact support-6 residue permanent improves the per-relation constant but is
not enough by itself.  Compared with the blunt `64*c1^6 = 7.35714264`, the
maximum normalized constants are:

```text
ambient d=6: 0.643084862  ratio 0.0874
ambient d=7: 0.321542431  ratio 0.0437
ambient d=8: 0.168426988  ratio 0.0229
ambient d=9: 0.091869266  ratio 0.0125
```

This explains why the next lemma must keep signed summation over relation
hyperplanes, not only replace `c1` by a sharper absolute permanent constant.

## Reframing

The phrase "execute the Minkowski count" is now too coarse.  There are three
distinct layers:

- **Counting layer:** the support relation equation couples coefficients and
  removes the free harmonic divergence.
- **Permanent layer:** exact support 6 is a sector permanent depending only on
  residues mod 7 after normalization.
- **Oscillatory layer:** the signed sum over each support hyperplane is far
  smaller than the absolute envelope, especially for large-coordinate supports.

The likely proof target is a relative signed theta estimate for each support
hyperplane, uniform after finitely many low-height subset-sum walls are removed.

## THM-410 Analogy

THM-410 says that reversing a long interval edge creates directed triangles
according to the number of vertices inside the interval.  The analogous object
here is:

```text
large offset e_L is dangerous only when bounded core coefficients can sum to e_L
```

Low-height walls such as `1+2+3+4+5+7=22` are support-6 subset-sum intervals.
They are finite resonance walls, not an infinite harmonic tail.  The long
offset is not dangerous because it is long; it is dangerous only by how many
bounded witnesses lie inside its subset-sum interval.

## Tournament Analysis Note

This session challenged the default vertex choice.  The useful tournament
vertices were not runners or arcs, but proof obligations:

- six-support hyperplanes for the count itself;
- offsets for a participation fingerprint;
- subset-sum walls for the finite resonance split.

The observable used in the script is support-6 envelope participation.  The
switch orients an edge toward the larger participant, and ties follow the
Hamiltonian path induced by increasing offset.  The resulting fingerprints were
mostly transitive, with informative edge flips only when a far stranger became
less participatory than the dense core.  This preserves the LRC predicate
"support-6 tail below the cap margin" but destroys phase-location data, so it
is a counting quotient rather than a witness-time quotient.

## Status

Open.  S11 did not prove LRC(14), and one farther-out exploratory signed-tail
probe timed out before producing data.  The saved computation does give a better
next theorem statement: prove a relative signed support-6 permanent/hyperplane
tail bound, not a bare full-lattice absolute Minkowski bound.
