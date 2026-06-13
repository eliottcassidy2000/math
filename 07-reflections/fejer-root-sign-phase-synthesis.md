# Fejer Kernels As Root-Sign Character Shadows

**Session:** codex-2026-05-30-S359
**Status:** exploratory synthesis plus finite cyclic probe
**Prompt:** integrate Fejer kernels and root-sign, then explore new tangents.

## Executive Thesis

The root-sign lens and the Fejer-kernel lens are not two analogies.  In the
circulant prime case they are literally two coordinates of the same object.

Start with a tournament on `Z/pZ`.  Translation invariance collapses the type-A
root sign cube to one binary choice for each cyclic root pair:

```text
{d,-d},  d = 1,...,m,  m=(p-1)/2.
```

Let

```text
sigma_d = +1  if the tournament chooses step d,
sigma_d = -1  if it chooses step -d.
```

Then every nonzero Fourier character `k` reads this root-sign vector by a sine
projection:

```text
lambda_k = sum_{d=1}^m exp(2*pi*i*k*sigma_d*d/p)
         = -1/2 + i * sum_{d=1}^m sigma_d sin(2*pi*k*d/p),
```

so

```text
|lambda_k|^2 = 1/4 + (sum_d sigma_d sin(2*pi*k*d/p))^2.
```

The interval tournament is the all-one chamber sign vector:

```text
sigma = (1,1,...,1).
```

For that vector, the same formula becomes the Fejer kernel:

```text
|lambda_k|^2 = sin^2(pi*m*k/p) / sin^2(pi*k/p).
```

So:

```text
Fejer kernel = Fourier character shadow of the all-one cyclic root-sign chamber.
```

This is the cleanest bridge between the representation-lens note and the
circulant maximizer/phase-profile program.

## What The Probe Did

I added:

```text
04-computation/fejer_root_sign_phase_probe_s359.py
05-knowledge/results/fejer_root_sign_phase_probe_s359.out
```

The script enumerates all circulant tournaments for:

```text
p = 7, 11, 13, 17, 19, 23, 29, 31
```

and computes:

```text
root_sign_formula_error
interval_fejer_error
top spectral fraction
low-frequency pair fraction
inverse participation ratio
Fejer alignment
root-sign change count
chamber bias
additive energy
interval-unit-orbit membership
Paley membership where p == 3 mod 4
Hamiltonian path count H for p <= 13
```

The two identity checks are numerical-zero:

```text
p= 7 max_root_formula_err=1.60e-14 interval_fejer_err=1.05e-15
p=11 max_root_formula_err=3.91e-14 interval_fejer_err=2.89e-15
p=13 max_root_formula_err=8.17e-14 interval_fejer_err=6.75e-14
p=17 max_root_formula_err=1.42e-13 interval_fejer_err=4.97e-14
p=19 max_root_formula_err=9.38e-13 interval_fejer_err=2.84e-14
p=23 max_root_formula_err=3.13e-13 interval_fejer_err=7.82e-14
p=29 max_root_formula_err=3.91e-12 interval_fejer_err=1.85e-13
p=31 max_root_formula_err=8.95e-13 interval_fejer_err=4.26e-13
```

The probe also recovered the old Fejer extremal pattern in a more
representation-native language: top spectral fraction maximizers are exactly
in the unit orbit of the interval chamber in the tested primes:

```text
p=7:  top_fraction maximizers=6,  in interval orbit=6
p=11: top_fraction maximizers=10, in interval orbit=10
p=13: top_fraction maximizers=12, in interval orbit=12
p=17: top_fraction maximizers=16, in interval orbit=16
```

That suggests the correct invariant is not "the literal interval sign vector"
but "the Weyl/chamber object modulo the cyclic unit action."

## The Three Layers That Must Not Be Collapsed

This session separated three things that had been too casually identified.

### 1. Chamber Coordinate

The root-sign vector in the ordered basis

```text
(1,2,...,m)
```

has a combinatorial roughness:

```text
sign_changes(sigma).
```

This is a chamber-coordinate statistic.  It depends on the chosen ordered
positive-root half.

It is useful, but it is not the Fourier phase profile.  The probe found
`corr(H,sign_changes)=0` for `p=7,11,13` in the full circulant samples.  That
does not make sign roughness useless.  It means roughness is a coordinate
feature, not the character channel itself.

### 2. Character Projection

The Fourier character sees

```text
S_k(sigma) = sum_d sigma_d sin(2*pi*k*d/p).
```

This is the phase channel.  The adjacency spectrum forgets the sign of
`S_k` and keeps:

```text
1/4 + S_k^2.
```

The collapse from signed sine projection to squared magnitude is already a
loss of phase.  It is still much richer than score or deletion residue, but it
is not the whole root-sign object.

### 3. Orbit-Invariant Concentration

Fejer alignment with the literal interval vector is coordinate-dependent.
Multiplying the connection set by a unit moves the interval to another chamber
with the same spectral concentration and same `H`, but its fixed Fejer
alignment can drop sharply.

The output makes this visible at `p=13`:

```text
H=3711175 interval_orbit top=0.4097 fejer=1.0000 changes=0
H=3711175 interval_orbit top=0.4097 fejer=0.1362 changes=2
H=3711175 interval_orbit top=0.4097 fejer=0.0930 changes=4
```

All three are interval-orbit copies.  Literal Fejer alignment sees only one
coordinate representative.  Top spectral concentration sees the orbit.

This is the practical lesson for `phase_profile`:

```text
record orbit-invariant concentration, not only fixed-coordinate Fejer alignment.
```

## Paley Versus Interval, Re-read

HYP-1801 says the Paley/interval transition is phase-dominant.  This probe
sharpens the word "phase".

Paley is a multiplicative-character root sign:

```text
sigma_d = Legendre(d).
```

for `p == 3 mod 4`.  Its spectrum is flat.  Interval is an additive/chamber
root sign.  Its spectrum is Fejer-peaked.

The old story was:

```text
Paley = flat spectrum
Interval = Fejer concentration
```

The root-sign story refines this to:

```text
Paley    = multiplicative character channel on cyclic root pairs
Interval = additive chamber channel on cyclic root pairs
```

The competition is not "random versus structured."  It is two kinds of
structured signs:

```text
multiplicative phase flatness  vs.  additive chamber concentration.
```

The small-prime `H` data remains deliciously non-monotone:

```text
p=7:  Paley H=189, Interval H=175
p=11: Paley H=95095, Interval H=93027
p=13: Interval-orbit H=3711175 wins the sample
```

The correlations say the same thing:

```text
p=7:  corr(H,top_fraction)=-1.000000, corr(H,ipr)=-1.000000
p=11: corr(H,top_fraction)=-0.106581, corr(H,ipr)=-0.087860
p=13: corr(H,top_fraction)=+0.444395, corr(H,ipr)=+0.509255
```

So Fejer concentration is not a universal monotone.  It becomes the right
direction after the Paley small-prime advantage stops dominating.

This is exactly what "phase transition" should mean here.

## Why Low-Frequency Pair Fraction Failed

I expected the low-frequency pair

```text
|lambda_1|^2 + |lambda_{p-1}|^2
```

to be a good feature.  It was not:

```text
corr(H,low_pair_fraction)=0 for p=7,11,13.
```

The reason is hindsight-obvious.  For tournaments, conjugate characters occur
in pairs.  Many unit rotations move the Fejer spike away from `k=1` while
preserving the spectrum as a multiset.  A fixed low-frequency window is not
unit-invariant.

The replacement feature should be:

```text
sorted spectrum profile
top_fraction
ipr
entropy of |lambda_k|^2
possibly the orbit of peak locations under units.
```

In other words, "low mode" is a coordinate choice; "concentration" is the
representation fact.

## How This Fits The Root-Sign Representation Note

The representation-lens note proposed:

```text
Root sign cube -> S_n quotient -> residue/phase/incidence probes -> scalar invariants.
```

The circulant Fejer story is a worked example:

```text
type-A root signs
  -> quotient by cyclic translations
  -> signs on root-pair orbits {d,-d}
  -> C_p character projections
  -> sine channel S_k(sigma)
  -> squared spectral magnitudes
  -> concentration statistics
  -> H / OCF / path-homology consequences.
```

This gives a usable model for a future representation-refined OCF.  Before the
scalar `H`, there should be a character-indexed packet count:

```text
odd-cycle packets by root support, decomposed by character channel.
```

For circulants, those channels are literal Fourier modes.  For general
tournaments, the analogue may be Krawtchouk/S_n isotypic channels:

```text
trivial
standard / score / Cartan
two-row residual
higher packet channels
```

The Fejer calculation is therefore a toy model for how a chamber sign can
become a sharply concentrated phase profile after quotienting.

## New Tangents

### Tension 1: Fejer Is An Extremal Shadow, Not A Proof Of H-Maximality

The interval chamber extremizes spectral concentration among tested
circulants, but `H` does not follow monotonically at small primes.  This
matters because previous notes occasionally sounded like:

```text
Fejer extremal -> interval maximizes H.
```

The correct version is weaker and more interesting:

```text
Fejer extremal supplies the additive chamber endpoint of the phase landscape.
H chooses between this endpoint and the multiplicative Paley endpoint.
```

At small `p`, Paley's flatness can generate more total cycle mass.  Later,
interval concentration appears to pack cycles better in the OCF sense.

### Tension 2: Unit Orbits Are The Hidden Weyl Group In The Circulant Quotient

In the full type-A story, vertex relabeling is the Weyl group `S_p`.  In the
circulant quotient, only affine transformations preserve circulant form:

```text
x -> ux + a.
```

Translation is already quotiented out.  The surviving action is the unit group
`(Z/pZ)^*`.

Thus "interval" is not one sign vector.  It is a unit orbit of chamber
half-sets:

```text
u * {1,2,...,m}.
```

This explains why literal Fejer alignment is fragile while spectral
concentration is stable.

### Tension 3: Multiplicative Characters Are Another Chamber System

Paley is often treated as the opposite of interval: flat versus concentrated.
But root-sign language suggests Paley is a chamber too, just not an additive
one.  It is a chamber cut by the multiplicative quadratic character.

So the real dichotomy is:

```text
additive order chamber        interval / Fejer
multiplicative character cut  Paley / Gauss flatness
```

This hints at a broader family: replace Legendre by other low-complexity
characters or by Bohr-set cuts, then ask which phase channel controls `H`.

### Tension 4: Root Roughness Is Not Phase

The sign-change count is tempting because it feels like smoothness.  The probe
says it is too coordinate-bound.  It should live in a separate feature block:

```text
root_chamber_features:
    sign_changes
    chamber_bias
    run lengths

phase_character_features:
    spectrum concentration
    phase entropy
    trace signs
    additive energy
```

That separation mirrors the larger residue/phase/incidence split.

## Near-Term Proof Targets

### Target A: Prove The Root-Sign Spectral Formula Cleanly

This is basically done, but should be written as a lemma:

For any prime `p`, any circulant tournament with root-pair signs `sigma_d`,
and any `k != 0`,

```text
lambda_k = -1/2 + i * sum_d sigma_d sin(2*pi*k*d/p).
```

This is the exact bridge between root signs and character channels.

### Target B: Classify Spectral-Concentration Maximizers

The probe suggests:

```text
max top_fraction among circulant tournaments on Z/pZ
is achieved exactly by the unit orbit of intervals.
```

This is stronger and cleaner than "the interval maximizes IPR" because it is
orbit-aware.  It sounds like a finite cyclic rearrangement theorem:

```text
Among sign choices sigma_d, maximize max_k |sum_d sigma_d sin(2*pi*k*d/p)|.
The maximizers are threshold/chamber signs after multiplying k by a unit.
```

This is plausibly accessible by a discrete rearrangement inequality.

### Target C: Build An Orbit-Invariant Phase Profile

The next script should emit:

```text
spectrum_sorted
top_fraction
ipr
entropy
peak_orbit
additive_energy
trace_sign_vector
root_chamber_run_profile
```

Then join this with:

```text
H
Omega alpha
deletion residue
path-homology ranks
```

for circulant `p=7,11,13`, and sampled/structured `p=17,19,23`.

### Target D: Character-Resolved OCF

The dream object is:

```text
H(T) = sum independent odd-cycle packets 2^k
```

refined by character channel.  For circulants, every odd cycle has a total
difference residue and therefore a Fourier character footprint.  That should
produce a character-resolved OCF ledger:

```text
alpha_j(T; character data)
```

If interval and Paley differ by cycle packing in phase space, this ledger
should show where.

## Summary

The clean mathematical sentence from this session is:

```text
The Fejer kernel is the squared character shadow of the all-one cyclic
root-sign chamber.
```

The clean caution is:

```text
Fejer concentration is an extremal phase coordinate, not automatically an
H-maximization theorem.
```

The clean next move is:

```text
build orbit-invariant phase profiles and then try a character-resolved OCF.
```

That feels like the right synthesis of Fejer, root signs, Paley/interval, and
the broader residue/phase/incidence program.
