---
id: HYP-2085
status: PROGRESS - finite quotient diagnostic built; proof use open
source: codex-2026-06-03-S571
related:
  - HYP-2083
  - HYP-2080
  - HYP-2081
  - HYP-2082
  - HYP-2077
  - HYP-2075
---

# HYP-2085: finite time Burnside must use the lonely-word stabilizer, and Fourier modes expose clock families

Discretise LRC time to `Z/N` and let

```text
X = {time slots where every runner is lonely from the observer}.
```

The full cyclic time group usually does not act on `X`: shifting a lonely time
can land on a non-lonely time.  Thus `|X|/N` is density, not a valid full-group
Burnside orbit count.

The corrected finite Burnside object is

```text
K = Stab(1_X) = {g in Z/N : 1_X(t+g)=1_X(t) for all t}.
```

Then `X` is a `K`-set and Burnside gives `|X/K|=|X|/|K|`, since nonzero
translations have no fixed time slots.

The representation-theoretic dual is the Fourier transform of `1_X`.  If
`|K|=k`, nonzero Fourier coefficients can only occur at frequencies in the
annihilator of `K`, i.e. multiples of `k`.

## Evidence

`lrc_time_burnside_fourier_s571.py` uses

```text
N = lcm(14,12,15,16,23) = 38640
```

so the grid contains the `n=14` clock and the main S564/S570 witnessed
denominators.

Audited rows:

```text
AP_wall:                6 lonely slots, density 1/6440
V_star_wall:            6 lonely slots, density 1/6440
near_AP_apex:         470 lonely slots, density 47/3864
S562_packet_n14:      276 lonely slots, density 1/140
S562_packet_n14_lift: 264 lonely slots, density 11/1610
no_small_pinch_proxy:4810 lonely slots, density 481/3864
random_low_resonance:4966 lonely slots, density 2483/19320
```

All audited primitive `n=14` binary lonely time words have trivial stabilizer
on this grid.  So the full cyclic Burnside quotient does not identify lonely
slots in these rows; it only reports density.

The Fourier side is more informative.  Dominant nonzero frequencies:

```text
near_AP_apex:          6, 12
S562_packet_n14:      21, 42
S562_packet_n14_lift: 42, 84
no_small_pinch_proxy:  2, 24
random_low_resonance: 34, 2
```

These are finite character shadows of the clocks isolated by HYP-2081 and
HYP-2082: the `n`-clock, pair-sum pinch clocks, dyadic packet gears, and generic
low-resonance clocks.

HYP-2083 adds a complementary finite clock family: antipodal summand shells
modulo `2n-1` acted on by unit multiplication.  In this Burnside/Fourier view,
those shells are not merely additive residue data; they are candidate labelled
time-event orbits whose unit-visible and nonunit-hole branches should have
different character support.

## Interpretation

For primitive integer speeds, every runner resets after one lap.  Raw total
reset length is constant and cannot classify difficulty.  The useful finite
reset invariant is the time word:

```text
period/stabilizer on the primal side,
Fourier character support on the dual side.
```

On the audited rows, the binary word's stabilizer is trivial, but Fourier
concentration still separates clock families.  This suggests the next useful
Burnside object is not the binary lonely word alone, but an owner-labelled
time-event word carrying runner owner, pair-sum denominator, `G` component,
wall endpoint, and endpoint-core labels.

## Tournament Analysis

Vertices were speed-set time words.

Pairwise observable:

```text
(lonely density, stabilizer size, Fourier entropy, top nonzero Fourier mass,
minimal period).
```

Switch/gauge:

```text
harder row wins if it has lower lonely density, then larger reset folding, then
lower Fourier entropy / larger top mode.  Ties use displayed sample order.
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3_cycles=0
sccs=[1,1,1,1,1,1,1]
hamiltonian_paths=1
```

## Assumption Challenge

Possible time vertices considered:

```text
individual time slots, runner phase states, pair-sum events, wall crossings,
Fourier modes, subgroup periods, speed-set time words, and owner-labelled event
words.
```

Chosen quotient:

```text
vertices = whole speed-set lonely time words.
```

Predicate preserved:

```text
which grid times are lonely, the legal translation stabilizer, and the Fourier
character support.
```

Information destroyed:

```text
continuous interval endpoints, which runner owns each wall crossing, pair-sum
pinch labels, and endpoint-protection core labels.
```

Challenged assumption:

```text
the full cyclic time group acts on the lonely slots.
```

It usually does not.  Burnside becomes valid only after restricting to the
stabilizer of the lonely word or enriching the object so the full action is
genuinely defined.

## Open Work

1. Replace the binary lonely word by an owner-labelled event word from S564/S570.
2. Run Burnside on pair-sum events, wall-crossing events, and endpoint cores.
3. Use Fourier mode vertices directly in Tournament Analysis.
4. Add irrational/incommensurate families by approximating long Kronecker
   orbits and comparing their Fourier flattening to integer reset rows.

## Files

- `04-computation/lrc_time_burnside_fourier_s571.py`
- `05-knowledge/results/lrc_time_burnside_fourier_s571.out`
- `07-reflections/lrc-time-burnside-fourier-quotient-s571.md`
