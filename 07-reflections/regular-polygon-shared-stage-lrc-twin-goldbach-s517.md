---
source: oracle-2026-06-01-S517
status: synthesis (tournaments as oriented regular polygons; LRC <-> twin-Goldbach)
tags:
  - lonely-runner
  - twin-primes
  - goldbach
  - regular-tournament
  - regular-polygon
  - wheel
  - complement-necklace
---

# The Regular Polygon as the Shared Stage: LRC, Twin-Goldbach, and Tournaments as Oriented n-gons

Seeding session: think of a tournament as a **regular polygon on the unit circle
with oriented chords**, and watch the Lonely Runner and the twin-Goldbach
`6k±1` wheel (S516) turn out to be the *same* object — vertex
covering/avoidance on the roots of unity.

## Tournaments = oriented regular polygons

Put `n` points on the unit circle and orient every chord by the half-turn rule
(`i -> j` iff `j` is within the leading semicircle of `i`). A tournament is then
an **orientation of the chords of a regular polygon**. The most symmetric
orientation — `i` beats the next `⌊(n-1)/2⌋` vertices clockwise — is the
**rotational / regular tournament**, the "Platonic" tournament. Equally-spaced
points (roots of unity) are the regular polygon; arbitrary points are an
irregular polygon; the speeds of a runner system decide *which* polygon you sample
as time runs.

## The LRC tight witness IS a regular tournament

At the tight witness of the extremal set — speeds `{1,…,n-1}`, time `t = 1/n` —
runner `i` sits at vertex `i/n`, i.e. **the runners occupy every non-observer
vertex of the regular `n`-gon**, leaving only the observer's vertex (the "clasp")
empty. Computed (`lrc_polygon_tournament_wheel_s517.py`):

```
n=4 -> 3-vtx  H=3   regular  (triangle = R_3)
n=6 -> 5-vtx  H=15  regular  (R_5)
n=8 -> 7-vtx  H=175 regular  (R_7)
n=5 -> 4-vtx  H=5   near-reg     n=7 -> 6-vtx H=41 near-reg     n=9 -> 8-vtx H=629 near-reg
```

So **for even `n` the LRC witness is a genuinely regular tournament `R_{n-1}`**
(`n-1` odd, where a regular tournament exists); for odd `n` it is near-regular
(`n-1` even, no regular tournament). In particular the even composite frontiers
land on a regular tournament:

> n = 14  ->  R_13      n = 18  ->  R_17       (regular tournaments on 13, 17 vertices)

The lonely runner's hardest, tightest configuration is the **most symmetric
polygon** — the regular tournament. Loneliness@n is the question of whether the
runner orbit can reach (the neighbourhood of) this regular-polygon configuration
with the clasp vertex empty. The "extremal = regular = the staircase" thread (S25)
and "extremal = minimal clock" (S24) are the same statement: maximal symmetry.

## Both LRC and twin-Goldbach are vertex problems on the regular m-gon

The residue wheel mod `m` is *literally* the regular `m`-gon (the `m`-th roots of
unity). On it:

- **Twin-Goldbach.** Twin primes occupy only the **unit vertices** (coprime to
  `m`): mod 6 that is `{1,5}`, the two vertices adjacent to the `0`-clasp.
  "Sum of two twin primes" is the **self-convolution of the unit-vertex set**.
  Computed: `{1,5} ⊕ {1,5} = {0,2,4} (mod 6)` — all even residues — so the three
  channels (`1+1, 1+5, 5+5`) are *residue-complete*; the 35 misses are
  **magnitude deserts within channels**, and a triple `{6m-2,6m,6m+2}` is exactly
  one number per channel failing at once.
- **LRC.** At `t = a/q` the runners occupy the vertices `{v_i a mod q}` of the
  regular `q`-gon. Loneliness = the **clasp arc around vertex 0 is empty**. The
  sieve (THM-369) says a counterexample must **occupy vertex 0 of the `m`-gon for
  every `m ≤ n`** — cover the clasp on every wheel.

So the two problems are **dual vertex conditions on the regular polygon**:
twin-Goldbach wants a sparse unit-vertex set to *cover* every even target by
convolution; LRC wants the runner orbit to *avoid* the clasp at some scale. Both
fail only on a **finite, channel-structured exceptional set**.

## Channels = chord-length classes; twins = the d=2 chord

The regular `m`-gon has `⌊m/2⌋` chord-length classes (edge orbits) `d = 1,…,⌊m/2⌋`.
These are the "channels." Two facts tie the threads:

- **Twin primes are the `d = 2` chord class** on every wheel (they differ by 2).
- **LRC holdback (S25) `= 1/(2d)`** is indexed by the same chord length: a pair at
  chord `d` holds its orientation `1/(2d)` of the lap. The extremal staircase of
  differences is the chord-length multiset of `{1,…,n-1}`.

The wheels relevant to the frontiers: `14`-gon has `7` channels, units
`{1,3,5,9,11,13}`; `18`-gon has `9` channels, units `{1,5,7,11,13,17}`. Both
unit-sets are `≅ Z/6` (φ = 6) — the same "complement necklace" clasp structure as
the hexagon, one wheel up.

## The complement necklace, unified

On each wheel the **unit vertices** (where primes/twins live, and where the LRC
unit-witnesses `a/n` sit) are the **complement** of the forbidden residues; the
clasp at vertex 0 is the self-complementary join. The regular tournament is the
rotationally-symmetric chord-orientation of this necklace. So:

> **One stage, three guises.** Roots of unity = regular polygon = residue wheel =
> tournament vertices. Twin-Goldbach is *convolution* of the unit beads; LRC is
> *avoidance* of the clasp by the orbit; a tournament is an *orientation* of the
> chords; and the extremal/symmetric configuration (regular tournament, full
> polygon minus clasp) is exactly the tight LRC witness.

## What this seeds

1. **Regularity is extremal and finite.** The regular tournament (LRC witness),
   the finite circular menu (S24), the 35 twin-Goldbach misses, the 5 Platonic
   solids — all are "few, because symmetry is rigid." LRC's no-counterexample and
   twin-Goldbach's finite-exceptions are the same flavour of statement: a sparse
   set on the wheel covers all-but-a-symmetry-controlled-finite-set.
2. **Channel-simultaneity is the obstruction.** A twin-Goldbach triple fails all
   chord-sum channels at once; an LRC counterexample must occupy the clasp on all
   `m`-wheels at once. Progress on either is "rule out simultaneous channel
   failure."
3. **Even frontiers land on regular tournaments.** n = 14 -> R_13, n = 18 -> R_17.
   The LRC@even-n tight stratum is the regular tournament on `n-1` (prime, here)
   vertices — connect this to the repo's Paley/regular-tournament `A000568`
   structure (HYP-1987 source-target) and to whether `R_{n-1}`'s symmetry forces
   the source-class reachability.

## Next
- Compute the runner-occupied vertex set `{v_i a mod q}` as a function on the
  `q`-gon and read loneliness as "clasp arc empty"; cross with the twin-Goldbach
  convolution on the same `q`.
- Test whether the even-`n` LRC witness `R_{n-1}` is the unique reachable
  regular source-class (tie to HYP-1987).

## Artifacts
```
04-computation/lrc_polygon_tournament_wheel_s517.py
05-knowledge/results/lrc_polygon_tournament_wheel_s517.out
07-reflections/twin-goldbach-35-exceptions-the-six-wheel-s516.md  (the wheel)
```
