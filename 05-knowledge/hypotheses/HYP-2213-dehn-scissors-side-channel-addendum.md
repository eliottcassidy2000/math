# HYP-2213: Dehn and Scissors Side Channels Extend the Rational Shadow Carrier

**Status:** OPEN addendum to HYP-2211.  HYP-2211 owns the corrected
`e+pi`/`e*pi` two-shadow carrier obstruction; this hypothesis extends that
thread through Dehn invariants, cuboid/simplex scissors congruence, and
unit-distance impairment tests.

## Claim

Rational versus irrational should not be treated like odd versus even.  Parity
is a genuine quotient `Z -> Z/2`; irrationality is not closed under addition or
multiplication.  HYP-2211's usable replacement is a retained side channel:

- for `e+pi` and `e*pi`, the side channel is algebraic dependence over the
  algebraic numbers;
- for cuboids and simplices, the side channel is the Dehn invariant beyond
  volume;
- for tournament `H`, the side channel is the strong-component/fiber packet
  beyond the scalar Hamiltonian-path count;
- for unit-distance constructions, the side channel is the spine/bulk,
  traceability, direction-support, and embedding packet beyond edge count;
- for LRC, the side channel is the `C=2n-1` unit/nonunit/lift/carry/owner
  packet beyond the raw pair-sum modulus.

The corrected mathematical seed is:

```text
At least one of e+pi and e*pi is transcendental.
```

Indeed, if `s=e+pi` and `p=e*pi` were both algebraic, then `e` and `pi` would
be roots of `x^2 - s*x + p`, hence algebraic over the algebraic numbers, a
contradiction.  It is not known that exactly one of `e+pi` and `e*pi` is
rational; if either one were algebraic or rational, the other would have to be
transcendental.

## Evidence

S637 adds `04-computation/dehn_scissors_shadow_s637.py`, which records the
scissors addendum:

| channel | scalar quotient | retained side channel |
|---------|-----------------|-----------------------|
| `e+pi`, `e*pi` | two symmetric coordinates | algebraic-dependence trap |
| cuboid/simplex | volume | Dehn invariant in `R/(pi Q)` |
| tournament `H` | Hamiltonian-path count | SCC multiset, `beta/c3`, OCF fibers |
| unit-distance `n=21` | `57` edges | `20`-edge spine plus `37=C_hex(3)` bulk |
| LRC pair-sum sieve | `C=2n-1` | unit/nonunit/lift/carry/owner labels |
| SC tournaments | fixed count | `Aut/Anti` transporter and perspective flip |

The Dehn row uses the regular tetrahedron: its dihedral angle has
`cos(theta)=1/3`, so `theta/pi` is irrational by Niven's rational-angle
theorem.  Hence the regular tetrahedron has nonzero Dehn side data, while a
cuboid has zero Dehn invariant.  Equal volume or literal packing is therefore
not enough for scissors congruence.

S637's proof-lens Tournament Analysis uses proof lenses as vertices.  Across
algebraic-certainty, scissors-retention, and repo-portability gauges, the
majority tournament is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}
directed_3cycles=0
H=1
```

The top lenses are `Dehn invariant`, `carrier compression`,
`algebraic-dependence channel`, and `symmetric-polynomial trap`; the raw
`exactly one rational` claim ranks last.

## Technique Program

For unit-distance construction/proof work, deliberately build small impaired
pairs that agree in a scalar and disagree in a side channel:

1. Same edge count, different traceability or number of unit Hamiltonian
   spines.
2. Same edge count, different direction-support mask or frontier-gain price.
3. Same `n=21`/`57` scalar, different `20+37` spine/bulk realization.
4. Same graph carrier, different embeddability or totally-unfaithful
   obstruction label.
5. Same tournament `H`, different `beta/c3` packet or strong-component
   scissors class.

The small-size goal is not just to find bigger examples.  It is to learn which
side channel changes first when a proof method is damaged.

## Next Tests

1. Add a unit-distance side-channel jackknife that fixes edge count but varies
   traceability, direction support, and endpoint-compatible ears on small
   Moser/triangular carriers.
2. Build an `H`-fiber atlas for low tournament values that records
   strong-component packet, `beta1`, `c3`, and OCF independence data, then
   compare it to the phantom-volume gaps `{7,21}`.
3. For LRC, treat the `C=2n-1` modulus identity as the volume row and measure
   which unit/nonunit/lift/carry labels act most like Dehn side data under
   quotienting.
4. For hard sequences, extend S633 by adding "same scalar, different side
   channel" shadow pairs rather than only fixed/merged/nonfixed companions.

## See Also

`04-computation/dehn_scissors_shadow_s637.py`;
`05-knowledge/results/dehn_scissors_shadow_s637.out`;
`07-reflections/dehn-scissors-side-channel-addendum-s637.md`;
HYP-2211; HYP-2210; HYP-2209; HYP-2208; HYP-2207; HYP-2206; HYP-2187; HYP-2186.
