# LRC14 gK8 and q6 Arithmetic Kernels

S79 started on the HYP-2812 concentration-extremality route and pulled a
concurrent mac-mini S23 result that landed directly on the same gap.  The new
q6 certificate is the missing shape of the small-f guard:

```text
q6(B u {f}) = q6(B)/7 + endpoint_periodic_error(B,f)/f.
```

For fixed bounded `B`, the numerator `f*(q6(B u {f}) - q6(B)/7)` is a
THM-563-style sawtooth sum over endpoints of the all-missed set `Q6(B)`, hence
periodic in `f` with period `7*lcm(B)`.  This turns the vague "q6 suppresses by
about 1/7" mechanism into a finite endpoint-period obligation.

Exact incoming arithmetic:

| row | q6 base | period-max | f>=15 ratio bound |
| --- | --- | --- | --- |
| consec k=9 base `1..8` | `1/56` | `6/49` | `3/5` |
| consec k=10 base `1..9` | `1/63` | `6/49` | `23/35` |
| 15-base scout worst | reported | reported | `33/35` |

I added `TournamentH7.LRCQ6Contraction` to record exactly these rational
reductions.  This is intentionally only an arithmetic kernel: the sawtooth
identity and the period scans remain in the Python certificates.

I also extended `TournamentH7.LRCFactorialAtom` with the full all-binding-row
`gK8` finite-check table for k=8..13.  The row-independent Delsarte domination
`10*q0 <= 10q0+q3+10q6` was already formalized; the new theorem
`capClear_gK8_all_binding_rows` packages the scalar cap comparisons:

```text
k=8:  2633/735   <= 2243/588
k=9:  3259/735   <= 9895/2002
k=10: 37/7       <= 550/91
k=11: 26603/4410 <= 660/91
k=12: 29287/4410 <= 60/7
k=13: 61529/8820 <= 10
```

Focused builds succeeded:

```text
lake build TournamentH7.LRCQ6Contraction
lake build TournamentH7.LRCFactorialAtom
```

I attempted `lake build TournamentH7.Verify`, but it expanded into broad
unrelated Mathlib/category-theory compilation and was stopped.  No aggregate
Verify result is claimed.

## Proof State

The arithmetic ledger is now mostly clean:

```text
per-shape Delsarte domination
  + all-row bounded gK8 table
  + q6 strict-contraction arithmetic
```

What remains is not another rational comparison.  The live theorem is a
generated-profile smoothing or Krawtchouk-majorization statement:

```text
wide far insertion loses enough weighted q6 mass,
and admissible q0/q3 movement cannot overtake that loss,
so L_yK8(wide) <= max_bounded L_yK8.
```

The q6 endpoint-period result proves one extreme atom contracts, including at
small `f>=15`.  It does not by itself control q0 gains or q3 redistribution.
That is the remaining structural lemma.  If it resists, the fallback remains
the generalized-doublet address atlas plus the closed-form `12*zeta(3)`
Tornheim R-tail.

## Tournament Analysis

Vertices considered this session:

- all-row gK8 cap-clear kernels;
- q6 endpoint-period obligations;
- q6 contraction rows;
- generated-profile smoothing lemmas;
- generalized-doublet/R-tail fallback atlases;
- raw runner configurations.

Observable: proof leverage toward `L_yK8 <= 10cap`.

Switch: orient toward the vertex that preserves more of the active proof
predicate with less unproved geometry.  Ties use the Hamiltonian path below.

Hamiltonian path:

```text
generated_profile_smoothing
> q6_endpoint_period_bound
> all_row_gK8_Lean_table
> bounded_LyK8_certificate
> generalized_doublet_Rtail_fallback
> scalar_p0_doublet_max
> raw_runner_vertices
```

The challenged assumption is that runner-level vertices are the natural proof
objects.  They are not.  The useful vertices are proof obligations and profile
coordinates; raw runners obscure the difference between a proved q6
contraction and the still-open q0/q3 profile movement.

## Post-Fetch Addendum: Moment Region And Relation Depth

After pulling HYP-2823 and HYP-2828, the clean target changed from "use gK8
somehow" to a precise feasible-moment inequality:

```text
10 - 10*S1 + 10*S2 - 9*S3 + 6*S4 <= 10*cap_k.
```

I added Lean names for both faces of this identity.  In
`TournamentH7.LRCFactorialAtom`, `LyK8_moment_form` and
`LyK8_probability_moment_form` package the degree-4 factorial-moment side,
while `LyK8_extremeMass_readout` and
`LyK8_moment_extremeMass_identity` package the `10*(q0+q6)+q3` side.

HYP-2828 then becomes an exception taxonomy for the same moment-region problem.
Depth-2/two-peel bounded rows are not a new proof route; they are the finite
resonant lane to generalized-doublet/R-tail if the global moment certificate
needs local discharge.  Depth>=3 is the branch where the desired proof should
show decorrelation pushes mass out of the extremes and into the middle
miss-count atoms.
