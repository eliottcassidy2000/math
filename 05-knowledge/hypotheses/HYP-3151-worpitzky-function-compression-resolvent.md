---
id: HYP-3151
title: Worpitzky and ordered functions execute the function-compression wall scout
status: SYNTHESIS / exact finite scout; not an LRC14 proof
source: codex-2026-06-28-S278
tangent: T1216
technique: LTI-277
tournament_technique: LTT-175
script: 04-computation/lrc14_worpitzky_function_compression_resolvent_codex_s278.py
result: 05-knowledge/results/lrc14_worpitzky_function_compression_resolvent_codex_s278.out
reflection: 07-reflections/worpitzky-function-compression-resolvent-codex-s278.md
related:
  - HYP-3150
  - HYP-3152
  - HYP-3153
  - HYP-3161
  - HYP-3149
  - HYP-3148
  - HYP-3147
  - HYP-3146
  - HYP-3145
  - HYP-3144
  - HYP-3143
  - HYP-3142
  - HYP-3141
  - HYP-3140
  - HYP-3139
  - HYP-3138
  - HYP-3137
  - HYP-3132
  - HYP-3129
  - HYP-3124
  - OPEN-Q-108
external:
  - https://github.com/davidturturean/erdos-870
---

# HYP-3151: Worpitzky Function-Compression Resolvent Bridge

This is the executable continuation of HYP-3150.  HYP-3150 reserved the
factor-through wall:

```text
compression q: X -> Y
observable  f: X -> Z
legal iff   f factors through q, or the lost coordinate is sidecarred
```

HYP-3151 runs the exact finite scout and adds two extra diagnostics: no affine
`GF(2)^3 -> GF(2)^2` map can replace the n=4 OR compression, and the
proof-carrier tournament over legality packets is transitive.

## Claim

The proof-facing object is not a tournament class, a scalar value, or a raw
generating function.  It is a function together with the compression under
which that function is being used.  A compression is legal only when the target
function is constant on its fibers, or when an ordered/canary sidecar
reconstructs the destroyed coordinate.

This rule merges the user's four function examples:

```text
a+b, a*b     symmetric shadows; safe under unordered pair quotient
a^b, b^a     ordered channels; require orientation sidecar
```

with the recent tournament lanes:

```text
n=3: C/T edge-flip kernel + Worpitzky descent word
n=4: fixed-path cube compressed by nonlinear OR canary map
k=8: quartic dual centered to the biquadratic u^4 - 5u^2 + 4
```

## Evidence

The scout `04-computation/lrc14_worpitzky_function_compression_resolvent_codex_s278.py`
verifies the exact finite data.

For the n=3 edge perspective, the two class kernel is:

```text
P(C -> T) = 1
P(T -> C) = 1/3
P(T -> T) = 2/3
stationary(C,T) = (1/4, 3/4)
eigenvalues = 1, -1/3
```

Worpitzky refines the transitive fiber by descent words.  The degree-3 and
degree-4 rows are verified exactly:

```text
x^3: 1,4,1
x^4: 1,11,11,1
```

The function audit over ordered pairs `1..7` records that `a+b` and `a*b`
are invariant on every oriented pair, while `a^b` is invariant only at the
accidental equalities such as `(2,4)` and `(4,2)`.  Thus symmetric functions
are useful shadows, but ordered exponentials are orientation payloads.

For the n=4 prompt tables, the class-preserving map from the fixed-path cube
to the two-bit table is:

```text
x = a OR c
y = b OR c
```

This nonlinear OR compression is class-preserving on all eight states.  The
scout also exhausts all affine maps `GF(2)^3 -> GF(2)^2` and finds:

```text
affine_class_preserving_count = 0
```

So the canary/filler move is not a harmless linear quotient.  It is a
function-specific nonlinear compression that is legal for the four-class
target only because the whole `S` bulk fiber is sent to `(x,y)=(1,1)`.

For the bounded-core hard node, the k=8 dual remains below the
Abel-Ruffini wall:

```text
g(t) = (t-1)(t-2)(t-4)(t-5)
     = t^4 - 12t^3 + 49t^2 - 78t + 40

u = t - 3:
g(u+3) = u^4 - 5u^2 + 4
odd centered coefficients = (0,0)
discriminant in u^2 = 9
```

The degree ceiling is `4`, and the centered quartic is biquadratic.  This
supports the meta-point: the LRC14 bounded core reaches the deepest quartic
case at k=8 but never crosses into the generic quintic obstruction.

## Incoming Mainline Integration

After rebasing over incoming main, HYP-3151 should be read as the executable
continuation of two new signals.

The S71 arc-cube parity scout reports that score-sequence compression is
bijective on tournament isomorphism classes for `n<=4` but fails at `n=5`.
The cap dip turns on at the same boundary: `|P|=5` for k=8, with large dip
`1081/76440`, while `|P|<=3` rows have no dip and k=9 at `|P|=4` has only
the tiny `1/4004` dip.  It also splits the k=8 correction:

```text
even side: +6 S4   symmetric / biquadratic / a+b,a*b
odd side:  -9 S3   antisymmetric / Worpitzky / a^b,b^a
```

with `|odd|/|even| = 3.15` on `consec_8`.  This is exactly the HYP-3151
message in sharper form: the solvable even resolvent is not enough unless the
ordered/Worpitzky sidecar is retained.

The KPS Lee-Yang lambda scout reports that AP/consec minimizes the coverage
PGF off-circle variance `lambda(phi4)` at k=8 and k=9, and that correlation
between `lambda` and the coverage gap grows from about `+0.695` at k=8 to
`+0.851` at k=10.  Thus the function-compression packet should not measure
only `q0` or one scalar; it should retain the root-curve / off-circle
`phi4` sidecar that tells whether the even biquadratic page is legally
controlling the odd correction.

## LRC14 Transfer

Add these packet fields before scalarizing edge witnesses, fiber PGFs, n=4
tournament quotients, or k=8 resolvent pages:

```text
target_function_id
function_swap_parity
symmetric_shadow_status
ordered_sidecar_status
edge_flip_kernel_mode
minority_edge_gate
worpitzky_descent_word
compression_map_kind
compression_fiber_function_constancy
canary_or_restoration_sidecar
resolvent_degree
centered_odd_coefficient_status
abel_ruffini_wall_status
terminal_exit_or_named_debt
```

The immediate proof-frontier use is HYP-3141/HYP-3124 edge rows and
HYP-3140/HYP-3142 coefficient/moment packets.  A packet can forget a coordinate
only after its target function is fiber-constant, symmetric-safe, ordered-sidecar
reconstructed, canary-restored, descended to the biquadratic resolvent, or
routed to named debt.

## Assumption Challenged

Tournament vertices need not be runners, arcs, or isomorphism classes.  Here
the live vertices are functions, compression fibers, ordered sidecars,
Worpitzky descent words, canary coordinates, and resolvent-degree obligations.
The preserved LRC predicate is target-function legality under quotienting; the
destroyed data is orientation, canary state, odd centered coefficients, and
fiber membership when no sidecar is retained.
