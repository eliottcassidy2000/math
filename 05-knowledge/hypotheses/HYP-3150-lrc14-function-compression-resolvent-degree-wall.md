---
id: HYP-3150
title: LRC14 function-compression packets may keep the hard core below the Abel-Ruffini wall
status: EVIDENCE / executable factor-through scout and proof-packet proposal; not a proof
source: codex-2026-06-28
tangent: T1215
technique: LTI-276
tournament_technique: LTT-174
script: 04-computation/lrc14_function_compression_resolvent_wall_codex_20260628.py
result: 05-knowledge/results/lrc14_function_compression_resolvent_wall_codex_20260628.out
reflection: 07-reflections/lrc14-function-compression-resolvent-wall-codex-20260628.md
related:
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
  - HYP-3135
  - HYP-3134
  - HYP-3133
  - HYP-3132
  - HYP-3129
  - HYP-3122
  - THM-084
  - THM-577
  - OPEN-Q-108
external_sources:
  - https://github.com/davidturturean/erdos-870
  - https://www.erdosproblems.com/870
mainline_integrations:
  - 04-computation/lrc_arc_cube_compression_parity_macmini_S71.py
  - 05-knowledge/results/lrc_arc_cube_compression_parity_macmini_S71.out
  - 04-computation/lrc_coverage_lee_yang_lambda_kps.py
  - 05-knowledge/results/lrc_coverage_lee_yang_lambda_kps.out
---

# HYP-3150: Function-Compression Below the Resolvent-Degree Wall

## Claim

The user's Worpitzky/function/K3/K4/resolvent prompt is usefully formalized as
a single legality calculus:

```text
compression q: X -> Y
observable  f: X -> Z
legal iff   f factors through q,
or the missing part is carried by a named sidecar.
```

The executable scout verifies this across three adjacent scales:

- HYP-3147/HYP-3144: K3 edge flips collapse to a two-class kernel, but the
  minority edge, Worpitzky descent word, state-level PGF curve, and ordered
  pair functions may not factor through the score-class quotient.
- HYP-3149/HYP-3146/HYP-3148: the K4 fixed-path cube factors to a legal
  two-bit class table only through a non-linear compression
  `x = a OR c`, `y = b OR c`, or by fixing the canary/filler coordinate.
- HYP-3142/HYP-3132/HYP-3139: the k=8 bounded-core dual is a quartic
  resolvent, but after the reflection coordinate `u=t-3` it is even:
  `u^4 - 5u^2 + 4 = (u^2-1)(u^2-4)`, so the proof-facing variable is
  the quadratic `v=u^2` together with resurrection sidecars for odd data.

The speculative meta-point to test is:

```text
LRC14 stays solvable because every hard compression seen so far has
effective degree <= 4, and the deepest k=8 node is a biquadratic below the
generic quintic / Abel-Ruffini wall.
```

The evidence supports this as a proof-packet heuristic, not as a theorem by
itself.  The exact verified content is factor-through legality plus named
sidecar debt.

## Evidence

The scout
`04-computation/lrc14_function_compression_resolvent_wall_codex_20260628.py`
stores its output at
`05-knowledge/results/lrc14_function_compression_resolvent_wall_codex_20260628.out`.

Pair-function audit:

```text
a+b      FACTORS through unordered-pair quotient
a*b      FACTORS through unordered-pair quotient
a^b      SIDE_CAR_REQUIRED
b^a      SIDE_CAR_REQUIRED
{a^b,b^a} FACTORS as an unordered multiset, but loses orientation
```

This exactly matches the prompt's order distinction: sum/product are
commutative shadows, while the exponentials are ordered channels.

K3 audit:

```text
class_transition_signature FACTORS through C/T
exit_edge_set              SIDE_CAR_REQUIRED
state_forward_edge_PGF     SIDE_CAR_REQUIRED
Worpitzky_descent_word     SIDE_CAR_REQUIRED
```

The exact kernel is the prior one:

```text
class_counts = {C:2, T:6}
raw_counts   = C->{T:6}, T->{T:12,C:6}
normalized rows C,T = [[0,1],[1/3,2/3]]
eigenvalues = (1,-1/3)
```

The strongest warning remains the PGF curve:

```text
C aggregate F = (1,4,1), state curves (0,2,1), (1,2,0)
T aggregate F = (1,4,1), state curves (0,0,1), 4*(0,1,0), (1,0,0)
```

So even a whole aggregate class PGF can hide state-level curve structure.

K4 audit:

```text
score_class              FACTORS through x=a OR c, y=b OR c
flip_weight              SIDE_CAR_REQUIRED
flip_word                SIDE_CAR_REQUIRED
c_canary_status          SIDE_CAR_REQUIRED
delete_one_stable_status SIDE_CAR_REQUIRED
```

The OR fibers are:

```text
(0,0): T, {E}
(1,0): +, {a}
(0,1): -, {b}
(1,1): S, {c, ab, ac, bc, abc}, weight PGF z+3z^2+z^3
```

Thus the two-bit scaffold is a legal class quotient only after the
fiber-PGF, canary, and deletion sidecars are named.

k=8 resolvent audit:

```text
g(t)=(t-1)(t-2)(t-4)(t-5)
coeffs in t: [1,-12,49,-78,40]
u=t-3 roots: (-2,-1,1,2)
g(u+3)=u^4-5u^2+4
v=u^2 gives v^2-5v+4, discriminant 9
```

The even value `g_even(u)` factors through `v=u^2`, but `u` and `u^3` do not.
This is the exact sidecar lesson from HYP-3138/HYP-3139/HYP-3142: the hard
node is not a generic quartic, but the odd coordinate is not erased for free.

## Two Maps

Compression tower:

```text
ordered pair values
  -> unordered pair quotient
     safe: a+b, a*b
     sidecar: a^b, b^a

K3 edge cube
  -> C/T score quotient
     safe: transition kernel
     sidecars: minority edge, Worpitzky descent, state PGF curve

K4 fixed-path cube
  -> OR compression x=a OR c, y=b OR c
     safe: score class
     sidecars: fiber PGF, canary c, deletion stability

k=8 resolvent
  -> even fold v=u^2
     safe: biquadratic value
     sidecars: odd coordinate, reflection resurrection, boundary leakage
```

Sidecar taxonomy:

```text
order sidecar       = endpoint orientation, a^b vs b^a
edge-role sidecar   = K3 minority edge / exit edge
curve sidecar       = state-level PGF/root curve, not one value
canary sidecar      = fixed filler coordinate and collapse slice
deletion sidecar    = which representatives survive coordinate deletion
odd sidecar         = sign/odd coordinate after even resolvent fold
degree sidecar      = effective algebraic degree and wall alarm
```

## Tournament Analysis

The scout uses proof carriers as vertices, not runners or raw arcs.

Pairwise observable:

```text
majority over factor-through exactness, sidecar hygiene, LRC transfer,
GF/root retention, degree control, tournament kernel signal,
canary/deletion signal, and raw-scalar warning.
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,6:2,7:1,8:2,10:1}
directed_3cycles=3
hamiltonian_path_count=9
edge_flips_against_raw_scalar_warning_gauge=13
largest SCC =
  {k8_even_resolvent_v_u2_fold,
   k4_OR_class_compression,
   k3_edge_flip_kernel,
   fiber_PGF_curve_sidecar,
   c_canary_deletion_sidecar}
```

Selected path by score:

```text
factor_through_compression_audit
-> k8_even_resolvent_v_u2_fold
-> k4_OR_class_compression
-> fiber_PGF_curve_sidecar
-> k3_edge_flip_kernel
-> c_canary_deletion_sidecar
-> worpitzky_descent_function_sidecar
-> ordered_pair_exponent_sidecar
-> generic_quintic_wall_alarm
-> symmetric_sum_product_shadow
-> raw_scalar_or_class_value
```

The nontrivial SCC is the useful surprise.  The small tournament kernels,
PGF curves, canary/deletion fields, and k=8 even fold are not a linear ladder.
They are a coupled packet: choosing one sidecar changes which compression is
legal next.

## Mainline Integration

While this note was being completed, incoming mainline work added two directly
relevant scouts.

The S71 arc-cube compression scout reports:

```text
score sequence -> tournament isomorphism class:
n=3: #iso=2,  #score=2, BIJECTIVE=True
n=4: #iso=4,  #score=4, BIJECTIVE=True
n=5: #iso=12, #score=9, BIJECTIVE=False
```

This is exactly the same factor-through lesson in a larger tournament cube:
the commutative score/`a+b` face determines the quotient through n=4 and fails
at n=5, where order-sensitive data is needed.

S71 also aligns that failure with the LRC cap dip:

```text
k=8:  |P|=5, dip=1081/76440, large binding row
k=9:  |P|=4, dip=1/4004, tiny edge row
k>=10 |P|<=3, dip=0
```

So the k=8 obstruction appears exactly at the `n=5` compression boundary.
Even more sharply, S71 splits the k=8 correction

```text
L_yK8 = 10S0 - 10S1 + 10S2 - 9S3 + 6S4
```

into:

```text
ODD  -9S3 = -12.135  orientation / Worpitzky / ordered side
EVEN +6S4 =  +3.853  symmetric / biquadratic / sum-product side
|odd|/|even| = 3.15
```

This upgrades HYP-3150's heuristic: the k=8 dip is not only
`even biquadratic plus sidecar debt`; the odd/Worpitzky face dominates the
local correction magnitude.

The kps Lee-Yang lambda scout adds the root-curve version of the same lesson.
For k=8,9,10 it measures off-circle variance
`lambda = Var(|roots|/R)` of the coverage PGF.  Consecutive/AP rows are the
most circular rows for k=8 and k=9 in the tested bank, and lambda correlates
with the coverage gap:

```text
k=8:  corr(lambda, cap-q0)=+0.695
k=9:  corr(lambda, cap-q0)=+0.752
k=10: corr(lambda, cap-q0)=+0.851
```

This reinforces the HYP-3140/HYP-3109 rule: the whole PGF/root curve is proof
payload.  A scalar value can be right for the wrong reason if it has forgotten
root circularity or off-circle variance.

## LRC14 Transfer

The proof packet should add:

```text
compression_map
observable_factors_through
fiber_collision_class
ordered_sidecar_required
fiber_pgf_curve_status
canary_deletion_status
even_resolvent_variable
effective_degree
quintic_wall_alarm
terminal_exit_or_named_debt
```

Immediate uses:

1. HYP-3140 fiber-PGF rows should be checked at curve/fiber level before a
   coefficient or first moment is scalarized.
2. HYP-3141 edge packets should declare which observables factor through the
   edge quotient and which are tip/tail/order sidecars.
3. HYP-3138/HYP-3139/HYP-3142 should state the k=8 fold as a degree-two
   visible variable plus an odd-coordinate resurrection obligation.
4. HYP-3143/HYP-3146/HYP-3149 K4 quotients should treat the OR compression as
   class-legal but fiber-illegal unless the `S` PGF and canary/deletion debt
   are retained.

## Assumption Challenge

The vertex set for the scout is not fixed in advance.  It should explicitly
compare at least these carriers:

```text
functions, quotient fibers, edge roles, fixed-path arc words,
canary/filler coordinates, generating-function curves, resolvent variables,
Fourier/SPEC modes, and proof obligations.
```

The quotient preserves the LRC predicate only when the target observable is
constant on visible fibers or reconstructible from sidecars.  It destroys
orientation, order, deletion robustness, odd coordinates, or full PGF curves
unless those are named.

## Guardrail

Do not promote the Abel-Ruffini wall phrase into proof.  The exact content is
smaller and stronger: every compression used in the current LRC14 endgame must
either be a finite factor-through map of bounded degree or carry the missing
coordinate as an explicit sidecar.  If a future step produces a genuinely
generic five-branch obstruction with no reflection, filler, or sidecar fold,
that should be treated as a route alarm, not as solved structure.
