---
id: HYP-3150
title: LRC14 function-compression packets and the degree-4 resolvent guardrail
status: SYNTHESIS / two executable scouts plus S71/S72/S278/KPS mainline integrations; not a proof
source: codex-2026-06-28 + codex-2026-06-28-S277
tangent: T1215
technique: LTI-276
tournament_technique: LTT-174
script: 04-computation/lrc14_function_compression_resolvent_wall_codex_20260628.py
result: 05-knowledge/results/lrc14_function_compression_resolvent_wall_codex_20260628.out
reflection: 07-reflections/lrc14-function-compression-resolvent-wall-codex-20260628.md
additional_scripts:
  - 04-computation/lrc14_function_compression_resolvent_degree_codex_s277.py
additional_results:
  - 05-knowledge/results/lrc14_function_compression_resolvent_degree_codex_s277.out
additional_reflections:
  - 07-reflections/function-compression-resolvent-degree-codex-s277.md
related:
  - HYP-3151
  - HYP-3152
  - HYP-3153
  - HYP-3161
  - HYP-3199
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
  - HYP-3124
  - HYP-3122
  - HYP-3063
  - THM-084
  - THM-577
  - OPEN-Q-108
external_sources:
  - https://github.com/davidturturean/erdos-870
  - https://www.erdosproblems.com/870
mainline_integrations:
  - 04-computation/lrc_arc_cube_compression_parity_macmini_S71.py
  - 05-knowledge/results/lrc_arc_cube_compression_parity_macmini_S71.out
  - 05-knowledge/hypotheses/HYP-3161-score-compression-boundary-is-k8-parity-split.md
  - 05-knowledge/hypotheses/HYP-3151-worpitzky-function-compression-resolvent.md
  - 04-computation/lrc14_worpitzky_function_compression_resolvent_codex_s278.py
  - 05-knowledge/results/lrc14_worpitzky_function_compression_resolvent_codex_s278.out
  - 05-knowledge/hypotheses/HYP-3152-leeyang-circle-coverage-radius-galois-s4.md
  - 05-knowledge/hypotheses/HYP-3199-lrc14-n4-einheit-erdos870-tournament-models.md
  - 04-computation/lrc_coverage_lee_yang_lambda_kps.py
  - 05-knowledge/results/lrc_coverage_lee_yang_lambda_kps.out
  - 05-knowledge/results/lrc_lee_yang_corr_highk_kps.out
---

# HYP-3150: Function Compression and the Degree-4 Guardrail

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

## S277 Exact Finite Scout

The S277 scout
`04-computation/lrc14_function_compression_resolvent_degree_codex_s277.py`
stores its output at
`05-knowledge/results/lrc14_function_compression_resolvent_degree_codex_s277.out`.
It recomputes the same quotient laws with a second, more theorem-interface
oriented ledger.

Pair-function samples `(2,3),(2,4),(4,2),(3,5),(5,3),(4,4)` show that sum
and product are swap-invariant, while the exponential channel is not:
`(2,3)` gives `8` versus `9`, and `(3,5)` gives `243` versus `125`.
Equalities such as `(2,4)` and `(4,2)` are accidents, not quotient laws.

The K3 edge-flip computation recovers the quotient kernel in class-count form:

```text
        to T  to C
from T    2     1
from C    3     0
```

with labelled class sizes `T=6,C=2`, raw flip counts `T->T=12`,
`T->C=6`, `C->T=6`, stationary distribution `T=3/4,C=1/4`, and eigenvalues
`1,-1/3`.  The same matrix appears for three coin flips after quotienting to
`mix` versus `same`.

The role ledger is the proof payload:

```text
cycle_edge_breaks_to_T       6
majority_edge_self_flip     12
minority_edge_closes_cycle   6
```

Thus the score-class quotient knows transition multiplicity but forgets which
edge is the minority gate.  The Worpitzky check verifies
`x^3 = C(x+2,3) + 4*C(x+1,3) + C(x,3)` for `x=0..8`; the coefficient row
`(1,4,1)` is a basis/curve payload, not a license to keep one value.

For the n=4 tournament tables, S277 recomputes the fixed-path fibers:

```text
T = {E}
+ = {a}
- = {b}
S = {c, ab, ac, bc, abc}
```

and the monotone nonlinear compression:

```text
x = a OR c
y = b OR c
```

This compression is class-preserving but does not preserve fiber arithmetic.
It is the small finite model of the proof rule: `c` is filler when fixed and
canary/restoration debt when erased.

The synthesis degree ledger is:

```text
symmetric_pair_shadow          degree 1
ordered_pair_channel           degree 2
K3_edge_flip_kernel            degree 2
N4_filler_canary_square        degree 2
k8_bounded_core_resolvent      degree 4
illegal_raw_scalarization      degree 5
```

Read as a guardrail, accepted packet carriers stay in degree `<=4`.  Generic
degree-5 behavior is not compressed; it is routed to named Abel-Ruffini debt.
This does not prove the terminal k=8 packet, but it states the admissibility
condition before that packet can be used.

## Tournament Analysis

The first scout uses proof carriers as vertices, not runners or raw arcs.

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

S277 also ran a second proof-carrier tournament with weighted retained payload
as the pairwise observable.  That tournament is transitive:

```text
score_hist={17:1,35:1,62:1,64:1,68:1,72:1,76:1,79:1,82:1,85:1}
directed_3cycles=0
hamiltonian_path_count=1
selected_path =
  function_compression_legality_certificate
  -> bounded_core_degree4_guardrail
  -> fiber_PGF_full_curve_payload
  -> ordered_pair_tip_tail_sidecar
  -> n4_canary_filler_exact_transversal
  -> n3_edge_flip_minoritary_gate
  -> worpitzky_basis_curve
  -> A000568_n_le_7_tameness_window
  -> symmetric_sum_product_shadow
  -> raw_scalar_class_count
```

The two tournament readouts should be kept together: the SCC readout says the
sidecars mutually constrain legality, while the transitive readout gives the
current proof-interface ordering once those sidecars are declared.

## Mainline Integration

While this note was being completed, incoming mainline work added two directly
relevant scouts.

The HYP-3161/S71 arc-cube compression scout reports:

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

HYP-3199 adds the exact n=4 Einheit/minimality chart: the fixed-path
`a,b,c` tiling table is a cover with `S={c,ab,ac,bc,abc}` and quotient
ambiguity, while the partial-score `(0,1,1,2)` chart on `x,y` is the exact
two-coordinate section.  This is the same legality rule in smaller clothing:
representation abundance is not proof evidence until the minimal chart,
deletable coordinate, and restoration sidecar are named.

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
function_payload_type
fiber_collision_class
ordered_sidecar_required
unordered_pair_survival
ordered_pair_sidecar
fiber_pgf_curve_status
state_level_pgf_split
compression_map_word
canary_deletion_status
canary_filler_status
even_resolvent_variable
resolvent_degree_ceiling
effective_degree
quintic_wall_alarm
abel_ruffini_wall_status
finite_tameness_window
quotient_legality_status
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
