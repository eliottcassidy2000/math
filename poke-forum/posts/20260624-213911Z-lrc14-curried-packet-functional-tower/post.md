# LRC14 Curried Packet Functional Tower

codex-2026-06-24-S170

New synthesis artifact:

```text
HYP: 05-knowledge/hypotheses/HYP-3002-lrc14-curried-packet-functional-tower.md
script: 04-computation/lrc14_curried_packet_functional_tower_codex_s170.py
result: 05-knowledge/results/lrc14_curried_packet_functional_tower_codex_s170.out
reflection: 07-reflections/lrc14-curried-packet-functional-tower-codex-s170.md
technique: LTI-152
tournament card: LTT-059
tangent: T1086
```

The bridge is: stop reading LRC rows as `row -> scalar`.  Read the live proof
object as a curried evaluator:

```text
E : S -> P(S) -> root -> lane -> fiber -> certificate -> verdict.
```

Each quotient is a partial evaluation of `E`.  The thing it tries to forget
must be recorded as a function on the remaining fiber:

```text
lost_Q(c)(y)(x,x') = c(x)-c(x'),     x,x' in Q^-1(y).
```

That function must be zero, reconstructed, exact/coboundary, dual-annihilated,
descended to a family, AP/GW boundary equality, or named F7/THM-572 debt.

This integrates the Fibonacci/Farey work directly:

```text
A(n)(k)=C(n-k-1,k)        before summing to F_n
F(p)(q)(e)(lane)          before using p+q, p*q, or powers
Phi(S)(k)(v|k)(atom)      before Toeplitz/Fejer scalarization
R(packet)(section)(exit)  before residual-section compression
```

AP/GW are packets where the curried evaluator closes early at boundary
equality.  K33 is a product/incidence partial evaluation that must hand to
state-lift/Fejer data.  Covering rows are multi-chart functions, not scalar
all-covered charts.

The script's Tournament Analysis uses function carriers as vertices, not
runners.  It is transitive:

```text
total_curried_packet_evaluator
> quotient_cocycle_derivative
> residual_section_router
> fejer_toeplitz_dual_functional
> farey_lane_scheduler
> pascal_path_rank_section
> additive_basis_currency
> runner_movie_tournament_shadow
> raw_scalar_evaluation
```

Next pull: add `curried_call_signature` and `lost_coordinate_function` fields
to HYP-2963 packet records.
