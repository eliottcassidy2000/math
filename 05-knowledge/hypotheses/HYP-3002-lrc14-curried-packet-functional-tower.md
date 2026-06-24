---
id: HYP-3002
title: LRC14 curried packet functional tower
status: SYNTHESIS / proof-interface guardrail; not a proof
source: codex-2026-06-24-S170
script: 04-computation/lrc14_curried_packet_functional_tower_codex_s170.py
result: 05-knowledge/results/lrc14_curried_packet_functional_tower_codex_s170.out
related:
  - HYP-3000
  - HYP-2999
  - HYP-2998
  - HYP-2997
  - HYP-2996
  - HYP-2995
  - HYP-2990
  - HYP-2974
  - HYP-2963
  - HYP-2954
  - THM-572
  - LTI-152
  - LTT-059
---

# HYP-3002: LRC14 Curried Packet Functional Tower

## Claim

The clean way to integrate HYP-2998/HYP-2999/HYP-3000 with the LRC14 proof
stack is to treat every proof object as a curried function.  The central object
is not a scalar invariant but a total evaluator:

```text
E : S -> P(S) -> root -> lane -> fiber -> certificate -> verdict
```

where

```text
S           = runner row / candidate row
P(S)        = labelled packet record
root        = exact M/Farey node, q, excess e=14p-q
lane        = root, p+q, p*q, q^p, p^q, additive-basis, harmonic, section
fiber       = the remaining packet fiber after current quotient choices
certificate = Haar, Fejer, Ramanujan, section, boundary moment, state lift
verdict     = strict-safe, AP/GW boundary equality, or named residual debt
```

Currying matters because a quotient is exactly a partial application of this
evaluator.  The proof obligation is to record what the partial application has
fixed and what coordinate it has tried to forget.

## Quotient Derivative

For a quotient

```text
Q : P -> Y
```

and a coordinate `c` that `Q` forgets, attach the curried lost-coordinate
function

```text
lost_Q(c)(y)(x,x') = c(x) - c(x'),        x,x' in Q^{-1}(y).
```

The quotient is legal only if this function is:

```text
zero on each fiber,
reconstructible from retained labels,
a coboundary / exact derivative on the fiber,
annihilated by a dual certificate,
descended to a named family,
an AP/GW boundary equality atom,
or routed to F7/THM-572 residual debt.
```

This is the HYP-2990 no-free-slider rule, HYP-2995 packet-cocycle rule, and
HYP-2997 cocycle normal form written as an ordinary higher-order function.

## Fibonacci / Farey Currying

HYP-3000 says the user's row

```text
1, 1, 1+1, 1+2, 1+3+1, 1+4+3, 1+5+6+1, ...
```

is the partially evaluated path-rank function

```text
A(n)(k) = C(n-k-1,k),        F_n = sum_k A(n)(k).
```

So the scalar Fibonacci number is the final evaluation `sum_k`.  Before that
evaluation, the row still remembers independent-set size, carry width, and
no-adjacent support.  Zeckendorf is a section of the resulting evaluation map:
it chooses a canonical legal support rather than only counting supports.

The Farey side should be read the same way:

```text
F(p)(q)(e)(lane)
```

with `e=14p-q`.  For the LRC14 unit-excess chain `p/q=p/(14p-1)`:

```text
root:    F(p)(14p-1)(1)(root)    = (p,14p-1,1)
sum:     F(p)(14p-1)(1)(sum)     = 15p-1
product: F(p)(14p-1)(1)(product) = 14p^2-p
power:   F(p)(14p-1)(1)(power)   = ((14p-1)^p,p^(14p-1)).
```

The additive lane is a safe partial evaluation only with exact root retained.
The product lane is an incidence-side function, not an order scalar.  The
power lane is a stress-test function for magnitude-blind quotients.

## LRC Integration

The LRC predicate itself is already curried:

```text
danger(S)(t)(v) = 1_{||vt|| < 1/14}
C(S)(t)         = sum_v danger(S)(t)(v)
safe(S)(t)      = 1_{C(S)(t)=0}.
```

The Fourier/Fejer route curries the same object differently:

```text
Phi(S)(k)(v | k)(packet_atom) -> coefficient / quadratic form.
```

That is why the divisor-curried Toeplitz formula in HYP-2974 was a real bridge:
low harmonics query speed-divisor fibers rather than raw speed counts.  The
current HYP-3002 addition says every other route should expose an equally
typed partial evaluation:

```text
Residual sections:     R(packet)(section)(grid_exit)
Farey scheduler:       F(p)(q)(e)(lane)
Pascal path ranks:     A(n)(k)
Additive basis:        B(atom_system)(N)(budget)
Packet quotient audit: lost_Q(c)(y)(x,x')
```

## Named-Row Readout

```text
AP:
  E(AP)(boundary_packet)(e=0)(root/sum)(same_tile)(no-open)
  -> AP/GW boundary equality atom.

GW:
  E(GW)(C27_boundary_packet)(e=0)(root/sum)(same_tile)(no-open)
  -> AP/GW boundary equality atom with C27 transfer label.

K33 12->36:
  E(row)(K33_packet)(p=3)(product)(cross_handoff)(Fejer degree 159)
  -> strict-safe certificate or named K33/state-lift debt.

P10+GW:
  E(row)(unit_petal_packet)(root)(sum/product)(owner_strip)(Fejer degree 280)
  -> strict-safe certificate after petal/GW labels are retained.

covering 12->84:
  E(row)(covering_packet)(root)(sum)(nested_refinement)(boundary_moment)
  -> positive covering/boundary-moment packet, not scalar all-covered.
```

The pattern is the same in every case: a row is dangerous only if the curried
evaluator has no remaining legal strict-safe return, no AP/GW boundary return,
and no named residual return.

## Computation

Script:

```text
04-computation/lrc14_curried_packet_functional_tower_codex_s170.py
```

Output:

```text
05-knowledge/results/lrc14_curried_packet_functional_tower_codex_s170.out
```

The script records the function signatures, prints the path-rank rows

```text
A(5)(k)=1+3+1,     sum=5
A(6)(k)=1+4+3,     sum=8
A(7)(k)=1+5+6+1,   sum=13
A(8)(k)=1+6+10+4,  sum=21
```

and lists unit-excess Farey partial evaluations `F(p)(14p-1)(lane)` for
`p=1..6`.

## Tournament Analysis

Vertices are curried function carriers / proof obligations, not runners.

Pairwise observable:

```text
predicate retention
fiber labels
argument order
dual certificate availability
named residual
Farey scale
additive carry
anti-scalar guardrail
```

Gauge: majority of strictly larger retention coordinates.  Tie Hamiltonian path:

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

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1,1]
hamiltonian_path_count=1
```

## Proof Target

The proposed proof target is a curried packet theorem:

```text
For every primitive LRC14 packet, the total evaluator E returns one of:
  strict-safe certificate,
  AP/GW boundary equality,
  named F7/THM-572 residual debt.
```

Equivalently: no primitive packet may reach `raw_scalar_evaluation` while still
having an untyped forgotten coordinate.  The final scalar is allowed only after
all earlier partial applications have declared their lost-coordinate functions.

## Next Hooks

1. Add a `curried_call_signature` field to HYP-2963 packet records.
2. For each existing route, emit `lost_Q(c)(fiber)` and one of the allowed
   exits.
3. Build a family template where Fejer certificates are functions of packet
   family, degree, and atom bank rather than per-row floating evaluations.
4. Test whether AP/GW are exactly the packets for which all current curried
   functions close at boundary equality before the strict-safe certificate
   argument is supplied.
5. Test whether K33 is the first product/incidence lane where the product
   partial evaluation must hand to state-lift rather than additive scale.
