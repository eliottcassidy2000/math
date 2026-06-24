---
id: HYP-2938
title: LRC14 pi approximants, flower-family aliasing, and unital semantics
status: PROOF-GUARDRAIL / normalization and terminology atlas; not a proof
source: codex-2026-06-23-S137
related:
  - HYP-2937
  - HYP-2936
  - HYP-2934
  - HYP-2932
  - HYP-2931
  - HYP-2920
  - HYP-2894
  - HYP-2892
  - HYP-2899
  - HYP-2900
  - OPEN-Q-108
results:
  - 04-computation/lrc14_pi_unital_flower_alias_codex_s137.py
  - 05-knowledge/results/lrc14_pi_unital_flower_alias_codex_s137.out
---

# HYP-2938: pi approximants, flower-family aliasing, and unital semantics

S137 adds a normalization guardrail to the S130-S136 LRC14 carrier chain.  The
prompt supplies three useful signals:

```text
pi ~= 22/7,
pi ~= cuberoot(31),
flower petals show 22 families when rotated by 1/pi,
```

and asks to fold in the several meanings of "unital."  The conclusion is that
all four signals are useful only if their carrier type stays labelled.

## Numerical facts

The diagnostic script
`04-computation/lrc14_pi_unital_flower_alias_codex_s137.py` stores output at
`05-knowledge/results/lrc14_pi_unital_flower_alias_codex_s137.out`.

It records:

```text
22/7          = 3.142857142857143, error +1.264489267349678e-3.
cuberoot(31) = 3.141380652391393, error -2.120011984003689e-4.
```

Thus `cuberoot(31)` is about `5.96` times better than `22/7` by absolute
error, with

```text
pi^3 - 31 = 6.276680299816206e-3.
```

The proof-role split is important:

- `22/7` is a rational/Farey approximant.  It gives denominator and numerator
  labels, especially `1/pi ~= 7/22`.
- `cuberoot(31)` is a cubic algebraic approximant.  It gives a `31` scalar
  or volume tag, not a denominator-31 Farey branch.

So `31` may be a useful echo for the existing Q31/fiber language, but only as
a cubic carrier unless a separate rational-denominator statement is proved.

## Flower normalization

The phrase "turn `1/pi` radians" has two different mathematical readings.

Literal radians:

```text
theta = 1/pi radians.
turn_fraction = 1/(2*pi^2) ~= 0.050660591821169.
steps per turn = 2*pi^2 ~= 19.739208802178716.
best q<=100 closure is q=79, four turns, residual 0.01373977980112073 rad.
q=22 residual is 0.7196321888638089 rad.
```

Therefore literal `1/pi` radians does not produce an exact or even especially
good 22-period flower.

Full-turn normalization:

```text
step = 1/pi of a full turn = 2 radians.
22/7 gives 1/pi ~= 7/22.
22 steps ~= 7 full turns, residual 0.01770284974289638 rad.
```

This is the clean source of the "22 families" observation.  The flower count is
an aliasing witness for the full-turn normalization and the rational
approximant `7/22`, not an invariant of the literal-radian rotation.

## Unital semantics

The prompt's unital paragraph should be attached as a terminology guardrail.

Block-design unital:

```text
v = q^3 + 1, block size = q + 1, lambda = 1.
For q=3: v=28, block size=4, blocks=63, point replication=9.
```

This is exactly the HYP-2892/HYP-2894 pair-frame numerology:
`28 = C(8,2)`.  It remains a secondary weighted frame after a category-1
AP/Goddyn-Wong labelling is chosen; HYP-2894 already refutes the naive
canonical `S8`-invariant design identification.

Algebraic/functional-analytic unital:

```text
an object has an identity, or a map preserves identity: f(1_A)=1_B.
```

This belongs in the repo as a quotient-map metaphor: source/observer,
identity, apex, and denominator labels should be preserved by any proof map
that claims to be conservative.  It is not the same object as the block design.

Residue unit shell:

```text
gcd-class visibility such as g1 versus g3/g9 in the C=27 shell quotient.
```

This is the primary finite label in HYP-2937.  "Unit-visible shell hole" means
residue-unit visibility, not a unital algebra or a unital design block.

Stable unitality is a still looser analogy: persistence of an identity after
stabilization resembles a carrier surviving lift/dilation, but no theorem
should use that phrase unless the stabilizing operation and preserved identity
are named.

## Tournament Analysis

Vertices:

```text
angle_normalization,
exact_M_Farey_branch,
C27_unit_nonunit_shell_transfer,
pi_approximant_labels,
q3_unital_pair_frame,
algebraic_unital_identity_metaphor,
raw_flower_family_count.
```

Pairwise observable:

```text
normalization exactness,
LRC predicate faithfulness,
label retention,
finite certifiability,
scalarization risk.
```

The role order is transitive:

```text
angle normalization
> exact M/Farey branch
> C=27 unit/nonunit shell transfer
> 22/7 and cuberoot(31) pi-approximant labels
> q=3 unital pair-frame after category-1 labelling
> algebraic unital identity-preservation metaphor
> raw flower/family visual count.
```

Assumption challenged:

```text
A visible 22-family flower is not an exact literal-radian invariant.
The useful 22 comes from the full-turn normalization and 7/22 ~= 1/pi.
```

Second assumption challenged:

```text
"unit", "unital design", and "unital identity map" are three different
carriers.  They can be compared only after their preserved labels are named.
```

## Proof Target

Any future visual/flower carrier should first be normalized into one of:

```text
(A) rational-turn approximant p/q with a q-family alias and residual bound,
(B) literal-radian irrational rotation with closure denominators computed
    from theta/(2*pi),
(C) algebraic pi approximant such as cuberoot(31), stored as a cubic scalar tag.
```

Only case (A) should feed petal-family counts.  Case (B) belongs to
irrational-rotation/continued-fraction diagnostics.  Case (C) can be compared
with Q31-style fiber language, but not promoted to a denominator or Farey
branch without a separate rational statement.

For the LRC14 chain, the practical rule is:

```text
exact M/q/Farey labels first;
C=27 unit/nonunit shell transfer second;
pi-approximant and flower-family labels only after angle normalization;
unital block designs only after AP/GW category-1 labelling;
algebraic unitality only as an identity-preservation guardrail for maps.
```

This does not prove LRC14.  It prevents a promising set of numerological and
terminological carriers from being used at the wrong layer of the proof.
