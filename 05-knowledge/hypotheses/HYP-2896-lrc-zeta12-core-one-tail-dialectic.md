---
id: HYP-2896
status: PROOF FRAGMENT / one-tail disproof branch closed
source: codex-2026-06-22-S109
tags: [lrc, lrc14, zeta, bernoulli, q-witness, covering, binding-pair, resonance-killing, tournament-analysis]
related:
  - HYP-2895
  - HYP-2893
  - HYP-2890
  - HYP-2856
  - THM-523
  - THM-524
  - THM-560
  - THM-563
  - THM-566
  - HYP-+2878
  - OPEN-Q-108
results:
  - 04-computation/lrc_zeta12_core_dialectic_codex_s109.py
  - 05-knowledge/results/lrc_zeta12_core_dialectic_codex_s109.out
  - 04-computation/lrc_resonance_hierarchy_kps.py
  - 04-computation/lrc_disproof_killer_attack_kps.py
  - 04-computation/lrc_sharpest_disproof_kps.py
  - 07-reflections/the-resonance-killing-game-and-the-zeta-duality-of-the-lonely-runner.md
---

# HYP-2896: the zeta -1/12 core closes the one-tail disproof branch

The owner asked to hold together the exact LRC fact

```text
M({1,2,...,11,13}) = 1/12
```

and the regularized identity

```text
1 + 2 + 3 + ... = zeta(-1) = -1/12.
```

Incoming KPS S31p gives the broad resonance-killing narrative and the sharper
killing-budget branch.  The AP kills the initial resonance ladder by using the
least-spread divisibility killers `1,2,...,n`; a counterexample attempt tries
to kill more small resonances with highly divisible speeds, but that spreading
opens larger lonely gaps and pushes `M` upward.  In the sibling row
`{1,...,12,v}`, one speed can kill `13` or `14` cheaply but killing both forces
`v` divisible by `182`, and the exact scan has minimum `1/14` at `v=13`.  S109
adds the complementary precise theorem at the first nontrivial missing-killer
core.

Let

```text
C = {1,2,...,11,13}.
```

Then `M(C)=1/12`, witnessed at `t=5/12` with active speeds `(5,7)`.  For every
single added positive integer speed `w`, the row `C union {w}` is safe for
LRC14.  More precisely, exactly one of the following holds:

```text
12 does not divide w:
  the q=12 witness survives, so M(C union {w}) >= 1/12.

12 divides w and 14 does not divide w:
  the q=14 witness survives, so M(C union {w}) >= 1/14.

84 divides w:
  write w=84m.  Then the covering row has the exact witness
    t_m = (35m+2)/(84m+5)
  and the exact value
    M(C union {84m}) = 7m/(84m+5) > 1/14.
```

The covering branch is the useful new fragment.  At `D=84m+5`, the numerator
distances at `t_m` are

```text
speed  1: 35m+2
speed  2: 14m+1
speed  3: 21m+1
speed  4: 28m+2
speed  5: 7m
speed  6: 42m+2
speed  7: 7m+1
speed  8: 28m+1
speed  9: 21m+2
speed 10: 14m
speed 11: 35m+3
speed 13: 35m+1
speed 84m: 7m
```

so every distance is at least `7m/D`, equality occurs at the binding pair
`(5,84m)`, and

```text
7m/(84m+5) > 1/14  iff  14m > 5.
```

This is the exact finite version of the zeta hint.  The formal negative
constant `-1/12` is a Bernoulli boundary term in

```text
sum_{n>=1} n exp(-epsilon n)
  = epsilon^-2 - 1/12 + epsilon^2/240 - epsilon^4/6048 + ...
```

but a finite row cannot contain a negative runner.  It has to realize
regularized cancellation through divisibility, q-witnesses, and binding pairs.
For the first missing-killer core `C`, that finite realization is completely
self-defeating: adding one integer speed either leaves a q-witness alive or,
when it kills both `12` and `14`, becomes the positive covering family
`84m` with an explicit margin.

## Disproof lesson

The sharp disproof attempt says:

```text
kill the 1/12 witness;
then kill the 1/14 witness;
hope the first survivor is below 1/14.
```

For one-tail rows this attempt is exhausted.  Killing `12` forces
`12 | w`; killing `14` as well forces `84 | w`; and that very divisibility
creates the binding-pair formula above.  Thus the one-tail counterexample
branch is closed.

This also corrects an overly simple reading of the resonance hierarchy.  In
non-covering rows, the largest witness is often well described by the first
surviving small denominator.  In the covering one-tail row `C union {84m}`, the
relevant value is not just `1/b` for the next small denominator; it is the
binding-pair value `7m/(84m+5)` at a denominator depending affinely on the
large speed.  The proof object is therefore a q-witness/binding-pair ledger,
not a scalar killed-denominator list.

## Proof lesson

S109 strengthens HYP-2895's decomposition:

```text
non-covering rows:
  closed by THM-523 q-witnesses.

AP budget sibling:
  {1,...,12,v} is minimized at v=13; killing 13 and 14 together forces
  the large lcm killer 182 and moves back toward a 1/13 scale.

one-tail covering over the 1/12 core:
  closed by the exact 84m binding-pair formula.

remaining battlefield:
  multi-large / moderate resonant covering rows, where several tails can
  interact before equidistribution fully takes over.
```

So the next LRC14 proof obligation should not be another search for a single
finite denominator basis or a single extra speed over `{1..11,13}`.  It should
target the multi-tail resonance middle: exact-period packets, support-six
residual leak (HYP-2890), HYP-+2878 atom over-determination, and the
Clebsch/Bruhat/octahedral signed-carrier route.

## Tournament Analysis

Vertices are proof/disproof mechanisms:

```text
formal_zeta_negative_constant
finite_divisibility_realization
q12_witness
q14_witness
covering_tail_84m
binding_pair_formula
multi_large_resonant_middle
```

Hamiltonian path:

```text
formal_zeta_negative_constant
  > finite_divisibility_realization
  > q12_witness
  > q14_witness
  > covering_tail_84m
  > binding_pair_formula
  > multi_large_resonant_middle
```

Observable: whether a formal cancellation survives as a finite integer
counterexample.  Switch/gauge: the q-covering predicate.  Tie path: the
proof/disproof pressure order above.

Preserved predicate: whether the row is still LRC14-safe by a q-witness or by a
binding-pair formula.  Destroyed information: the geometry of arbitrary
multi-tail resonance interactions.  Challenged assumption: the zeta
regularized `-1/12` can be realized by one finite positive speed; the finite
realization turns into divisibility and positive margin instead.
