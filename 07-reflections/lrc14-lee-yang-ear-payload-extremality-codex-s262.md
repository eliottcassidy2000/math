# LRC14 Lee-Yang ear-payload extremality

*codex-2026-06-27-S262b.  Owner prompt: work on Lee-Yang extremality, the
whole PGF curve and root structure, phi4-style distribution thinking, and ear
decomposition facts as inspiration for the LRC14 proof and tournament
connections.*

## The scalar was never the curve

HYP-3103 already made the key shift: the LRC inclusion-exclusion is not just a
number.  It is the evaluation at `z=0` of a miss-count PGF

```text
G_E(z) = sum_t q_t z^t.
```

This session adds the missing dynamic piece.  If a new runner is appended, the
root motion of `G_E` is controlled by an exact hidden payload:

```text
A_t(E,a) = P(N_E=t and runner a hits a sector empty for E)
q_{E+a}[t] = q_E[t] - A_t + A_{t+1}.
```

That identity is the LRC version of an observer-extension/cut payload.  It is
also the right place to import the ear-decomposition intuition: adding a runner
is an ear only if the proof keeps the payload that says how the new ear changes
the partition function.

## What the scout saw

The exact S262b scout computes rational miss-count distributions, roots, and
one-runner payloads.

The AP/consecutive and even-AP rows have the same PGF:

```text
q = [481/1470, 359/1470, 25/147, 26/245, 17/210, 5/98, 1/49]
p0 = 0.327211
q6 = 0.020408
p0+q6 = 0.347619
real roots = 0/6
nearest root = 1.4886
dist(roots,[-1,0]) = 0.9119
axis gap = 25.84 deg
```

The spread and broken rows look different:

```text
break_mid:     real roots=2/6, nearest=0.1212, dist([-1,0])=0
random_spread: real roots=4/6, nearest=0.0472, dist([-1,0])=0
```

The far-resonant row keeps the same `#real=0` stratum but sits much closer to
the danger interval:

```text
single_far_21: real roots=0/6, nearest=1.1467, dist([-1,0])=0.2786.
```

So the signal is sharper than "complex roots good, real roots bad."  The
Lee-Yang distance to `[-1,0]` separates AP/consecutive from far-resonant
complex-root packets.

## Ear payloads are the missing sidecar

The final AP/consecutive ear and the single-far ear have different payload
levels:

```text
nested AP/consec final +7:
  A_total=0.362585
  A_mean=1.965291
  A_even/A_odd=0.110544/0.252041
  dist(roots,[-1,0])=0.912

single-far final +21:
  A_total=0.313605
  A_mean=2.993492
  A_even/A_odd=0.130612/0.182993
  dist(roots,[-1,0])=0.279
```

This suggests a concrete proof heuristic:

```text
low-level nested ears stabilize roots away from [-1,0];
high-level far ears leave roots close to [-1,0];
nonnested/broken ears create real-root contacts or named proof debt.
```

That is a useful replacement for a vague Lee-Yang analogy.  It says exactly
what must be measured next: `A_t`, not just the final root signature.

## Ear decompositions as grammar, not imported theorems

The prompt's graph facts are the right shape:

```text
strong directed graph <-> directed ear decomposition
factor-critical graph <-> odd ear decomposition
2-vertex-connected series-parallel graph <-> nested ear decomposition
```

For LRC14 they should be used as controlled-forgetting grammar:

```text
directed ear        -> one-runner extension with retained A_t
odd ear             -> parity split in the payload
nested ear          -> AP/consecutive legal refinement
nonnested ear       -> real-root collision or named sidecar debt
```

HYP-2879 already gives the tournament analogue: strong tournaments have a
one-vertex-ear calculus.  HYP-3112 gives the LRC analogue: the one-runner ear
has an exact polynomial payload, and this payload reconstructs the next PGF.

## The phi4 / Lee-Yang lesson

The density `exp(-lambda S^4 - b S^2) dS` is useful here as a warning, not as
a proof import.  It says the object of interest is an entire distribution and
its transform/root structure.  In the LRC setting, `p0` is just one evaluation
of `G_E`; moments are projections; root movement and zero-free regions are
global data.

The proof target should therefore be ledger-shaped:

```text
row -> G_E -> roots -> legal ear A_t -> root motion -> terminal exit/debt.
```

That ledger can talk to HYP-3104's maximizer-signal atlas, HYP-3106's
perspective functors, HYP-3101's component packets, HYP-3102's first
obstructions, and HYP-3105's obstruction-transfer guardrails.

## Best next theorem candidate

The most promising exact statement is not yet a cap inequality.  It is a
separation theorem:

```text
If an LRC14 packet has roots approaching or meeting [-1,0], then it has
one of the known debts:
  high-mean far ear payload,
  nonnested ear,
  component-bound debt,
  first-obstruction syndrome,
  K33/THM-572 state-lift debt,
  or AP/GW boundary status.
```

The contrapositive would make AP/consecutive extremality a zero-free-region
statement for packets whose sidecars are legally discharged.

Related: HYP-3112, HYP-3111, HYP-3109, HYP-3108, HYP-3107, HYP-3106,
HYP-3105, HYP-3104, HYP-3103, HYP-3102, HYP-3101, HYP-3085, HYP-2879,
THM-577, THM-576, THM-573, OPEN-Q-108.
