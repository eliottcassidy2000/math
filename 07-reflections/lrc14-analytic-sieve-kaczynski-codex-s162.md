# LRC14 Analytic Sieve / Kaczynski Synthesis - Codex S162

This session merged the user's analytic-number-theory prompt into the current
LRC14 proof program.  The useful import is not "use primes."  It is the proof
architecture around primes:

```text
main packet,
minor-arc exponential-sum control,
explicit smoothing and boundary errors,
finite exceptional verification.
```

The LRC version of the main packet is labelled:

```text
qdiv,
exact M/Farey branch,
Ramanujan/totient exact-period packet,
Haar/Baire boundary status,
endpoint owners,
C27/K33/state-lift labels,
dual-certificate fields.
```

## What The Prime-Sum Motif Should Mean Here

The repo already learned that the LRC singular-series analogy is archimedean,
not an Euler product.  That guardrail stays.  Prime sums and Goldbach enter as
discipline, not as literal variables.

In Helfgott/Vinogradov, the proof is not one global estimate.  It is:

```text
major arcs,
minor arcs,
careful smoothing,
explicit constants,
finite exceptional checks.
```

For LRC14 this becomes:

```text
q/Farey/Ramanujan major packets,
off-resonance Fourier modes,
smoothing kernels for danger arcs,
exact endpoint-owner boundary defects,
finite AP/GW/K33/petal/covering packet atlas.
```

The missing theorem should have that shape.

## The Mobius/Totient Lesson

S162 computed:

```text
A(N) = sum_{n<=N} mu(n)/n
B(N) = sum_{n<=N} mu(n)^2/phi(n)
C(N) = sum_{n<=N} mu(n)^2/n
D(N) = average phi(n)/n.
```

The important reading is:

```text
sum mu/n       -> cancellation diagnostic
sum mu^2/phi   -> primitive squarefree packet capacity
sum mu^2/n     -> squarefree mass
avg phi/n      -> coprime-density floor
```

The second line is the most LRC-native.  It says exact-period primitive
capacity should be attached to Ramanujan/Farey packets before scalarization.
This connects directly to HYP-2978/HYP-2979.

## Kaczynski Is The Boundary Word We Needed

THM-548 and HYP-2679 already made the true-wide branch look like a boundary
function problem.  S162 makes that language more operational:

```text
Fatou-style ordinary approach      -> decorrelated discharge
Kaczynski/Bagemihl ambiguous arc   -> resonance exception
resonance exception                -> finite packet atlas
smoothing disagreement             -> boundary defect term
```

AP and Goddyn-Wong are then not just tight examples.  They are boundary atoms:
zero strict-open mass, zero-credit endpoint owner pairs, and no strict
primitive exact-period witness in the bounded Ramanujan audit.

The next proof should try to show that no other packet can be ambiguous in all
these ways at once.

## Why Large Sieve Must Be Labelled

A large-sieve or upper-bound sieve estimate can be useful only after the
predicate-carrying packet has been named.  Otherwise it repeats the old
failure mode:

```text
bound a scalar,
forget endpoint owners,
lose AP/GW boundary debt,
misclassify K33 or covering rows.
```

The safe use is:

```text
large-sieve bound on a labelled packet family
+ explicit residual term for the labels it cannot see.
```

The same applies to smoothing.  Different smoothing functions are not
cosmetic.  If a packet's classification changes when the kernel changes, that
difference is a boundary defect and belongs in the proof ledger.

## Tournament Readout

The S162 tournament uses proof modules as vertices.  The pairwise observable is
retention of:

```text
LRC predicate,
exact arithmetic packet,
local/global split,
exponential-sum control,
boundary/resonance exceptions,
smoothing adaptability,
auditability.
```

The nontrivial SCC is:

```text
mobius_totient_packet_ledger
smoothing_saddle_explicit_formula
circle_large_sieve_decomposition
```

That is a real signal.  The analytic middle layer is coupled.  We should not
try to first prove a raw Mobius/totient statement, then separately add a circle
method split, then separately choose smoothing.  These choices interact before
the LRC predicate is visible.

## Candidate Next Lemma

The most useful theorem shape now is:

```text
For every primitive LRC14 source packet, a chosen smoothing family gives an
explicit finite Fourier/Ramanujan expansion.  Either the off-resonance modes
are bounded strongly enough to produce a strict safe interval, or the resonant
modes force the packet into AP/GW equality, unit-petal/two-block discharge,
K33/state-lift debt, or a previously classified covering/lift packet.
```

The proof burden is to make "off-resonance" exact and to write the smoothing
boundary error as a named Kaczynski defect, not as an untracked epsilon.
