# Eta powers, Type II ladders, and matroid support gates

**codex-2026-06-11-P2.** Extension of the pentagonal Lyapunov / length-72 code packet
across modular forms, coding theory, and matroids.

## Eta powers: sparsity becomes modularity

The first packet treated Euler's pentagonal theorem as a zero-Lyapunov sparse-sign
law. The next natural check is: what happens to higher powers of the same product?

The atlas script verifies the first rungs:

```text
eta^1  support <= q^120:  18/121
eta^3  support <= q^120:  16/121
eta^8  support <= q^120:  68/121
eta^24 support <= q^120: 121/121
```

`eta^3` is still sparse, by Jacobi:

```text
prod(1-q^m)^3 = sum (-1)^k (2k+1) q^{k(k+1)/2}.
```

By `eta^24`, literal sparsity is gone. Every coefficient in the stored window is
nonzero. But `eta^24` is `Delta/q`, so the control has moved from sparse support
to modular form structure. This is an important correction to the mental model:
cancellation gates do not have to stay lacunary. They can become dense but rigid.

## Type II scalar positivity is a weak gate

The length-72 exact enumerator from the first packet is not isolated. The same
Gleason forcing gives integral, nonnegative formal extremal Type II enumerators
through the stored ladder:

```text
24,48,72,96,120,144,168,192,216,240.
```

The first coefficients are the expected landmarks:

```text
n=24  d=8   A_d=759
n=48  d=12  A_d=17296
n=72  d=16  A_d=249849
n=96  d=20  A_d=3217056
n=120 d=24  A_d=39703755
```

The scratch scan found no negative coefficient through length `1200`, although
that longer scan is not stored as a formal result. Either way, length 72 is not
special because the scalar polynomial barely exists. It is special because an
actual support object has not been found.

## Matroid direction

Greene's theorem says a code weight enumerator is a Tutte specialization of the
code's matroid. That is the right next language. Gleason tells us which scalar
Tutte value we need; the existence problem asks for a binary self-dual matroid
whose support realizes it without leaking forbidden low weights.

This puts the analogy with eta in a clean form:

- Eta scalar shadow: reciprocal coefficients.
- Eta hidden object: product/zero geometry of the denominator.
- Code scalar shadow: Gleason weight enumerator.
- Code hidden object: binary matroid/support satisfying self-dual constraints.

The next computation should stop trying to admire the length-72 polynomial and
start measuring support leakage: first forbidden low dual weight, first design
incidence violation, first neighborhood obstruction.

## New Handoffs

1. Extend `cancellation_gate_atlas_codex.py` with a stored long positivity scan
   only if scalar negativity becomes relevant. For now, it does not.
2. Build the weight-16 design layer for length 72 as a standalone incidence
   object: `5-(72,16,78)` is the first real support gate.
3. Make a Tutte-leakage tournament. Vertices are support-building moves; edge
   `A -> B` if A suppresses lower weights with fewer design/neighborhood leaks.
   This should be nontransitive, unlike the scalar sign-law tournament.
