# Lonely Runner Hypothesis Noise S385

**Session:** `codex-2026-05-31-S385`
**Script:** `04-computation/lonely_runner_hypothesis_noise_s385.py`
**Stored output:** `05-knowledge/results/lonely_runner_hypothesis_noise_s385.out`

This was a deliberately noisy exploration pass.  I fed recent repo threads
back into the Lonely Runner frontier: endpoint transfer, Fejer/Riesz pressure,
quotient transport, natural operation shadows, product-sum gates, and the
fourteen-runner composite-denominator anomaly.

The useful outcome was not one new candidate counterexample.  It was three
more precise models for why the existing near-counterexamples keep failing.

## Endpoint Matroid

S385 treats the endpoint-protection relation as a binary matrix over `GF(2)`.
Rows are forbidden endpoints; columns are pulled-back forbidden intervals; a
row records which intervals strictly protect that endpoint.

The striking table entry is `coreE=0` across all sampled examples:

```text
initial n=14:       coreE=0, privateE=128
n14 seven-ladder:  coreE=0, privateE=710
n14 14-ladder:     coreE=0, privateE=1420
n15 3x5 ladder:    coreE=0, privateE=204
n15 mixed gates:   coreE=0, privateE=232
```

Even the high-rank-looking systems have private pivots and peel away.  This is
now HYP-1842: the proof may be a private-pivot matroid theorem for
integer-realizable endpoint systems.

## Gap Width Is the Wrong Danger Metric

The tiny-gap examples remain bad candidates when viewed through the max-min
landscape.  The best composite `n=14` debt ladder has:

```text
gap/th = 0.002706
gap-probe critical/th = 1.208333
unprotected endpoints = 168
coreE = 0
```

The visible open gap is tiny, but the safe valley is steep.  A true disproof
candidate should have critical-radius surplus close to zero and a protected
endpoint core.  This became HYP-1843.

## Debt Export

S380 showed that the `14`-multiple ladder improves the seven-ladder gap ratio.
S385 explains the cost:

```text
seven-ladder:      unprotected=84,  low-layer debt=36, privateE=710
14-ladder debt:    unprotected=168, low-layer debt=72, privateE=1420
```

The composite gates protect the designed low layers, but the endpoint
obligation reappears on descendant denominators.  The one-swap provocations
around both ladder families preserve the same positive-gap/unprotected pattern.
This became HYP-1844.

## Reframes Worth Keeping

1. **Endpoint-incidence matroid:** a counterexample is a leafless, private-
   pivot-free endpoint-protection matrix that is also integer-realizable.
2. **Morse landscape:** tiny threshold gaps can be steep valleys; use critical
   radius and endpoint core size as the main danger coordinates.
3. **Debt export:** divisor gates move leakage through quotient layers rather
   than eliminating it.
4. **Two operation shadows:** addition supplies the scalar Dirichlet spine;
   multiplication supplies divisor gates and endpoint-protection channels.
5. **Abstract-arc mirage:** all-protected circular-arc systems might exist
   topologically but fail integer-speed realizability.

## Next Experiments

1. Search abstract circular-arc endpoint systems for nonzero protection cores,
   then impose integer-speed realizability.
2. Add exact or bounded critical-radius estimates to the S373-S385 candidate
   database and rerank candidates by `(critical surplus, coreE, debt)`.
3. Build a denominator-layer transfer matrix for the `14`-multiple ladder:
   protected low-layer targets in, descendant unprotected endpoints out.
4. Repeat the debt-export ledger for `n=15` order-3 CRT families.
