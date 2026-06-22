---
id: HYP-2875
status: OPEN sharp proof target
source: codex-2026-06-22-S97b
tags: [LRC14, Beurling-Selberg, spectrum-sum, boundary-functions, decorrelation, complement-even, OPEN-Q-108]
related:
  - HYP-2857
  - HYP-2861
  - HYP-2868
  - HYP-2871
  - HYP-2873
  - THM-537
  - THM-548
  - OPEN-Q-108
---

# HYP-2875: the LRC14 floor should use an admissible bandlimited boundary certificate, not a literal pointwise Beurling-Selberg majorant

S33's Beurling-Selberg prompt is live, but only after the repo's earlier
guardrails are enforced.  The literal pointwise route is blocked: THM-537 and
the Angle-A sessions show that sharp indicator majorants/minorants are either
too weak or impossible in the required signed inclusion-exclusion setting.

The admissible version is the spectrum-sum boundary certificate already
emerging in HYP-2857, HYP-2861, and HYP-2871.

Let

```text
A = G_P,
B = GOOD_E = {x : maxgap(frac(e*x), e in E) > 1/7}.
```

Write

```text
meas(A cap B) = meas(A) meas(B) + SPEC(P,E),
SPEC(P,E) = sum_{n != 0} chat_A(n) conj(ghat_B(n)).
```

The proof target is to choose a fixed cutoff `H` and prove uniformly

```text
sum_{0<|n|<=H} chat_A(n) conj(ghat_B(n))
  - TailBound_H(P,E)
    >= -(1-c) meas(A) meas(B)
```

for some explicit `c>0`, where `TailBound_H` is the better of the L2
Cauchy-Schwarz tail and the Abel/Dirichlet signed tail.  HYP-2871 reports
an exact finite certificate with `R' >= 0.434` at `H=42` on the binding bank;
this hypothesis is the uniform proof program for making that certificate
structural rather than bank-specific.

## Why this is the corrected Beurling-Selberg move

The Beurling-Selberg content is not a pointwise finite-degree surrogate for
the event indicator.  It is a bandlimited signed boundary decomposition:

- keep the finitely many low resonances exactly;
- use signed L2/Dirichlet cancellation for the high tail;
- route coherent low packets to AP/Freiman/GAP finite atlases;
- route incoherent packets to Parseval/energy cancellation.

This matches HYP-2873's energy lesson.  Additive energy is the spectral fourth
moment, and intervals are Fejer-concentrated.  In this route, high additive
energy is not a scalar proof by itself; it is the flag that the low packet has
Freiman structure and should be finite-atlased.

## Boundary-function interpretation

This is the LRC analogue of the S97 admissible-boundary lesson.  The allowed
approach is a fixed low-frequency boundary value plus a controlled signed
tail.  Arbitrary pointwise majorant/minorant attempts are wild curvilinear
approaches: they may preserve some visual shadow of the event, but they do not
preserve the signed quantity that has to be bounded.

So the next proof obligation is not "find a universal Selberg majorant."  It is
"prove the fixed `H` low-mode certificate plus the signed tail bound uniformly,
after the HYP-2871 complement-even reduction."

## Concrete next checks

1. Re-run the HYP-2871 `H=42` certificate on non-bank complement-even stress
   rows and record the exact low-mode deficit distribution.
2. Express the low-mode finite certificate in terms of complement-pair
   autocorrelations and low missed-sector profiles, not raw runner labels.
3. Split low packets by additive energy / Freiman dimension: coherent packets
   get bounded atlases, incoherent packets get the L2 tail.
4. Formalize the tail inequality as the concrete analytic input to the existing
   `LRCWitnessFloorConcrete` / `LRCEventMeasureBridge` route.

