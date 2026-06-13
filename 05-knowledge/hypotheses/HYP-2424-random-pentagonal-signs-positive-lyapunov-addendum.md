# HYP-2424 - Random pentagonal sign families have positive reciprocal Lyapunov addenda

**Status:** OPEN addendum. Euler side is theorem-level; THM-485/HYP-2416-2417
already give the stronger paired-sign Lyapunov/sign-rigidity mainline. This
file records the reciprocal-denominator zero-locus formulation, independent
term-sign scout, and deterministic zero-Lyapunov classification target.
**Source:** codex-2026-06-11-P1.
**Companions:** THM-485, HYP-2416, HYP-2417, OPEN-Q-064, T783.
**Artifacts:** `04-computation/pentagonal_lyapunov_code72_codex.py`,
`05-knowledge/results/pentagonal_lyapunov_code72_codex.out`.

## Statement

Let

```text
D_eps(q) = 1 + sum_{k>=1} eps_{k,-} q^{k(3k-1)/2}
             + eps_{k,+} q^{k(3k+1)/2}
```

with signs either paired by `k` or independent by pentagonal exponent, and let
`A_eps(q)=1/D_eps(q)=sum a_eps(n) q^n`. Euler's paired signs
`eps_{k,-}=eps_{k,+}=(-1)^k` are a zero ordinary Lyapunov point:

```text
lim_n log p(n)/n = 0,
log p(n) ~ pi sqrt(2n/3).
```

For random signs, conjecturally

```text
lambda(eps) = limsup_n log |a_eps(n)| / n > 0
```

almost surely in the reciprocal-denominator model. Equivalently, `D_eps` almost
surely has at least one zero in the open unit disk, and `lambda(eps)` is
controlled by the nearest such zero. THM-485 attacks the paired-sign recurrence
temperature directly; this addendum asks for the analytic zero theorem and the
larger independent-term sign law.

## Evidence

The script computes exact reciprocal coefficients through `n=650`, using a
model complementary to the THM-485 transfer-operator computation.

- Euler: finite regression still sees the Hardy-Ramanujan sqrt scale; this is not
  a limiting Lyapunov estimate.
- Paired-random signs: all `160/160` samples had positive finite-window slopes,
  with median about `0.244645`.
- Independent term-random signs: all `160/160` samples had positive finite-window
  slopes, with median about `0.261340`.
- Single early paired flips away from Euler create visible exponential growth;
  late flips beyond the sampled window become indistinguishable from Euler at
  this cutoff.

The all-plus deterministic control has low finite-window slope, while THM-485
separately identifies a large fixed-sign rate. Thus the claim is not "every
non-Euler law is obviously exponential." The honest target is a random
analytic-function/interior-zero theorem and a deterministic zero-Lyapunov
classification, not a regression theorem.

## Proof Route

Prove that a lacunary random series on generalized pentagonal exponents,
normalized with constant term `1`, has an interior zero almost surely. Possible
tools: Jensen formula for random analytic functions, small-ball estimates on
circles `|q|=r`, or a two-radius Rouché argument using the first few random
pentagonal terms against the deterministic tail.

The Euler exception is explained by product structure:

```text
D_Euler(q) = prod_{m>=1} (1-q^m),
```

so its zeros sit on roots of unity, not inside the disk. The cancellation is a
factorization theorem, not a generic sparse-sign phenomenon.

## Mistake Guards

- Finite linear slopes do not estimate a zero limiting exponent for partition
  numbers; they still see `sqrt(n)` growth.
- A positive finite-window slope is evidence for an interior zero, not a proof.
- The right object is the denominator's zero set. Coefficient recurrences are
  diagnostics.
