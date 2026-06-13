# HYP-2439 - Midpoint Faulhaber and tournament reversal are parity-dual centered cancellation gates

**Status:** COMPUTATION CONFIRMED on the stored scout; theorem-level bridge OPEN.
**Source:** codex-2026-06-12-P5.
**Companions:** HYP-2128, HYP-2426, THM-091, THM-092.
**Artifacts:** `04-computation/parity_dual_cancellation_gate_codex.py`,
`05-knowledge/results/parity_dual_cancellation_gate_codex.out`.

## Statement

Two repo mechanisms share one centered-cancellation shape, but with opposite
parity:

1. **Power-anchor / midpoint gate.** In

   ```text
   F_p(c,n) = c^p + sum_{j=1}^n ((c-j)^p - (c+j)^p),
   ```

   the midpoint antisymmetrization kills every even `j`-power, so only odd
   Faulhaber moments `S_1,S_3,S_5,...` survive. The root `a_p(n)` is therefore
   controlled by an odd-only Bernoulli/Faulhaber channel.

2. **Tournament / reversal gate.** For forward-edge count,

   ```text
   fwd(sigma) + fwd(sigma^rev) = n-1,
   ```

   so the centered distribution of `fwd-(n-1)/2` kills every odd centered
   moment/cumulant. The tournament side is therefore controlled by an even-only
   cumulant hierarchy.

Conjectural bridge: these are **parity-dual centered cancellation gates**. The
same “center first, then see what survives” operation appears in both settings,
with the midpoint scalar problem retaining odd channels and the reversal
distribution problem retaining even channels.

## Evidence

The stored scout checks:

- for `p=1..10`, the midpoint difference `(c-j)^p-(c+j)^p` has support only on
  odd powers of `j`;
- direct midpoint balances match the odd-moment expansion exactly on sample
  instances;
- the anchor approximation

  ```text
  a_p(n) = p n^2 + (p-1)n + alpha_p + beta_p/(n(n+1)) + O(n^-4)
  ```

  with

  ```text
  alpha_p = (p-1)(p-2)/(12p),
  beta_p = -(p-1)(p-2)(2p^2-4p-1)/(180 p^3)
  ```

  leaves a residual that stabilizes at `u^-2`, `u=n(n+1)`;
- exhaustive `n=5` tournament scan (`1024` tournaments) finds zero failures of
  vanishing odd centered forward moments up to order `5`, while even centered
  moments vary over multiple values.

## Use

This is not an OCF identity. It is a proof-routing hypothesis:

- when a scalar/tower/LRC quantity is naturally centered at a midpoint, rewrite
  it before expanding;
- when a tournament observable is naturally paired by reversal, center it before
  tracking moments;
- keep the surviving parity channel, not the raw moment list, as the stable
  invariant.

That lens may sharpen:

- HYP-2128's triangular carrier;
- HYP-2426's cancellation-gate program;
- midpoint-style scalar excision routes in LRC, where the useful object may be
  the centered residual rather than the raw endpoint sum.
