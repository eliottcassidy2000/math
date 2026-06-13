# Parity-dual centered cancellation gates

**codex-2026-06-12-P5.** Prompt: spend a longer exploratory session in the math
repo, search for nearby places where the power-anchor / odd-Faulhaber idea could
 make real progress, refine the hypotheses, and push the result.

## Where I looked

Three routes were live in the repo scan.

1. **Triangular/LRC carrier route.** `HYP-2128` already had the midpoint anchor
   asymptotics and the odd-square/triangular bridge. The risk here was writing
   only another inspirational synthesis without changing the proof surface.
2. **Tournament moment route.** THM-091/092 already isolate reversal symmetry:
   all odd centered forward-edge cumulants vanish, and the even cumulants form a
   cycle hierarchy.
3. **Generic cancellation-gate route.** `HYP-2426` was broad enough to absorb
   the new anchor identity, but too broad to sharpen unless the new mechanism
   could be stated in the same exact “what survives after centering?” language.

The productive move was to combine 1 and 2 first, then feed the result back
into 3.

## What the scout found

The midpoint power-anchor equation and the tournament reversal identity are the
same operation with opposite parity.

On the anchor side,

```text
F_p(c,n) = c^p + sum_{j=1}^n ((c-j)^p - (c+j)^p)
```

is a center-antisymmetrization. Expanding in `j`, every even `j`-power cancels,
so only odd Faulhaber moments remain:

```text
F_p(c,n) = c^p - 2 sum_{r odd} binom(p,r) c^(p-r) S_r(n).
```

On the tournament side, reversal gives

```text
fwd(sigma) + fwd(sigma^rev) = n-1.
```

So after centering at `(n-1)/2`, every odd centered moment vanishes. The
surviving hierarchy is even.

That is the bridge:

```text
midpoint scalar gate    -> odd channels survive
reversal tournament gate -> even channels survive
```

The stored scout makes this finite and explicit:

- for `p=1..10`, midpoint support lands only on odd powers of `j`;
- direct midpoint balances match the odd-moment expansion exactly on test cases;
- the two-term anchor asymptotic leaves a clean `u^-2`, `u=n(n+1)`, residual;
- exhaustive `n=5` tournaments (`1024`) have zero failures of vanishing odd
  centered forward moments up to order `5`, while even centered moments vary.

## Why this is better than a loose analogy

This is not saying “OCF secretly equals Faulhaber.” That would be too coarse and
likely false. The real common structure is narrower:

- center the object at its natural involution;
- expand only after centering;
- the involution kills one parity sector completely;
- the surviving parity sector carries the correction hierarchy.

That is a reusable proof instruction. It says the stable object is not the raw
sum, raw endpoint formula, or raw moment sequence. It is the **surviving parity
channel after centering**.

## Best nearby uses

1. **HYP-2128 / triangular carriers.** The anchor family should now be read as
   an odd-channel carrier, not just a triangular packing story. That gives a
   cleaner explanation for why `p=1,2` are exact and `p>=3` deform only through
   Bernoulli corrections.
2. **HYP-2426 / cancellation gates.** Eta, Gleason, OCF, and now midpoint
   Faulhaber are all “kill a dangerous parity/layer by centering or symmetry”
   mechanisms. The new one is finite and scalar, which makes it a good toy model
   for the larger gate program.
3. **LRC midpoint/scalar-excision routes.** The scalar-ramp midpoint identity
   work around S371 is already phrased in centered variables. The next honest
   question is whether any LRC residual is controlled by an odd-only centered
   Bernoulli channel rather than by raw endpoint data.

## What I am not claiming

- no new theorem in canon;
- no direct proof of an LRC statement;
- no exact OCF/Faulhaber equivalence.

The refined claim is methodological but exact: the midpoint anchor and reversal
moment constructions are parity-dual centered cancellation gates.

That is now recorded as `HYP-2439`, with script/result pair:

- `04-computation/parity_dual_cancellation_gate_codex.py`
- `05-knowledge/results/parity_dual_cancellation_gate_codex.out`
