# LRC14 Analytic Sieve Packet Weights

The useful merge from the prompt is not "use the large sieve" in the abstract.
It is more surgical:

```text
analytic estimate = kernel + smoothing + exceptional set + retained packet labels
```

Helfgott's ternary Goldbach proof is valuable here because it is explicit about
the proof architecture: circle method, large sieve, exponential sums, explicit
formulas, and smoothing choices.  The repo's LRC work has independently arrived
at the same moral through HYP-2978 and HYP-2981: a quotient is only allowed when
its forgotten coordinates are reconstructed, annihilated, or recorded as
residual debt.

## The Computed Split

The atlas computes two quantities that should never be confused:

```text
Phi(Q) = sum_{q<=Q} phi(q)
G(Q)   = sum_{d<=Q} mu(d)^2 / phi(d)
```

`Phi(Q)` is primitive denominator capacity.  It grows quadratically and is the
right language for "how many unit packets are killed or survive."

`G(Q)` is the inverse-unit normalizer that appears in dimension-one
Selberg/large-sieve thinking.  It grows like `log Q + C`, with the computation
already showing `G(200000)-log(200000) ~= 1.3325`.

So LRC has a new guardrail:

```text
Phi-side statements count primitive packets.
G-side statements normalize upper-bound sieve weights.
Moving between them is a quotient and needs a certificate.
```

## Kaczynski And Kaczorowski

The repo's "Kaczynski" lane is boundary functions and curvilinear convergence:
an approach path and a boundary value are separate data.  The Goldbach-adjacent
thread is Kaczorowski-Perelli-Pintz: exceptional sets in short intervals.

These are not the same theorem, but they give the same LRC instruction:

```text
exceptional set = boundary approach class that must be labelled
```

This fits HYP-2679's far-runner boundary-function view and HYP-2901's lcm wall.
The hard part is not the generic bulk; it is the labelled exceptional approach
where scalar averaging has forgotten which boundary it approached.

## Toward LRC14

The analytic-sieve packet should sit between HYP-2901 and HYP-2981:

```text
fixed denominator basis fails
  -> delete/killed denominator packets
  -> apply G(z)-type inverse-unit normalization with declared smoothing
  -> produce twist/Fejer/Ramanujan witness or route to state-lift debt
```

That is a plausible high-leverage lane because it avoids pretending one finite
denominator ladder proves everything, but it also avoids loose equidistribution
language.  The packet must still carry endpoint owners, route labels, safe
measure, exact-period Ramanujan data, and Fejer interval certificate fields.

The final slogan is simple enough to be dangerous:

```text
large sieve without packet labels is just another scalar quotient.
large sieve with packet labels is a possible LRC14 proof bridge.
```

