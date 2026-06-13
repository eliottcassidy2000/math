---
source: codex-2026-06-03-S601
status: conceptual synthesis plus deterministic Helly entropy audit
tags: [LRC, Helly, iterated-log, determinant, certificate-entropy, CRT]
---

# Helly Scale And Logarithms

The useful move this session was to stop asking, vaguely, "where does the log
come from?" and instead ask:

```text
What is the search space whose logarithm the proof is paying?
```

For a Helly proof, the search space is not automatically denominator height.  It
can be the family of possible local certificates.  If there are `M` live
component languages and the proof only needs a witness of rank at most `H`, the
natural currency is

```text
Lambda_H(M) = log sum_{h<=H} binom(M,h).
```

For bounded `H`, this is `H log M + O(1)`.  That little formula is the whole
point.  If the live component count is already compressed to `M ~ log N`, the
Helly search pays `loglog N`.  If it is compressed to `M ~ loglog N`, it pays
`logloglog N`.  The iterated log is then a receipt for local certificate
entropy, not for the global CRT modulus.

## What The S601 Audit Says

I reran the S599 two-block determinant classifier through this lens.  Across
the deterministic sample:

```text
certified_empty = 1113
high_order_empty = 0
bounded_live = 0
h=1 certificates = 1084
h=2 certificates = 29
h>=3 certificates = 0
```

So the current determinant branch is not paying a large Helly tax.  It is
earning a Helly dividend: singleton walls and pair incompatibilities appear
before a global CRT intersection is needed.

That matters for n=14.  In the `BC_only` regime the fresh sample saw

```text
n=14: h1=25, h2=0, high=0, live=0.
```

In the full stack, n=14 was entirely preempted.  This is exactly the shape we
want if HYP-2144 is going to become a human proof: local determinant rows, not
large residue products.

## The Helly Log Template

Here is the template I would keep:

```text
If every live residual with M component languages has an empty subfamily
of size at most H, then the local certificate entropy is Lambda_H(M).
```

Then:

```text
M ~ N^a        gives ordinary log N
M ~ log N      gives loglog N
M ~ loglog N   gives logloglog N
```

And if the proof removes residuals by ranks harmonically,

```text
R_H <= R_1 exp(-alpha sum_{h<=H} 1/h) <= R_1 H^{-alpha},
```

so a logarithmic rank cutoff `H` becomes another source of loglog behavior.

## Relation To Tao And Collatz

The incoming Collatz tower work reads logs as successive floors of orbit
linearization.  This Helly session reads logs as successive floors of proof
certificate entropy.  Those are different ladders, but the same meta-question:

```text
What did the proof compress before it found averageable mass?
```

For Collatz, the compressed object is an orbit/excursion scale.  For this LRC
determinant branch, the compressed object is the family of component-language
subsets.

## Next Proof Target

The concrete LRC route now looks sharper:

1. Formalize singleton determinant walls.
2. Formalize pair determinant incompatibility.
3. Bound the number `M` of component languages surviving cap/owner gates.
4. Only then worry about high-rank Helly or CRT/ZDD machinery.

If `H<=2` survives formalization, the logarithm is mostly bookkeeping.  If some
family forces `H` to grow, HYP-2152 tells us which log is being paid and why.
