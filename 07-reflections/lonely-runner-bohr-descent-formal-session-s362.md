---
source: codex-2026-05-30-S362
status: formalization session
tags:
  - lonely-runner
  - diophantine-approximation
  - endpoint-protection
  - bohr-boundary-descent
  - peeling
---

# Lonely Runner Bohr-Descent Formal Session

## What Became Formal

The S361 slogan was:

```text
Lonely Runner = finite anti-Bohr boundary theorem.
```

This session formalized three pieces of that slogan.

## 1. Initial Segments Are Dirichlet Equality Cases

THM-358 proves the unit-boundary skeleton for the initial segment

```text
V_n = {1,2,...,n-1}.
```

The proof is exactly Dirichlet's pigeonhole argument, but read at equality.
Look at

```text
0, t, 2t, ..., (n-1)t
```

on the circle.  If two points are closer than `1/n`, their difference gives a
forbidden speed `v` with `||v t|| < 1/n`.  If this never happens, all `n` gaps
are at least `1/n`, hence exactly `1/n`, so the points form the regular
`n`-gon.  Therefore

```text
t = a/n,  gcd(a,n)=1.
```

These are the safe boundary points.  They are also forbidden endpoints, because
for a unit `a` there is a speed `v` with `av == +/-1 mod n`.

So the initial-segment tight examples are not mysterious.  They are the rigid
equality case of the most elementary simultaneous approximation lemma.

This also clarifies a small language bug: `0` belongs to the pigeonhole orbit,
but it is a forbidden center, not a forbidden endpoint.  The endpoint skeleton
is the nonzero unit group:

```text
(Z/nZ)^* / n.
```

## 2. Peeling Has A Core Theorem

THM-359 formalizes the endpoint/interval peeling operation.

Given endpoints `E`, intervals `I`, interval boundaries `B(i)`, and endpoint
protectors `P(e)`, a protection core is a pair `(E',I')` such that:

```text
each e in E' has P(e) cap I' nonempty
each i in I' has B(i) subset E'
```

The peeling algorithm removes unprotected endpoints, then removes intervals
whose boundary was removed.  THM-359 proves this process computes the largest
core.  Empty terminal core is therefore a real finite certificate for this
incidence notion, not merely a diagnostic.

## 3. Unit Protection Forces Divisibility

THM-360 proves the first exact quotient filter.

For `n=k+1`, a unit endpoint

```text
t = a/n,  gcd(a,n)=1
```

can only be strictly protected by a speed divisible by `n`.  If `v` is not
divisible by `n`, then `va mod n` is nonzero, so

```text
||v a/n|| >= 1/n.
```

Thus any full-open-cover counterexample must contain at least one speed
divisible by `k+1`.  Conversely, one speed divisible by `n` protects every unit
endpoint at once.

## Computed Evidence

`04-computation/lonely_runner_bohr_descent_s362.py` implements the exact
endpoint/interval system over rational arithmetic.

Main output:

```text
THM-358 audit: initial segments n=2..36, failures=0
full-measure inherited boxes: no nonempty cores
```

The full-measure primitive-box audit is:

```text
k=3, max_speed=24: 1 full-measure case, 0 nonempty cores
k=4, max_speed=24: 2 full-measure cases, 0 nonempty cores
k=5, max_speed=20: 2 full-measure cases, 0 nonempty cores
k=6, max_speed=16: 1 full-measure case, 0 nonempty cores
k=7, max_speed=14: 3 full-measure cases, 0 nonempty cores
```

In every full-measure case, the first peel layer is the unit quotient
`unit_mod_n`.  Near-tight positive-gap examples begin at higher denominator
layers instead.

## What This Suggests

The speculative descent picture now has a hard base case and a finite core
certificate:

```text
initial segment       -> Dirichlet equality/unit quotient
tight sporadic cases  -> unit quotient first, then finite peel cascade
positive-gap cases    -> higher denominator residue layers
counterexample        -> would need a nonempty terminal protection core
```

If an LRC counterexample exists, it cannot merely imitate the known tight
families.  It would need to protect the unit layer and also prevent every later
endpoint/interval peel from exposing a quotient layer.  That is a much sharper
target than "find a lonely time."

## Next Formal Targets

1. Extend the peeling audit beyond full-measure cases by ranking positive-gap
   examples by peel depth and first quotient layer.
2. Add protector-speed ledgers for the unit layer: which speeds can protect
   `a/n`, and what additional endpoints do those speeds expose?
3. Compare peel layers with Fejer/Riesz safe-product kernels from HYP-1812.
4. Search artificially generated interval systems for nonempty cores, to learn
   what a genuine all-protected core must look like.
