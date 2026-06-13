---
source: codex-2026-05-31-S389
status: proof attempt / structural synthesis
tags:
  - lonely-runner
  - n16
  - cayley-dickson
  - dyadic
  - endpoint-debt
---

# Lonely Runner n=16 as a Cayley-Dickson Gate

The prompt was to attempt an `n=16` proof while thinking Cayley-Dickson-ly.
That turns out to be a good instruction, if interpreted literally enough.

The Cayley-Dickson tower is not just "powers of two are mystical."  It is a
ledger of what survives under doubling and what starts leaking.  At `16`, the
sedenion layer still exists as a doubled algebra, but zero divisors appear.  A
formal operation that looked like it should preserve a norm now exposes a new
kernel.

The `n=16` Lonely Runner gate has the same shape.

## The Real Lemma

There is one clean proof branch already.

At denominator `16`, the odd unit endpoints are

```text
1/16, 3/16, 5/16, ..., 15/16.
```

For odd `a`, speed `v` strictly protects `a/16` exactly when

```text
||v*a/16|| < 1/16.
```

Since `v*a mod 16` is an integer residue, this inequality forces

```text
v*a == 0 mod 16.
```

Because `a` is odd, that means `16 | v`.  So:

```text
No 16-multiple -> no open cover.
```

That is not analogy.  That is a usable first branch of an `n=16` proof.

## The Gate Is a Zero-Divisor Move

The tempting thought is that adding a `16`-multiple closes the problem: it
protects all odd unit endpoints at once.  But the endpoint audit says the
opposite.  The `16`-gate closes the old unit layer by opening a descendant
layer.

The dyadic ladders are the cleanest witness:

```text
d=2  ladder: gap/th=1/33,  unprotected=34
d=4  ladder: gap/th=1/66,  unprotected=68
d=8  ladder: gap/th=1/132, unprotected=140
d=16 ladder: gap/th=1/264, unprotected=280
```

The visible gap halves.  The endpoint debt doubles.  The product

```text
unprotected * (gap / threshold)
```

stays essentially constant (`34/33`, `34/33`, `35/33`, `35/33`).  This is the
Cayley-Dickson feeling: a norm-like quantity survives, but what it measures has
changed.  You did not annihilate the obstruction; you moved it into a new
kernel.

The endpoint layers make the movement literal:

```text
d=2  leaks at v2(denominator)=5
d=4  leaks at v2(denominator)=6
d=8  leaks at v2(denominator)=7
d=16 leaks at v2(denominator)=8
```

So the `16`-gate is not the proof.  It is the zero-divisor branch: it kills the
old witnesses and creates new ones with doubled denominators.

## The Scalar Moat

The normalized scalar quotient gives the second half of the attempted proof.

For `n=16`, the exact cell system has:

```text
1152 alpha patterns
18432 candidate cells
```

All scalar ramps cover this finite system.  But moving away from a scalar ramp
has a moat:

```text
one-coordinate defect: best missed = 128
two-coordinate defect: best missed = 160
```

The best one-defect is a half-turn in the last coordinate.  The best two-defect
uses half-turns at coordinates `10` and `15`.  This is again a pure dyadic
phenomenon: the only near-dangerous defects are half-turns, and even they leak
many cells.

## The Attempted Proof Shape

The proof I would now try to write is:

```text
Assume a primitive n=16 counterexample.

1. By the unit-gate lemma it has a 16-multiple.
2. Quotient by the scalar ramp/gauge direction.
3. If the residue is close to scalar, the scalar moat gives a positive gap.
4. If the residue is dyadic-ladder-like, the debt norm gives a positive gap
   or an exposed endpoint.
5. Otherwise, the labelled endpoint-protection graph has a private dyadic leaf,
   so THM-365's endpoint cycle cannot be realized.
```

The missing lemma is exactly step 5:

```text
Every primitive n=16 all-protected endpoint system with a 16-gate
either descends to the dyadic debt ledger, or has a private dyadic leaf.
```

This is not yet a proof.  But it is now a smaller and much more shaped problem.
The analogy did useful work: it told us not to celebrate the `16`-gate, but to
ask what norm it fails to preserve and where the zero divisors leak.

That feels like the right `n=16` question.
