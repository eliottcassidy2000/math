---
source: codex-2026-05-31-S397
status: exploratory reflection
tags:
  - lonely-runner
  - n16
  - egyptian-fractions
  - q-counterterm
  - reciprocal-residue
---

# LRC, Egyptian Fractions, and q

The S394 `q` idea says:

```text
p(x)=min(x,1/x)=1/x+q(x),
q(x)=x-1/x on (0,1].
```

That is not a positive decomposition.  It is a finite-part decomposition: the
infinite logarithmic tail of `1/x` is cancelled by an anti-reciprocal
counterterm.

Egyptian fractions are the discrete version of the same gesture.  They ask
whether a reciprocal mass can be split into simpler reciprocal atoms:

```text
k/N = 1/b + 1/c.
```

The classical divisor transform is:

```text
(kb-N)(kc-N)=N^2.
```

So the split exists exactly when some divisor `d|N^2` has the right residue:

```text
d == -N mod k.
```

That is already close to the n=16 Lonely Runner endpoint test.  For an endpoint

```text
t=(16m+eps)/(16u),
```

protector `p` works exactly when

```text
|p*(16m+eps)-16*a*u| < u
```

for some integer `a`.  Egyptian fractions have an exact divisor congruence.
LRC has a strict small-residue inequality.  The same grammar is hiding in both:
split reciprocal mass without leaving residue.

## The Coordinate Ledger

S397 compares the coordinates `1..15` by:

1. whether `16/c` splits into two unit fractions;
2. whether `c` is a product-sum target;
3. how many `u=16` endpoints residue `c` protects;
4. how many scalar-quotient cells a half-turn defect at `c` misses;
5. whether `c` lies in the exact nine-residue lower cover of the `16`-gate.

The table found the expected old fact in a new language:

```text
c=15 is best one-coordinate half-turn defect:
  missed=128
  16/15 = 1 + 1/15
  15 = -1 mod 16
```

This is the anti-reciprocal coordinate.  It is the finite residue version of
`q0=x-1/x`.

The best support-2 half-turn pair is:

```text
(10,15), missed=160.
```

That pair is more interesting than either coordinate alone.  Coordinate `10`
is unsplit at numerator `16`; coordinate `15` is the anti-residue with the
trivial split.  This looks exactly like a q-counterterm pair:

```text
unsplit debt + anti-reciprocal repair.
```

It may be the two-coordinate shadow of the endpoint-debt proof.

## The Harmonic Shell

The Egyptian-fraction approximation

```text
ln 2 = lim_m sum_{m<j<=2m} 1/j
```

is the rectangle-sum version of the S394 hyperbola area.  At `m=16`, the shell
`17..32` is the first denominator octave beyond the `16`-gate:

```text
sum_{17}^{32} 1/j = 13981692518567/20629078984800
error = -0.015380978...
```

This error is small but not zero.  In the LRC analogy, that error is endpoint
debt: passing from the continuous log area to finite reciprocal atoms leaves a
residue that must be paid.

## Breaker Speeds

The exact local lower cover for the endpoints of a `16`-gate is:

```text
(1,3,5,7,8,9,11,13,15).
```

That leaves six even non-cover residues:

```text
(2,4,6,10,12,14).
```

Adding the `16`-gate gives ten speeds, so a 15-speed system can include only
five of those six simple even breakers.  S397 tested every such choice.  All
stay positive-gap.  The closest simple row omits `14` and still has:

```text
gap/th = 0.060606
unprotected endpoints = 12.
```

The high-speed probes tied to the `(10,15)` half-turn pair also stay
positive-gap and expose many endpoints.

So the Egyptian angle did not produce a counterexample.  It produced a new
proof target.

## Proof Target

Try to prove:

```text
Any n=16 all-protected endpoint cycle induces an Egyptian split ledger.
Unsplit residues require even breaker speeds.
Primitive even breakers create positive q-defect.
Positive q-defect is witnessed as a gap or private endpoint leaf.
```

This would merge four recent n=16 routes:

1. THM-366/367 dyadic endpoint sieve and capacities.
2. S392 antipodal fan incompatibility.
3. S393 folded danger forcing endpoint leaves.
4. S394 q-counterterms as finite reciprocal residues.

The aesthetic version is:

```text
LRC n=16 is an Egyptian fraction problem with inequalities instead of equality.
```

The practical version is:

```text
Add Egyptian splitability and q-defect to the branch-and-bound feature vector.
```

That may give the search a much better nose for which residues are genuinely
dangerous and which are only large by dyadic valuation.
