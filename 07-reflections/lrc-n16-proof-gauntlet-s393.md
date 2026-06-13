---
source: codex-2026-05-31-S393
status: proof-search reflection
tags:
  - lonely-runner
  - n16
  - folded-antipodal
  - dyadic-debt
  - endpoint-leaf
---

# LRC n=16 Proof Gauntlet

This was the session where one of the pretty ideas broke in a useful way.

The pretty idea was the antipodal fold.  At `n=16`, set `s=2t`.  Even speeds
descend cleanly: if `v=2w`, then the two antipodal preimages `t=s/2` and
`t=s/2+1/2` have the same status for `v`.  Odd speeds are one-sided: they can
kill the left preimage or the right preimage, but only by putting `u*s/2` very
near `0` or very near `1/2`.

So I hoped for a theorem of this flavor:

```text
Find s safe for all even speeds.
Then at least one side of the antipodal pair is safe for all odd speeds.
```

That would prove the `16` case almost too beautifully.

It is false as stated.

The initial segment and the single `16`-gate row already have zero folded
witness mass.  Later, the sieve-aware beam found high-coverage rows with
`folded_witness=0` as well.  The fold can lose boundary information, and some
sieve-complete systems genuinely arrange odd lows and odd highs over the
even-safe quotient.

But the failure was not a dead end.  Every folded-danger row the beam found was
still positive-gap and endpoint-core-empty.

That changes the proof target:

```text
Folded danger should force an endpoint leaf.
```

## The Three Surviving Routes

### 1. Folded parity plus endpoint leaves

The fold partitions the problem into:

```text
even_safe(s)
odd_low(s):  side t=s/2 killed
odd_high(s): side t=s/2+1/2 killed
```

If `even_safe` is not contained in `odd_low cap odd_high`, we get an actual
lonely antipodal preimage.  If it is contained, the system is folded-dangerous.
The computations suggest folded-dangerous systems have exposed dyadic endpoint
leaves.

That is a more realistic theorem than the pure fold.

### 2. Local cover rigidity

The endpoints of a `16`-gate have a rigid lower-residue cover.  To cover all
`32` endpoints of `v=16` using lower residues, the exact minimum is:

```text
(1, 3, 5, 7, 8, 9, 11, 13, 15)
```

That is already nine speeds.  Add the `16`-gate itself and only five slots
remain.  I forced that seed and let a beam search choose five more speeds from
a structured pool.  The best completion was sieve-complete:

```text
(1,2,3,5,7,8,9,10,11,12,13,15,16,42,120)
```

but still had exact gap ratio `0.041667`, `36` unprotected endpoints, and
`coreE=0`.

This suggests a local rigidity lemma: the exact local cover of the `16`-gate
spends too much of the speed budget to also close the global endpoint core.

### 3. The half-turn moat

The full normalized two-torsion cube at `n=16` is finite and exact.  Among all
half-turn defect vectors, the best one still misses `128` cells:

```text
support 1: coordinate 15, missed 128
support 2: coordinates 10 and 15, missed 160
```

This is the cleanest scalar-quotient obstruction.  It says the most dyadic
non-scalar moves cannot even cover the finite quotient cell model.  A proof
could try to compress arbitrary residues toward half-turn data without losing
coverage; if that compression exists, the `128`-cell moat becomes a certificate.

## What I Would Try Next

The honest proof route now looks like a three-lock mechanism:

```text
1. THM-366 forces small-denominator gates.
2. The antipodal fold says either we get a quotient witness or enter folded danger.
3. Folded danger forces a private dyadic endpoint leaf.
```

The endpoint leaf is the key.  Folded danger means odd speeds are arranged in
opposite extreme bands over the even-safe quotient.  That arrangement seems to
spend exactly the labels needed to protect the `16`-gate endpoints, leaving
descendant dyadic endpoints private.

So the next lemma should not say "folded danger impossible."  It should say:

```text
folded danger has positive endpoint divergence.
```

That is a better problem.  It has a shape.
