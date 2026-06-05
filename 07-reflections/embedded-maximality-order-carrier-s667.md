# Embedded Maximality Order Carrier (S667)

Prompt: merge in the concept of embedded maximality; think of the rational
numbers with the ordering less than; extend beyond that.

The clean formulation is:

```text
maximal(object, ambient embedding, allowed extensions).
```

That is the whole trick.  "Maximal" is not unary.  A finite search can have a
perfectly real maximum, but if the problem actually lives in a denser ambient
space, that maximum may be only a coordinate artifact.

## The Order Toy

In `(Q,<)`, a finite chain `[0,1,2]` has maximum `2`.  But the dense ambient can
insert `3` above it, or insert rationals into any of the internal cuts.  The
finite maximum is internal, not ambient.

The Dedekind-cut version is sharper.  Among rationals below `sqrt(2)` with
denominator at most `12`, the best lower approximation is `7/5`.  But `24/17`
is better and still below `sqrt(2)`.  The finite search maximum is destroyed by
enlarging the denominator bound; the completed cut has no rational maximum at
all.

The same point can also change status by changing the ambient boundary.  The
rational `1` is maximal in `{0,1/2,1}` and in the closed interval `[0,1]_Q`, but
not in `(0,1)_Q` or in all of `Q`.

So there are at least three distinct claims:

- internal maximum in a finite image;
- boundary maximum in a named closed ambient;
- embedded maximum against the allowed extension relation.

## Why This Hits LRC14

S666 found that AP, Vstar, and `2AP` are floor atoms in the visible `Res_27`
shell, but that visible shell alone is a bad ambient.  It has three huge mixed
buckets: each bucket contains one floor atom and `377` strict local carry rows.

That is exactly the finite-chain phenomenon.  The visible quotient has a local
maximum/floor label, but the ambient carry extension still moves through it.

The private-owner deletion flag acts like a cut address.  It does not merely
add more scalar information; it records which residue privately owns a proof
obligation, so deletion of that residue exposes the boundary.  In S666, adding
that one address kills all mixed floor/strict fibers while staying carry-free.

The proof target becomes:

```text
visible shell floor + embedded owner/carry address
=> AP, Vstar, 2AP, or strict looseness.
```

This is better than asking whether the visible quotient is maximal.  The
visible quotient is not the ambient.

## What Clocks Matter

Embedded maximality tells us which clocks are worth keeping.  A clock matters
when it indexes an allowed extension cut that can destroy local maximality.

For LRC14, that points to:

- the `+27` carry cocycle;
- owner-private deletion bits;
- endpoint protector owners;
- gcd-shell and lift addresses over `C=27`;
- HYP-2165 owner-route labels;
- coherent scalar lifts, because they are legitimate floor-preserving ambient
  moves rather than local leaks.

Clocks that only vary inside a separated embedded fiber can be ignored for this
proof stage.  They may still be good explanatory coordinates, but they are not
the obstruction boundary.

## Tournament Analysis

The S667 computation treats embedded-maximality lanes as tournament vertices.
The observable is exact evidence, ambient extendability, address need,
derivative leverage, LRC transfer, and actionability.

The resulting tournament is transitive:

- `directed_3cycles=0`;
- singleton SCCs;
- `hamiltonian_paths=1`;
- leader: `LRC14 owner-private stability`.

That transitivity is informative.  This session is not uncovering a wild
nontransitive analogy cluster.  It is tightening a definition: embedded
maximality names exactly why HYP-2241 felt stronger than another scalar repair.

## Next Move

The next local theorem should quantify allowed extensions, not arbitrary
perturbations.  For each extension family inside the visible `Res_27` fiber, ask
whether it:

1. changes the owner/deletion/carry address;
2. pays positive maximin tax;
3. is a globally coherent scalar floor lift.

That is the Dedekind-completion analogue for the LRC quotient: list the cuts,
attach addresses, and prove no unnamed cut remains.
