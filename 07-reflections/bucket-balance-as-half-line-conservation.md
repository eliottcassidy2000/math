# Bucket Balance as Half-Line Conservation

**Session:** kind-pasteur-2026-05-30-S1  
**Status:** reflection after Lean formalization  
**Related:** THM-346, THM-348, INV-194, HYP-1775

The merged bucket balance theorem looked like a tiling-cube fact:

```text
2*self_b(M) + cross_b(M) = |bucket_b| * |M|.
```

After formalizing the first layer in Lean, the theorem splits into two ideas
with different characters.

The first idea is pure conservation.  For any finite bucket map and any finite
move set, every incident half-line from a bucket either stays in the bucket or
leaves it:

```text
|selfHalf_b(M)| + |crossHalf_b(M)| = |fiber_b| * |M|.
```

This part knows nothing about tournaments, masks, hypercubes, or quotient
graphs.  It is just a partition of a finite product.  That is why it belongs
in a generic Lean module (`TournamentH7.BucketBalance`) and why it is useful
as an engineering checksum: any finite-state feature extractor with buckets
and sampled moves should satisfy the same row identity.

The second idea is where the geometry enters.  In the tiling cube, the map

```text
(x,u) <-> (x xor u,u)
```

pairs internal half-lines with no fixed points.  That pairing is what turns
oriented internal half-lines into unordered internal lines counted twice.
Cross-lines, viewed from one bucket, are counted once.  So THM-346 is not a
single counting trick; it is:

```text
finite conservation + fixed-point-free involution = unordered quotient balance.
```

This factorization is useful mathematically because it tells the next Lean
step exactly what is missing.  We do not need to rediscover the quotient
formula inside the tournament model; we need an abstract involution-pairing
lemma and then a Boolean-mask instance.

It is also useful conceptually.  The same conservation law sits under several
objects we keep returning to: merged tournament buckets, even-graph buckets,
good-cut levels, projection defects, and future TDA perturbation features.
The interesting information is not the conservation itself.  The interesting
information is the excess and its routing: which half-lines escape, which stay
silent, and whether the escaping mass goes through spine, ribs, or sea.

That makes a clean division of labor:

- Lean should own the conservation and pairing skeletons.
- Computation should measure the transport excess by geometry.
- Engineering should use the identity as a cheap invariant/checksum before
  trusting any quotient feature matrix.

The theorem is small, but it has the right smell: a local row identity that
can be reused every time we turn a huge tiling cube into a finite quotient.
