# Equinumerosity Versus Outer-Extension Usability, S669

The useful compression from this session is:

```text
equinumerosity is a shadow;
usability is a profile that survives extension.
```

The user's examples are exactly the right stress test.  An open interval, the
real line, and finite-dimensional Euclidean spaces all live at continuum
cardinality, but cardinality does not remember dimension.  Likewise, tournaments
and Royle-even graphs can be count-equivalent while the repo still needs a
predicate-preserving fiber map before any theorem can transfer.

The finite toy made that failure visible.  A path on 27 vertices, a 3x9 grid,
and a 3x3x3 grid all have the same finite cardinal shadow.  If the target
predicate is dimension, cardinality gives one mixed bucket.  Deletion profiles
are better: puncturing the path can disconnect it, while puncturing the grid or
cube does not.  But deletion alone still mixes 2D and 3D.  Degree, local growth,
and small embedding profiles separate all three.

That is a nice miniature of the topological fact: cardinality is not
homeomorphism.  The discrete deletion profile is the first homology-like
warning, and local growth/embedding is the finite dimension address.

The outer-extension check sharpened it.  Adding an outer layer changes size and
some boundary counts, but the smaller dimension signature stays stable:

```text
branching threshold, C4 presence, puncture splitting
```

This is the word "usable" doing work.  A profile is usable when it is stable
enough under the allowed extension to keep the predicate pure.

S617 now reads cleaner too.  Royle-even/tournament equinumerosity is not wrong;
it is just at the wrong proof layer unless we have the fiber functor.  Degree
even and naive Burnside-even were false substitutes, and even raw tournament
`H` splits by `beta1` and packet data.  So the map we want is not:

```text
same count -> same theorem.
```

It is:

```text
same count + retained H/beta1/packet/embedding profile -> transfer candidate.
```

This folds back into HYP-2243.  Friedman's homeomorphic-embedding world is the
high-proof-strength warning that the uniform theorem lives in the embedding
profile, not in the local count.  S668 saw that directly: raw outer address
leaked, embedding/downset profiles repaired.

For LRC14, this says not to ask whether `Res_27` or a raw cover count is the
right object by itself.  Ask what stable embedding profile sits over it:
owner-private deletion bits, carry cocycles, endpoint protector fragments,
pair-pinch fragments, and D/U/N owner routes.  The count is the address label on
the envelope; the proof is inside the retained fiber.

The next build should be an LRC14 proof-obligation embedding profile over
coherent carry subspaces.  If the analogy is doing real work, the profile should
separate floor/strict rows the way local growth separates equal-cardinality
dimensions, and the remaining failures should be named scalar floor lifts.
