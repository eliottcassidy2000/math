# LRC14 Farey Mutation Certificate Scheduler

The product mutation looked useless in the S130 scalar reading because `p*q` created too many inversions. The literal reading is better: if the numerator is multiplied by the denominator inside the fraction, `(p*q)/q` collapses to `p`.

That collapse is dangerous globally, but it is exactly the right destructive quotient after the unit-excess gate. Once `e=14p-q=1`, the numerator is no longer just a numerator. It is the route index: `p=1` is the right-neighbor parent, `p=2` is the C27/petal/two-block strip, and `p>=3` is the K33/state-lift lane.

This is a small but useful repair to the Farey thread. The right sequence is not "choose q or p or p*q." It is:

```text
retain M=p/q -> compute e -> only then collapse to p if e=1.
```

The sum mutation is almost boring in the best way: `(p+q)/q=1+M`, so it is an affine check that a quotient preserved the actual Farey value. The power mutations go the other direction. They are too violent to prove the local inequality, but they are good stress tests for any invariant that claims to be robust while forgetting magnitude.

The proof stack now has a cleaner front/back split. Farey mutations schedule the packet. Fejer, Ramanujan, Kaczynski, endpoint-owner, C27, and K33/state-lift labels certify it.
