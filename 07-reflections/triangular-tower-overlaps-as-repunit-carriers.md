# Triangular Tower Overlaps As Repunit Carriers

The tempting story was that `[21,24]` should be the master bridge.
It is the unique exact side shared by both towers, so it feels like the place
where everything lines up.

But once the overlap is turned into a polynomial carrier, that instinct is
wrong in a useful way.

A consecutive block of length `L` is not just an interval.  In the fixed-path
coefficient-row grammar it is the support row

```text
R_L(x)=1+x+...+x^(L-1),
```

up to a shift by `x^a`.  Then the exact hinge has `L=4`, so its row is

```text
1+x+x^2+x^3=(x+1)(x^2+1),
```

reducible at once.  Exactness of overlap is not atomness of the polynomial
carrier.

That pushes the eye to the Pell families instead of the singleton hinge.
The whole-equation family `T_n=2T_m` has lengths

```text
5,29,169,985,5741,33461,...
```

and already in the stored window it supplies four prime lengths.  Those give
exact cyclotomic irreducibility without any heuristic step, and `5` and `29`
also land on prime Mersenne values at base `2`.  So the richest
irreducibility lane is not the exact hinge but the looser symmetric
containment family.

This feels very repo-native.  The scalar-looking invariant is not the end of
the story, but neither is the most literal structural coincidence.  The right
carrier sits one step sideways: not “which overlap is exact?” but “which
overlap lengths are atomic enough to survive as irreducible coefficient rows?”

The family tournament makes that clash visible.  Support/exactness and
prime-length supply do not rank the families the same way; the prime-length
observable flips eight edges relative to the support-only comparison and
produces a genuinely nontransitive tournament.  That is the useful outcome.
If the ranking had stayed transitive, the polynomial bridge would have been a
cosmetic relabeling.  Instead it changes the geometry.

The next question is cleaner than the first one.  Forget the unique hinge for
a moment and ask:

```text
which Pell overlap families contain infinitely many prime lengths?
```

A positive answer would give infinitely many exact cyclotomic carriers, and
any extra Mersenne/Cohn hits would be even stronger.  That is a sharper
bridge from the triangular-tower story to irreducible polynomials than the
original exact-overlap fascination.
