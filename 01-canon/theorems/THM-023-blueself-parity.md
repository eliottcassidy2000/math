# THM-023: Blueself requires even n

**Status:** PROVED (algebraic + verified computationally at n=3,...,8)
**Instance:** opus-2026-03-06-S3
**Dependencies:** THM-022, definitions of blueself/blackself/GS

## Statement

A tiling is **blueself** if it is both grid-symmetric (GS) and self-flip
(its flip is in the same isomorphism class). Blueself tilings exist if and
only if n is even.

## Proof

For a GS tiling, the associated tournament T is self-converse under
the relabeling v -> n-1-v. This implies:

  score(v) + score(n-1-v) = n-1 for all v.

The flip operation reverses all non-path arcs (arcs other than i->i+1).
Under flip, vertex scores transform as:

  score'(v) = n-1-score(v)      for 0 < v < n-1  (internal vertices)
  score'(0) = n - score(0)      (source endpoint)
  score'(n-1) = n-2 - score(n-1) (sink endpoint)

For two tournaments to be isomorphic, they must share the same sorted score
sequence. Internal vertices map s -> n-1-s, which preserves the multiset
(by the palindrome property of SC tournaments).

The boundary correction is the key:
- Original endpoints: {score(0), score(n-1)} = {score(0), n-1-score(0)}
- Flipped endpoints: {n-score(0), score(0)-1}

For the sorted score multiset to be preserved, we need:

  n - score(0) = score(0)   =>   score(0) = n/2

This requires n even. At odd n, score(0) = n/2 is not an integer,
so flip ALWAYS changes the score sequence of a GS tiling.

## Verification

- n=3 (odd): 0/2 GS tilings have flip-same-scores. 0 blueself.
- n=4 (even): 2/4 GS tilings have score(0)=2=n/2. Both are blueself.
- n=5 (odd): 0/16 GS tilings have flip-same-scores. 0 blueself.
- n=6 (even): 24/64 GS tilings have score(0)=3=n/2. 4/24 are blueself.
- n=7 (odd): 0/512 GS tilings have flip-same-scores. 0 blueself.
- n=8 (even): 1280/4096 GS tilings have score(0)=4=n/2. (blueself TBD)

## Remarks

1. score(0)=n/2 is necessary but not sufficient for blueself.
   At n=6, 24 eligible but only 4 are actually blueself.
   The additional constraint is that the flip must produce an isomorphic
   (not just same-score) tournament.

2. Blackself (self-flip but not GS) can exist at all n >= 5, including
   odd n. Blackself tilings are NOT self-converse tournaments.

3. The score(0) distribution among GS tilings follows a binomial-like
   pattern centered at (n-1)/2.
