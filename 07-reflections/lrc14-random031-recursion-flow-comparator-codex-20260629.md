# LRC14 Random031 Recursion Flow Comparator, codex-2026-06-29

HYP-3481 already made the topology visible: four mirror-paired punctures and a
bypassed saddle seam.  This pass asked which older recursion language actually
predicts the bypass.

The answer is a span, not a winner.

The `n*2` / two-adic picture explains the moving flow.  The `12` lower-delta
bypass hits are not just a count; they are two exact six-point blocks in
`u=2t`:

```text
branch0: phases 7,8,9,10,11,12 on component 54
branch1: phases 2,3,4,5,6,7 on component 43
```

They pair by mirror as `a <-> q-a`, flip branch, and swap component.  That is
too structured to be a generic bypass count.

The `n+2` picture explains the stationary boundary.  The seam carries all seven
owners

```text
(23,45,93,113,147,169,173)
```

while the bypass sees only `(23,93,113)`.  So the seam is not fake; it is the
boundary-owner debt.  It just is not the phase-flow route.

The proof shape should therefore be:

```text
seven-owner seam boundary
  + two-adic mirror bypass blocks
  + four puncture island current
  => discharge or named owner/two-adic/SPEC debt
```

Do not collapse the two bypass blocks to the scalar `12`.  The useful datum is
that they are ordered six-hit mirror blocks in the doubled phase coordinate.
