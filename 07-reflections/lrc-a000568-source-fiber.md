# LRC A000568 source fiber

The useful correction is small but sharp: A000568 is not the answer to LRC, but
it is hiding in the answer.

Raw unlabelled tournament classes mix lonely and non-lonely states.  That was
already visible in S509/S512/S535, and it remains true here: moving-runner
classes can contain both good and bad states.  The new computation isolates
the clean layer underneath.  If the stationary observer is attached by
threshold edges, then loneliness is exactly the observer being a source.  Once
the observer is a source, deleting it leaves an ordinary unlabelled tournament
class on the moving runners.  Adding the source back is faithful.

So the picture is:

```text
raw A000568 class          too coarse
observer-rooted phase      still loses threshold data
observer-source lift       exact LRC predicate
source-deleted class       A000568 base that survives at a witness
```

This reconciles the old A000568 enthusiasm with the projection-defect warnings.
The target is not "which raw class is good?" but "can the arithmetic clock reach
a source cone over some deleted class?"

For LRC14 this matters because THM-497 says cardinality will not save us:
thirteen bands can have enough total mass to cover the unit group.  The missing
proof ingredient must be alignment and fiber structure.  The source-deleted
A000568 class is one structural coordinate; Q27/Q31, divisor fibers,
7-ideal/13-clock debt, owner pressure, and blocking-height support are the
multiplicative-cover coordinates.

The creative proof route now feels less foggy.  A no-witness row is not merely
"avoiding all good times."  It is avoiding all source-cone entries in a
decorated quotient stack.  If it avoids them for a long dilated-band interval,
then the avoidance itself should impose balanced-cover congruences.  That is
exactly where HYP-2471/HYP-2480/THM-497 already point: the Q31 portal, the
7-ideal packet, and the 13-clock escape.

The next computation should pair each hard-row shell with two labels at once:
the deleted A000568 fingerprint of the runner phase and the multiplicative
cover certificate in `(Z/q)^*`.  Repeated deleted fingerprints with different
first witnesses would prove that the A000568 base is necessary but not
sufficient; repeated cover certificates with different deleted fingerprints
would show the opposite.  The proof likely lives in the intersection.
