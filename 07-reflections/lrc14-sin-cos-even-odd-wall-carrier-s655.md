# LRC14 Sin/Cos Even/Odd Wall Carrier S655

The even/odd thread connects to LRC `n=14`, but not in the blunt way.  This is
the wall-pair companion to HYP-2230/S654: that packet identifies the carry
coordinate linking parity and the `mod 14` obstruction, while this one asks
which odd walls and even slack labels the carry proof must preserve.

The tempting scalar picture is the S642/S643 duplicate branch:

```text
7 + 7 = 14
7 + 2*7 = 21.
```

That is real as a pair-projection shadow.  It is not the active LRC floor wall.
At `n=14`, AP and `Vstar` hit the floor at times

```text
1/14, 3/14, 5/14, 9/14, 11/14, 13/14.
```

The active odd complement pairs are

```text
(1,13), (5,9), (3,11), (3,11), (5,9), (1,13).
```

So the wall is off-diagonal and symmetric around `7`.  The even number `14` is
the pair sum.  The odd number `21` is the diagonal scalar shadow.

The trig prompt gives the correct typing:

```text
sin  -> odd boundary coordinate
cos  -> derivative/slack beside that boundary
cot  -> log-derivative scalar, useful but lossy
```

For AP at `t=1/14`, the odd layer hits `1/14`, while the even depth strata sit
at `2*delta`, `2*delta`, and `6*delta`.  That is HYP-2116 in local-slope form:
the odd speeds bind; the doubled prime-7 even layer supplies slack.  The active
slopes are also odd-pair data: `[1,-13]`, `[5,-9]`, `[3,-11]`, then the mirrored
signs.

This makes the S653 Basel lesson useful for LRC.  A clean scalar identity is
downstream of a retained packet object.  For Basel, the sine product retains
the elementary packets before zeta values appear.  For LRC `n=14`, the sine/cos
carrier retains:

```text
odd boundary pairs
even derivative slack
pair-sum pinch clocks
C=27 gcd shells
owner/carry labels
```

The `C=27` carrier is the load-bearing composite part.  AP and `Vstar` have the
same gcd-shell counts:

```text
{1: 9, 3: 3, 9: 1}.
```

That is exactly the HYP-2222 fixed-clock pocket.  The new information from S655
is that this pocket also has the same odd complement wall.  In the AP
single-swap census through `n=14`, no row falls below the floor, and the only
tight non-AP row is the known `12 -> 24` move.

The cotangent check gives the caution:

```text
prod_{k=1}^{(C-1)/2} cot(pi*k/C)^2 = 1/C
```

for the odd clocks tested.  Lovely scalar, bad proof object.  It cannot tell a
prime clock from the composite `C=27` carrier unless the gcd-shell labels are
kept separately.

The progress target is therefore not "explain 14 and 21."  It is:

```text
odd wall fixed + C=27 gcd carrier fixed + carry/owner visible
  => AP, Vstar, or strict looseness.
```

That is a much smaller proof obligation than a fresh global LRC `n=14` search.
It asks for a no-leak lemma over retained side channels.  The quotient can
forget raw runner order for the first pass, but it cannot forget which odd
boundary pair owns the floor and which `C=27` gcd shell owns the lift.
