# Two-seed mixed Collatz descent with a hole-safe reverse compiler

**Proved:** Every positive odd integer outside `{1,5,7,21}` has a strictly
decreasing mixed-sign odd Collatz block of at most 10 edges, each with
exact valuation at most 5. Each odd endpoint and every count in the
corresponding reverse dyadic ladder avoids `{7,21}`. Thus these arithmetic
moves have precisely two descending roots, 1 and 5. This is neither
ordinary Collatz nor a proof that the moves lift to tournaments.

## 1. The ladder and the arithmetic guard

Write an odd step as (3n+s=2^a m), where (s=\pm1),
(a=v_2(3n+s)) and (m) is positive odd. Necessarily (3\nmid m).
Define (D_t(x)=2x+t). Reversal can be written
`m → 2m−s → 4m−s → ... → 2^a m−s = 3n → n`:
first (D_{-s}), then (a-1) applications of (D_s), finally division
by three. Call the edge *hole-safe* if (m) and every intermediate odd
count (2^i m-s), (1\le i\le a), avoid 7 and 21.

Solving (2^i m-s=7) or (21) with (m) odd gives exactly
((m,s,i)=(1,+1,3),(11,+1,1),(3,-1,1),(5,-1,2)).
Discard the third because (3\nmid m). For (a\le5), a bad
ladder in the remaining three cases occurs exactly at the sources
(5,27,29,117): substitute (n=(2^a m-s)/3), requiring integral (n)
and (a\ge i). The case (a=i) is already impossible since
(3n\notin\{7,21\}) for allowed sources. The only bad odd endpoint
is 7, because 21 is divisible by three.

## 2. The uniform descent

For an allowed odd (n>1), let (s\in\{\pm1\}) be congruent to
(n\pmod4), and write (3n+s=2^a k), so (a\ge2).
If (a\le5), the immediate endpoint (k) is smaller than (n)
since (3n+s<4n). It is hole-safe except for the four ladder exceptions
above and for (k=7), whose sources are exactly (9,19,37,75).

If (a\ge6), use the other sign (-s) three times. Exact valuations
are ((1,2,2)), since all the powers of two displayed below have
exponent at least one:
[
 y=2^{a-1}k-s,qquad z=3\cdot2^{a-3}k-s,qquad
 w=9\cdot2^{a-5}k-s=(27n-23s)/32.
]
The outputs satisfy (y\ge31,z\ge23,w\ge17), so none equals 7
or can cause one of the three admissible bad ladders. No odd output
is divisible by 3. Also (w<n): only the sign (s=-1) needs checking,
and then (n>23/5), as already follows from (a\ge6).

Use the following small exception table (sign strings correspond
left to right to the edges):

| source | odd itinerary | signs |
|---|---|---|
| 9 | 9,13,5 | −+ |
| 19 | 19,29,43,65,49,37,55,41,31,23,17 | +−+++−−+−− |
| 27 | 27,41,31,23 | ++− |
| 29 | 29,43,65,49,37,55,41,31,23 | −+++−−+− |
| 37 | 37,55,41,31 | −−+ |
| 75 | 75,113,85,127,95,71 | ++−−− |
| 117 | 117,175,131,49 | −−− |

The length 10 is sharp with valuation cap 5: the exact minima in
the first nine endpoint layers from 19 are
`29,43,65,49,37,55,41,31,23`. The next layer reaches 17. The
verifier exhausts all safe sign choices at each layer, not just
the exhibited route. The expanded peak (including the dyadic ladders)
is at most (195n/19), also sharp: the first three edges at 19 force
29,43,65, and leaving 65 forces the ladder endpoint 195. The exhibited
route never exceeds it. In general each ladder peaks at (3x) for its
source (x); direct cases peak at (3n), high-valuation reroutes peak
at most ((9n+3)/2), and the displayed exceptions give the rest.

At 5, the minus edge directly hits 7; the plus edge goes to 1 but
its reverse dyadic ladder (1\to1\to3\to7\to15) hits 7.
So no admissible first edge exists at 5 regardless of the valuation cap.
At 1 no strictly smaller positive integer exists. Iterating the theorem
gives termination at one of the two seeds by strict descent at block
boundaries. An arbitrary sequence of sign choices is not asserted to
terminate.

## 3. Why keep both signs, and why include valuations ≥3?

One fixed sign does not admit a uniform bounded stopping-time argument:
for (3n+1), positive (n\equiv-1\pmod{2^{L+2}}) have the first
(L) odd steps at exact valuation 1, by shadowing the rational affine
fixed point (-1) and losing one binary digit per step. For (3n-1)
the mirror point is (+1); use arbitrarily large positive
(n\equiv1\pmod{2^{L+2}}). Each of those (L) steps strictly
increases positive (n>1). These are finite shadow obstructions, not
statements about those integer orbits afterwards.

Even with both signs, restricting to exact valuations 1 and 2 prevents
any uniform bound on block length. In the 2-adic odds, the rational
cycle (-1/5\to1/5\to-1/5) uses plus then minus, each at exact
valuation 1. At each point the other sign produces numerator
(\pm8/5), of exact valuation 3. Positive integers arbitrarily
far out in the residue class (-1/5\pmod{2^{L+3}}) shadow (L)
steps of this cycle: the valuation-1 move loses just one bit of
precision, leaving enough to enforce the same choice and verify that
the other has valuation 3. Every allowed step in that prefix is
strictly increasing. Thus valuation cap 2 cannot support bounded
descent, even ignoring holes. The result does not decide caps 3 or 4.

## 4. What would a tournament lift have to prove?

Arithmetic seeds are not readiness seeds automatically. One sufficient
interface would attach to every realized odd target a *ready* tournament
representative, containing ready 1 and 5 seeds, so that any of the
short **reverse blocks** above can be executed from a ready output
through the indicated dyadic ladders and certified divisions by 3,
ending in a ready representative of the source. Only reentry at a
block boundary is required; internal interface states may differ.
Reverse whole decreasing blocks in well-founded order; internal
odd or ladder counts can exceed the target, and cannot be supplied by
an induction hypothesis. One may alternatively synthesize each block
as a macro and avoid scalar lifting of its individual edges.

The integer result does not assert either a universal (D_\pm) move
or arbitrary division by 3 for tournament counts. Session 303's
guarded triangle peel is a candidate only where a particular constructed
tournament actually has its eligible triangle. In fact an (H=15)
tournament can have a peel, yet the arithmetic plus edge at 5 is still
unusable by this *factorization* from 1 because it asks for 7 in the
middle. A physical macro need not share that obstruction.
