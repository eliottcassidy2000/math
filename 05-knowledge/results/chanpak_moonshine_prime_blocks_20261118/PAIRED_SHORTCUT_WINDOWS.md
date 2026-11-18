# Three signed shrinking steps with physical witnesses across a 33-count window

The [signed arithmetic selector](CARRY_PARENT.md) chooses the sign s∈{−1,+1} for which 4 divides 3x+s and strips all powers of two. The earlier individual physical examples extend to an entire window:

> For every m≥16 and every odd beta in [−33,31], start at x=2^m+beta. The first three signed shrinking steps have literal tournament witnesses at both ends, on both the ordinary and the sevenfold axes.

The upper statements concern 7 times each visited value, not the ordinary Collatz trajectory of the multiplied value. Neither axis claims one fixed attachment turns every parent graph into every child. Rather the target chooses a context from the table at each step. This is a concrete instance of the freedom to re-encode the graph while retaining an exact arithmetic itinerary.

## Exact construction, not a bounded count scan

As long as there are sufficient powers of two in the large term, applying the selector to 3^j·2^w+b is exact: compute s from b modulo 4, put v=v2(3b+s), and replace (j,w,b) by (j+1,w−v,(3b+s)/2^v). For all the rows here m≥16 supplies w>v at each required step, so there is no carry into the valuation. Each signed step is strictly decreasing because v≥2 and x>1.

For every required (level,offset,axis) the table gives a literal marked wrapper from the one- or two-run construction. Its scalar read is G·2^n+intercept, where G=scale·3^level·2^a and intercept=scale·offset, with scale 1 or 7. The signed marked counts verify these coefficients and the read applies for any positive transitive width n. The target itself determines n=w−a; all rows have n≥1 at the common bound m≥16. The replay checks all 136 rows and 272 literal small graph controls independently of the arithmetic shortcut. In particular it includes the edge for beta=−33, whose second ordinary parent has offset −19 at mantissa 9. For the upper target 7(27·2^w+1), a different contact arrangement is needed. Put a transitive spine of length w+2 before the deletion D=(0,1,3,23,103,79,143,63), and let the apex beat the two first spine vertices and the D-mask 79. Here (h,h+R,kappa−2h)=(13,63,7), where R is the nonempty first-slot ending read into the apex and kappa counts D plus the apex. Odd-cycle packing gives 3·2^w(h+R)+kappa−2h: cycles through the apex confined to D contribute the whole core base kappa; a cycle confined to the spine adds h(3·2^w−2); and cycles meeting both add 3·2^w R because either chosen first spine position has a free parity bit. This gives exactly 7(27·2^w+1), rather than asserting the two-run wrapper already has that row.

Example: beta=−17, m=16 gives three prime counts 65,519→49,139→18,427, with actual sevenfold companions 458,633→343,973→128,989. Primality is checked by trial division; the general theorem makes no prime-distribution claim. The third step is also physical for this example, though we do not claim that its fourth count is prime.

## Where further iteration simplifies, and where it does not

On each offset alone, the signed stripping strictly decreases its magnitude above 1 and preserves its sign. Thus offsets eventually become ±1. For all later steps before the width is exhausted, the transition is just
```
3^j·2^w ± 1  →  3^(j+1)·2^(w−2) ± 1.
```
This identifies the natural infinite-family obligation rather than demanding a new offset table at each depth: supply signed central rails with unbounded powers of three and enough width, then connect their terminal small-width states to actual smaller witnesses. The archived all-order positive ordinary rail (the commuting clone and spine construction in “Both 81 rails”) already supplies 3^j·2^w+1 for j,w≥3 by commuting extensions. That can carry the ordinary positive lane for a number of steps growing linearly with initial width once its offset is +1; it does not give the corresponding upper or minus lane or a fixed terminal root. This is why a bounded window is useful but still not a proof of universal placement.

Replay: `python3 paired_collatz33_replay.py paired_collatz33.tsv`.
