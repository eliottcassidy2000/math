# One height, several physical ports: exact disjoint inverse bands

The [carry-parent construction](CARRY_PARENT.md) used the sink or penultimate vertex after gluing a transitive tail. Here every interior tail position gives a different reading while the ordinary height stays the same. This supplies distinct parent choices on a larger state space.

Start from any rooted tournament with state (h,e,t,d) in the [paired recursion](PAIRED_RECURSION.md). Attach a transitive tail of k+1 vertices by reversing only the cross-arc from its sink back to the old root. Mark any tail vertex with i preceding and j following tail vertices, i,j≥1 and k=i+j. The new state is
```
N = h+2^k t,     E = 2^i t,
T = 2^j h+2^(i+1)3^(j−1)t,
D = h+2^(k−1)t.
```
Here N is the ordinary height and T the distinguished-start optional-two-slot read, not a second graph. For E, an ordinary path ending at the mark must use the reversed edge; the i earlier tail vertices can lie before that edge or after the returned old-block suffix. For T, without the reverse edge the j later vertices can be placed in either slot while the old block stays in the ordinary slot. With the reverse edge either slot can use it; for each choice the i earlier vertices have two possible pieces and the j−1 vertices between mark and sink have three. The old vertices split into the distinguished suffix starting at the old root and an optional path, giving t. Deleting the mark shortens the tail by one. For the first tail mark i=0 with j≥2 the same proof works but subtracts e from T: the old optional path can use the lost forward arc precisely if all tail vertices other than the sink lie in the marked path.

## Which port can the target have used?

Put C_j=(3/2)^(j−1) and w=(N−h)/h. For i≥1,
```
T/N = C_j + (2^j−C_j)/(1+w),       w≥2^(j+1),
```
since t≥h. Hence C_j<T/N<C_(j+1), and these are disjoint bands. A proposed (N,T) picks **one** j and the exact smaller parent height h=(T−C_j N)/(2^j−C_j). It leaves finitely many scales k≥j+1 with t=(N−h)/2^k integral and at least h. No rounding or asymptotic estimate is involved. For the first mark i=0, the correction term is −e/N, w≥2^j, and j≥2. The same strict band persists: 2^j−C_j>1, e≤h, and (2^j−C_j)/(1+2^j)<C_j/2. But recovering h then uses the additional coordinate D=(N+h)/2 or a supplied e. Thus this boundary can be a genuine alternative port rather than an ambiguous band's worth of choices.

If the two inputs scale in just (h,t), all interior output reads scale. The first-mark branch also reads e. This is why the interior operation can synchronize independently supplied old states before a new sink reset. At a fixed seed (h,t), a two-stage exact inverse tries i,j≥1 within N<p and tests whether (p−N)/T is a positive power of two at least two. Successful choices determine the second tail width.

From the 5/35 entry roots in the replay, the target **421** selects (i,j,r)=(1,2,2), and **613** selects (1,3,1). In the first case the intermediate (N,T) is (53,92); in the second it is (101,256). A common sink tail with r+1 vertices gives respectively 421 and 2,947, or 613 and 4,291. All graphs were literally counted. The replay also verifies 174 small four-state transfers including the first-mark boundary. This gives target-selected changes in the actual continuation while keeping the intermediate height fixed. It does not prove a feasible marked T exists for each target.

Replay: `python3 interior_ports.py`.
