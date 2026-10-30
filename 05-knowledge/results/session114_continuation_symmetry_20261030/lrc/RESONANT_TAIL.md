# Exact resonant-tail transfer and a complete thirteen-speed family

Session 114, 2026-10-30. **PROVED:** the transfer theorem and the complete family below. **OPEN:** unrestricted Lonely Runner Conjecture, including fourteen runners. No literature-priority claim.

## Result

For a finite set of positive integer speeds V, put

M(V) = max over t in R/Z of min over v in V of ||vt||.

Let C={1,2,...,11,13}, and append any positive integer speed v. Then

- If 12 does not divide v, M(C union {v})=1/12.
- If v=12 or v=24, M(C union {v})=1/14.
- If v=12m with m>=3, M(C union {v})=m/(12m+5).

Thus this entire one-parameter thirteen-speed family satisfies the fourteen-runner bound. The case of a repeated speed is interpreted as a set, and obeys the same formula. The last line closes the all-m obligation left conjectural beyond m=10 in [THM-1230](../../../../01-canon/theorems/THM-1230-n14-gap-is-populated.md). Its m=1,2 exceptions are essential: substituting them in m/(12m+5) gives false values 1/17 and 2/29.

For m>=3 the maximum is attained at exactly two times:

 t_m=(5m+2)/(12m+5), and 1-t_m.

The active speeds there are 5 and 12m. The complete last-cell mass slope is

 B(C union {12m}, M)=2(1/5+1/(12m)).

Here B means mu S(M-epsilon)=epsilon B for sufficiently small positive epsilon, with closed endpoints retained. This is a concrete use of the marked endpoint information developed in Session 113: the relevant continuation statistic is the smallest **directional binding owner**, not an unweighted count of symmetry-related witnesses or just the sum of their reciprocal-owner weights.

## 1. A general zero-phase tail theorem

Let a core W have maximum clearance rho<1/2, with all its maximizing times in (1/Q)Z/Z. Suppose that for every r in (r0,rho], its survivor set is exactly the disjoint windows

 [c_j-(rho-r)/p_j, c_j+(rho-r)/n_j],

where Q c_j is integral and p_j,n_j are positive binding-owner speeds. Ties have already been resolved by choosing the **maximum** speed among owners on the same side, as required by physical radius relaxation. Let

 s=min_j min(p_j,n_j).

For a tail v divisible by Q, define r_v=v rho/(v+s). If r_v>r0, then

 M(W union {v})=r_v.

Moreover, at r_v every maximizing point is an outer endpoint on a side whose binding owner is s. These are all the maximizers (identify duplicate endpoints if needed).

**Proof.** The core bounds the maximum by rho. At any point t=c_j+x in a surviving window at radius r>r0,

 ||vt||=||vx|| <= v |x| <= v(rho-r)/s,

because vc_j is integral. Hence safety of the tail forces r<=v(rho-r)/s, or r<=r_v. Since r_v>r0, this excludes every radius above r_v, regardless of what happens below r0. Conversely choose a side with owner s and its outer endpoint at radius r_v. There |vx|=v(rho-r_v)/s=r_v<1/2, so the tail distance is exactly r_v; core safety follows from the window hypothesis. Equality in the preceding chain also proves the maximizer classification. No equidistribution assumption is used. Q need not be the denominator of rho.

When all core maximizers are unit fractions a/Q and rho=1/Q, every tail not divisible by Q preserves rho: its phase at each such center is a nonzero multiple of 1/Q. The exceptional multiples are precisely the tails that vanish simultaneously on the full unit orbit.

The theorem is conditional on an actual complete top-window description. It does not assert that every core has the required common zero-phase structure, nor that a symmetry orbit alone transports the directional owner data.

## 2. An explicit complete top-window certificate

The AP subcore {1,...,11} has maximum 1/12. Indeed, among the twelve points 0,t,...,11t one cyclic adjacent gap is at most 1/12, giving ||jt||<=1/12 for some 1<=j<=11. Equality requires all twelve gaps to equal 1/12, hence t=a/12 with gcd(a,12)=1. The extra speed 13 preserves these four times. Therefore

 M(C)=1/12, and argmax C={1/12,5/12,7/12,11/12}.

For every 1/14<r<=1/12, the complete closed survivor is exactly

 I_1(r)=[(1+r)/13,(1-r)/11],
 I_5(r)=[(2+r)/5,(3-r)/7],
 1-I_5(r), and 1-I_1(r).

These are respectively the windows with centers 1/12,5/12,7/12,11/12 and directional owner pairs

 (13,11), (5,7), (7,5), (11,13).

In particular s=5, not 1: at the center 1/12 the old AP owner 1 is tied with 13, and the faster speed 13 binds under relaxation.

Here is a finite exact certificate of completeness. At r0=1/14, intersecting the twelve one-speed safe sets gives exactly the four intervals

 [15/182,13/154], [29/70,41/98],
 [57/98,41/70], [141/154,167/182],

and the six singleton points a/14 for a=1,3,5,9,11,13. The checker independently verifies this equality both by rational interval intersection and by inspecting every wall and every intervening open cell; these exhaust the circle, not a sampled grid.

Each singleton has core clearance exactly r0, so none survives above r0. By nesting, every surviving point above r0 lies in one of the four listed intervals. In the first interval the positive 13-wall and negative 11-wall impose precisely I_1(r); in the second the positive 5-wall and negative 7-wall impose I_5(r); reflection gives the other two. Conversely the two endpoints of each proposed interval vary affinely between safe endpoints at r0 and the appropriate safe center at rho=1/12. For each core speed they stay in one fixed safe phase branch. The linear inequalities r<=vt-k<=1-r are convex in (t,r), so interpolation proves safety of both endpoints and the whole interval. This proves the formula throughout the entire top chamber, not merely in an unspecified neighborhood of rho.

## 3. Proof of the complete tail classification

If 12 does not divide v, each unit-twelfth center has nonzero phase modulo 1 under v. Its tail distance is at least 1/12. Thus the core maximum 1/12 is retained.

For v=12m the transfer theorem gives

 r_v=(12m)(1/12)/(12m+5)=m/(12m+5).

It lies above 1/14 exactly when 14m>12m+5, equivalently m>=3. The unique two directions with owner 5 are the left side about 5/12 and its reflected right side about 7/12. Their endpoints are exactly t_m and 1-t_m. At t_m the phase of speed 5 is r_v and the phase of speed 12m is 1-r_v. All other core speeds are strictly inside their safe arcs, since r_v is strictly inside the certified chamber. Hence the terminal slope is 2(1/5+1/(12m)).

For m=1,2, at any radius r>1/14 the same necessary inequality r<=m/(12m+5) contradicts r>1/14. At radius 1/14, all six unit-fourteenth points survive both the core and each tail, since neither 12 nor 24 is divisible by 14. The four positive-length core windows contain no tail-safe point: their largest possible displacement phase is at most (12m)(1/12-1/14)/5=m/35<1/14 for m=1,2. Thus their complete maximizer set is exactly the six unit-fourteenth points. This proves both exceptional values and excludes extra maximizers.

## 4. The whole family's terminal derivative

For a nonmultiple v of 12, the four unit-twelfth centers are still exactly the maximizers. Define the core residue owners

 c_1=13, c_5=5, c_7=7, c_11=11.

Then

 B(C union {v},1/12)=2 sum over u in {1,5,7,11} of 1/d_u,

where d_u=max(c_u,v) if v is congruent to u modulo 12, and d_u=c_u otherwise. A nonunit residue of v creates no new active wall at these centers. A unit residue changes two reflected boundary sides, and the maximum speed resolves any tie. This derives the formula directly from physical phases rather than assuming a modular symmetry preserves off-grid derivatives.

For v=12 or24, the derivative is the common AP14/GW value

 B=23324/6435.

For v=12m,m>=3, it is 2/5+1/(6m), as above. Thus the resonance at v=0 modulo12 changes not only the maximum clearance but also its number and owner structure of maximizing points. The limiting clearance 1/12 does not make these finite-v boundary structures interchangeable.

## 5. Verification and scope

`check_resonant_tail.py` uses exact Fraction arithmetic and explicit exception gates (active under python -O). It checks the full core certificate by two different finite algorithms, all pair-opposite-slope maximum candidates for a range of tails, the maximizing-time and derivative formulas, and large symbolic m witnesses. `results.json` records the actual tested ranges and check count. Independent root checks also reconstruct maxima from all pair-wall candidate clocks.

The proof applies for arbitrarily large integer v, whereas the finite checks are audits of its formulas. This closes the specific THM-1230 ladder and classifies all tails appended to this fixed core. It neither supplies an arbitrary twelve-speed core nor proves a universal lower bound for Session113's guarded relaxed mass target. Unrestricted LRC(14) remains open.
