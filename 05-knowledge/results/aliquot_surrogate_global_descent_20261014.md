# Global descent for the even prime cofactor surrogate

Research theorem, 14 October 2026. Every orbit of the affine surrogate defined below eventually reaches a fixed point or one of two explicit nonconstant cycles. In particular, every finite aliquot cycle whose members are all even with prime odd part is constant.

This proves both conjectures C1 and C2 from the earlier [all-length frontier](aliquot_all_length_surrogate_frontier_20261014.md). That predecessor is preserved as a record of the previously open questions. The proof here is unbounded and elementary; the companion computations check its identities and branches rather than supplying a larger cycle census. Global novelty has not been established.

## 1 The map and the exact arithmetic scope

For a positive even integer `n`, write uniquely

`n=x p`, where `x=2^v2(n)` and `p` is positive and odd.

Define

`F(n)=(x-1)p+2x-1=n-p+2x-1`.

It maps positive even integers to positive even integers. Its fixed points are exactly

`n=x(2x-1)` for powers of two `x>=2`.

The odd part need not be prime for this surrogate. Its discrepancy from the genuine aliquot map `s(n)=sigma(n)-n` is

`s(xp)-F(xp)=(2x-1)(sigma(p)-p-1)`.

Thus `F=s` precisely when the positive odd part `p` is prime. Composite odd parts and the value `p=1` must retain their actual divisor sums when discussing aliquot dynamics.

**Theorem.** Every positive even integer has an `F` orbit that eventually reaches either a fixed point or one of the two cycles

`8 -> 22 -> 14 -> 10 -> 8`,

`16 -> 46 -> 26 -> 16`.

These are the only nonconstant primitive cycles, up to cyclic rotation. Both contain a power of two. Consequently, all integral dyadic affine template closures are repetitions of fixed points or these two cycles. This last formulation uses the exact template-to-`F` equivalence proved in the predecessor: integrality forces the prescribed cofactors to be odd, so the dyadic cores are genuine.

## 2 The bounded descent statement

The key statement is stronger than a restriction on periodic orbits:

**Descent lemma.** Every positive even `n` is a fixed point, belongs to

`E={2,4,8,12,16}`,

or satisfies `F^j(n)<n` for some integer `1<=j<=6`.

The exceptional set is explicit, not an unclassified remainder: `2 -> 4 -> 10` enters the first cycle, `4 -> 10` enters it, `12 -> 16` enters the second, and `8,16` are already on the displayed cycles.

Write `n=xp` as above. Since

`F(n)-n=2x-1-p`,

the cases `p>2x-1` and `p=2x-1` give an immediate descent and a fixed point, respectively. It remains to treat `p<2x-1`, where the first step increases. The proof controls the next few exact dyadic cores rather than assuming that size or valuation decreases on every edge.

## 3 Powers of two and the first quadratic exception

For a power of two `M>=8`, direct calculation gives

`F(M)=3M-2`,

`F^2(M)=3M/2+2`,

`F^3(M)=3M/4+4`.                                      (3.1)

The first two displayed successors have exact dyadic core `2`: their odd parts are `3M/2-1` and `3M/4+1`, respectively. In particular, `M>=32` gives `F^3(M)<M`. The smaller powers `2,4,8,16` lie in `E`.

Now suppose `p>1` and `p=x-1`. Then

`n=x(x-1) -> x^2`.

This increasing edge has increased the dyadic core, but its target is a power of two. For `x>=8`, (3.1) yields

`F^4(n)=3x^2/4+4`,

`4(n-F^4(n))=x^2-4x-16=(x-8)(x+4)+16>0`.

The remaining case is `x=4`, namely `n=12` in `E`. This explicitly handles a power-of-two hit that is larger than the initial state; such a hit alone would not justify strong induction.

## 4 Two steps force a Mersenne-shaped odd part at a cycle minimum

In the remaining increasing cases, factor

`p+1=yq`, where `y=2^v2(p+1)` and `q` is odd.

Because `p+1<2x` and `p+1!=x`, one has `y<x`. Let `r=x/y`, a power of two at least two. Then

`n_0=x(yq-1)`,

`n_1=F(n_0)=y[(x-1)q+r]`.

The bracket is odd, so the exact dyadic core of `n_1` is `y`. Applying the corresponding branch gives the exact identity

`n_2-n_0=-(q-2)(x+y-1)+1-r`.                          (4.1)

If `q>=3`, the right-hand side is negative. Thus descent occurs at step two.

The only unresolved possibility is `q=1`, giving

`n_0=x(y-1)`, with `y>=4` because `p>1`.

In particular, any hypothetical nonconstant cycle that avoids cofactor one would have a minimum with this Mersenne-shaped odd part. This is the global entry point that the local size/core tradeoff by itself lacked.

## 5 The four return branches

Keep `q=1`, `y>=4`, and `r=x/y>=2`. The first two iterates simplify to

`n_0=r y(y-1)`,

`n_1=r y^2+r y-y`,

`n_2=r y^2-r+y`.                                      (5.1)

There are four dyadic possibilities.

### The quotient is smaller than the second core

Suppose `r<y`. Then `r<=y/2`, and `y/r` is an even integer. The quotient

`n_2/r=y^2+y/r-1`

is odd, so the exact core of `n_2` is `r`. Consequently

`n_3-n_0=-y^2+ry+y-y/r+r`.

It satisfies

`n_3-n_0 <= -(y^2-3y+4)/2 < 0`.

For completeness, the difference between the upper bound and the preceding exact expression is

`(y-2r)[(y+1)/2+1/r]>=0`,

and `y^2-3y+4=(y-4)(y+1)+8>0` for `y>=4`.

### The quotient is at least four times the second core

Suppose `r>=4y`, and write `R=r/y>=4`. Since

`n_2/y=ry-R+1`

is odd, `n_2` has exact core `y`. Therefore

`n_3-n_0=-R(y-1)+3y-2`

`          =2-y-(R-4)(y-1)<0`.

### The quotient is twice the second core

Suppose `r=2y`. The same exact-core calculation gives

`n_3=n_0+y=y(2y^2-2y+1)`.

Its bracket is odd, so `n_3` also has exact core `y`. One further step gives

`n_4-n_0=-2y^2+5y-2=-(2y-1)(y-2)<0`.

### The quotient equals the second core

Suppose `r=y`. Equation (5.1) gives the power of two

`n_2=y^3`, while `n_0=y^3-y^2`.

For `y>=8`, applying (3.1) to this power yields

`n_5=3y^3/4+4`,

`4(n_0-n_5)=y^3-4y^2-16`

`            =(y-8)(y^2+4y+32)+240>0`.

The sole smaller possibility, `y=4`, gives the exact six-step descent

`48 -> 76 -> 64 -> 190 -> 98 -> 52 -> 46 < 48`.

Because `r` and `y` are powers of two, the cases `r<y`, `r=y`, `r=2y`, and `r>=4y` exhaust all possibilities. This completes the descent lemma.

## 6 Global conclusion and coarse quantitative bounds

Apply strong induction on the positive even integer `n`. A fixed point or a member of `E` has the claimed eventual behavior. Every other state reaches a smaller positive even state within at most six steps, and that smaller state has the claimed behavior by induction. The theorem follows. In particular, a periodic orbit must be one of the fixed points or the two explicitly displayed cycles.

The same induction gives a coarse effective bound: first entry into a fixed point or one of the two cycles occurs within at most `3n` iterates. Indeed, a strict decrease between descent anchors is at least two, so a block of at most six steps costs at most `3(n-m)` when it lands at `m<n`. The explicit exceptional starts satisfy the bound directly.

Also `F(m)<=3m-2<3m`. Between successive descent anchors there are at most six steps, all anchors are at most the initial `n`, and the final nonconstant cycles have maximum `46`. Thus every orbit value is at most `3^6 n=729n`. These bounds are deliberately loose; their role is to make the global termination statement explicit.

## 7 Consequence for genuine aliquot dynamics

**Corollary.** Every finite aliquot cycle confined to even integers with prime odd part is constant. Such a constant member has the form

`2^a(2^(a+1)-1)`, with `a>=1` and `2^(a+1)-1` prime.

For these members the true aliquot map agrees with `F` on every edge. The two nonconstant `F` cycles are unavailable because they contain odd cofactor one. In fact, any genuine aliquot trajectory that stays in the stated prime-odd-part class forever must eventually reach a perfect number.

The displayed surrogate cycles are not genuine aliquot cycles: `s(8)=7` and `s(16)=15`. The theorem does not establish convergence or a cycle classification for general aliquot dynamics, for odd members, or for members with composite odd part. It does settle the precise all-length prime-cofactor chamber target, including the earlier length-three target.

## 8 Exact certificates

The dependency-free companion script checks polynomial identities in an exact integer polynomial ring, verifies every displayed exceptional orbit, and evaluates the descent classifier on even inputs through `4096`. That finite branch check lies inside the predecessor's existing window and is not a new census. The unbounded proof is the case partition and the factored inequalities above.

Checks use explicit exceptions and remain active under Python `-O`. The normal and optimized JSON outputs are identical. Earlier notes, scripts, and data are preserved unchanged.
