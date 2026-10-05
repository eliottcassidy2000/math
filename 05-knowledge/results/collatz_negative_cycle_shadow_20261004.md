# The -5 shadow defeats every bounded-depth checkpoint bank

2026-10-04. **PROVED:** a whole mixed congruence cell excludes every
bounded-depth smaller-child join whose affine intercept is greater than -1;
the resulting infinite progression of missed powers of 3; and the exact
valuation budget for repeated `(1,2)` blocks. **FINITE-EXACT:** the controls
below. Universal coverage remains **OPEN**. This is a scoped obstruction,
not a claim about all affine joins or a literature-priority claim.

Artifacts: [program](../../04-computation/experiments/collatz_negative_cycle_shadow_20261004.py)
and [saved output](collatz_negative_cycle_shadow_20261004.out).

## 1. Recovery and the extra coordinate

The [earlier whole-cell obstruction](collatz_partition_cover_20261004.md)
uses the all-ones shadow at -1 and imposes no intercept restriction. Its
deep binary cells do not contain odd powers of 3. The
[later-checkpoint compiler](collatz_checkpoint_reroute_20261004.md) supplies
the extra fact used here: every candidate child has

    h = lambda*n+b,      b=beta-1 > -1.

That compiler successfully handles an infinite power-of-three family; the
question here is whether bounded source and child depths could handle them
all. The [signed rational-anchor compiler](rational_anchor_returns_20261004.md)
already keeps actual words and denominator guards on negative shadows.
[THM-4512, coefficient-descent classes](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md)
supplies the inherited distinction between a formal integral endpoint and
an exact odd-word cylinder.

The canonical hostile to the first attempted proof is
`-11/3 --(1)--> -5`: a negative rational strictly between -5 and zero can
enter the -5 cycle. One cannot silently make the formal child an integer.
The repaired proof uses both the ternary source guard and the intercept
bound. The least-used sidecar is this affine intercept, not only the slope.

The concept board is **signed shadow / exact words / two depths / ternary
denominator / intercept / smaller obligation**. Anchor: a limitation on a
precise coverage grammar. Niche: its arithmetic intersection with powers of
3. Wildcard: connect the shadow's repeated word to an exact valuation budget.

The source is a checked positive common-future diagram; the target is its
evaluation at the signed source -5. The map preserves its affine identity
and, under the binary guard, its exact valuation words. It loses positive
sign and may lose integrality. The ternary guard and the bound `b>-1`
restore enough information for the obstruction. No convergence assertion
is transported from the negative cycle.

## 2. The whole-cell theorem

Use `U(x)=oddpart(3x+1)` on positive odd integers, retaining `U(1)=1` for
arithmetic bookkeeping. For a positive valuation word w of length r and
total cost A, write

    F_w(x)=(3^r*x+B_w)/2^A.

The empty word has length and cost zero and carry zero. A word is *actual*
only when each listed exponent is the exact valuation at that step.

**PROVED.** Fix integers `R,S>=0` and put `K=floor(3R/2)+1`. Suppose a
positive odd integer n satisfies

    n = -5 mod2^K,            n = 0 mod3^S.                    (1)

There is no positive odd h<n with actual valuation words w,v of lengths
`r<=R, s<=S`, respectively, such that `F_w(n)=F_v(h)` and the induced affine
join has intercept `b>-1`:

    h=lambda*n+b,
    lambda=3^(r-s)*2^(D-A),
    b=(B_w*2^D-B_v*2^A)/(2^A*3^s).                           (2)

Here D is the total cost of v. Individual child exponents and D are
unrestricted. Both empty-word cases are included. No source-size threshold,
root certificate, or assumption about other negative basins is needed.

The later-checkpoint compiler satisfies `b>-1` by its positive beta formula.
For its common-future diagram, the source depth includes the next source
edge: `r_checkpoint+1`; the child depth is `ell+1`. The bounds in this
theorem refer to these full words, not just the displayed checkpoint.

### Proof

The signed cycle is `-5 --(1)--> -7 --(2)--> -5`. Its first r letters have
cost `A_r=floor(3r/2)`. Condition (1) fixes the first R actual source letters
to this alternating word. Thus `A=A_r` and

    J=F_w(-5) is either -5 or -7.

Because `3^S|n`, the number `lambda*n` is integral at the prime 3. Since h
is integral, b has no factor 3 in its reduced denominator. Formula (2)
then gives `b in Z[1/2]`. Set

    y=-5*lambda+b.

The difference `h-y=lambda*(n+5)` has 2-adic valuation at least `D+1`.
Therefore y is an odd rational with the same exact child word v as h, and
`F_v(y)=J`. Every positive valuation map preserves positivity, including
from zero after its first step; hence `y<0`. If `lambda>=1`, then
`h-n=(lambda-1)n+b>-1`, impossible for a negative integer. Consequently

    lambda<1,             -6<y<0.                            (3)

**Case r>=s.** Both lambda and b are dyadic, so y is dyadic. Since y is
odd 2-adically, it is an odd integer. The only possibilities in (3) are
-5,-3,-1. The latter two enter the -1 fixed point, not J; hence y=-5.
Now both words start at -5 and finish at the same cycle phase. Their
length difference is nonnegative and even, and

    lambda=(9/8)^((r-s)/2)>=1,

contradicting (3). This includes `r=s=0` and the equality slope one.

**Case r<s.** Put `d=s-r>0`, and let p be the first d letters of v, of
total cost C. The number y has exact denominator `3^d` at the prime 3:
its term `-5*lambda` does, whereas b is dyadic. Thus `z=F_p(y)` is dyadic.
It is also odd 2-adically by the actual-word guard, so z is an odd integer.
The remaining r child steps take it to J. It must be negative, and it
cannot be -1 or -3. Therefore `z<=-5`.

Write

    j=D-A-C,       z=-5*2^j+zeta,       zeta=F_p(b).            (4)

Each map `(3x+1)/2^a`, a>=1, preserves the interval `(-1,infinity)`;
therefore `zeta>-1`. If j<0 then `z>-7/2`, a contradiction. If j=0,
the odd integer z lies above -6, so z=-5 and zeta=0. This gives
`b=-B_p/3^d<0`. Since b is also dyadic it is an integer, contradicting
`b>-1`. Hence j>=1.

Now `zeta=z+5*2^j` is integral, so
`b=(2^C*zeta-B_p)/3^d` is integral too. As y is odd and
`v2(lambda)=C+j>=1`, b is odd. From (3) and `b>-1` we get

    1<=b<5*lambda<5;       thus b is 1 or 3.

Moreover `v2(y-b)=C+j>=C+1`. The actual first d child letters at y are
therefore also the actual first d letters at the positive integer b.
For b=1 their total cost is `C=2d`. For b=3 it is `C=1` when d=1,
and `C=2d+1` when d>=2, from `3 --(1)-->5 --(4)-->1` followed by root
loops. In all cases

    lambda=2^(C+j)/3^d>1,

the final contradiction. Root padding here only enlarges the allowed
arithmetic words; it does not export a padded first-hit certificate.

The two moduli in (1) are coprime. This obstruction cell has exact natural
density `1/(2^K*3^S)` among integers and `1/(2^(K-1)*3^S)` among odd
integers. These are densities of failure of the stated join model, not
of nonconvergence or failure of every possible proof method.

## 3. Every finite checkpoint bank misses infinitely many residual powers

For K>=5, the elementary identity

    v2(3^(2^j)-1)=j+2,       j>=1,

shows that `<3^8>` modulo `2^K` is exactly the subgroup of residues one
modulo32. Since `3^3=-5 mod32`, there is a unique exponent class

    a=a_K mod2^(K-2),       a_K=3 mod8,
    3^a=-5 mod2^K.                                           (5)

Given uniform full-word depth bounds R,S, take
`K=max(5,floor(3R/2)+1)`. Every exponent in (5) with a>=S makes `n=3^a`
satisfy (1), so every such power escapes the entire bounded-depth bank.
This applies even to infinitely many valuation-cost choices or additional
source filters when their two word depths stay bounded and `b>-1`.
A finite bank of fixed rows is a special case.

Some exact phases are:

| K | exponent a | exponent period |
|---:|---:|---:|
|5|3|8|
|7|11|32|
|10|11|256|
|12|267|1024|
|16|15627|16384|
|18|15627|65536|

The successful checkpoint family `a=483 mod1024` has a different phase;
there is no conflict. This result leaves depth-growing controllers and
changed join grammars available. It does not rule out rows with `b<=-1`;
no sharpness claim at that boundary is made. An infinite progression of
exponents is not a positive-density set of integers `n=3^a`.

## 4. Exact repeated-block budget and the first failed shortcut

The affine block identity is

    F_(1,2)(n)+5=(9/8)*(n+5).                                (6)

**PROVED.** For every positive odd n, put `H=v2(n+5)>=1`. The number of
complete actual `(1,2)` blocks before the first mismatch is exactly

    floor((H-1)/3).                                         (7)

A complete actual block requires `n=-5 mod16`, or H>=4; after it, (6)
decreases H by exactly three. Once H is 1,2, or3, that guard fails, proving
(7). Merely requiring a formal integral endpoint would incorrectly give
`floor(H/3)`. At n=3, H=3 and the formal block endpoint is4, but the actual
word begins `(1,4)`, not `(1,2)`.

The map `n -> n+5` makes the block multiplicative and exposes its consumed
binary precision. It does not prove that the orbit returns to this chart
after a mismatch or give a global decreasing rank. The signed-anchor
compiler and this obstruction use the same exact shadow for different
purposes: a controlled exit can certify a family, while insufficient depth
cannot provide every smaller-child join inside its shadow.

## 5. Reproducible checks and limits

Run either command; their standard output agrees exactly:

```text
python 04-computation/experiments/collatz_negative_cycle_shadow_20261004.py
python -O 04-computation/experiments/collatz_negative_cycle_shadow_20261004.py
```

The script uses explicit exceptions, not disabled-under-optimization
assertions, and imports no existing implementation or home oracle.

* For all `r,s=0..5`, it enumerates every positive child word satisfying
  the necessary slope inequality `2^D*3^r<2^A*3^s`. This proves the finite
  total-cost cutoff; there is no independent cap on a valuation letter.
  The 190 words comprise162 with r<s,13 with r=s,15 with r>s. All signed
  endpoint words replay exactly. None has a dyadic intercept, so this
  particular finite universe does not test sharpness of `b>-1`.
* An independent literal graph enumerates odd starts up to21627 through
  four steps. It tests74 source instances: the first three positive cell
  members for each `R,S=0..4`, omitting1. It finds no smaller joins at all
  in that finite universe. The stronger finite observation is not promoted
  to an unrestricted-intercept theorem.
* Binary exponent lifting is checked at K=5..18, with56 progression
  controls and independent complete phase searches at K=5..12. No huge
  integer `3^a` is materialized. Twenty literal sources verify (7).
* Sixteen small-positive-prefix controls check the b=1 and b=3 cost contradiction.
  Twenty members of the positive family `155+4096t ->111+2916t` verify its
  actual common future, and direct modular exponentiation verifies its
  power phase483 as well as the original discovery phase27107.
* Boundary controls retain the rational hostile `-11/3 -> -5`, the
  illegal formal `(1,2)` endpoint at3, the legal predecessor `3 ->5`
  outside the ternary guard, and the descent `9 ->7` outside the binary
  guard. Six invalid exact domains are rejected.

Independent proof audit: the repaired r<s argument, the equality case
y=-5 for r>=s, and the finite-bank scope were checked by a second agent
before recording this result. The mathematical theorem is proved above;
finite replay is additional evidence with its own stated universe.
