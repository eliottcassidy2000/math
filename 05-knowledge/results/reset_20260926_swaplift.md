# Recreated precision can certify a selected descent

**PROVED elementary lifting and clipping theorems / FINITE-EXACT bank
q odd1..341 / OPEN Collatz.** No novelty claim for finite-word descent
certificates. The result applies uniformly in the unbounded actual source,
not merely to a tested numerical range.

For every hard swap whose coefficient q is odd and at most341, recreated
precision R>=64 guarantees an actual smaller odd iterate within41 U
steps, including the swap itself. The per-q thresholds are sharper.
The [full certificate table](reset_20260926_swaplift.out) records them.
The table also certifies65 complete dyadic residue classes, after1431
exact checks of their small-source exceptions.

## 1. Inheritance and the remaining hard swap

The closest mechanism is the exact consume/swap/collision controller in
[the precision note](nextforest_20260926_precision.md). Its hostile is
precision recreation; the corrected near miss is treating valuation
subtraction as a Euclidean rank. The least-used sidecar is the actual
small core3q created by a swap. The board is q, the recreated precision,
an exact core prefix, its final affine coefficient, collision clipping,
and the original source to which descent is compared.

Write U(n)=oddpart(3n+1). For a state `y=q*2^r+u`, with q,u positive
odd and `a=v2(3u+1)>r`, the next iterate divides by exactly2^r.
If r>=2 then `U(y)=(3y+1)/2^r<y`; this case already descends.
The hard case is r=1. Write a=R+1 and `3u+1=2^(R+1)v`. Then

    y=(2^(R+1)v+6q-1)/3,
    z=U(y)=v*2^R+3q,                                (S1)

where q,v are positive odd, R>=1, and

    2^(R+1)v = 1 mod3.                              (S2)

These conditions are sufficient as well as necessary for this typed
swap: u is a positive odd integer, a=R+1, and the odd number in(S1)
is the actual successor of y. The coefficient v is unbounded.

## 2. A verified core crossing is enough

Take an actual core prefix

    c_0=3q, c_i=U(c_(i-1)), i=1,...,J,
    A_i=sum_(h=1)^i v2(3c_(h-1)+1),
    A=A_J, P=A_(J-1), b=c_J,

and require only

    b<2q.                                           (S3)

No convergence theorem is assumed for arbitrary q. A finite verified
prefix is the input certificate. Reaching1 supplies one, but often uses
unnecessarily many steps.

Composition gives the integer affine identity

    2^A b=3^(J+1)q+C,   C>0.                        (S4)

Hence(S3) automatically implies

    lambda=3^(J+1)/2^(A+1)<1.                        (S5)

In particular a core prefix already reaching1 needs no appended1-loops:
then lambda<1/(2q). This removes a redundant continuation from the
initial proposed proof target.

If R>A, the core valuations are preserved throughout the first J steps
after the swap, so

    U^J(z)=3^J v*2^(R-A)+b.                         (S6)

The preservation proof is inductive: before each core division, the
high summand still has strictly larger2-adic valuation than the core
numerator. The smaller valuation therefore belongs to the actual source,
not just the formally substituted core.

Comparing(S6) with the original y, rather than with z, gives exactly

    U^J(z)<y
      iff 2^(R+1)v(1-lambda)>3b-6q+1.               (S7)

The left side is positive and the right side is at most-2. Thus every
admissible v gives descent after J+1 actual U steps.

**The equality boundary R=A is also valid.** All preceding core steps
still preserve their valuations because A_i<A for i<J. At the final
step the two odd terms meet, giving

    U^J(z)=oddpart(3^J v+b) < 3^J v+b.              (S8)

The formal sum on the right already lies below y by(S7) at R=A.
Extra division at the collision only strengthens descent. For example,
q=1,J=2,A=5,R=5,v=1 gives23->35->53->5; its formal final bound is10.

**Simple theorem.** A certificate(S3) proves selected descent within
J+1 steps for every hard swap(S1)--(S2) with R>=A.

## 3. Clip only the final core division

Suppose instead

    P<R<A.

All but the last core step still shadow exactly. The last one swaps
the high and low roles, yielding the exact odd endpoint

    U^J(z)=3^J v+2^(A-R)b.                          (S9)

If `2^(A-R)b<2q`, applying(S4) after multiplying by2^(-R) proves
`3^(J+1)<2^(R+1)`. Both the high coefficient and the constant term in
the comparison with y therefore have the required sign. A sufficient
uniform threshold is

    R>=max(P+1, A-floor(log2(floor((2q-1)/b)))).      (S10)

The floor is defined because(S3) makes its argument at least1. This
includes the simple theorem when R>=A.

The final comparison can be optimized without requiring each term to
have its own sign. For P<R<=A define exact integers

    D_R=2^(R+1)-3^(J+1),
    B_R=3*2^(A-R)b-6q+1.                            (S11)

For R<A, equation(S9) gives descent iff `D_R v>B_R`.
The sufficient test

    D_R>max(0,B_R)                                  (S12)

works for every v>=1, so it respects the stronger admissibility(S2)
without needing a real approximation. D_R increases and B_R decreases
with R. A passing R therefore certifies every larger R through A;
(S6)--(S8) certify the rest. There is always a passing value at A by
(S3)--(S5).

Let R_* be the smallest integer in[P+1,A] passing(S12). This is the
optimized threshold in this particular sufficient scheme. No claim says
that smaller R cannot descend by another word or certificate. The four
extra improvements beyond(S10) in the finite bank are:

| q | J | A | P | b | threshold(S10) | R_* |
|---:|---:|---:|---:|---:|---:|---:|
|1|2|5|1|1|5|4|
|5|4|8|3|5|8|7|
|9|40|66|61|5|65|64|
|61|28|46|42|61|46|45|

The full-output affine carries and chronology, not a finite colour, make
these inequalities valid at one fixed source.

## 4. The finite certificate bank and its all-height consequence

For each odd q in1..341, the experiment directly follows3q until it
first falls below2q. Every one of the171 scans succeeds. The largest
J is40 and the largest A is66, both at q=9, whose core27 reaches5.
The largest R_* is64, again at q=9. Thus the headline R>=64 and
J+1<=41 result follows uniformly for all admissible v.

The first-crossing choice minimizes R_* among later prefixes satisfying
(S3): every later prefix has P at least the old A, hence its permitted
threshold is at least old A+1, above the old R_*. This is optimal only
within the final-division clipping scheme; it is not optimality among
all possible Collatz certificates.

For comparison, following every core all the way to1 gives largest
J=56,A=99, at q=339. This extra work is unnecessary for the selected
descent. The source set is not bounded by341: q is the coefficient,
while v and R, and therefore y, have no upper bound.

For a given row of the certificate table, put K=R_*+1 and

    rho_q=(6q-1)*3^(-1) mod2^K.

Every positive integer n>2q in the residue class

    n=rho_q mod2^K                                 (S13)

has `a=v2(3n-6q+1)>=K` and R=a-1>=R_*. Its positive odd v and u=n-2q
give precisely(S1), so the theorem proves `U^(J+1)(n)<n`.
The finitely many positive class members n<=2q are checked separately.
All1431 such(q,n) pairs pass the same J+1-step inequality. The result
therefore covers the *entire positive residue classes*, not merely tails.

Dyadic residue cylinders are either nested or disjoint. Removing classes
contained in earlier shorter classes leaves65 pairwise disjoint classes.
Their exact natural density is

    sum 2^(-K)
      =6985206796614369409/36893488147419103232
      =0.18933440960374526... .                      (S14)

All are odd, so their relative density among odd integers is twice this
value, about37.8669%. The density statement concerns this explicit finite
union only. It supplies no conclusion for an uncovered fixed source and
is not a claim to improve classical global density results.

## 5. An unbounded coefficient family with fixed required precision

There is a useful infinite positive control beyond the finite bank:

    q_j=(64^j-1)/9, j>=1.

It is a positive odd integer and its core3q_j reaches1 in one U step,
with A=6j. Nevertheless(S10) needs only R>=3: indeed
`floor(log2(2q_j-1))=6j-3`. Every admissible hard swap in this family
therefore descends in two U steps at precision at least3, uniformly in
both j and v.

The thresholds1 and2 actually fail two-step descent here. When R<6j,
the exact endpoint is `3v+2^(6j-R)`. Three times its difference from y is

    (9-2^(R+1))v+3*2^(6j-R)-6q_j+1.                (S15)

For R=1,2 both the v coefficient and constant term are positive.
For R=3 they are negative (j>=1), and the already proved threshold
argument handles all larger R, including the final collision.

Since q_j=7 mod8, the corresponding precision>=3 cylinder is simply
`n=3 mod16`, apart from the initially excluded finite n<=2q_j.
Using q_1=7 and checking n=3 recovers that entire familiar two-step
descent class. The unbounded q-family explains the mechanism; it does
not add new coverage beyond that same residue class.

## 6. Fixed small precision still allows arbitrarily long actual growth

Fix q=1,R=1. This covers exactly all positive n=3 mod8, by
`n=(4v+5)/3` with positive odd v=1 mod3. A short core proof for3
does not give a uniform descent horizon here:27 is one member and
`U^3(27)=47>27`, even though the core3 reaches1 in two steps.

There is an all-height obstruction to *any* such uniform horizon, now
including a family that starts at27. For c in{2,4}, put

    n_(c,k)=c*8^k-5, k>=1.

These all have q=1,R=1, since
`v2(3(n_(c,k)-2)+1)=2`. Direct induction proves, for0<=j<=k,

    U^(2j)(n_(c,k))=c*9^j*8^(k-j)-5,                (S16)

and for0<=j<k,

    U^(2j+1)(n_(c,k))=(3c/2)*9^j*8^(k-j)-7.         (S17)

At each pair the valuations are exactly1 and2: the first value is3
mod8, and the intermediate is1 mod8 with its numerator divisible by4
but not8. The pair multiplies n+5 by9/8. Hence every one of the first
2k iterates is strictly greater than n_(c,k). With first-descent time
allowed to be infinity, this proves `tau(n_(c,k))>2k`; it does not prove
that tau is finite for every k. No bounded first-descent horizon exists
over this single fixed(q,R) fibre.

The c=4 family is

    27,251,2043,16379,131067,1048571,8388603,67108859,...,
    n_(4,k+1)=8*n_(4,k)+35.                          (S18)

Its first-descent times are **provably unbounded**, unlike the earlier
constant-defect family `(4^k+17)/3`, whose large members have first
descent6. This is a stronger answer to the search for families beyond27:
the mechanism is the exact two-step expanding word, not an extrapolation
from a few long computed orbits.

Exact finite computations for the first eight members give:

| k | n_(4,k) | first-descent odd steps | first smaller endpoint |
|---:|---:|---:|---:|
|1|27|37|23|
|2|251|17|61|
|3|2043|11|691|
|4|16379|34|7583|
|5|131067|31|35953|
|6|1048571|27|909041|
|7|8388603|22|3830695|
|8|67108859|39|29486219|

The times are not monotone. The theorem is the lower bound greater
than2k, not a proposed exact stopping-time formula. The companion output
also retains the first eight times for c=2, beginning11,123,1019,... .
For X>=27, the exact number of c=4 family members at most X is

    floor(log_8((X+5)/4)).                           (S19)

This is a logarithmic, density-zero family. Its recursion magnifies the
same negative-cycle address; it supplies no natural-density estimate
for all long Collatz trajectories.

These are positive finite shadows of the negative odd cycle
`-5->-7->-5`. The positive source varies with k; no divergent positive
orbit is constructed. The unbounded coefficient v remains essential.
This is why large recreated precision can be paid by the theorem while
fixed small precision cannot be dismissed using the same finite core.

## 7. Verification, transfer and remaining target

The source is one hard-swap integer(S1); the target is its actual
selected U iterate. The map is the exact affine lift of a *verified*
core prefix, retaining both strict valuation guards and its last possible
collision. It preserves the endpoint and the comparison with the original
source. A core convergence statement alone loses the available precision;
the threshold R_* repairs precisely that coordinate. The fixed(q,R) family
(S16) is the decisive hostile below the certified threshold.

Run [the independent experiment](../../04-computation/experiments/reset_20260926_swaplift.py)
normally and under `-O`. The [full output](reset_20260926_swaplift.out)
contains all171 rows `(q,J,A,P,b,clipped_R,R_*,residue,K)`, the65 selected
class labels,21,546 elementary hard-swap checks,10,955 exact selected-source
checks including2311 final collisions and1716 clipped swaps,1431 small
exception pairs, and18,932 directly checked class members below100000.
It also checks the unbounded-q formulas through j=64, both positive-shadow
families through k=128, and the first eight descent times of each family
with an explicit cap10000.
Universal scope comes from the proofs above, not these ranges.

The remaining target is an effective selected descent across unbounded q
and insufficient precision, or a proof that an actual unresolved source
must enter a certified region. The finite bank and its positive density
do not prove that entry statement. Collatz remains open.
