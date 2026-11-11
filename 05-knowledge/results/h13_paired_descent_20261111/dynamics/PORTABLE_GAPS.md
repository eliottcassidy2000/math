# Portable rooted switches of every positive even width

**PROVED IN NOTE.** The existing two-flank H13 count difference is a pure-H
rooted difference for every pair of addresses, not only on the A=0 unit-gap
axis. Its complete source sidecars are explicit. The compiler and the H13
signatures are inherited from the
[golden switch note](../../five_core_web_20261111/history/continuation/PROOF.md).

Use core Q=(18,8,3,5,14), lower cut6 and upper cut3. Let α=A mod2,
β=B mod2, b=bitlength(B), and r=rev_b(B), with b=r=0 at B=0. The lower
and upper total counts are

    H− = 13(g(A)+g(B)−1)+104AB+47A+9α+32B+6β+18,
    H+ = H−+3A−α+2.

Their rooted states have exactly the same three sidecars:

    E = 7 if B=0, and 13 if B>0,
    U = 2^b(104A+45)−13r,
    d = 13.                                           (1)

**Proof.** Deleting the root leaves the ordinal sum of a transitive left
flank, Q, and a transitive right flank; its count is13. If the right flank
is nonempty, its last vertex beats the root, so all13 full deletion paths
can end there and then enter the root. If it is empty, the full paths end
in Q, whose selected ending sum is6 for either cut, giving E=13−6=7.

To compute U, the root-alone path contributes13. For a marked path starting
in the left flank, the weighted core splitting sum is d+k+=13+91=104.
The left selected-start weights sum to A, and each right vertex has two
path-bin choices, giving104 A 2^b. The selected core-start sum is32 for
either cut and contributes32·2^b. Starts in the right flank contribute
13(2^b−1−r), using the compiler's complementary reversed right-word
convention. Adding the four disjoint cases gives (1). QED.

Consequently Δ=(H+−H−,0,0,0) with width

    w(A)=3A−α+2 = 6q+2 when A=2q,
                   6q+4 when A=2q+1.                 (2)

As A ranges over the nonnegative integers, this is a bijection onto the
positive even integers not divisible by6. Every width in that set is
available at every right address B and survives every common O/Z word.
The width is an actual rooted response difference, not merely a difference
of two marginal count formulas.

**Corollary: every positive even width has actual pure-H paired sources.**
Write the desired width as `w=3^a v`, where v is even and not divisible
by3. Choose A by (2) so w(A)=v, and take any right address B. Prepend the
ordinal sum of a directed triangles to both twins, retaining the old root
in the right factor. Every rooted coordinate is multiplied by3^a. For U,
the path starting at the root cannot enter the preceding factor, so that
entire factor lies in the other path and contributes exactly its Hamilton
count3^a. The same product rule for H,E,d is immediate. The resulting
actual pair has rooted difference `(w,0,0,0)`. The case a=0 uses no prefix.
In particular divisibility by6 is a boundary of the unscaled H13 menu,
not a boundary of realizable portable differences.

## Root changes turn the difference into a real construction response

After one common J, the difference becomes

    w(A)·(1,1,1,1)

in `(H,E,U,d)` coordinates. This is w(A) times the state of a rooted
singleton K1, with H(empty)=1. Thus for every subsequent O/Z/J word W,

    H(W J T+)−H(W J T−) = w(A) H(W K1).              (3)

Here WJ denotes execution of J first, then W; both sides are counts of
literal graph constructions. Equation (3) follows from the inherited
linear physical maps. It is not an assertion that every tournament is
constructible from K1 by O/Z/J.

This produces the same abstract gap multiplication as ordinal join with
an arbitrary tournament, but in a different physical setting: the two
sources evolve under a common program, and their baseline is its actual
rooted response rather than the ordinal product of two scalar counts.

For repeated J, write a_n=H(J^n K1). The exact response matrix on
`(h,u,t)` in the common frame is

    (h,u,t) -> (t,h+u,h+2u).

Its characteristic polynomial is λ³−λ²−λ−1. Therefore

    a_0=1, a_1=1, a_2=3,
    a_(n+3)=a_(n+2)+a_(n+1)+a_n.

The sequence begins1,1,3,5,9,17,31,57. Thus the tribonacci recurrence
reappears here as the exact transport of a selectable H13 difference.
It is an operator calculation, not a spectral analogy alone. The existing
golden P4 comparison still concerns triangle-overlap spectra; the rooted
response supplies the information needed for these count recurrences.

## The count-restoration gate

For the new H13-to-H5 suffix rewrite, the arbitrary-address debt is
`8f(B)+2β`, with f=B+1−β−g(B). Its bit recursion is

    f(2q)=f(q),
    f(2q+1)=f(q)+q+(q mod2).

The resulting odd-step debt increment is `8q+6(q mod2)+2`. It belongs to
the width set (2) exactly when q is0 or1 modulo3. At q=2 modulo3 it is
a multiple of6. The corollary realizes that width too, by an actual ordinal
prefix, so no independent-switch addition is needed merely to realize the
required difference somewhere. But the resulting pair has its own absolute
baseline and sidecars. There is still no demonstrated operation installing
that difference into the specific source whose debt is being repaired.
That compatibility is the remaining construction question exposed by the
binary recursion. Realizing every possible gap does not prove that every
possible absolute count occurs.

Controls in `portable_probe.py` use literal subset DP for both rooted
states at A=0..7,B=0..3 and compare the independent J iteration with the
tribonacci recurrence. The all-address proof is the splitting argument.
