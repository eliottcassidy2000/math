# A Frobenius boundary: more symmetry, larger recurrence order

**PROVED elementary formulas below; FINITE-EXACT literal tournament actions
in frobenius.py.** Groups are specified subgroups of automorphisms; no full
automorphism census at 27 vertices is claimed.

Work in F_27=F_3[t]/(t^3-t-1). The cubic has no root in F_3 and is irreducible.
Let H be the 13 nonzero squares. The Paley tournament has x->y iff y-x is
in H. Since -1 is a nonsquare this orients every distinct pair. Consider

    G_aff={x -> ax+b : a in H, b in F_27},                |G_aff|=351,
    G_sem={x -> a x^(3^j)+b : a in H, b in F_27, j=0,1,2}, |G_sem|=1053.

The Frobenius x->x^3 preserves H and has order three. frobenius.py constructs
all these permutations and checks every tournament adjacency, not merely
their predicted cycle counts.

## S1. Root laws and all-depth recurrences

The fixed-point histograms are

| number fixed | affine multiplicity | semilinear multiplicity |
|---|---:|---:|
| 0 | 26 | 494 |
| 1 | 324 | 324 |
| 3 | 0 | 234 |
| 27 | 1 | 1 |

For j=1,2, the F_3-linear map x -> a x^(3^j)-x has kernel size three:
its nonzero roots solve x^2=a^(-1) or x^8=a^(-1), respectively. The
exponents have gcd two with 26, and a is a square. Its image has size nine,
so nine translations give three fixed points and eighteen give zero.
This proves the table without relying on enumeration.

For the iterated wreath action of G_sem, write F_k(z)=E exp(z X_k), with
F_0=exp(z). Then

    F_(k+1)=(494+324F_k+234F_k^3+F_k^27)/1053.            (S1)

This is another exact full recursive relation on actual tournament actions.
The root has mean one, variance two and third factorial moment eighteen.
Thus its ordered port counts are

    a_(k,1)=1,
    a_(k,2)=2k+1,
    a_(k,3)=6k^2+18k+1.

The affine subgroup instead has offspring polynomial
(26+324s+s^27)/351 and third port count 6k^2+50k+1. At depth one the
three-port counts are 25 and 57, although both first and second moments
agree. Criticality and variance alone miss the Frobenius action.

## S2. Three-generator gaps and exact linear memory

Root fixed support is {0,1,3,27}, with positive masses. As in M2, at depth k
the support is the sum of N=27^(k-1) copies of this set: smaller root sums
can always be padded by zero. Its reflection at 27^k is the set of sums
of at most N generators from {24,26,27}.

For exactly l nonzero generators, all attainable sums are

    [24l,27l] minus {24l+1},                l>=1.

Indeed subtracting 24l leaves a sum of l terms from {0,2,3}, which realizes
0 and every integer from 2 through 3l. For N>=9 the union over l<=N is
<24,26,27> intersect [0,27N]. All integers from 218 to 27N occur: successive
intervals overlap, and from l=10 each interval's exceptional point was
already covered by the previous interval. The last semigroup gap is 217.

There are 109 gaps: the eight between-interval gaps have lengths
23,20,17,14,11,8,5,2, totaling 100; the remaining nine are 24l+1 for
l=1,...,9. M1 therefore gives the all-depth minimal arity orders

    d_sem(1)=4,
    d_sem(k)=27^k+1-109,                   k>=2.          (S2)

The affine orders are 3 at depth one and 27^k+1-325 thereafter. At depth
two these are respectively 621 and 405. Larger symmetry makes each orbit
count no larger, but introduces additional fixed-count frequencies into
the moment sequence and increases its minimal linear recurrence order.
This comparison concerns the entire sequence in arity, not the number of
orbits at one arity. The script checks both support recursions through
depth three and the independent coin-reachability gaps.

## S3. Why the prime multiplier compression cannot simply be reused

The negation centralizer at the semilinear root is

    H semidirect <Frob>,           of order 39,
    Frob M_a Frob^(-1)=M_(a^3).

It is nonabelian: multiplication by a and Frobenius do not commute for
a!=1 in H. Their composition carries the actual phase rule

    (a,j)(b,l)=(a b^(3^j), j+l mod 3).

An unmodified product of local multipliers is therefore no longer a
homomorphism. The cyclic group-algebra charges of FOURIER do not encode
this root action. Its fixed histogram is 12 elements fixing one point,
26 fixing three, and the identity fixing 27, again checked literally.

At depth one the square-root marking connection is still flat. Negation
conjugates translations to their inverses and fixes the entire semilinear
linear group. Every reversal involution consequently has the form
x -> -x+b; the product of two such markings is a translation, and its unique
square root is translation by half the displacement. These translations
compose additively around every triangle, with zero error. Thus Hol_1=1,
despite the nonabelian centralizer of order 39. At higher depth a nonabelian
charge or representation carrier is an OPEN target, not supplied by F1.

This is a literal finite-field transport connection. It is separate from
binary normal-basis rotation on F_64 and F_256: those earlier controls remain
valid in their own decoration action, but they do not determine this odd
field's semilinear group or its port recurrence.
