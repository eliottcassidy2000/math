# Fixed-gap centers are equal-area triangle hypotenuses and doubled elliptic x-coordinates

Status: PROVED elementary equivalence; FINITE-EXACT checks; no literature novelty claim and no solution or global exclusion for the magic-square-of-squares problem.

## 1. Exact center and magic-square interfaces

Fix a positive rational d and set

C_d = {c in Q_{>0}: c-d, c, c+d are rational squares}.

Every c in C_d actually satisfies c>d: equality would make both d and 2d nonzero rational squares, which is impossible. For rational C,e with e>0, the condition

C-e, C, C+e in C_d

says exactly that all nine values C+i d+j e, i,j in {-1,0,1}, are rational squares. Their arrangement

    C+d       C-d-e       C+e
    C-d+e     C           C+d-e
    C-e       C+d+e       C-d

has all eight line sums 3C. Conversely any 3-by-3 rational magic square has opposite entries summing to twice the central entry C. Calling its top-left offset d and top-right offset e, the remaining entries are forced to be the displayed form. Independent sign reversals of d,e and interchange do not change the nine-value set. For a distinct square, d and e are nonzero, so both may be chosen positive when discussing existence.

The nine values are distinct if and only if

    e not in {d,2d,d/2}.

Indeed a collision is k d+l e=0 with k,l in {-2,-1,0,1,2}, not both zero; for d,e>0 its possible ratios are exactly 1,2,1/2. Positivity is C>d+e. This is automatic from the three C_d memberships, but worth retaining as an explicit certificate guard.

This is the fixed-gap transpose of the fixed-center approach D_C={t: C-t,C+t are squares}: instead of seeking offsets d,e,d+e,d-e in one D_C, seek a three-term AP of centers in one C_d. It is a change of slicing of the same nine square equations, not a new relaxation. Clearing a common denominator L of the nine square roots multiplies every entry, C,d,e by L² and gives an integer magic square with the same distinctness and line sums. Conversely an integer solution is rational. Individual rescalings are not permitted.

## 2. Equal-area triangle equivalence

Let N=d/4. For a positive rational right triangle (a,b,h) of area N,

    a²+b²=h²,     ab/2=N,
    (a+b)²=h²+d,  (a-b)²=h²-d.

Thus c=h² belongs to C_d. Conversely let u,h,w be the positive rational square roots of c-d,c,c+d. Then

    a=(w+u)/2, b=(w-u)/2

are positive rational legs, since w>u, and

    a²+b²=(w²+u²)/2=c=h²,
    ab/2=(w²-u²)/8=d/4=N.

The unordered positive-leg triangle is recovered uniquely from c, after forgetting root signs and leg order; the canonical order above has a>b. Consequently a guarded three-center AP in C_d is exactly three positive rational triangles of one common area N whose hypotenuse squares are in that AP, with the full nine-square collision guard retained. Existence of one equal-area triangle, or even infinitely many, does not prove that such a triple exists.

## 3. Doubled elliptic x-coordinate, with a proved inverse

On E_N: Y²=X³-N²X, take a finite rational point P=(X,Y) with X>N and Y nonzero. Define positive lengths

    a=|(X²-N²)/Y|,
    b=|2NX/Y|,
    h=|(X²+N²)/Y|.

Substitution verifies the right-triangle identity and area N. Conversely any positive rational triangle (a,b,h) of area N gives a positive-branch point

    X=N(a+h)/b,   Y=2N²(a+h)/b².

Here X>N and Y>0. The identities b²=(h-a)(h+a), ab=2N show Y²=X(X²-N²), and the forward formulas recover the ordered legs exactly. No external surjectivity theorem or rank computation is needed.

The chord-tangent law gives

    x(2P)=(X²+N²)²/(4Y²)
          =(X²+N²)²/[4X(X²-N²)].

Therefore

    h²=4x(2P).

The point-to-center map is surjective, not one-to-one: signs and 2-torsion translates can yield the same doubled x-coordinate. Combining the proved inverse with Section 2 yields the exact image identity

    C_d = {4x(2P): P in E_{d/4}(Q), X>d/4, Y nonzero}.

Thus the existence problem is equivalent to three such eligible points P_-,P_0,P_+ for which the three doubled x-coordinates form a nonconstant AP. Write C=4x(2P_0) and e=4(x(2P_0)-x(2P_-)); sort the coordinates to make e>0. Require e not in {d,2d,d/2}. Membership provides C>d+e automatically.

For an independent direct-coordinate check, roots u,h,w of c-d,c,c+d give Q=(c/4,huw/8) on E_N. Its three quantities x(Q)-N, x(Q), x(Q)+N are respectively (u/2)²,(h/2)²,(w/2)². The explicit triangle inverse above supplies a rational half P of this Q (with the canonical a>b ordering), so this is not merely a point on E_N but a doubled point.

One may instead allow all rational points with Y nonzero and use absolute sides; the positive-branch formulation avoids unnecessary sign cases and is already surjective. Merely asking for three arbitrary rational x-coordinates in AP on E_N is not this condition: the image must be doubled points, and the gap/distinctness guards remain. The torsion triple (-N,0),(0,0),(N,0) is an immediate hostile for an unqualified arbitrary-x-AP claim. Likewise the group AP P,2P,3P does not imply an AP of the rational functions x(2P),x(4P),x(6P). The canonical-height parallelogram identity is not an identity for these x-coordinates.

Elliptic translation by a 2-torsion point T satisfies 2(P+T)=2P. It therefore changes the point presentation without producing a new center. Point counts must be deduplicated by h², not by point address or sign. The checker verifies these collisions directly on E_6.

## 4. Fixed integer gap is a finite divisor problem

For an integer gap d>0, integer centers in C_d correspond exactly to positive integer pairs r>s with

    d=2rs,   c=r²+s² a square.

Proof: since c and c±d are integers and rational squares, their roots u,h,w are integers. From w²+u²=2h², reduction modulo 4 forces w and u to have the same parity. Hence r=(w+u)/2 and s=(w-u)/2 are integers, r>s>0, with d=2rs and c=r²+s². Conversely such a pair gives c-d=(r-s)² and c+d=(r+s)². A pair r=s would require 2r² to be a square and cannot occur. Thus enumerate factor pairs of d/2 and test the hypotenuse square; no bound on c is needed for each fixed d. In particular C_24 intersect Z={25}, coming from (r,s,h)=(4,3,5).

This is not a finite enumeration of rational centers at fixed d. A rational triangle may have arbitrarily large denominators, and clearing those denominators replaces d by L²d. Confusing the finite integral slice with the rational curve would produce a false global exclusion.

Under common side scaling by L, every triangle area becomes L²N and every hypotenuse square becomes L²c. The curve changes to E_{L²N} under (X,Y) -> (L²X,L³Y). Scaling each triangle independently generally destroys both common area and the center AP. Taking individual primitive representatives is therefore not a valid normalization of a prospective three-triangle witness.

## 5. Reproducible finite certificates and hostiles

Run `python3 check.py`. Only Python standard-library exact integers and Fraction arithmetic are used. `certificate.json` contains every sampled point and triangle and the full integer-gap catalog.

- All nP for 1<=n<=24, with P=(12,36) on E_6, pass curve membership, area, square-center, inverse-map, and duplication identities. Their 24 centers have no distinct three-term AP.
- Two maximal-cardinality finite hostiles are C_3360 intersect Z={3364,5476,12769} and C_13440 intersect Z={13456,21904,51076}; neither triple is an AP.
- All integer gaps 1<=d<=20000 are exhaustively enumerated by the proved factor-pair formula. There are 105 integer centers across 96 nonempty gaps, with at most three centers per gap. No guarded center AP occurs in that universe.
- The integer catalog is independently cross-checked by direct square-root enumeration for hypotenuse h<=200 and d<=20000.
- The ordinary magic linear equations and exact collision ratios are checked for every integer pair 1<=d,e<=14.
- The first three elliptic multiples have nonzero center-AP defect, exactly -2257995926870999701/4268239692012450.
- The three nonzero rational 2-torsion translates of P all have the same double, demonstrating repeated-center contamination.
- One common side denominator for the first three sample triangles is 92393070; common scaling preserves equal area as L²N.

The finite negatives say nothing about unsampled rational points, larger integer gaps, other elliptic curves, or the global magic-square problem. No claim that P generates E_6(Q) is required or made.

## 6. Inheritance and primary-source priority boundary

Read in the existing repository:

- `05-knowledge/results/ap_square_pythagorean_20261019.md`: elementary three-squares/Pythagorean maps, denominator and positivity discipline, with no novelty claim.
- `05-knowledge/results/congruent_square_gap_classes_20261019.md`, Sections 2.3 and 3: the same congruent-curve triangle map and duplication formula; its target is AP square-sum triangles and gap squareclasses, not the present simultaneous nine-square condition.
- Session55 fixed-center integer search and diamond interface: the present fixed-gap construction transposes that interface rather than replacing its constraints.

Primary sources inspected on 2026-10-25:

1. William Stein, *The Congruent Number Problem*, Proposition 6.5 (three squares and triangle/curve correspondence): https://wstein.org/edu/2007/spring/ent/ent-html/node92.html . Its displayed signed-point convention differs from our chosen positive branch; the equivalence itself is classical.
2. Onno M. Cain, *Gaussian Integers, Rings, Finite Fields, and the Magic Square of Squares*, arXiv:1908.03236v2, Section 1: https://arxiv.org/html/1908.03236v2 . Cain explicitly identifies the earlier connections to congruent numbers, elliptic curves, and rational right triangles, citing John P. Robertson, *Magic Squares of Squares*, Mathematics Magazine 69(4), 1996, pp. 289–293, https://www.jstor.org/stable/2690537 . Robertson is identified through Cain's bibliography; its full text was not inspected here.

Accordingly this package claims a verified, guard-preserving synthesis and concrete bounded certificates, not discovery of the magic-square/congruent-number connection. The useful frontier is the simultaneous doubled-x AP constraint on one fixed curve; the triangle or elliptic reformulation alone does not solve it.
