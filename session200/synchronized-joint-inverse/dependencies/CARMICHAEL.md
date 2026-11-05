# Source-fiber twins, dyadic splitting, and the rigid nonbalanced chamber

Session 200. **PROVED:** every nontrivial balanced target has two distinct
physical source profiles with the same C1 readout. Their exact continuations
split into two explicit dyadic families. Any represented target outside the
balanced image has a rigid source-word boundary, and its two-run subfamily
is exactly a shifted product of Mersenne numbers.

**FINITE-EXACT:** through 50,001 the only singleton fibers of the full C1
readout on physical numerical profiles are 3 and 11. **OPEN:** the analogous
all-order statement, and the {7,21} coverage problem. None of these claims
counts graph isomorphism classes by counting profiles.

## 1. Inherited physical domain and the fiber question

Write a nontrivial rooted start-one source as

    s=(D,E,U,d)=(H(V),e_c(V),R_c(V),H(V-c)).

The complete physical grammar starts at t=(1,0,2,1) and uses chronological
letters

    O(s)=(D,0,2U,d),
    Z(s)=(D-E+U,d,2U-d,d),
    J(s)=(D-E+U,D-E,D+U,D).

All these moves are literal vertex births; only J changes the root. The
inherited physical inequalities include

    D-E>=d>=1,  D<=U,  U>=2d.

They follow from the positive coordinates in the frozen adjoint dependency.
In particular U>d. Put

    F(s)=D-E+U,   G(s)=D+2U.

The full strong C1 spectrum is F(S). Closing a source by J produces its
C1 graph. The balanced spectrum is {3} union G(S), where the separate 3
comes from the singleton precursor; closing by OJ realizes G.

The multiplicity considered here is

    m(n)=#{s in S:F(s)=n},

where identical numerical quadruples are counted once. This is finite:
F(s)>=U+1, and the complete U-bounded source inventory is finite. It is
neither the number of words nor the number of unrooted graph classes.

## 2. The elementary Carmichael import, with its exact scope

The totient-multiplicity conjecture asks whether every nonempty fiber of
Euler's phi has at least two positive integers. Its elementary first move
pairs odd m with 2m: phi(m)=phi(2m). Thus a singleton preimage, if any,
must be divisible by 4. This uses phi(2)=1 and coprimality, not an unguarded
formal multiplication rule.

The import is to seek a physical partner operation and to isolate its
uncovered chamber. It is not a map taking totients to path counts. Earlier
repository Carmichael-number work concerns Korselt conditions, a different
question. Its retained-guard and fiber-occupancy methods are useful, but its
number-theoretic conclusions are not premises here.

## 3. A guaranteed physical twin and an exact splitting law

**Theorem 1.** Every s in S has two distinct profiles O(s),J(s) with

    F(O(s))=F(J(s))=G(s).                         (1)

Their cover/deletion pairs are

    (U_O,d_O)=(2U,d),  (U_J,d_J)=(D+U,D).        (2)

In particular the J choice weakly lowers U while weakly raising d, at the
cost of one more root change in the supplied recipe.

**Proof.** Substitute into F. The O output has endpoint coordinate zero;
the J output has endpoint D-E>=d>=1, so the profiles are distinct.
Equation (2) and d<=D<=U prove the remaining assertions. QED.

The swap is an involution on *tagged presentations* (s,O) and (s,J). It is
not asserted to be a globally defined involution on profiles without the
parent recipe: O forgets the parent's E, and J forgets its d.

This guarantees profile multiplicity but does not imply any constant
graph-isomorphism multiplicity.
For example tO and tJ are distinct profiles with F=5; closing them by J
gives isomorphic copies of the unique strong order-four tournament.

**Theorem 2.** For k>=0,

    F(Z^k O(s))=G(s)+(2^(k+1)-2)(2U-d),
    F(Z^k J(s))=G(s)+(2^(k+1)-2)U
              =D+2^(k+1)U.                     (3)

Thus for k>=1 the two continuations differ by exactly

    (2^(k+1)-2)(U-d)>0.                         (4)

Every value in (3) is physically C1-realized. The second branch is also
balanced: it equals G(O^k(s)).

**Proof.** Direct induction in Z gives, on any source q,

    U(Z^k q)=2^k(U(q)-d(q))+d(q),
    F(Z^k q)=F(q)+(2^(k+1)-2)(U(q)-d(q)).        (5)

The F increment in one Z step is 2(U-d), and the new U-d is twice the old
one. Apply (5) to the twin profiles in (2). QED.

This is a constructive output: a given physical source compiles two exact
infinite target families. It also exhibits the precise missing continuation
coordinate. Equality of scalar F does not license future source substitution.
A single Z distinguishes every pair in Theorem 1. The needed sidecar is at
least U-d for this continuation, and the full four-coordinate response for
arbitrary O/Z/J continuations.

At t the twin value is 5. The first split gives 11 versus 9, then 23 versus
17, then 47 versus 33. The 11 source is not a mysterious numerical exception:
it is the first nonbalanced branch of this smallest physical split.

## 4. The chamber left after removing automatic twins

The inherited physical adjoint sends a word w to w-dagger: reverse it and
exchange O with Z while fixing J. It proves

    F(t.w)=F(t.w-dagger).                       (6)

This preserves actual source recipes, not merely a relaxed cone.

**Theorem 3.** Suppose n lies in F(S) but not in {3} union G(S). Every word w
from t representing n has the following properties:

1. w is nonempty, begins with O and ends with Z;
2. its first letter not equal to O is Z;
3. its last letter not equal to Z is O.

Consequently such a word either has exactly two runs O^a Z^b, a,b>=1, or
has the form O^a Z [arbitrary middle] O Z^b. Any J is inside these two gates.

**Proof.** Ending in O or J makes the target a G value by (1). Beginning
in Z or J has the same consequence after (6). The empty word gives 3.
Moreover a suffix J Z^k gives a G value by (3), and a prefix O^k J gives
such a suffix after dagger. These exclude the stated first/last alternatives.
Removing the initial O run and final Z run gives the final description.
QED.

Thus, apart from 3, every hypothetical singleton m(n)=1 is a nonbalanced
target in this chamber. All its representing words satisfy these guards.
This is a necessary restriction, not a descent proof eliminating the chamber.
The dagger can supply a second profile, but may fix a word or yield the same
profile. Counting two words without checking their profiles would be invalid.

## 5. A complete arithmetic compiler for the two-run chamber

**Theorem 4.** For a,b>=1 the physical source word O^a Z^b gives

    D=2+(2^b-1)(2^(a+1)-1), E=d=1,
    U=1+2^b(2^(a+1)-1),
    F=2+(2^(a+1)-1)(2^(b+1)-1).                (7)

Conversely every shifted product in (7), with a,b>=1, has that physical
recipe. Exchanging a and b preserves F and changes U by 2^a-2^b. Therefore
all the off-diagonal pairs a!=b have m(F)>=2, using deletion-one sources
alone.

**Proof.** O^a t has D=d=1,E=0,U=2^(a+1). Iterating Z gives the displayed
coordinates, or use the exact binary-block compiler. The factorization and
the cover difference follow by expansion. The dagger exchanges a and b.
QED.

The diagonal a=b is not disposed of by swapping. It gives

    n=(2^(a+1)-1)^2+2,

beginning 11,51,227,963,3971,... . Only 11 is a singleton in the exact
finite census below. No all-order statement that these later diagonal values
are balanced is proved here. This diagonal is a sharply specified hostile
for an attempted universal partner proof.

## 6. Finite-exact census and verification

The author census exhausts the complete physical inventory through U<=50,000,
not words of a guessed maximum length. It visits 3,816,297 profiles. By
F>=U+1 this includes every source with F<=50,001. Among the represented
values, precisely two fibers are singletons:

    n=3:  (1,0,2,1), word empty;
    n=11: (5,1,7,1), word OZ.

The census uses the raw source grammar with predecessor words. The verifier
included here replays this census when run with --cap 50001; its default cap
10,001 is smaller and fast. It separately validates the formulas on a full
formal coordinate basis, checks all source words through length 6, and counts
literal graph examples. Normal and python -O runs must agree; all assertions
use explicit exceptions.

An independent auditor used a separately written source enumeration through
U<=5,000: 126,910 profiles, 2,471 represented values F<=5,001, again only
singleton fibers 3 and 11. The auditor also confirmed the order-four graph
collision, preventing the wrong interpretation of m(n).

**Open source-fiber conjecture.** For every n in F(S) other than 3 and 11,
m(n)>=2.

This is deliberately conditional on membership in F(S). Even a proof would
not show that an absent target has a source. The {7,21} coverage target asks
for nonempty fibers in the relevant C1/C2/C3 union; a singleton prohibition
and a positive multiplicity average do not answer that question. The useful
transfer is the explicit partner, its exact continuation defect, and the
smaller guarded chamber that an eventual descent must handle.

## 7. Canonical physical descent to a balanced parent below half-height

This is the target-existence continuation of the fiber argument. Retain the
joint physical image

    P={(m,h)=(D+2U,2U-d):s=(D,E,U,d) in S}.           (8)

It is essential that one source realizes both coordinates. All D,d are odd,
and E and U have the same parity, by induction in the source grammar.
Consequently m,h are odd. The physical inequalities prove

    3<=h<m<2h.                                      (9)

Indeed h>=3d>=3, m-h=D+d>0, and
2h-m=2U-D-2d>=U-2d>=0. The last difference is odd and hence positive.

**Theorem 5 (exact existence compiler).** With B={3} union G(S),

    F(S)=B union {m+2(2^k-1)h:(m,h) in P, k>=1}.     (10)

Every value on the right has a literal source recipe and a C1 closing.

**Proof.** The second family is exactly s.O.Z^k followed by the closing J,
by Theorem 2. Conversely take any source word w. The empty word gives 3.
If w consists solely of r>=1 Z letters, its value is 2^(r+1)+1, also
G(O^(r-1)t), and is balanced. Otherwise remove the terminal Z run, which
may have length zero. The last remaining letter is O or J. At length zero
both give G. At positive length a last J gives a balanced value by (3),
while a last O gives the second family in (10). QED.

**Theorem 6 (canonical descent on each nonbalanced profile).** Let q be an
actual physical profile with n=F(q) outside B. There is a uniquely determined
pair of integers

    k=v_2(U_q-d_q)>=1,
    h=(U_q-d_q)/2^k,                             (11)

and a uniquely determined scalar parent

    m=n-2(2^k-1)h,                               (12)

such that (m,h) lies in P. The parent is balanced and satisfies m<n/2.
Every source-word presentation of q ends in exactly k terminal Z letters,
preceded by O. Removing that O and these Z letters gives an actual source
s with G(s)=m and 2U_s-d_s=h.

The *parent profile* is not asserted unique: the O operation forgets E.
Different profiles q in the same scalar fiber may also have different k,h,m.

**Proof.** Theorem 5's cases show that any presentation of a nonbalanced
q must end O Z^k with k>=1. If s precedes this suffix, then

    d_q=d_s, U_q-d_q=2^k(2U_s-d_s).

The parent d_s is odd, so the factor 2U_s-d_s is odd. Thus its exact 2-adic
valuation is k, proving (11) for every presentation and not just a chosen
word. Equation (3) gives (12), and the original predecessor is the physical
witness for (m,h) in P. By (9), h>m/2, so
n=m+2(2^k-1)h>2m for k>=1. QED.

This descends the numerical target, not merely graph order. It does not
assert that every odd target has a descending parent. The remaining joint
condition is explicit: for a prescribed n and k, choose an odd h satisfying

    n/2^(k+1)<h<n/(2^(k+1)-1),
    m=n-2(2^k-1)h,

and require (m,h) in P. One may enumerate the actual sources for this test;
separate membership of m in B is insufficient. The bound h>=3 gives
2^(k+1)<=(n+1)/3, so there are only finitely many k.

For a fixed candidate (m,h), the source has

    D=m-2U, d=2U-h,

and necessary physical bounds

    max((h+1)/2,ceil(m/3)) <= U
       <= min(floor((m+h)/4),floor(2h/3)).        (13)

These are pruning bounds, not a physical membership theorem. The predecessor
recipe remains the constructive certificate.

### A decisive hostile to the one-step descent shortcut

The target 443 belongs to F(S) but not to B. Its entire numerical source
fiber, found by a complete U<=442 enumeration, is

    (219,1,225,1), word OOZZZZZ,
    (191,1,253,1), word OOOOOZZ.

Their canonical data (k,h,m) are respectively

    (5,7,9), (2,63,65).

There is no k=1 parent, even though many candidate scalar parents are
balanced. This is not an incomplete recipe search: any k=1 presentation
would produce a third profile or force v_2(U-d)=1 in one of these two,
contradicting the complete fiber and Theorem 6. The independent verifier
also directly checks every possible k=1 physical parent using n>=5U+1.
The first target missed by restricting (10) to k=1, while keeping the whole
balanced image, is 443. At cap50,001 the only such misses are 443 and617.
Those last two statements are FINITE-EXACT, not a universal two-step theorem.

The physical adjoint can rebalance an existing scalar witness to two
square-root-height sources, but it cannot create a k=1 presentation of443.
Its scalar height bound does not preserve h by itself. To transport the
joint constraint exactly, h(s)=K(z,s) uses the auxiliary vector
z=(1,0,2,0), while G(s)=K(b,s) uses b=(3,1,3,1). The adjoint can transport
both boundaries along the same word, but z has deletion coordinate zero
and is not an actual source. Retain it as an algebraic sidecar, never as an
extra physical member of S. This states precisely what a midpoint-based
existence proof would still have to synchronize.
