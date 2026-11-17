# Root migration makes a neutral profile usable in the reflected bank

This combines the [Farey reflected flank](../../chanpak_free_seam_second_move_20261117/extra/PROOF.md)
with a new root-moving operation. Five vertices already give two candidate
flanks indistinguishable to every unmarked slot polynomial, but not to a
common rooted edit. Both choices and their +2 reflection partners are
actual tournaments.

## The three-coordinate root transfer

Let A have a distinguished vertex r, at which exactly one full Hamilton
path ends. Record h=H(A), R=count of complementary path pairs with one
nonempty path ending r (weight each such distinguished path by H of its
complement), and S=count of full paths starting r. Adjoin a transitive
block B on k vertices with its increasing order, and a new root s
with a k-bit contact word u (bit 1 means s beats the corresponding
vertex). Orient all of B and s toward A, except for r→s. Regard s as
the new distinguished vertex. For k>=1 the transfer is

    h' = F(u)h+(u+1)R-[u even]S,
    R' = 2^k R + (2^k-rev_k(u))h,
    S' = [first bit is 1] h.

For k=0 the transfer is (h+R−S,h+R,h−S). The new root has exactly
one full path end at every step: after deleting s, B→A, so a full
path ending s must traverse B and then the unique full path in A
ending r. For k=0 the same holds with B absent.

For the first read, delete s and use the adjacent marked-interval
law on B→A. Its marked start count is u, signed start −[u odd],
and the unmarked end on A has weight R, signed end −S. The contribution
within A+s has count h+R−S by the same signed path-end formula, so
h'=F(u)h+(u+1)R−[u even]S. Equivalently split the possible directed
odd cycle using r→s: it has a path starting s in B+s and a path
ending r in A. Their length parities must be opposite. The signed
weighted end in A equals minus the full start count S, as proved
by deletion of the distinguished vertex in the Farey note.

For R', paths ending s either use r→s, with 2^k arbitrary choices
of earlier B vertices, or stay within B+s. In the latter case s
alone contributes h, and a last B-vertex at position j can be
preceded by any subset of its j predecessors, provided its contact
bit is zero. This gives (1+RC_k(u))h=(2^k−rev_k(u))h. For S' when
k>=1, after s the old deletion tournament B→A must start at the
first vertex of B; when k=0 the first old vertex can be anything
but r. These formulas work even if the new A is not strong. In the
reflected construction, positive outer prefix on the left and the
ordinary flank on the right give the assembled strong tournament.

## A small neutral pair, including every slot count

Start from a directed triangle with any root (h,R,S)=(3,3,1).
There are two 5-vertex descendants:

| migration words (k,u) | outgoing masks | root | (h,R,S) | P2 |
|---|---|---:|---|---:|
| (1,1) | (18,4,1,7,14) | 4 | (9,9,3) | 84 |
| (0,0), (0,0) | (10,4,1,22,7) | 4 | (9,11,3) | 84 |

Call them X0,X1. Each has one full path end at its root. In fact
they have identical *entire* ordered-slot polynomial, not merely
the displayed two values:

    P_q(X0)=P_q(X1)=q^5+6q^3+2q.

Indeed a 5-vertex odd-cycle packing has at most one nontrivial
cycle. The reads h=9,P2=84 determine the two coefficients
uniquely. Thus equality survives every composition seeing only
the unmarked P_q of this slot, with no loss of a tagged degree
of freedom. A common singleton migration has counts 15 and 17
and still preserves the unit end. This neutral pair is a literal
change of rooted topology rather than a high-address collision.

One can retain that difference as a whole polynomial. Prepend a
transitive dominating block of t vertices to each X before the
same final singleton migration. For the original rooted X define
E_q,O_q by selecting a path ending r of positive even or odd
length, weighting its complement by P_q. With differences X1−X0,
we have Delta E_q=q and Delta O_q=q². For E_q the possible degrees
are 3 and 1; the first cancels because both roots have equal
indegree (check their masks); the value at q=1 gives the second.
For O_q the degree-4 singleton contributions cancel, and the
constant full-path terms cancel by the common unit end, leaving
degree 2 determined at q=1. Prepending the transitive block
convolves the parity split with (q+1)^t and (q−1)^t: each added
vertex goes in the distinguished path or in the residual slots.
The new singleton can lie in an odd mixed cycle only with an
even selected path ending r, so its mixed-cycle contribution is
2q E_q. The full difference after migration is therefore

    Delta P_q = q²((q+1)^(t+1) - (q-1)^(t+1)).

In particular the count shift is 2^(t+1), and the two-cover shift
is 4(3^(t+1)−1). The reader can separate the 2- and 3-powered
responses physically, after a starting interface where *all*
unmarked slot counts were equal. This is an exact profile
fusion rule, not a claim that arbitrary independent powers are
available in one physical source.

## Fuse the neutral pair into the actual +2 bank

Use Xi as the nontransitive flank, marked at its distinguished
root, after an arbitrary transitive contact prefix u; take the
other flank to be an even positive ordinary address w. The
Farey determinant is still one for both i=0,1, so each physical
lower count has its exact +2 reflection mate. Their marked
profiles (h,P,L,R,Lminus,Rminus,K) are

    (9,84,21,66,-1,-6,29) and (9,84,23,64,-1,-6,31).

The weighted starts differ by 2 and the weighted unmarked ends
by −2, while all signed coordinates agree. Therefore in the
lower Farey read the differences in the left marked start
and internal apex correction are 2 and 2−2u, respectively.
The two reflected pairs differ by

    5(2−2u)+142·2+256·2w = 294−10u+512w.

For any prescribed even displacement D, choose positive even
w in the class 512w+294−D=0 mod10 and put
u=(512w+294−D)/10. Infinitely many positive choices exist.
For D=4 write w=10j,u=512j+29. At j=10 the four realized
counts are 11136264313, 11136264315, 11136264317, 11136264319;
the last two are primes (elementary trial division in replay).
No infinite-prime assertion is made. This branch still chooses
a skew line and cannot select an arbitrary absolute starting
count. Its usefulness is the preservation of the unmarked
slot interface and a second exact physical control, not
horizontal density by itself.

## A reset boundary in the padded migration language

With a leading zero in each nonempty word, the new S is zero,
and a preceding S=0 state uses the 2 by 2 matrix

    [ F(u)      u+1 ]
    [ 2^k-rev   2^k ].

This matrix has rank one exactly when u=2^j−1 (including j=0),
and strictly positive determinant otherwise. Let rho(u) denote
the bit-reversed dyadic fraction rev_k(u)/2^k, independent of
leading zeros, and f(u)=F(u)−(u+1)(1−rho(u)). Directly,
f(2v+1)=f(v), and f(2v)=f(v)+[v odd]−rho(v)/2.
Across a run of j>=1 zeros appended to a positive odd v the
gain is 1−rho(v)(1−2^(-j))>0. Initial leading zeros contribute
zero. Thus the determinant, 2^k f(u), vanishes exactly when
the effective word has no zero run. In a rank-one migration
from S=0, R'/h'=2^(k−j) is integral and S'=0. Once this reset
has occurred at h'>1, all three active coordinates are multiples
of h', so every subsequent count under any of the integer
root-transfer rules is divisible by h'. Further strict growth
can therefore never create a prime. Prime-target controllers
must retain the determinant rather than reset the hidden read.
