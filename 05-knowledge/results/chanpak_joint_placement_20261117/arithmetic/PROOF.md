# An infinite exact collision fiber with distinguishable continuations

**Status: PROVED**, conditional only on the already established physical free-seam reader. This is an all-order family, not an extrapolated finite collision census. It does not prove prescribed-count coverage.

The inherited reader, for positive odd u and positive even w, is

    H(u,w)=81+142(u+w)+256uw+5(F(u)+F(w)),
    F(0)=1, F(2a+1)=F(a), F(2a)=F(a)+a+[a odd].

Each address pair represents a strong tournament, with a reflected strong tournament of the same order and count H+2. The three-contact template of `free_second_move/three_contacts/PROOF.md` preserves this reflected pair.

## 1. Exact fibers, with an unbounded free parameter

Fix G=2^k, k>=2. Set r=3 and s=G-4, so that r+s=G-1. For k-bit reverse complements C_r,C_s put

    A=16384G^2-5,
    B=5120G^2-1024G(s-r)+10-20(C_r-C_s).

For any fixed even delta, choose the unique z0 in {1,...,15} with

    delta+A z0+B=0 mod15,
    y0=(delta+A z0+B)/15.

This is possible since A is a unit modulo15. For all sufficiently large integers t>=0 set

    z=z0+15t, y=y0+A t, x=y+z,
    u=G(16x+5)+3, w=G(16y+4)+G-4,
    v=u+4G, q=w-4G.

Then x is even: delta+B is even and A is odd, hence y and z have the same parity. The addresses are positive for sufficiently large t, have unchanged leading bits under the transformation, and the three-contact formula gives

    H(v,q)-H(u,w)=15y-Az-B=delta.

For delta=0 this yields infinitely many **exact count collisions**. The two presentations also have exactly the same sum u+w=v+q, the same complementary low-word sum, and the same order. Thus retaining both absolute count and address sum does not remove this ambiguity.

## 2. A literal one-arc continuation separates the collisions

Both w and q are divisible by4. The operation w->w+2, and respectively q->q+2, toggles bit1, a single actual apex contact. It preserves leading bits and parity, and is compatible with the reflected +2 partner.

For any positive odd u and w=4a,

    H(u,w+2)-H(u,w)=512u+294-5(a+[a odd]).

Under (u,w)->(v,q) we have u increased by4G and a decreased byG. Since G is even, the parity term is identical. Consequently

    (H(v,q+2)-H(v,q))-(H(u,w+2)-H(u,w))
        =2048G+5G=2053G.

Combining with the preceding exact displacement,

    H(v,q+2)-H(u,w+2)=delta+2053G.

In particular, **identical starting counts and identical address sums have different counts after the same named low-bit continuation**, with exact separation2053G. This supplies a certified nontrivial continuation coordinate on a scalar count fiber. Reflection duplicates the entire diagram with every count increased by2.

For G=4,t=0,delta=0 the two address pairs are

    (3757847,3757648), (3757863,3757632).

Their common lower count is3614891650465537. The continued lower counts are3614893569786435 and3614893569794647, differing by8212. Numerical controls use the recurrence independently of the affine displacement formula.

## 3. Interpretation and limits

The diagram is a physical mixed-difference law: an exact zero-count reconfiguration changes the response to a subsequent contact. It explains why a quotient by count, even enriched by the total address sum and order, does not carry a well-defined deterministic low-bit dynamics. An exact-count search can profit by retaining several representatives of one count rather than collapsing them prematurely.

This is not an independently repeatable2053G adjustment at an arbitrary count. Its starting values are constrained by the collision line. For fixed G and delta, as t increases,

    H(u,w)=65536 G^2 A(A+15)t^2+O_G,delta(t).

The bound follows from F(n)<=n+1. Therefore the family is sparse, with order sqrt(X) parameter values below X; its usefulness is a reusable structural interface distinction, not numerical density.

The next constructive obligation is to place such a collision-bearing interface at a prescribed target, or to compose collision fibers while controlling their baseline and response. Neither the exact zero displacement nor the nonzero commutator solves that obligation alone.

## Reproduction

Run `python joint_placement/arithmetic/collisions.py`. It performs5408 exact controls across13 dyadic scales,104 parameter choices (including200-bit parameters), and four fixed shifts. It checks the two exact laws and unchanged binary lengths. `RECEIPT.json` records the output. These checks supplement, rather than replace, the algebraic proof.

## 4. A closed recurrence subfamily, not just a digit recurrence

The collision fiber contains an exact quadratic dyadic subsequence. Fix G and delta, and choose the positive initial solution above (or translate t by a fixed integer first). Write

    u(t)=U1 t+U0, w(t)=W1 t+W0.

Here U1=16G(A+15), W1=16GA are even. Choose b with U0,W0<2^b and put D=2^b. On t=D*2^j, the suffix identity gives, for c=U0 and p=U1 (and likewise for w),

    F(pD*2^j+c)
      =p(D-rev_b(c))*2^j +F(p)-p+F(c)-1.

Indeed reversal in b+j bits is2^j rev_b(c); the possible leading-zero parity term vanishes because p is even. Hence the physical count is **exactly**

    H_j=alpha*4^j+beta*2^j+gamma,
    alpha=256 U1 W1 D^2,
    beta=256D(U1 W0+W1 U0)+142D(U1+W1)
         +5[U1(D-rev_b(U0))+W1(D-rev_b(W0))],
    gamma=256U0W0+142(U0+W0)+81
          +5[F(U1)-U1+F(U0)-1+F(W1)-W1+F(W0)-1].

Thus H_(j+3)=7H_(j+2)-14H_(j+1)+8H_j at every j>=0. Every term carries the exact collision and the distinct continuation responses. This separates three levels cleanly: the whole collision line is governed by a binary digit reader; a dyadic subsequence has a finite constant-coefficient recurrence; neither statement implies all odd targets occur.
