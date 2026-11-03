# A complete projective decoder and an all length return obstruction

**PROVED for the eight specified source macros and the supplied seed ray.**
This does not classify all source words, all tournaments, or arbitrary native
complement and re-marking programs. Proof review is by the same agent;
the finite controls use independent algorithms.

## P1 The actual source family

Use the inherited coordinates (v=(h,x,y)), correction

    tau=(2,-1,0),  beta(v)=(Y,X)=(y,h+2x-y),  g=h-y.

Here h is H, y-x counts paths ending at the supplied port, and x+y counts
covers with a distinguished path starting at that port, including the
empty second path. These are invariants of the rooted tournament. The
integer coordinates (Y,X,g) satisfy X=g mod 2.

The three literal seed states are

    V0=(31,14,25), V1=(95,41,75), V2=(97,40,75),
    Vj=3 V0+j tau (j=1,2),
    beta(V0)=(25,34), beta(V1)=beta(V2)=3(25,34).

Their source orders are 7,8,9. The seeds and source matrices are inherited
from [the exact inverse packet](../golden_inverse_defects_20261103/TRANSPORT.md).
Words below are in **processing order**: append a letter to the processing
word by right multiplying its matrix. Its actual source component is
prepended, so the component order reverses this chronological word.

There is an exact explanation of the pair 96-1,96+1 in these coordinates.
V1 and V2 have the same bulk (75,102), with charges g=20,22. Holding the
bulk fixed, integer states require g=X mod 2 and differ by integer multiples
of tau. They are consecutive points on that integer fibre. Their midpoint
is (96,81/2,75), which is not an integer state. This is a literal lattice
statement; it does not turn the numerical midpoint into a tournament.

Use B=b, Q=q2, S=s78 and the five inherited error loops:

| Macro | Processing word | Added vertices | tau multiplier |
|---|---|---:|---:|
|B|b|3|-1|
|Q|q2|4|1|
|S|s78|6|1|
|L2|r24 a a a|8|-2|
|L6|r24 q4|9|6|
|L45|r19 r32 r48 q4|19|45|
|L13|r19 r19 a a r18|17|13|
|Lp2|r24 a a r21|12|2|

The first four loop words and gains are [inherited](../return_channels_20261103/THEORY.md),
and S is [inherited](../ternary_repair_clock_20261103/CARRIES.md).
Lp2 arrived in [Session155's mixed-return certificate](../session155_spectrum_20261103/returns/MIXED.md)
during this session and is explicitly included in the final alphabet.
Their literal 3 by 3 matrices are recomputed from the source atlas in the
receipt. Define (V_j(w)=V_jP_w). Each is realized by the ordered source
join followed by its actual spine complement. The retained spine and the
last-to-first edge give a directed Hamilton cycle, so all these targets
are strong. No scalar reset is an operation in this definition.

## P2 Eight disjoint ratio images

Because tau is an eigenrow of every macro, each has an autonomous bulk
matrix Q satisfying beta(vP)=beta(v)Q. They are

    B   [[5,8],[2,3]],                 Q   [[11,18],[4,7]],
    S   [[77,152],[28,55]],             L2  [[240,339],[90,127]],
    L6  [[720,1101],[270,413]],
    L45 [[1020432,1560769],[366879,561148]],
    L13 [[211965,378780],[76220,136205]],
    Lp2 [[4773,8628],[1789,3234]].

For a row (Y,X) with Y>0 put r=X/Y. A matrix [[a,b],[c,d]] acts by

    f(r)=(b+dr)/(a+cr).

Every matrix is positive and nonsingular. Its exact image of J=[1,2] is:

| Macro | Lower endpoint | Upper endpoint |
|---|---:|---:|
|L2|593/420|233/165|
|L6|757/495|1927/1260|
|L45|2121917/1387311|178871/116946|
|B|14/9|11/7|
|Q|5/3|32/19|
|L13|102997/57637|130238/72881|
|Lp2|5931/3281|15096/8351|
|S|262/133|69/35|

These closed intervals are pairwise disjoint and lie strictly inside J.
The seed ratio 34/25 lies in J below all eight. These are rational strict
inequalities, checked by cross multiplication, not numerical separation.

**Projective code theorem.** If two words w,z satisfy

    (25,34) Q_w = c (25,34) Q_z  for any c>0,

then w=z and c=1.

Proof: a nonempty word's ratio belongs to exactly the interval of its
last letter. Equal ratios force equal last letters. Cancel their invertible
fractional linear maps and repeat. If only one word becomes empty, its
ratio is 34/25, outside every last-letter image, a contradiction. Thus the
whole words coincide. Their nonzero vectors then force c=1.

This proves freeness on the seed ray, not merely absence of short relations.
It includes arbitrary mixed uses of the nonunit loops. In particular
the incoming D=L2 B L2 loop and every repetition of S are covered.

An additional metric statement is immediate: on J,

    |f'(r)|=|det Q|/(a+cr)^2 <= 1/45

for every macro, with maximum attained at the Q boundary r=1. The nested
images give an injective coding of infinite outer-to-inner words into a
compact subset of J. This assertion uses nested cylinder order, not an
unqualified convergence claim for arbitrary chronological forward orbits.
The n-level cylinders have total length at most (8/45)^n, so the limit
set has Lebesgue measure zero. This geometry is separate from density of
integer Hamilton counts.

## P3 Exact source membership and counting

Given an integer candidate v with positive Y and r in J, identify its last
macro from the table, invert its 2 by 2 bulk matrix, and reject unless the
predecessor is a positive integer bulk with ratio in J. Repeat until the
ratio equals 34/25. Each forward matrix multiplies Y by at least 7 on J,
so every accepted backward step strictly decreases positive integer Y.
This is a terminating algorithm; no arbitrary word-length cap is needed.
Finally test the resulting word against all three actual seed states using
the full 3 by 3 matrices. That last test retains the scale and charge and
rejects a rational projective address with the wrong integer state.

The projective theorem proves uniqueness of the decoded word. Scaling
distinguishes V0 from V1,V2; the remaining charge distinguishes V1 from V2.
Every accepted state has an explicit literal source word, not just a
positive solution of matrix equations. Reconstructing that word reconstructs
its rooted tournament within this family.

All states V_j(w) are distinct. Therefore the generating function for
distinct states, equivalently for the rooted isomorphism classes in this
particular constructed family, is exactly

    (z^7+z^8+z^9)/(1-z^3-z^4-z^6-z^8-z^9-z^12-z^17-z^19).       (P1)

Indeed different v values distinguish rooted classes, and each word has
one supplied rooted graph. This is stronger than counting presentations
with possible state collisions. It is not a series for distinct scalar H
values, unrooted classes, or A000568.

## P4 Complete classification of ternary affine edges inside the family

Suppose two family states satisfy

    V_k(z)=3 V_j(w)+d tau,  d in {0,1,2}.              (P2)

Their bulk directions are equal. P2's projective code theorem forces z=w;
the seed bulk scales then require k in {1,2} and j=0. For that common word,

    V_k(w)-3V0(w)=k lambda(w) tau,
    lambda(w)=(-1)^(#B+#L2) 2^(#L2+#Lp2) 6^#L6 45^#L45 13^#L13.

Consequently:

- Both positive offsets occur exactly when lambda(w)=1: no L2,L6,L45,L13,Lp2
  occurs, and the number of B letters is even; Q and S are arbitrary.
- The +2 offset uses V1(w). The +4 offset uses V2(w) when lambda=1,
  or V1(w) when lambda=2. The latter means exactly one L2 and odd B parity, or exactly one Lp2 and even B parity,
  with no other nonneutral macro.
- No 3H edge occurs inside this three-seed family.
- No V1(w) or V2(w) has either next affine offset inside this family,
  at any size or macro depth.

Thus the first 31-to-95/97 construction cannot iterate inside this entire
infinite macro bank. This is independent of a cutoff on source orders and
does not assume that return determinants are units at three.

More generally let a finite seed bank have bulks d_i(25,34), with positive
integers d_i, and use these same eight macros. Every affine ternary edge
must retain its word and satisfy d_target=3 d_source. The seed graph is
acyclic; a chain has length at most max v3(d_i)-min v3(d_i). This extends
the *method* of the incoming unit-return valuation obstruction to these
nonunit loops by a different hypothesis: disjoint projective images.

This does not exclude source families with other seed directions, other
macros, non-affine carry changes, or actual complement/re-marking resets.

## P5 Three Hamilton probes reconstruct this family

The inherited Q observation matrix has determinant 96:

    [h,H(QT),H(Q^2 T)] = v [[1,5,73],[0,8,144],[0,6,120]].

Writing the observations as (h,h1,h2), its explicit inverse is

    y=(h2-18h1+17h)/12,
    x=(20h1-27h-h2)/16.                              (P3)

Integrality requires the displayed numerators to vanish modulo 12 and 16,
respectively. The observation Smith factors are (1,2,48), giving image
index 96. Exact observations recover v, after which the decoder in P3
decides family membership and reconstructs the word, seed and vertex clock.
Three probes are sufficient; no minimality of their number on this sparse
family is asserted.

The seed Q observations are (31,417,7279). Changing only 7279 to 7281
retains positive odd numbers but fails the integer-image guard. The
independent fact that the two seed outputs 95 and 97 surround 96 comes
from their midpoint 3*31+3, not from the determinant computation. Equal
numerals here are not an identification of these two operations.

All five Q4 phases remain available in the ambient source language. Their
correction rows are respectively

    q0: (-6,-1,-2), q1: (-6,0,-3), q2: (2,-1,0),
    q4: (0,4,3),    q5: (0,3,2).

Only q2 preserves this correction line. The others are explicit ways to
leave the theorem's hypothesis; they cannot be merged into the eight-letter
projective code without retaining the changed correction channel.
