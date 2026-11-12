# Compensated carries make a four-point physical interval

**Status: PROVED using the inherited physical sandwich compiler.**
The constructions below give infinitely many four-point consecutive-odd
intervals. Their starting points are not arbitrary, their widths do not grow,
and Chan–Pak remains OPEN.

Use the primitive C9 central bank15381/15383 from the digits lane. Put

    g(0)=0,
    g(2a+1)=g(a),
    g(2a)=g(a)+a+(a mod2),

and write its lower output as

    F(A,B)=15381+3357(g(A)+g(B))+50004AB
                      +25915A+20732B+1865(A mod2)+1492(B mod2).  (1)

Every A,B≥0 is physically constructible, as proved in `../digits/PROOF.md`.
The central alternative supplies F(A,B)+2 with the same rooted interface.
No integer subtraction below is interpreted as a graph operation.

## 1. Primitive content does not make neighboring digit banks join

For all A,B≥0, the two-point banks at (A,B) and (A+1,B) never overlap or
abut as odd integer intervals. The same holds between (A,B) and (A,B+1).
In fact their lower endpoints differ in absolute value by at least6.

Proof. Let A=2a. Then

    F(A+1,B)−F(A,B)=50004B+25915+1865−3357(a+(a mod2)).

Modulo gcd(50004,3357)=9 this is6. It is even, so it cannot belong to
{-4,-2,0,2,4}. For odd A=2a+1 the increment is

    50004B+25915−1865+3357(g(2a+2)−g(2a+1)).

The recursion proves0≤g(x)≤x by induction. Hence

    g(2a+2)−g(2a+1)=g(a+1)+a+1+((a+1) mod2)−g(a)≥1.

Both g values are even, so the difference is at least2. The entire increment
is therefore at least30764, much greater than4. Interchanging A,B gives an
even-to-odd increment congruent3 modulo9 and an odd-to-even increment at
least25954. This proves the claim. QED.

This is an obstruction to a particular *parameter-grid overlap argument*,
not a metagraph disconnection theorem and not an omission theorem for (1).
Parameter adjacency is not assumed to be a single tournament tile move.
The point is that even a primitive bank with every odd residue modulo6714
cannot be filled by simply joining neighboring integer-argument banks.

## 2. Retain both carries instead of moving one argument

For k≥1 and d≥0 define

    D_k(d)=(2^(k+1)−3)d+2^k−(d mod2).

Then

    g(2^k(2d+1))−g(2^k(2d+1)−1)=D_k(d).                 (2)

Proof. The lower word consists of the even prefix2d followed by k ones,
so its g value is g(2d). The upper word consists of the odd prefix2d+1
followed by k zeros, hence its g value is

    g(2d+1)+(2^k−1)(2d+1)+1.

Since g(2d)=g(d)+d+(d mod2) and g(2d+1)=g(d), subtraction gives(2).
The proof includes d=0 and k=1.

Take A odd and B positive even, and write uniquely

    A+1=2^k(2d+1), B=2^l(2e+1).

The compensated diagonal move (A,B)→(A+1,B−1) has exact response

    Delta=50004(B−A−1)+4810+3357(D_k(d)−D_l(e)).          (3)

This pays the changed product term and both parity fees. Omitting the
−50004 in B−A−1 would invalidate the construction.

More generally, for a sandwich
c+h(g(A)+g(B))+qAB+uA+vB+p(A mod2)+s(B mod2), put

    a_k=(h−q)2^(k+1)−3h,
    K=u−v−p+s.

Fix parities d=2r+alpha, e=2t+beta. Formula(3) becomes

    Delta=2a_k r−2a_l t+C,
    C=(q−h)(2^l−2^k)+K+a_k alpha−a_l beta+h(beta−alpha).  (4)

When q>h>0, both a_k,a_l are negative. Therefore a desired even response
Delta has a nonnegative integer solution r,t iff

    2 gcd(a_k,a_l) divides Delta−C.                       (5)

Indeed the divisibility is necessary and sufficient for an integer solution.
The homogeneous positive direction
(|a_l|/gcd(a_k,a_l), |a_k|/gcd(a_k,a_l)) moves any solution into the
nonnegative quadrant and then yields infinitely many. All resulting A,B
remain nonnegative actual physical inputs. This is an exact carry-stratum
existence criterion, not merely a modular heuristic: the same Diophantine
solution gives the integer operands and constructs the tournament.
It does not say an arbitrary absolute F value is attainable.

## 3. An explicit infinite four-count bank

For the C9 packet, choose k=3,l=2,alpha=1,beta=0, and Delta=4.
One solution is d=1845,e=3642, hence A=29527,B=29140.
Here a_3=−756423, a_2=−383247 and gcd=9. The positive solution direction
in the half-variables r,t is(42583,84047). Thus for every integer j≥0 set

    A_j=29527+1362656j,
    B_j=29140+1344752j.

Their carry lengths and parities stay fixed, and(4) proves exactly

    F(A_j+1,B_j−1)=F(A_j,B_j)+4.                          (6)

Consequently four actual tournaments have counts

    N_j, N_j+2, N_j+4, N_j+6,  where N_j=F(A_j,B_j).

At j=0, N_0 is recorded by the replay. Both operand pairs are produced by
finite binary words; the two choices of the central cut produce the+2
inside each pair. The two operand pairs generally have different rooted
interfaces: these four outputs are NOT claimed to form a single central
interface bank, and a continuation must retain which pair supplied its input.
There is no need to infer existence from a residue class,
normalize away a common factor, or assume independent physical subtraction.

This is a genuine exact overlap/abutment ingredient: two central intervals
of width2 join to an interval of width6 by jointly compensating two carries.
It improves on isolated+2 siblings. It neither joins the N_j across j nor
supplies a widening return map. N_j grows at least quadratically with j
because its positive product term does, so these constant-width intervals
alone cannot be claimed to fill a tail.

## 4. The four-point bank is locally optimal in this parameter geometry

Define a graph on all pairs(A,B) of nonnegative integers. Join two pairs
only when both coordinate changes have absolute value at most1, the pairs
are different, and their two-point central banks overlap or abut in the
odd-integer lattice; equivalently their lower outputs differ by at most4.

**Theorem.** This entire infinite graph is a matching together with isolated
vertices. Every edge is (A,B)→(A+1,B−1), with A odd and B positive even,
and its directed count change is exactly4. The family in Section3 supplies
infinitely many such edges. Consequently every connected component of this
specific nearest-neighbor parameter-overlap graph supplies at most four
consecutive odd outputs, and that bound is attained.

Proof. Section1 rules out horizontal and vertical edges. For a same-sign
diagonal increment, the elementary lower bound

    g(x+1)−g(x)≥−x/2−1

and the paid product term give

    F(A+1,B+1)−F(A,B)
      ≥(50004−3357/2)(A+B)+86580>4.

This rules out both same-sign directions by reversal. For the remaining
direction(A,B)→(A+1,B−1), reduce its response modulo9. The product term
and both g terms vanish. According to the parities(A mod2,B mod2), the
four residuals are

    (0,0):8, (0,1):3, (1,0):4, (1,1):8.

Among the five allowed changes{-4,-2,0,2,4}, only+4 appears in that list,
and only in parity(1,0). Its endpoint has parity(0,1), so it cannot continue
in the same direction; its only possible reverse-direction neighbor is the
one it came from. No other directions survived. QED.

In particular the tempting two-step zigzag(+1,−1),(+1,+1) cannot thicken
this interval: its second move is excluded by the explicit same-sign invoice.
The opposite zigzag is excluded likewise. This is not a universal bound on
central intervals, and not a tournament-metagraph component theorem.
Nonlocal parameter jumps, another central core, another consumer, or a wider
initial menu are outside the graph just classified. The theorem identifies
precisely which operation must change if this overlap mechanism is to grow.

## 5. What the mechanism changes

The one-dimensional parameter move failed because it retained a local
mod9 obstruction that the global gcd had concealed. Keeping both actual
carry lengths converts the nonlinear digital corrections into a linear
Diophantine equation on each fixed stratum. The positive homogeneous
solution direction realizes infinitely many integer solutions, rather than
just demonstrating a lattice gcd. This is the distinction between a global
congruence observation and a constructive stratum theorem.

The earlier repo's Zeckendorf/carry theme and the digits lane's central
interface supply the mechanism-level analogy: a carry should be classified
and used, not erased. No external universality theorem is imported. The new
remaining task is to synchronize multiple such strata so their absolute
intervals join and their width grows. The independent product term qAB
must stay in every attempted synchronization. Section4 shows that nearest-
neighbor synchronization in this particular input plane cannot do it; a new
operation is genuinely required.

## Reproduction

Run `python3 -B .../carry/check.py`. It checks(2) on independent recursion
values, the all-input nearest-grid obstruction over a finite control box,
and arbitrary-precision instances of the all-j family. The all-order claims
rest on the displayed recurrences and Diophantine proof. `probe.py` records
exploratory discovery; no minimal-witness claim depends on its bounds.
