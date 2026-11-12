# Marked intervals, five recursive families, and the remaining coverage gap

**Status: PROVED IN THIS NOTE** for the displayed all-length families and
obstructions; **FINITE-EXACT** for the census; **OPEN** for its extrapolation.
These are same-agent proofs and independent computational formulations, not
an independent-agent audit or a literature-priority claim. Chan–Pak is open.

## 1. Inherited interface, with the endpoint information retained

Use the endpoint path-cover interface of
[seed-recursion births, PROOFS](../seed_recursion_births_20261111/PROOFS.md).
For a tournament T, H counts directed Hamilton paths. Let B_ij count ordered
two-path spanning covers whose first path ends at i and second starts at j.
Put W_ij=B_ij+B_ji, and U_j=s_j+sum_i B_ij, where s_j is the number of
Hamilton paths starting at j. For a nonempty selected vertex set O define

    delta(T,O)=H(T)−sum_{i<j in O} W_ij,
    U(T,O)=sum_{j in O} U_j.

The inherited tail/apex construction gives H_new=delta+2^t U for t>=1.
Thus delta=1 is a unit seed, and U=c2^j supplies the ray 1+c2^s only from
s=j+1 onward. Neither delta=1 alone nor a coefficient gcd fills lower heads
or supplies every odd c. The 95 strong-core choices and 49 stable pairs are
finite banks from the earlier packet, not uniform completeness theorems.

## 2. Complete endpoint classification in a three-interval chart

Let C3 point from block1 to2 to3 to1. Replace its blocks by transitive
tournaments of sizes (a+b+c+2,r,s), where a,b,c>=0 and (r,s) is (1,2) or
(2,1). In the long block select precisely the vertices at zero-based
positions a and a+b+1. Thus a is the unselected prefix length, b the gap,
and c the suffix length.

**Interval theorem.** This is a unit seed **if and only if b=0 or b=1**.
For b>=1 its exact defect is

    delta=6·2^(a+c)(2^(b+1)−3^b−1)+2^b−1,                 (1)

while b=0 has delta=1. For every b>=0 the slopes are

    U_(1,2)=6·2^a3^c(2·2^b+3·3^b)−[a=0],
    U_(2,1)=10·2^a3^c(2·2^b+3·3^b)−3[a=0].             (2)

Here [a=0] is an indicator; in particular that boundary correction is
constant in c. Omitting it gives wrong first-step recursion and wrong
aggregate generating functions.

Proof of the formulas: split each path at each maximal visit to a block.
If a free interval of length m is visited j times, assigning its ordered
vertices to nonempty visits contributes j! S(m,j). This is the inherited
walk-content substitution mechanism, **THM-3121 /
[path-cover-walk-content-substitution-kernel](../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md)**.
In this chart the fixed vertices bound every free interval's number of
visits. Enumerating these finite directed skeletons gives the complete
coefficient tensors in `c3_pair_family.json`, for all eight zero/positive
boundary strata and both orders of the small blocks. Substitution of
j!S(m,j)=sum_(d=1)^j (−1)^(j−d) binom(j,d)d^m and collection gives (1)–(2).
The companion compares the entire tensors, not a finite numerical fit, and
also checks their evaluations against vertex-level path counting.

For b=1, (1) is 1. For b>=2, set K_b=3^b+1−2^(b+1).
K_2=2 and K_(b+1)−2K_b=3^b−1>0, so K_b>=2^(b−1).
Consequently delta<=−6·2^(b−1)+2^b−1<1. This proves necessity
at every a,c, including the zero-length boundaries. Sufficiency at b=0
comes from the corresponding zero-gap tensor. ∎

For the two allowed gaps write U=K2^a3^c−r0[a=0]. The four rows are

| Small blocks | Gap | K | r0 |
|---|---:|---:|---:|
|1,2|0|30|1|
|1,2|1|78|1|
|2,1|0|50|3|
|2,1|1|130|3|

This supplies literal seed-preserving recursive moves:

* Append to the suffix: U'=3U+2r0 if a=0, and U'=3U otherwise.
* Add the first prefix vertex: U'=2U+2r0. Every further prefix vertex doubles U.
* Insert the first gap vertex: U'+r0[a=0]=(13/5)(U+r0[a=0]).
  A second gap vertex destroys delta=1.

The 5 and13 in the last line are 2+3 and2²+3² from the same path-cover
kernel. The next value 2³+3³=35 occurs in the slope formula but fails the
defect constraint. This is a structural connection with an exact boundary.
It does not identify these numbers with dimensions of elastic stresses.

The old U29(m)=10·3^m−1 and U47(m)=50·3^(m−1)−3 rays are
special cases of the zero-gap, zero-prefix rows. The present contribution
is the complete prefix/gap/suffix chart and its iff condition. At a=c=1,
the zero-gap chart is an order7 stable pair; inserting the allowed gap gives
an order8 stable pair with H=95 and U=468 or780. Thus a concrete move
transports an old small-core stable pair beyond its original core-size bank.

Together with the inherited book C3[TT_L,1,1], selected on the two singleton
blocks and having U=3·2^L, these are **five explicit infinite families**.
Transitive tournaments with a singleton mark form a sixth, elementary family.

## 3. A five-state recognition language and new sequences

The long-block marking words are exactly

    0* 1 (epsilon | 0) 1 0* .                            (3)

Their minimal complete deterministic automaton has five states: before the
first mark, just after it, after one intervening zero, after the second
mark, and a dead state. The four live states and the dead state have different
continuation languages; `paper_transfers.py` checks the minimal partition.
This is a precise additional species of five, generated by the gap theorem.

At order n>=5, the two small-block orientations and two allowed gaps supply
4n−18 distinct marked types. The book adds one and the transitive family n,
giving **5n−17**. The orientations are not identified by reversal: reversing
all arcs is not a tournament isomorphism in the census.

Let N(z) count this union of proved families, and M(z) sum U over it. Then

    N(z)=z/(1−z)^2 + z^3/(1−z) + 2z^5(1+z)/(1−z)^2,

    M(z)=z/((1−z)(1−2z)) + 6z^3/(1−2z)
         + z^5[(80+208z)/((1−2z)(1−3z))−4(1+z)/(1−z)].   (4)

These are proved all-size laws for the specified union. For M the recurrent
modes are 1,2,3, so its tail obeys

    M_n=6M_(n−1)−11M_(n−2)+6M_(n−3).

Here11 is the elementary symmetric sum 1·2+1·3+2·3. The same numeral in a
prime list or an automaton dimension is a different assertion until a map
is supplied. At n10 the family count is33 by the affine law, not by a claim
that33 universally controls metagraphs.

The complete finite round census agrees with these families plus eleven
exceptional marked types. Its seven C3 exceptions have (n,H,U):

    (4,5,16), (4,5,11), (5,11,37),
    (6,23,111), (6,23,99), (8,107,841), (8,107,596).

Four further C5 exceptions have

    (9,705,6288), (9,705,2718),
    (11,54123,358600), (11,54123,236287).

Actual adjacency rows, selected masks, automorphism quotients and first
proof controls are retained in the receipts. The finite exception count
polynomial and U-moment polynomial are respectively

    E(z)=2z4+z5+2z6+2z8+2z9+2z11,
    E_U(z)=27z4+37z5+210z6+1437z8+9006z9+594887z11.        (5)

**OPEN round-seed classification hypothesis:** the six proved families and
these eleven exceptions exhaust all locally transitive unit seeds. It
predicts the exact total generating functions N+E and M+E_U, and no
reduced cyclic core of order7 or larger. The census supports it only to its
declared finite order. It is stronger than the evidence; the next cheap
hostiles are interior and multiple marks in larger carousel blocks, then
multi-block inflations. Neither an eventual linear count nor a finite
exception theorem has been proved for the whole round universe.

## 4. All-length obstructions, not just failed small searches

Inflate one vertex of the regular carousel of order k=5,7,9,11. Select
any subset of the other singleton blocks and either no vertex of the long
block or exactly its first or last vertex. **No such chart contains a unit
seed at any positive length.** For k3 the book and two isolated endpoint
choices survive.

Proof certificates are in `carousel_clock_certificates.json.gz`. In each
chart the finite walk skeleton gives delta−1=sum_b c_b b^m. If nonzero,
the largest base dominates beyond an explicit integer cutoff; every smaller
m is tested exactly. The maximum cutoffs for k5,7,9,11 are6,14,24,35.
Every chart has literal length1 and length2 controls against vertex-level
counting. This is an all-length conclusion for the specified finite cores,
not an extrapolation in k. Interior and multiple marks are excluded.

The four C5 exceptions have no individually stable unselected port, and
none of their 24 first/last selected-clone gates is stable at all lengths.
Their complete defect tensors show failure already on the first extra clone.

There is also a complete two-clock classification in the C5 chart with
word0001010, free positions0,2 and selected position1 plus two of3,4,5,6.
For positive free lengths(a,b), the only unit solutions are

    O={1,3,5}: (a,b)=(1,3), U=6288;
    O={1,4,6}: (a,b)=(3,1), U=2718.                      (6)

Here is a short audit beyond the recorded40-by40 search. Set x=2^a,y=3^a,
u=2^b,v=3^b. The six defect polynomials in pair order34,35,36,45,46,56
are, with the last and fifth obtained by interchanging a,b in the first
and second,

    D34=−5−2u−6x+12xu−8xv−6yu,
    D35=−5−12u+4x+4xu−4xv−10yu+6yv,
    D36=−5−4u−4x+4xu−6xv−6yu+6yv,
    D45=−7−4u−4x+14xu−8xv−8yu.

D34,D45 are negative because v/u,y/x>=3/2. D35 at b1 equals
−29−2·3^a; at a1 it is3−34·2^b+10·3^b, equal to1 only at b3
(b1,b2 are below; b4 is above, and the positive region then increases).
For a,b>=2 it is at least D35(2,2)=9: the increment in a after subtracting
twice the preceding value is 3^a(6v−10u)+12u+5>0; at a2 the analogous
b increment is38·3^b−11>0. D36 at a1 is−13−14·2^b+6·3^b,
giving−23,−15,37 at b1,2,3 and increasing thereafter. In the quadrant
a,b>=2 its minimum is81, since its corresponding increments are positive.
Symmetry completes the six cases. Thus (6) consists of isolated resonances,
not two new freely extendable ports.

## 5. Four coefficient rays now have their entire lower heads

The four exceptional slopes have odd parts393,1359,44825,236287 and
2-adic valuations4,1,3,0. The inherited tail construction plus the explicit
eight heads in `heads.json` proves

    for c in {393,1359,44825,236287}, every 1+c2^s, s>=1,
    is the Hamilton-path count of an actual tournament.                 (7)

| c | First exponent supplied by the seed | Added finite head counts |
|---:|---:|---|
|393|5|787,1573,3145,6289|
|1359|2|2719|
|44825|4|89651,179301,358601|
|236287|1|none needed|

The heads use explicit one-apex transitive tournaments and ordered joins
of verified factors. Large one-apex counts are checked by a separate
two-endpoint formula; small factors also have vertex-level dynamic-program
controls. This does not enlarge a previously proved general bounded range
and does not assert that each c itself is an exact odd slope U.

Finally the odd parts supplied by the proved infinite families are contained
in

    {1,3}, {15,39,25,65}·3^c,
    30·3^c−1, 78·3^c−1, 50·3^c−3, 130·3^c−3, c>=0.

Only O(log X) of these are <=X. Hence their odd-coefficient supply has
density zero. This is a proved obstruction for these families, stronger
than saying that the finite search has gaps. A proof of every odd seed
must introduce another operation or leave this union. No unrestricted
nonexistence claim is made for coefficient1039 or any other missing target.
