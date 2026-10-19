# Gold partition and 2/5 balance under arbitrary upset/ideal collars

Date: 2026-10-19. Status: **PROVED**, with elementary hand counts and independent exact controls. The general Gold Partition and 1/3–2/3 conjectures remain **OPEN**.

## Main result and scope

An explicit nine-element core admits one fixed two-comparison Gold certificate under every nonnegative mixture of four boundary-flag measures. In particular the same comparisons work inside arbitrary exterior posets attached through a terminal upset on the left and an initial ideal on the right. These exterior sets can have arbitrarily many elements. The first comparison alone remains at least 2/5 balanced.

This strengthens the bounded single-port theorem in `gold_four_collar_20261019.md`. That companion gives a different strategy with the sharper uniform normalized Gold margin 109/1174 when each exterior exceptional set is a singleton. Here the exceptional sets have unbounded size; the Gold budget remains strict for every finite genuine context, but its normalized margin has infimum zero.

The final section proves that the entire unbounded seven-thin bridge-chain family from the companion avoids Peczarski's elementary majority/at-most-one-slave sufficient criterion. It is an all-length statement, with a finite explicit core certificate and an exact cross-bridge argument.

## 1. The conjecture and the linear certificate

For a finite nonchain poset P, let t_0=e(P) be its number of linear extensions. After one comparison and then a second comparison, let t_1,t_2 be the remaining extension counts. The Gold Partition Conjecture asks for two comparisons, fixed before the answers are known, such that

    t_0 >= t_1+t_2

for every possible answer path. If the first comparison completes the order, use t_2=1. We use the non-strict original formulation. The strict formulation printed in a later lexicographic-sum paper cannot hold even for a two-element antichain.

Primary source: Marcin Peczarski, *The Gold Partition Conjecture*, Order 23 (2006), 89–95, and his later account *Comments on the Golden Partition Conjecture*, Contributions to Discrete Mathematics 12 (2017), 106–109:
https://cdm.ucalgary.ca/article/download/62349/46785/178291

For a fixed pair of comparisons write the joint table as X=(x_ij), with the first answer selecting a row and the second a column. Then every realized path has deficit

    D_ij = t_0-t_1-t_2 = (opposite row total)-x_ij.

These are linear functions of extension weights. In particular a triangular table

    [[A,0],[B,C]]

is Gold exactly when

    max(B,C) <= A <= B+C.                                  (1)

The impossible upper-right path imposes no condition; its formal deficit B+C is nonnegative anyway.

## 2. Arbitrary terminal-upset and initial-ideal fibre theorem

Let B be a finite poset with a minimal element a and a distinct maximal element b. Let A,C be any finite posets. Choose an upset U of A and an ideal V of C; empty sets and empty exterior posets are allowed.

Start with the ordinal sum A<B<C. Make exactly the pairs u,a for u in U and b,v for v in V incomparable, retaining every other cross-relation. This is a poset: the upset condition prevents a retained relation x<a from being contradicted by an omitted predecessor relation, and the ideal condition gives the dual property. The extremality of a,b prevents any core vertex from restoring an omitted relation by transitivity.

For alpha in L(A), let L_U(alpha) be the length of its terminal consecutive run of elements belonging to U. For gamma in L(C), let L_V(gamma) be the length of its initial consecutive run belonging to V. Set

    r = E[L_U], s = E[L_V],
    I(sigma)=1{a first in sigma}, J(sigma)=1{b last in sigma}.

Both parameters are nonnegative. They are not bounded by one: choosing U=A to be a chain of length m gives r=m, and dually s can be arbitrarily large.

**Fibre theorem.** The number of extensions of the attached poset restricting to a particular sigma in L(B) is

    e(A)e(C) (1+r I(sigma))(1+s J(sigma)).                   (2)

**Proof.** Fix alpha,sigma,gamma. If a is not first in sigma, some other core vertex precedes it and must follow all of A; hence all of A precedes the entire core. If a is first, its only possible movement into alpha is across a terminal segment wholly contained in U. There are exactly 1+L_U(alpha) legal insertion gaps. Every other core element remains after all of A. The right side is dual: there are 1+L_V(gamma) choices precisely when b is last in sigma. These choices are independent because a and b are distinct and every element of A precedes every element of C. This describes each extension uniquely. Summing the product of the two choice counts over alpha and gamma proves (2). ∎

A singleton U={u} is an upset exactly when u is maximal, in which case r=Pr_A[u last]<=1. Thus the single-port four-corner kernel is a special case, rather than a separate assumption.

## 3. A certificate for the entire nonnegative flag cone

Let B have labels 0,...,8 and covers

    0<2<7<8,  3<4<5<8,  0<6,  1<6.

Its transitive predecessor masks are

    (0,0,1,0,8,24,3,5,189).

Take a=0 and b=6. Compare, in order,

    first: 5 versus 7;
    second: 4 versus 7.                                    (3)

The increasing answer is the first row or column. The upper-right cell is impossible, since 4<5.

For all extensions, extensions starting a, extensions ending b, and extensions satisfying both flags, respectively, the joint tables are

    F_0  = [[325,0],[201,138]], total 664;
    F_a  = [[112,0],[84,84]],   total 280;
    F_b  = [[80,0],[48,32]],    total 160;
    F_ab = [[28,0],[21,21]],    total 70.                    (4)

**Flag-cone theorem.** Give each core extension sigma the weight

    w(sigma)=alpha+beta I(sigma)+gamma J(sigma)+delta I(sigma)J(sigma),

where alpha,beta,gamma,delta are any nonnegative real numbers, not all zero. Then the fixed strategy (3) is Gold. Its first comparison has probability in [2/5,1/2], so is at least 2/5 balanced.

**Proof of Gold.** The table is alpha F_0+beta F_a+gamma F_b+delta F_ab. Its three nonzero entries are

    A=325alpha+112beta+80gamma+28delta,
    B=201alpha+84beta+48gamma+21delta,
    C=138alpha+84beta+32gamma+21delta.

The four formal deficits are

    D_00=B+C-A=14alpha+56beta+14delta,
    D_01=B+C=339alpha+168beta+80gamma+42delta,
    D_10=A-B=124alpha+28beta+32gamma+7delta,
    D_11=A-C=187alpha+28beta+48gamma+7delta.

Every coefficient is nonnegative, proving (1). This proof works for the whole four-dimensional nonnegative cone; product weights alpha,beta,gamma,delta=(1,r,s,rs) are only a subfamily. ∎

**Proof of balance.** The probability of the first increasing answer in the four basis measures is

    325/664, 112/280=2/5, 80/160=1/2, 28/70=2/5.

Normalization makes the mixture probability a convex combination of these four numbers. Thus it belongs to [2/5,1/2]. ∎

Applying (2) proves both claims for every arbitrary upset/ideal context. A genuine finite context has alpha=e(A)e(C)>0, so every realized Gold deficit is strictly positive and the first probability lies strictly between 2/5 and 1/2. There is no context-independent positive normalized Gold margin for this strategy: take r=0 and s tending to infinity. Then D_00/Z=14/(664+160s) tends to zero. The uniform 2/5 balance bound does survive.

### Elementary derivation of the four tables

Delete vertices 1 and 6. There are 20 shuffles of chains (0,2,7) and (3,4,5), followed by 8. Write t,q,h for the positions of 0,2,7 in the six-letter shuffle, so 1<=t<q<h<=6. Insert 1 in backbone gap j in {0,...,7}. The number of choices for 6, which must follow both 0 and 1, is

    W(t,j)=8-t if j<t, and 8-j if j>=t.

Summing over j gives K(t)=36-t(t+1)/2, with K(1),...,K(4)=35,33,30,26. The event 5<7 is h=6, and 4<7 is h>=5. Therefore the three cells are h=6, h=5, and h<=4. Their unflagged counts are

    A=4*35+3*33+2*30+26=325,
    B=3*35+2*33+30=201,
    C=35+2*35+33=138.

For a first, impose t=1 and j>=1. The slot sum is then 28 for each qualifying shuffle, and the three shuffle counts are 4,3,3. This gives 112,84,84. For b last, the slot weight becomes one, so the eight gaps times shuffle counts 10,6,4 give 80,48,32. Imposing both flags gives seven gaps times counts 4,3,3, hence 28,21,21. This proves (4) without a census.

## 4. All-length bridge chains avoid the elementary low-slave criterion

Take k copies B_1,...,B_k. Order every earlier block below every later block, except make b_i=6 in B_i incomparable with a_(i+1)=0 in B_(i+1). This is the bridge-chain family. Every block has singleton exceptional sets in its surrounding prefix and suffix, so both the present theorem and the sharper 109/1174 Gold-margin theorem apply for every k.

The companion proves width three, maximum incomparability degree seven, and ordinal indecomposability for every k. Here is an additional all-length obstruction to applying the elementary published sufficient condition.

For an ordered incomparable pair (x,y), a slave is an element z with either x<z and z incomparable with y, or z<y and z incomparable with x. Peczarski's sufficient condition requires an orientation with at most one slave and Pr(x<y)>=1/2. **No bridge chain in this family has such an orientation.**

**Internal pairs.** Any slave already in a block remains a slave in the full chain. Thus an internal orientation having at most one slave globally must be one of the following sixteen core orientations. The four numerators count x<y at context corners (r,s)=(0,0),(0,1),(1,0),(1,1); the corresponding denominators are 664,824,944,1174.

| Ordered pair | Core slaves | Four numerators |
|---|---:|---|
| (1,0) | 0 | 203,238,203,238 |
| (1,3) | 1 | 224,259,323,373 |
| (2,3) | 1 | 140,172,252,312 |
| (4,0) | 1 | 116,148,116,148 |
| (5,0) | 0 | 26,34,26,34 |
| (5,1) | 1 | 110,165,156,236 |
| (5,2) | 1 | 124,156,152,191 |
| (6,2) | 0 | 112,112,147,147 |
| (6,3) | 0 | 36,36,57,57 |
| (6,4) | 1 | 128,128,191,191 |
| (6,7) | 1 | 280,280,385,385 |
| (7,1) | 1 | 111,166,176,261 |
| (7,3) | 0 | 35,43,63,78 |
| (7,4) | 1 | 138,170,222,275 |
| (8,1) | 0 | 20,40,30,60 |
| (8,6) | 1 | 160,320,230,460 |

Every ratio is at most 35/83<1/2. This finite certificate is obtained directly from the displayed covers and the 20-shuffle insertion formula; `gold_four_collar.py` also checks it by independent actual-collar ideal counting. The list exhausts the core orientations with at most one slave. For singleton contexts r,s lie in [0,1], and each normalized probability is a convex combination of its four corner probabilities. Therefore none of these orientations is a majority in any bridge chain.

**Cross pairs.** The only cross-block incomparable pair is b_i,a_(i+1). Its forward orientation b_i<a_(i+1) has ten slaves: six vertices in B_i incomparable with b_i and below a_(i+1), and four vertices in B_(i+1) incomparable with a_(i+1) and above b_i. Every other vertex is comparable with both; no other slaves occur.

The reverse orientation a_(i+1)<b_i has zero slaves but is a strict minority. To see this, split the chain into the left prefix L and right suffix R at that bridge. All L<R except b=b_i and a=a_(i+1). Extensions with b<a number e(L)e(R). Extensions with a<b number e(L;b last)e(R;a first): they arise uniquely by swapping the adjacent b,a in a concatenation having the indicated boundary flags. The latter product is strictly smaller. Indeed b is not always last in L, since the same block has a distinct maximal element 8; and a is not always first in R, since the same block has another minimal element 1 (or 3). Consequently Pr(a<b)<1/2. This completes the all-length proof. ∎

The exclusion concerns this specific sufficient condition only. It is not a claim that every other published Gold criterion fails, nor that the family establishes Gold for all seven-thin posets.

## 5. Exact controls and continuation target

The mathematics above is the proof. The following controls use exact integers and fractions, no floating point and no disabled Python assertions:

- `gold_unbounded_collar.py` reproduces (4) independently by hand sums and core extension enumeration, then checks the fibre theorem against ideal counting on 400 actual arbitrary-upset/ideal contexts.
- It also checks 36 chain-collar pairs with r,s in {0,1,2,5,10,30}, reaching 69 vertices and explicitly leaving the bounded four-corner regime.
- `gold_four_collar.py` independently verifies all sixteen low-slave rows at the four actual corner posets, plus the sharper single-port certificate and every block of bridge chains through 12 copies.
- Normal and `python3 -O` runs produce identical outputs. Published JSON receipts are in `05-knowledge/results`; scripts print fresh receipts without modifying files.

The unresolved structural task is to find a broader class of cores whose four boundary-flag measures all lie in the triangular Gold cone (1), preferably via flag-preserving injections. The theorem here proves a fully arbitrary-context result for one explicit core; it does not assert that every nonchain poset contains such a core or admits such a boundary kernel.

Primary references:
- Peczarski's sufficient conditions and non-strict Gold formulation: https://cdm.ucalgary.ca/article/download/62349/46785/178291
- Equal-fibre lexicographic-sum context, useful for distinguishing the present unequal-fibre theorem: https://arxiv.org/html/2410.12494v1
