# A 92-step repair of the Collatz bank's original-height gap

19 October 2026. Status: PROVED from the stated affine-bank identities and a small, exhaustive exact offset certificate. The certificate covers 533 possible endpoint offsets, each on an infinite lawful source cylinder. It is not a census of integers or of bank words.

## Result in plain language

Take any fully lawful word in the existing expanding bank {K,V}, followed by a nonempty word in the reset bank {F,H}. Suppose its total affine multiplier is less than one. There is no restriction on its length or on the order within either bank.

Either the original source has an actual smaller iterate by at most 92 steps after this word ends, or the word is an exact return to its source. In the latter case the source is the minimum of a positive cycle, and the independently lawful adjacent companion proves that a different positive cycle contains the strictly smaller integer (x+1)/4.

Consequently the least positive integer failing to reach 1 cannot start any such contracting bank word. This conclusion is unconditional as an exclusion of a possible least counterexample. It does not assert that every integer starts or later enters the bank. Nor does it assert that an entry reached after growth repairs the height of the original, pre-entry source.

The previous endpoint zero-head question remains open: an endpoint could in principle be x+d with 1<=d<=533. The new result proves that every such positive overshoot is repaired by the actual subsequent orbit. It therefore bypasses the endpoint question for the first-descent goal.

## 1. The map, the banks, and their guards

Use the Terras map T(x)=x/2 for even x and T(x)=(3x+1)/2 for odd x. All lengths below count Terras steps, not bank blocks. Parity words are written in chronological order.

The expanding branch K has lower and upper words

    K:  110110110110110101011011101010
    K+: 001111111111111101010110000100.

Write

    P=3^19=1162261467, Q=2^30=1073741824,
    C0=4912087879, D=P-Q=88519643.

Its map and complete lower cylinder are

    K(x)=(P x+C0)/Q,
    x=815464443 (mod Q).

The upper word has carry C0+Q-P and cylinder 815464444 modulo Q. Thus the two independently guarded trajectories satisfy

    T^30(x+1)=K(x)+1.

The two reset branches are

    F: 10,       F+: 01,
       F(y)=(3y+1)/4,        y=1 (mod 4);
    H: 110101000, H+: 001101001,
       H(y)=(81y+101)/512,   y=43 (mod 512).

Their upper cylinders are one higher, and their upper carries are respectively 2 and 532. They likewise preserve the adjacent relation at their endpoints. Define V=K K K K F, with upper word K+ K+ K+ K+ F+. Its length is 122, and

    V(x)=(5474401089420219382077155933569751763 x
           +83829391903560839836419813590128166431)
         /5316911983139663491615228241121378304.

Its multiplier exceeds one. Every nonempty prefix of K and V has multiplier greater than one, so it strictly increases any positive source on which it is lawful. These finite prefix inequalities and all four independent companion guards are checked in the accompanying certificate.

A complete lawful word in this note consists of a nonempty expansion word over {K,V}, followed by a nonempty reset word over {F,H}. Actual source guards are imposed on the complete concatenated parity word. Applying the authenticated companion of each successive block propagates the adjacent relation. In particular, if the lower endpoint is y after N Terras steps, the actual upper endpoint is y+1 after the same N steps. This uses no coalescence theorem or probability argument.

## 2. The inherited carry envelope, re-derived

For an expanding affine map a x+b, put c=b/(a-1). The expanding bank has c<=C, where

    C=83829391903560839836419813590128166431
       /157489106280555890461927692448373459.

For a resetting map r x+d, put p=d/(1-r). The reset bank has p<=1. Composition inside the two banks therefore gives

    E(x)=A x+B,       A>1,   0<B<=C(A-1),
    R(x)=r x+d,       0<r<1, 0<=d<=1-r.

These inequalities follow by induction: appending an expander preserves B<=C(A-1), and appending a reset preserves d<=1-r. If the total multiplier s=rA is below one, the full affine map is

    T^N(x)=s x+e,
    0<e=rB+d <= C s+1-(C+1)r < C+1 <534.        (1)

Its actual endpoint and source are integers. Therefore any non-descending endpoint has the form

    T^N(x)=x+d,    d in {0,1,...,533}.           (2)

This finite displacement budget is the new use of the carry envelope. It holds regardless of how close s is to one, so no upper bound on N and no estimate for a near-resonant difference 2^N-3^h are needed.

## 3. A failed endpoint descent forces the K^4 prefix

An expansion word is either K, K^2, K^3, or begins with K^4: the first occurrence of V itself starts with K^4, and any earlier blocks are K's.

The first three cases already descend at their first reset. Indeed K increases positive real inputs, H(y)<F(y) for y>0, and

    F(K^l(x)) <= F(K^3(x)) < x   (l=1,2,3)

for all x>=815464443. The last inequality follows from the contracting slope of F K^3 and its exact inequality at x=815464443. Later complete F/H blocks do not increase positive integer height, so the final endpoint remains below x.

Thus any source satisfying (2) must start with K^4. Its exact lawful source class is

    x=r4 (mod 2^120),
    r4=588943162011782273059357097228697595.       (3)

In particular x>=r4>9000. To independently derive this residue, write the K^4 carry as C4=C0(P^4-Q^4)/(P-Q). Its cylinder is P^4 x+C4=0 modulo Q^4. Equivalently Dx+C0=0 modulo Q^4, since D and P are odd. Both equations authenticate (3) in the checker.

## 4. The finite offset lemma

**Lemma.** If x is any positive integer satisfying (3), then for every integer d with 1<=d<=533 there is a j<=92 such that

    T^j(x+d)<x.                                  (4)

Here is the entire certificate format. For each d, take the first j<=120 for which the first j parity bits of r4+d have multiplier numerator a=3^h less than m=2^j. Compute their carry b using

    (a,b,m)=(1,0,1) initially;
    append 0: (a,b,m) -> (a,b,2m);
    append 1: (a,b,m) -> (3a,3b+m,2m).

The exhaustive certificate verifies, for each of the 533 offsets,

    j<=92,
    a(r4+d)+b=0 (mod m),
    m divides 2^120,
    (a d+b)/(m-a)<9000.                           (5)

The congruence and parity-word bijection imply that every x in (3), not merely the anchor, follows these same j parity bits from x+d. Hence

    T^j(x+d)-x=[(a-m)x+a d+b]/m <0,

using x>=r4>9000 and (5). This proves (4).

For inspection, 526 offsets have j<=29. The seven longer certificates are

    d=4:   j=92
    d=76:  j=31
    d=172: j=67
    d=256: j=32
    d=288: j=40
    d=324: j=35
    d=500: j=59.

The largest threshold in (5) is exactly 116927/13<9000, attained at d=480, j=8. The certificate is compressed to 54 disjoint residue classes of offsets; their counts sum to 533. The complete words, carries, exponents, residue descriptions and class counts are in verification-receipt.json. The checker independently validates that every d in 1,...,533 occurs in exactly one class.

The first contracting prefix for d=4 really has length 92, so reducing this particular certificate's deadline to 91 would fail. The lemma makes no claim that 92 is the best bound obtainable with a different argument or additional source information.

Combining (2)-(5), every positive overshoot after the bank word gives an actual iterate below the original bank source by time N+92.

## 5. Zero displacement creates a smaller, distinct periodic obstruction

Suppose T^N(x)=x. The authenticated adjacent companion gives

    T^N(x+1)=x+1.                                 (6)

Its first two parity bits are 00, so

    z=T^2(x+1)=(x+1)/4<x,
    T^N(z)=T^2(T^N(x+1))=z.                        (7)

Since x>=r4, z>2. It is a periodic positive integer outside the trivial cycle {1,2}.

One can say more: the two cycles are distinct. Every nonempty proper lower prefix of the complete word is strictly above x. During the expansion segment this was proved in section 1. During the reset segment each complete block strictly decreases a value greater than 1. Also every proper F/H prefix, including its empty prefix, is above that block's final value whenever its input is greater than 1. For the latter assertion, subtract the full branch from each prefix branch: the slope difference is positive, and the difference at input 1 is nonnegative. There are only 2+9 such comparisons, checked exactly.

The completed reset boundaries decrease strictly to the final value x. Each intervening proper prefix lies above its next reset boundary. Thus no nonempty proper lower prefix is at or below x, namely T^t(x)>x for 1<=t<N. The lower trajectory is a cycle of least period N and minimum x. The companion's cycle contains z<x, so it cannot be that lower cycle.

This is an obstruction, not a proof that two nontrivial positive cycles cannot exist. It also does not turn the smaller periodic point into a smaller iterate of the lower trajectory.

## 6. Main theorem and exact scope

**Theorem.** Let x>0 follow any fully guarded word in the specified bank language, with nonempty expansion and reset segments and contracting total multiplier. Let N be its total number of Terras steps. Exactly one of the following conclusions applies:

1. Some T^t(x)<x with 1<=t<=N+92; or
2. T^N(x)=x, the lower trajectory is a positive cycle of least period N and minimum x, and a different nontrivial positive cycle contains z=(x+1)/4<x.

Proof: If the endpoint is below x, use t=N. The small expansion cases are covered by section 3. In the remaining cases apply (2). A positive d gives case 1 by the offset lemma. Zero d gives case 2. The strict minimum property in case 2 makes the cases disjoint.

**Corollary.** A least positive integer whose trajectory does not reach 1 cannot start a contracting word in this bank language.

If case 1 held, its smaller actual iterate would reach 1 by minimality, so x would too. If case 2 held, the smaller z>2 would be periodic and fail to reach 1, contradicting minimality. This is precisely where a least-counterexample hypothesis may be used. It cannot be transferred to an arbitrary later entry value y>x.

**Consequences and limits.**

- Every nonperiodic source in the specified contracting bank cylinders has finite stopping time, with the explicit upper bound N+92.
- Even the astronomically long bank words left outside the previous length obstruction are covered by this dichotomy.
- The result does not show that their endpoints themselves descend. The endpoint zero-head conjecture remains open.
- The result does not show that every positive integer enters these cylinders, or that every initially expanding bank path eventually has a lawful contracting reset.
- It does not exclude other positive cycles, prove global Collatz convergence, or replace the separate original-height problem when bank entry occurs after growth.

## 7. Reproduction, lineage, and literature boundary

From this directory run

    python3 verify_original_height_repair.py > /tmp/height-repair.json
    python3 -O verify_original_height_repair.py > /tmp/height-repair-O.json
    cmp /tmp/height-repair.json /tmp/height-repair-O.json
    cmp /tmp/height-repair.json verification-receipt.json

Both replays produce the same 5,539 explicit checks. They include direct lifted trajectory controls as an independent check on the symbolic certificates. Python assertions are not used, so optimization does not remove checks.

Repository inputs are the proved affine-bank envelopes in collatz_bank_reset_20261019.md; the forced-prefix argument in collatz_zero_head_length_obstruction_20261019.md; and the independent adjacent words in collatz_odd_adjacent_return_20261019.md, all under 05-knowledge/results. The present proof re-derives the specific interfaces it needs rather than invoking the old finite bank census.

The general question whether first coefficient contraction always equals first actual descent is Terras's coefficient-stopping-time conjecture. The scope and the existence of non-first-stopping paradoxical segments are discussed by Olivier Rozier and Claude Terracol, “Paradoxical behavior in Collatz sequences,” Discrete Mathematics 349 (2026), 115167, arXiv:2502.00948v5, https://arxiv.org/html/2502.00948v5 . The present theorem is a bank-specific repair with a periodic alternative; it is not a solution of that conjecture. No general novelty or priority claim is made for the elementary affine identities.
