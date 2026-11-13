# The exact two-run frame cone

**PROVED, independently audited by the root session; finite checks passed.**
This is a scoped research-note result, not a canon promotion. It concerns exact
Hamiltonian counts in uniform module frames, not arbitrary port operations.

Let (T,r) be a rooted tournament. Let a_j count legal words using every
vertex other than r once and r exactly j times, with no adjacent r copies.
Thus a_1=H(T), and uniform substitution has full count sum_j a_j K_j(M).
Call the frame two-run when a_j=0 for all j>=3. The matching-star sufficient
criterion and this run expansion are inherited from
[CORE_REPLACEMENTS](../../prime_ternary_lift_20261113/CORE_REPLACEMENTS.md).
The converse and weighted path comparison below are proved explicitly.

## C1. Matching-star equivalence

Put I=N^-(r), O=N^+(r). The frame is two-run if and only if the bipartite
graph of arcs O->I has matching number at most one.

If a word has three or more r copies, its first two internal excursions
from r to r contain two vertex-disjoint O->I arcs. Conversely, two disjoint
arcs o1->i1 and o2->i2 give the legal partial word

    r,o1,i1,r,o2,i2,r.

Insert every remaining exterior vertex by the ordinary path-insertion rule:
place it immediately before the first word entry it beats, or at the end
if it beats none. The preceding entry, when present, beats the inserted
vertex. This argument is valid even though r repeats, and insertion cannot
create adjacent r copies. The resulting full word witnesses a_3>0.

A nonempty bipartite graph with matching number one is a star. If r belongs
to a directed cycle, O->I is nonempty. Passing to the converse tournament
preserves every a_j by word reversal, so it suffices to handle a hub w in O:
every O->I arc has tail w. In particular all arcs I->(O minus {w}) point
forward, and w beats at least one vertex of I. The argument does not require
the whole tournament to be strong; cyclicity of r suffices.

## C2. Four exterior-word classes

Consider permutations of V(T) minus {r}. The following four classes are
disjoint; Q and D letters denote cardinalities.

* Q0: Hamilton paths of T-r with no O->I transition.
* Q1: Hamilton paths of T-r with one O->I transition.
* D0: words whose only illegal adjacency is i,w with i in I and w->i,
  and which have no O->I transition.
* D1: words with that unique illegal adjacency and one O->I transition.

There cannot be two legal O->I transitions: both would have to leave the
single vertex w. For a binary I/O word, the number of insertion slots for r
is

    #(I->O transitions) + [first entry in O] + [last entry in I]
      = 1 + #(O->I transitions).                         (11)

Deleting r from a Hamilton path leaves either a legal exterior path or a
word with the unique illegal adjacency i,w. Any other adjacency bridged by
r has the form I->O and is already legal unless its O endpoint is w. Thus
the four classes exhaust the relevant exterior words. In a D word the
illegal slot must receive a copy of r. Applying (11), for one or two copies,
gives the exact identities

    H(T)=Q0+2Q1+D0+D1,
    a_2=Q1+D1.                                         (12)

Let B=H(T[I])H(T[O]), with the usual empty convention. Here I and O are
nonempty because r is cyclic. Concatenate an arbitrary Hamilton path PI
of I and PO of O. Its boundary is legal unless PO starts at w and w beats
the last vertex of PI. Hence

    Q0+D0=B.                                           (13)

This also proves the reverse correspondence: a zero-transition exterior
word has all I vertices before all O vertices.

## C3. A recoverable injection into Q1

Fix a pair (PI,PO) counted by B. Split PI immediately before its first
vertex beaten by w:

    PI=Ipre,Ipost,       PO=Opre,w,Opost.

Ipost is nonempty, since w has an out-neighbour in I; each vertex of Ipre
beats w. Send the pair to

    Ipre,Opre,w,Ipost,Opost.                            (14)

Every adjacency is legal. If Ipre and Opre are nonempty, Ipre->Opre is
forward; if Opre is empty, the last vertex of Ipre beats w. The adjacency
Opre->w comes from PO, w beats the first vertex of Ipost, and Ipost->Opost
is forward when Opost is nonempty. Empty Ipre, Opre or Opost causes no
other boundary condition. The result has exactly one O->I transition,
namely w followed by Ipost, and therefore lies in Q1.

Projection onto the I vertices recovers PI and projection onto the O
vertices recovers PO. Thus (14) is injective and

    Q1>=B>0.                                           (15)

## C4. Cone and equality boundary

Combining (12)--(15),

    H(T)=a_2+B+Q1,
    a_2+B+Q1 <= a_2+2Q1 <= 3a_2.

Consequently every cyclic two-run frame satisfies

    a_2+2 <= H(T) <= 3a_2,                              (16)

in particular a_2<H(T)<=3a_2. The upper equality holds if and only if
Q1=B and D1=0. The stronger lower inequality uses B>=1 and Q1>=B.

For the coefficient notation H(T[r<-M])=a p+b h, with
a=a_2, b=H(T)-2a_2, this gives

    -a+2 <= b <= a.                                    (17)

The upper boundary is attained by a directed triangle rooted at any
vertex: H=3, a=1. For the lower boundary, take the directed three-cycle
with ordered blocks {r}, TT_k, {w}, k>=1, and inflate r by TT_m.
Its full count is

    (2^m-1)(2^k-1)+2,

so the frame has H=2^k+1 and a=2^k-1. One direct check uses the positive
odd-cycle formula: every directed cycle uses w and nonempty subsets of
both transitive blocks; it is odd exactly when their total chosen size
is even. There are ((2^m-1)(2^k-1)+1)/2 such cycles and none are disjoint.
Thus the additive lower bound in (16) is sharp for every k>=1.

If r is acyclic, every a_j for j>=2 is zero, and the response is simply
H(T) H(M). No positive-a exception is hidden in this case.

The cone is a necessary condition for frame existence. It is not asserted
that every lattice point satisfying (17) occurs. The previous empirical
order-eight cone is now justified by this independently audited all-order
argument. An exact-statement search in the inherited canon and
the routed current packets found no prior claim of (16); the run/matching
mechanism itself remains inherited.

## C5. An all-frame obstruction for the primitive small-debt family

Let Q=(2,4,9,3), M_s=Q->...->Q with s copies, s>=3, and put

    h=5^s, p=32^s, S=15p-h.

These are the actual modules and source from [PROOF, T3--T4](PROOF.md).
There is **no tournament containing M_s as a uniform module** with either
target count 3S+2 or 3S+4. Equivalently, no rooted tournament frame B of any
order can give either target by uniform substitution of this same M_s.
The exterior size is unrestricted: contracting the module supplies exactly
such a frame, and conversely every such inflation has the module.

**Proof.** Fix delta in {2,4}, and suppose a target frame exists. If its root
is acyclic, its output is a multiple of h and hence of five. But S=15p-h
is divisible by five, so 3S+delta is congruent to delta modulo five and is
not divisible by five. Thus the root must be cyclic.

The third ordered-cover response is

    K_3(M_s)=117^s-3*32^s+3*5^s.

For s>=3, (117/32)^s>48, so K_3(M_s)>45p>45p-3h+delta.
If B permits three or more visits, C1 gives a_3>=1; the nonnegative
run expansion then already exceeds the target in that single term.
Therefore B is necessarily two-run, even though no run bound or frame-size
bound was imposed on it in the statement.

Write its response as a p+b h. C4 gives -a+2<=b<=a. Since p>47h+4 for
s>=3, a<=44 would give

    a p+b h <=44p+44h <45p-3h+delta,

whereas a>=46 would give

    a p+b h >=46p-44h >45p-3h+delta.

Hence a=45, and equality with the target forces (b+3)h=delta, impossible
because h>=125. QED.

The first pair blocked this way is 1474187 and 1474189, from source
S=491395 at s=3; K_3(M_3)=1503684. This is an all-order restriction on
retaining M_s as one uniform module, not an assertion that either integer
is unattainable by an unrelated tournament. The cover content of every M_s
is one, so this obstruction is not the inherited fixed-content divisor.

This is the decisive outcome of the variable-correction-core probe for this
family: an arbitrarily large exterior alone cannot repair its small debt.
A successful construction must alter or split the module, or select another
source. The root session's shifted transitive-module construction is outside
this hypothesis, as it should be.

**Smaller parameters.** For s=2, the cone already excludes every two-run
target frame for the counts 46007 and 46009: its integer coefficient a
must lie between 44 and 46, whereas reduction modulo h=25 forces a=43
or 41 modulo 25, respectively. But K_3(M_2)=10692 is below these targets,
so the all-rank argument does not cover s=2. For s=1 the cone permits
three pairs (a,b) for each digit: (41,23),(46,-9),(51,-41) for delta=2,
and (42,17),(47,-15),(52,-47) for delta=4. These are necessary arithmetic
candidates, not claims that corresponding frames exist. No all-frame
absence claim is made for s=1 or s=2.

**Draft correction lineage.** The first C5 draft incorrectly reduced its
target to delta modulo h. That would require h to divide 45p and fails
already at s=3. The root reviewer caught this before session integration.
The correct and sufficient obstruction is modulo five, as written above:
five divides both h and S, but does not divide delta. The two-run
coefficient comparison and the final equation (b+3)h=delta were unaffected.

## Verification

The companion `cone_check.py` independently classifies the four exterior
word types, counts the r-repeated skeletons, verifies the injection, and
compares a_1 with literal Hamilton-path enumeration. Its universe is all
rooted labelled tournaments through order five satisfying the cyclic
two-run guard, plus the three named incoming frames. These checks corroborate
the universal proof; they do not substitute for it.
