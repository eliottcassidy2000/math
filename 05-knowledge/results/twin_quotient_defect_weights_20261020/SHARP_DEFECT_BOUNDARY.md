# A sharp order-twelve boundary for twin-packing defect cancellation

**Status: PROVED**, from the weighted false-twin formula in
`TWIN_QUOTIENT_DEFECTS.md`, with independent literal-packing and quotient
controls. 20 October 2026. The global defect sum is already zero by the
inherited prime-local theorem; the open problem is an explicit structural
cancellation, not its numerical existence.

Write delta(X)=b(X)-e(X). If X has a twin triple, delta=-D, with D the signed
cycle-index average in the companion note. Otherwise delta=0. A class is
called big if its false-twin size c is at least three, and k=floor(c/3).

## 1. Two useful general deductions

**Vanishing test.** If, for some size c with floor(c/3) odd, the number of
false-twin classes of exactly that size is odd, then D(X)=0.

Every weight-preserving quotient automorphism permutes those classes. An odd
number of them forces an odd-length cycle, whose factor F(k,ell) is zero.
Thus every summand in the cycle-index average vanishes separately.

**One-sided sign theorem.** Every globally sign-even graph satisfies
delta(X)<=0. Indeed chi=1 and all factors F(k,ell) are nonnegative, so D>=0.
The opposite implication is false: globally sign-odd graphs can have either
sign of defect, with unbounded magnitude, as proved below.

## 2. The sharp bound through eleven vertices

**Theorem.** Every graph on at most eleven vertices has

    delta(X) in {-1,0,1}.

Moreover a globally sign-odd graph in this range has delta(X) in {0,1}.

**Proof.** No twin triple means zero defect. Otherwise there are at most
three big classes. Small classes have k=0 and contribute factor one in every
cycle, even if their automorphism action is complicated.

- One big class is fixed by every weighted quotient automorphism. Its factor
  is zero or one, so D is either zero or the average of chi. Thus D=0 or 1,
  and for a globally odd graph D=0.
- With two big classes of different sizes both are fixed, and the same
  argument applies. With equal sizes c, the bound on n gives c<=5, so both
  have k=1. A nonzero summand must interchange these two classes, contributing
  a factor two. If no interchange exists, D=0. Otherwise the interchanges
  form one coset of an index-two kernel M. If chi is nontrivial on M their
  signed sum is zero. If chi is trivial on M, D=chi(s) for an interchange s,
  hence is +1 or -1. In a globally odd graph only the -1 case is possible.
- With three big classes, each size is at most five, so all three k-values
  are one. Every permutation of three classes has an odd-length cycle.
  Its factor vanishes, giving D=0.

These cases are exhaustive, and delta=-D. The argument does not assume a
bounded graph database or a bound on the size of a small-class automorphism
group. QED.

## 3. Complete first-failure classification at twelve

**Theorem.** At order twelve the only graph types with |delta|>1 are

    K_(6,6):                              delta=-2;
    P4 with each vertex replaced by 3
    independent vertices:                 delta=+2.

The only globally sign-odd type with negative defect at order twelve is

    K_(3,3) disjoint union K_(3,3):         delta=-1.

**Proof.** The earlier one- and two-big-class arguments still apply except
for two classes of size six. Those exhaust all vertices, must be adjacent
to be distinct false-twin classes, and give K_(6,6), with D=(1+3)/2=2.

For three big classes, every k is one unless the sizes are (6,3,3). In the
all-one case an odd cycle again forces zero. In the exceptional case the
six-class is fixed and contributes one, and the same two-class coset
argument bounds |D| by one. It cannot give positive D on a globally odd
graph: a nontrivial sign on the kernel cancels the sum; otherwise a nontrivial
character is negative on the interchange coset.

Four big classes must all have size three, with no other vertices. The
unweighted quotient has four vertices and no false twins. Up to isomorphism
the complete list is

    2K2, P4, K3 disjoint union K1, the paw, K4.

This list follows by edge count: at two edges only the matching qualifies;
at three edges the path and triangle-plus-isolated qualify but the star
does not; at four edges only the paw qualifies; at zero, one, and five
edges there are false twins; at six edges the complete graph qualifies.

For the triangle-plus-isolated and the paw a quotient vertex is fixed by
every automorphism, so its k=1 factor kills every term. For K4 every subset
of the four twin classes has a stabilizer containing a transposition, which
is odd; there are no sign-trivial packing-vector orbits, so D=0.

For P4 the quotient group consists of the identity and reversal. The identity
contributes zero. Reversal has two 2-cycles and negative graph sign, giving
D=-4/2=-2. The expanded graph is globally odd, so delta=2.

For 2K2 the expansion is 2K_(3,3). A packing is sign-trivial on its stabilizer
only if it chooses exactly one of the two twin triples in each component:
otherwise that component's part swap remains in the stabilizer and is odd.
All the valid packings form one orbit, with two triples and contribution -1
to b. A component swap is even, so this orbit is indeed eligible. The graph
is globally odd and e=0, giving delta=-1 and D=1.

All alternatives have now been covered. QED.

The independent audit exhausts all 64 labelled four-vertex quotients and
retains exactly the 32 with distinct open neighborhoods. Its D-values agree
with the five cases above. This finite check supports, rather than replaces,
the case proof.

## 4. Negative defects of odd graphs are unbounded

If U and V are nonisomorphic connected graphs, each with at least one edge,
then

    D(U disjoint union V)=D(U)D(V).                         (1)

No false twins can cross the two components, since every vertex has a
nonempty neighborhood in its own component. Automorphisms preserve the
nonisomorphic components, so both the quotient group and its character
factor. The cycle-index average is a product, proving (1).

For r>=1 let

    U=K_(3,3),  V_r=K_(6r+3,6r+3),  Z_r=U disjoint union V_r.

The balanced-bipartite formula gives D(U)=-1 and D(V_r)=-(r+1). Therefore

    delta(V_r)=r+1,       delta(Z_r)=-(r+1).

Both graphs are globally sign-odd. Thus even restricting to odd graphs,
positive and negative defect weights have unbounded magnitude. This refutes
not only a unit-weight split/join involution but also the proposed all-order
partition of defects into 'negative even graphs' and 'positive odd graphs'.
The original finite observations through eight vertices remain correct.

The smallest negative odd example, 2K_(3,3), has identical components and was
proved separately above; it is not an application of the nonisomorphic-factor
hypothesis in (1).

## 5. Exactly what remains to be constructed

The correct finite signed carrier at each n consists of:

1. Every pair (X,[P]), where [P] is an Aut(X)-orbit of nonempty disjoint twin
   triple packings with sign-trivial stabilizer, with sign (-1)^(|P|+1).
2. One additional negative marker for every globally even X having a twin
   triple.

Its signed cardinality is exactly sum_X delta(X)=0. Abstract equinumerosity
therefore guarantees some pairing after listing the positive and negative
terms. The genuinely open target is a direct, canonical, graph-structural
sign-reversing involution on this decorated carrier, preserving vertex count
and avoiding an arbitrary global ranking or an appeal to the already-proved
Walsh matching. The order-twelve triple of controls must be included in any
candidate's first hostile test.

This carrier retains a packing orbit and its stabilizer-sign condition.
Collapsing it to the underlying graph forgets multiplicity, and even remembering
whether the underlying graph is globally odd does not restore the sign of its
total defect. The branch-integrality obstruction has an analogous loss of
weights on another carrier; there is no claimed map identifying the two.

The odd-module sign lemma itself extends to every odd block size. A weighted
packing formula with floor(c/d) follows for independent d-blocks when d is odd.
The particular relation to a single-cycle tournament symmetry uses the unique
three-vertex cyclic tournament and is not automatically extended to d>3.
