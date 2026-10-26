# Free gluing is a kernel condition on named ports

2026-10-25. **CITED / PROVED specialization / FINITE-EXACT.** The general
kernel description is classical: Hidalgo, *Hyperelliptic quotients of
generalized Humbert curves*, arXiv:1705.09337v5, Theorem 1. The count and
five-port application below are derived explicitly. No novelty claim for
the classical description is made.

## A complete compatibility carrier

Let b>=5 distinct geometric branch points have inertia vectors v_1,...,v_b
in G=F_2^(b-1), with their only nonempty binary relation sum v_i=0.
These are the marked deck data of a generalized Humbert cover. A subgroup
H acts freely exactly when it contains none of the v_i. Equivalently, its
quotient map rho:G->F_2^c is encoded by named nonzero columns

    a_i=rho(v_i),   sum_i a_i=0,   span{a_i}=F_2^c.       (P1)

This is Hidalgo's kernel criterion in additive notation. It retains more
than pairwise compatibility: composing proposed identifications means
forming their generated subgroup and checking every marked inertia.

For an arbitrary proposed H with basis matrix, freeness is therefore a
finite list of b linear nonmembership tests. There is no need to enumerate
all 2^dim(H) group elements. The quotient's actual port columns a_i must
remain attached to the same branch points. Their multiset alone loses
which branch conditions and arithmetic twists belong to each port.

## Minimum quotient rank and an exact odd/even split

If b is even, c=1 is possible: assign a_i=1 to every port. It is the unique
rank-one free quotient kernel, of dimension b-2. A smaller quotient is
impossible because the trivial quotient kills all inertia.

If b is odd, rank one is impossible because sum a_i=1. Rank two is
possible: assign its three nonzero vectors multiplicities (b-2,1,1).
More generally, all three multiplicities must be odd. Indeed their sum
is zero precisely when their three parities agree; odd total forces
each parity to be one. In particular every valid assignment spans rank two.
Thus the maximum dimension of a free H is b-3.

The number of **distinct largest-dimensional free subgroups** for odd b is

    (3^(b-1)-1)/8.                                  (P2)

Proof: use the four additive characters of F_2^2 to count length-b words
in its three nonzero letters with XOR zero. The count is
(3^b+3*(-1)^b)/4=(3^b-3)/4. All such words have three odd multiplicities
and span; GL(2,2), of order six, acts freely by permuting their three
letters. Two surjections have the same kernel iff related by this action.
Dividing by six gives (P2). This counts largest dimension, not all subgroups
maximal under inclusion. Confusing those two predicates would overstate it.

Riemann--Hurwitz gives g(C)=1+2^(b-3)(b-4). A free H of dimension d has
g(C/H)=1+2^(b-3-d)(b-4). At maximum dimension this is (b-2)/2 for even b
and b-3 for odd b. These are geometric genera; no rational-point existence
is inferred.

For b=5, the multiplicities are (3,1,1), there are ten kernels, and their
quotients have genus two. For b=7 the two possible unordered multiplicity
patterns are (5,1,1) and (3,3,1); there are 91 kernels and genus-four
quotients. The recurring three-colour pattern now has an exact explanation:
the colours are the three nonzero elements of a two-dimensional binary
quotient, and the odd counts express the single port relation.

## Why the Petersen graph is exceptional at five ports

At b=5, the ten pair-sums exhaust the nonzero non-inertia elements. The sum
of two pair-sums on disjoint port pairs equals the remaining inertia.
At b=6, it instead equals the pair-sum on the complementary two ports.
At b=7, it is neither an inertia nor a pair-sum. Thus the same pairwise
disjointness graph no longer means the same geometric obstruction.

This is a cheap hostile to extending the Petersen rule by its shape alone.
The general rule is (P1), or directly H intersect {v_i}=empty. The special
five-port circuit makes disjointness detect that predicate at arity two;
the three-star hostile in FIXED_TWISTS.md already shows its arity-three loss.

## Arithmetic and algorithmic connections, with their boundaries

* **Non-UFD conductors:** a port column records ramification parity, but
  not the coefficient ideal or its codimension-two defect. Keep the ideals
  from CONDUCTOR_IDEALS.md as well as this binary carrier. THM-2685 supplies
  the separate unit/class-group sequence; neither parity nor a class-group
  coordinate is a replacement for actual conductor membership.
* **LCA:** a trie can compress a hierarchy of named ports, but generated
  subgroups form a lattice with crossing subspaces. Use the Session56
  partition-lattice hostile before replacing it by a single tree. A tree
  of a chosen computation is permissible; it is not the full compatibility
  relation.
* **Tournament repairs:** the incoming Session64 certificates preserve
  actual vertex identities even through a clone merger. Our kernel compiler
  similarly preserves actual branch identities, but no map from its genus
  or freeness predicate to black-move legality is claimed.

Source: [Hidalgo, Theorem 1 and Corollary 1](https://arxiv.org/html/1705.09337v5#S3.SS1).
Our independent finite control checks all nonzero column assignments at
b=5,7,9 and their GL(2,2) orbits, plus the b=5/6/7 disjoint-pair hostile.
