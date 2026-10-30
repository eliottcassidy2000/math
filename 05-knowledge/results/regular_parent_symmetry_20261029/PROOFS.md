# Regular-parent growth, random ports, and rigid forbidden-pattern examples

**Status: PROVED deductions in this packet; FINITE-EXACT native controls.
The inherited coloring recurrences are credited below. No claim of priority
in the literature, or resolution of Erdős–Hajnal, LRC, or Collatz, is made.**

Tournament automorphisms preserve the literal directed relation. A port is
a named vertex fixed pointwise unless a statement explicitly removes its mark.
A regular tournament of order 2m+1 has every outdegree m.

## R1. A regular parent makes the new vertex intrinsically fixed

Let T be regular of order 2m+1, let U be obtained by adding one vertex v,
and let G=Aut(T). Put S={u in V(T): v->u}. Then restriction gives

    Aut(U) = G_S,

with equality of permutation actions after extending G_S to fix v.

**Proof.** For any even tournament W of order 2m+2 define
R(W)={x: W-x is regular}. If x,y belong to R(W), then for every u other
than x,y,

    d_W^+(u)=m+1[u->x]=m+1[u->y].

Thus x,y have identical incidence to every outside vertex: they are clones
apart from their mutual oriented edge. If there were a third z in R(W),
then x,y would remain clones in the regular tournament W-z. Their
outdegrees there differ by one, an impossibility. Hence |R(W)|<=2.

A finite tournament automorphism group has odd order. Indeed, an
involution would exchange the ends of at least one arc and reverse that
arc; Cauchy's theorem then excludes even group order. The invariant set
R(W), having at most two elements, is therefore fixed pointwise by Aut(W).
Apply this with v in R(U). Restriction now yields exactly G_S. QED.

**Consequences and boundaries.** Every full child automorphism prime divides
|G|; no new prime can arise in this step. Regularity is essential: adding
one vertex to the transitive two-tournament can give C3. Even two successive
additions to a regular parent can introduce a prime: inside the carousel
R5, the induced sets {0,1,3}, {0,1,2,3}, {0,1,2,3,4} have automorphism
orders 3,1,5. The theorem concerns neither converse maps nor switching
equivalence. In particular the order-168 switching-class symmetry of a
source adjoined to Paley seven does not replace its order-21 literal
tournament automorphism group.

## R2. Exact unmarked extension count and its stabilizer refinement

Let e(T) count isomorphism classes of unmarked tournaments U with one
vertex-deleted card isomorphic to the regular tournament T. Then

    e(T) = |2^V(T)/G| - |V(T)/G|.                         (1)

Burnside equivalently writes this as

    e(T)=|G|^-1 sum_(g in G) [2^c(g)-fix(g)],

where c(g) counts all cycles and fix(g) counts fixed vertices.

More precisely, for every abstract group type K (or group order s), the
number of child classes with automorphism group K equals

    number of subset orbits with stabilizer K
      - number of vertex orbits with point stabilizer K. (2)

**Proof.** Rooted extensions are exactly subset orbits. R1 shows that an
unmarked child has either one or two possible regular-card roots, each
fixed by every child automorphism. Its regular cards are isomorphic:
if roots x,y are distinct, they are clones, so replacing y by x gives an
isomorphism W-x -> W-y.

Such a double-root child is obtained by cloning a vertex w of T, with
either orientation of the new mutual edge. The two rooted extensions
have new-vertex outdegrees m and m+1, so are distinct rooted orbits.
Forgetting the root merges precisely this pair. Clone children are
classified by the vertex orbit of w: an isomorphism must preserve the
two regular roots as a set, and their oriented edge distinguishes its
tail and head. The corresponding group is the point stabilizer G_w.
Subtract one for each vertex orbit, giving (1) and (2). QED.

In particular, extension families of nonisomorphic regular parents are
disjoint. Summing (2) over all regular parent classes therefore enumerates
exactly the even-order tournaments having a regular card. It does not
enumerate children with no regular card.

**An executable canonical augmentation rule.** Enumerate one subset
representative per G-orbit. Reject precisely those masks of the form
N_T^+(w) union {w}. Keep every other mask. For each clone collision this
keeps its new-root-outdegree-m representative and rejects the degree-m+1
representative; R2 proves there are no other collisions. This needs no
child isomorphism comparison or child automorphism computation. The full
permutation action on parent subsets is still required.

**Labelled identities.** If L_n counts labelled regular n-tournaments
(n odd), then the numbers of labelled (n+1)-tournaments with exactly one
and exactly two regular cards are respectively

    (n+1)(2^n-2n)L_n,       (n+1)n L_n.                  (2a)

Their sum is (n+1)(2^n-n)L_n. To prove this, count extensions with a
specified regular-card root: there are (n+1)2^n L_n. Each parent has
exactly 2n clone masks; a double-root child is counted twice. Thus there
are (n+1)n L_n double-root children, and subtraction gives (2a). The
clone masks are distinct, since equality for two different w would give
three regular roots, forbidden by R1. Equivalently, summing the refined
formula with weight 1/|Aut| gives (2^n-n)/|G| per parent.

**FINITE-EXACT audit.** All labelled even tournaments through order six,
all regular parent classes through order seven, every cut of each parent,
and all 512 cuts of C3[C3] were checked. The child groups and child
isomorphism classes were computed by native edge-consistent backtracking,
independently of (1). The even-order regular-card subclass has 1,3,7,75
classes at orders 2,4,6,8. At order eight the group-order inventory is

    59 trivial; 12 of order 3; 2 of order 7; 2 of order 21.

The 75 classes account for 2,555,520 labelled tournaments. At orders 2,4,6,
independently enumerated labelled regular-card counts agree with the
orbit-stabilizer mass of this formula.

## R3. An inherited rooted recurrence becomes an unmarked recurrence

Let T_0 be one vertex and T_{k+1}=C3[T_k,T_k,T_k], n=3^k. Its group is the
iterated imprimitive wreath product, of order 3^((n-1)/2). Let P_k(z)
count binary coloring orbits by the base-three logarithm of their
stabilizer order. The already proved recurrence from
[STABILIZERS.md](../stabilizer_necklace_kernels_20261028/STABILIZERS.md) is

    P_0(z)=2,
    P_{k+1}(z)=(P_k(z)^3-P_k(z^3))/3 + z P_k(z^3).        (3)

Nonconstant triples of child color-orbits have quotient orbit size three;
constant triples retain the top C3 and cube their internal stabilizers.
That is the mechanism of (3), not a new recurrence claimed here.

Since T_k is regular and vertex-transitive, R2 gives the following
**unmarked** child-group inventory:

    E_k(z)=P_k(z)-z^((3^k-1)/2-k).                        (4)

If F_k=P_k(1), then

    F_0=2,  F_{k+1}=(F_k^3+2F_k)/3,  e(T_k)=F_k-1.

At parent orders 1,3,9,27,81 the child counts are respectively
1,3,23,4623,32955830623. There are exactly two rigid unmarked children for
every k>=2; at k=0,1 there is exactly one. Indeed the constant coefficient
of P_k is always two, and the point stabilizer is nontrivial exactly when
k>=2. The old rigid-mask count 2|Aut(T_k)| now counts full unmarked rigid
extensions, because R1 removes the potential new-root ambiguity.

## R4. Growing palettes: a critical law for the residual group

Consider B=C3, R5 (the carousel five-tournament), or P7 (Paley seven).
Let b=|B|, h=|Aut(B)|, T_k=B^[k], n=b^k. Their groups are the full
iterated imprimitive wreath products, with order h^((n-1)/(b-1)).
This uses the prime substitution blocks of the canonical decomposition
recovered in Session68, not merely a visible subgroup.

Attach r named, pointwise fixed ports with independent fair incidence
bits. Each old vertex receives an independent uniform color from Q=2^r
possible bit vectors. Port-port edges may be arbitrary and fixed. The
group fixing all ports is exactly the color stabilizer on the old core.
For r>1 it need not equal the full group of the unmarked enlarged object.

For a permutation g of the b base vertices define its index
d(g)=b-c(g), where c(g) is its number of cycles including fixed points.
Let d=min_{g!=1}d(g). The complete base color counts are

| B | h | d | trivial stabilizer | other stabilizers, labelled colorings |
|---|---:|---:|---:|---|
| C3 | 3 | 2 | Q^3-Q | C3: Q |
| R5 | 5 | 4 | Q^5-Q | C5: Q |
| P7 | 21 | 4 | Q^7-7Q^3+6Q | C3: 7(Q^3-Q); C7 semidirect C3: Q |

For P7 there are seven order-three subgroups, each with orbit sizes
1,3,3. Nonconstant colorings fixed by distinct such subgroups are disjoint:
two of these subgroups generate the full group of order 21, whose
invariant colorings are constant. Nonidentity translations force a
constant coloring as well. These facts also follow from the complete
native base-group and coloring enumeration for Q=2,3,4.

The base polynomials, the iterated free-orbit formula, and fixed-Q
exponential rarity were already established in
[PALETTE_DYNAMICS.md](../exception_palette_dynamics_20261021/session24-palette-dynamics/PALETTE_DYNAMICS.md).
The growing-Q limit below is the extension pursued here.

**Why permutation index appears.** When Q is a prime power, the colors
may be viewed as F_Q. A permutation fixes the subspace defined by
x_u=x_(g(u)), of dimension c(g) and codimension d(g). For arbitrary Q
the same statement is a partition-equality count. Nontrivial symmetry
is a union of these fixed spaces; their intersections matter as well
as their individual dimensions. The cycle index counts fixed points
with multiplicity, whereas a regular-coloring count needs this union
or the subgroup/partition intersection data. Proposition R5 below is
an exact hostile to replacing the full action by cyclic-triple data.
This is a linear-constraint connection, not an identification with
every arithmetic squareclass matroid.

**Higher-level error bound.** At height j>=2 put m=b^(j-1) and
g_m=h^((m-1)/(b-1)). A particular child coloring-orbit has probability at
most g_m/Q^m. For any nonidentity base permutation sigma, equality of
child orbit types around its cycles requires at least d independent
equalities. If the child-orbit probabilities are p_i, then
sum_i p_i^l <= (max_i p_i)^(l-1). A union bound therefore gives

    eps(B,k,Q) =
      sum_{j=2}^k [n/b^j](h-1)(g_(b^(j-1))/Q^(b^(j-1)))^d.  (5)

Except on an event of probability at most eps, the residual group is
exactly the direct product of the stabilizers in the n/b bottom blocks.
No independence of the higher-level events is assumed.

If k->infinity, Q->infinity and n/Q^d->lambda<infinity, then eps->0:
with a=h^(d/(b-1))/Q^d the sum is bounded by a fixed multiple of
n sum_{m>=b} a^m, which is O(n Q^(-bd)). The independent bottom-block
events consequently give the following limits in distribution of
isomorphism types of finite groups:

    C3 towers:  C3^Z, Z ~ Poisson(lambda/3), rigidity -> exp(-lambda/3);
    R5 towers:  C5^Z, Z ~ Poisson(lambda/5), rigidity -> exp(-lambda/5);
    P7 towers:  C3^Z, Z ~ Poisson(lambda),   rigidity -> exp(-lambda).

For P7, full order-21 bottom stabilizers have total probability at most
(n/7)Q^-6, tending to zero in this regime. The order-three probability
is 7Q^-4-7Q^-6, giving the stated rate.

Thus r is critical around (log_2 n)/d. If r-(log_2 n)/d tends to positive
infinity the stabilizer is trivial with probability tending to one; if
it tends to negative infinity the independent bottom defects already
force that probability to zero. At bounded offsets Q is a discrete
power of two: the Poisson statements hold along subsequences with the
stated n/Q^d limit, not an asserted continuous-offset limit.

**Explicit triangle bound.** Put a=3/Q^2<1. At a node with m leaves
per child, the probability of all three child orbit types agreeing is
at most g_m^2/Q^(2m)=a^m/3. Summing the higher nodes gives

    eps <= n/[Q^6(1-3/Q^2)].                              (6)

The direct product of bottom groups is C3^Z with
Z~Binomial(n/3,Q^-2). Coupling Bernoulli(p) to Poisson(p) costs at most
2p^2 in total variation; tensoring and taking their sum gives

    TV(Law(Stab), Law(C3^Poisson(n/(3Q^2))))
       <= n/[Q^6(1-3/Q^2)] + 2n/(3Q^4).                  (7)

One elementary coupling bound follows by comparing probabilities of
0,1,>=2 for a single Poisson(p); all three errors are O(p^2), with the
displayed constant valid for 0<=p<=1.

**Scope.** This controls event probabilities and bounded observables.
It does not by itself control the expectation of the unbounded group
order, or a weighted return Gram matrix. Such a transfer would need
tail or uniform-integrability estimates. Arithmetic incidence vectors
are not automatically independent fair bits.

## R5. The 81/27 pair separates exact attachment counts

Let D be a source adjoined to C3. On three labelled blocks of size four
take the quotient C3, with all children D in A and one child D^op in A'.
These have the same 67 labelled cyclic triples, and full groups of
orders 81 and 27 respectively. Reversal changes a child's isomorphism
type and destroys the quotient C3 symmetry.

For Q named colors let t=Q(Q^3-Q)/3, the number of regular color-orbits
on D. Then the numbers of labelled colorings with trivial core
stabilizer are

    A:  27(t^3-t),       A': 27t^3.                      (8)

All three child colorings must be regular. In A their orbit types must
not all agree; in A' the mismatched child type has already destroyed
every nonidentity quotient rotation. Multiplying the corresponding
free-orbit counts by 81 and 27 gives (8).

The difference is 27t=9(Q^4-Q^2). At Q=2 the two counts are 1620 and
1728, independently exhausted over all 4096 masks on each labelled
core. These parents are not regular: (8) is a named-port/core-stabilizer
statement, not an application of R1 to an unmarked thirteen-vertex
tournament. Leading small-defect asymptotics alone cannot replace this
exact decorated action.

## R6. Explicit rigid P11-free examples with a controlled transitive size

Let B=P7 on residues 0..6, with x->y when y-x is in {1,2,4}.
Let T_k=B^[k]. The set S={0,1} has trivial setwise stabilizer in Aut(B).
For example, its stabilizer cannot swap two points in an odd-order
group, and the affine group x->a x+b (a in {1,2,4}) fixes two distinct
points only when it is the identity. Native enumeration confirms that
this is the full base group, of order 21.

Define binary coloring A_0=(0), and recursively place its complement
on the two blocks indexed by S and A_(k-1) on the other five blocks.
Its complement B_k is a second coloring. If a_k is the number of ones,

    a_0=0,  a_(k+1)=2*7^k+3a_k,
    a_k=(7^k-3^k)/2,  |B_k|_1=(7^k+3^k)/2.              (9)

Both colorings have trivial stabilizer by induction. Their different
weights distinguish their child orbits; a quotient automorphism
preserving A_k must preserve S and is therefore trivial, and every
child stabilizer is trivial. This is a constructive two-color
distinguishing family, not an appeal to a probabilistic existence
argument. The color-count transfer matrix [[5,2],[2,5]] has eigenvalues
7 and 3, retaining both total size and the signed imbalance.

Attach one vertex v that beats exactly the color-one vertices of A_k,
and call the resulting tournament U_k. T_k is regular; R1 gives

    |U_k|=7^k+1,   Aut(U_k)=1.                           (10)

Let tau be largest transitive subtournament size. For any substitution,

    tau(Q[T_i]) =
       max_{J inducing a transitive subtournament of Q} sum_{i in J} tau(T_i). (11)

Every transitive set projects to a transitive quotient set and is
transitive in each block; conversely such pieces concatenate.
This exact formula is inherited from
[tournament_open_problem_routes/PROOFS.md](../tournament_open_problem_routes_20261029/PROOFS.md).

Since tau(P7)=3, tau(T_k)=3^k. Moreover {2,3,4} is a transitive triple
outside S. Its k-fold product is an all-zero transitive set of size
3^k, all of which beats v. The general upper bound tau(U_k)<=tau(T_k)+1
is consequently attained:

    tau(U_k)=3^k+1.                                     (12)

**The forbidden-pattern guard.** A tournament H is prime here if it has
no nontrivial proper module. If H is prime and |H|>b, it cannot embed
in B^[k]: its intersection with each outer block would be a module
of H, so either it is contained in one block or it has at most one
vertex in each block. The latter is too small, and the former descends
inductively to a singleton.

The Paley tournament P11 is prime and every P11-v is prime of order
ten. This finite premise is checked by testing every nontrivial
subset of P11 and all eleven cards for the module condition, without
assuming a canonical labeling or a claimed group order. Both sizes
exceed seven. Thus P11 cannot occur in T_k; if an occurrence in U_k
used the new vertex, its deletion would embed a prime ten-tournament
in T_k, also impossible. Therefore

    U_k is rigid and P11-free for every k>=1.            (13)

This gives a precise Erdős–Hajnal obstruction test: even among rigid
P11-free tournaments, a guaranteed exponent cannot exceed log_7(3),
since log(tau(U_k))/log(|U_k|) tends to log_7(3).
It proves neither a matching lower bound nor the conjecture. Removing
all automorphisms has left the substitution extremal obstruction intact.

The general guard is reusable: if H and all H-v are prime with
|H|-1>b, every one-vertex extension of B^[k] remains H-free.
It is stronger than assuming that symmetry breaking preserves H-freeness.

**A general parity construction and a sharp guard test.** For any regular
prime base B of odd order b, choose a subset S of size s with trivial
setwise stabilizer. On words x=(x_1,...,x_k), color by

    A_k(x) = sum_i 1[x_i in S] mod 2.

The same induction gives two regular colorings with weights
(b^k-(b-2s)^k)/2 and (b^k+(b-2s)^k)/2. They are distinct in weight
because b is odd, so b-2s is nonzero. This proves the needed
distinguishing-orbit guard at every level rather than discarding it.
The signed imbalance is the product character sum (b-2s)^k.

For R5 take S={0}. The triple {1,2,3} is transitive outside S, giving
rigid one-vertex extensions of size 5^k+1 with transitive maximum 3^k+1.
Every P7-v is prime of size six, as the exact module audit verifies.
The same deletion guard makes this whole family P7-free. Thus the
rigid P7-free and P11-free examples have upper exponent barriers
log_5(3) and log_7(3), respectively.

The analogous guard fails for H=R5: R5-v has a nontrivial module.
This is a real obstruction, not just a failed sufficient test.
At depth two, color C3[C3] by 1[x_1=0] XOR 1[x_2=0], and attach v=9
to its color-one vertices. The child is rigid but its vertices
{0,2,4,6,9} induce R5, while the original ternary tower is R5-free.
Depth one has only four child vertices, so this is the first failure
within this explicit family. The exact verifier checks the five-vertex
isomorphism and the child's full group directly.
