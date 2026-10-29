# Weighted shared ports, general atoms, and the exact Paley holonomy deficit

2026-10-28. **PROVED, scoped / FINITE-EXACT; same-agent audit.** Local statements
P1-P6 belong to this file. The inherited ternary results and their charged
loop bound are credited to
[FULL_HOLONOMY](../lucas_prime_reflection_towers_20261028/FULL_HOLONOMY.md) and
[PORT_MOMENTS](../lucas_prime_reflection_towers_20261028/PORT_MOMENTS.md).

## P1. A full weighted recurrence, and a stricter hostile

For A_k, the automorphism group of the depth-k iterated triangle, define

    F_k(z)=E_(g uniform in A_k) product_(v fixed by g) z_v.

For C_k, the centralizer of digit negation t_k, define G_k analogously. With
m=3^(k-1) and z_i the variables of child i, the exact recursions are

    F_k(z)=2/3 + F_(k-1)(z_0) F_(k-1)(z_1) F_(k-1)(z_2)/3,
    G_k(z)=G_(k-1)(z_0) F_(k-1)((z_(1,u) z_(2,t(u)))_u).

The second product lies inside a single child polynomial: the paired branches
share the same automorphism. The bases are F_0=G_0=z. Conditioning on the
root rotation proves the first equation; the linked child action proves the
second. Arbitrary integer, rational or formal weights are allowed.

**Even this entire fixed-set polynomial does not classify transport.** On the
actual nine-vertex iterated triangle, let two C3 subgroups fix 0,1,2. One acts
on suffix coordinates in blocks {3,4,5} and {6,7,8} by (j,l)->(j+a,l+a); the
other by (j,l)->(j+a,l-a). Both have exactly the fixed-set distribution

    all nine vertices with probability 1/3;
    {0,1,2} with probability 2/3.

Thus every weighted fixed-set statistic, and every scalar fixed-count moment,
agrees. But the transporter (3,6)->(4,7) has size 1 for the first group and
0 for the second. Only the second subgroup centralizes digit negation.
Both are genuine subgroups of the tournament automorphism group; they are
not claimed to be two full centralizers of that same marking.

Their ordinary cycle indices also agree, namely
(p_1^9+2 p_1^3 p_3^2)/3. Thus every ordinary finite-palette orbit count agrees
as well. Incoming Session94 correctly distinguishes cycle inventories from
fixed-point moments for subset and palette queries; this example identifies
the next boundary, from complete inventory to realized labelled transport.

The missing data are the joint action and the realization equation gt=tg.
This is stronger than a failure of bounded moments: even all fixed-set
moments lose off-diagonal transport. `ports.py` checks this literally.

## P2. The joint transporter recursion

K_(A_k)(u,v) averages over the three possible root rotations. Discard a root
rotation unless it sends every source port's block to its target block.
For each surviving rotation, multiply the three child kernels on all tuple
positions belonging to each source block, then divide by three.

K_(C_k)(u,v) is zero unless every port preserves its root block. On block0
use one C_(k-1) kernel. Fold suffixes in block2 by t_(k-1), leaving block1
unchanged, and put **all** positions from blocks1 and2 into one A_(k-1)
kernel. Multiply these two kernels. Empty tuple kernels equal one.

This is a full all-order weighted transport recurrence. It follows by unique
wreath coordinates and uniform independent choices of the actual group
variables. Repeated ports are permitted. The implementation compares it with
literal permutation enumeration for every ordered source/target pair of ports
through depth2, and sampled triples at depth3.

LCA constraints are useful here because they organize the actual shared child
variable. They do not turn two marginal kernels into this diagonal kernel.
A Frobenius phase can represent a cyclic child action only after the field
coordinate map and its realization guard have been specified.

## P3. The general atom law: criticality is transitivity

Let Q be a finite tournament with specified automorphism group G. Iteratively
substitute identical copies of the previous tournament at every vertex, using
the iterated wreath action A_k. This action is always valid. Identifying it
with the *full* automorphism group requires the usual intrinsic-module
hypotheses; these hold for the indecomposable bases checked below.

Put f(s)=|G|^-1 sum_(g in G) s^fix(g). The fixed-count generating functions
satisfy F_0(s)=s and

    F_k(s)=f(F_(k-1)(s)).

Indeed only fixed root vertices can contain fixed leaves, and their child
actions are independent. This is exactly a Galton-Watson generation law with
offspring fix(g) for uniform g. Burnside gives

    E offspring = # G-orbits on V(Q).

Thus the law is critical exactly when G is transitive. If G is transitive
and rho is its number of orbits on ordered vertex pairs, then

    Var(offspring)=rho-1,
    E X_k=1,       Var(X_k)=(rho-1)k.

The variance measures the rank of the actual permutation action. It is not
specific to a three-vertex quotient and cannot classify strong atoms.

There is also a general shared-pair recursion. Suppose Q has an anti-involution
t. It fixes epsilon=0 or1 vertices, according to the parity of |Q|: it cannot
fix two distinct vertices because it reverses their arc. Write C=C_G(t).
For uniform c in C, let J(c) be the number of t-pairs whose two vertices c
fixes. An odd-order c cannot swap the two members of a pair. Set

    h_C(s)=E s^J(c).

For coordinatewise reversal on the iterated substitution, the centralizer
fixed-count generating function satisfies

    G_k(s)=G_(k-1)(s)^epsilon h_C(F_(k-1)(s^2)),   G_0(s)=s.

Proof: a fixed root vertex supplies a centralizer child; each fixed paired
root block supplies two conjugate copies of one uniform A_(k-1) action.
All unpaired choices are independent. Nonfixed root blocks contribute zero
fixed leaves. This proves the recurrence for even and odd bases with their
declared wreath actions.

For odd transitive Q, put mu=EJ, v=Var(J), sigma2=rho-1. Then

    E Y_k=1+2 mu k,
    Var(Y_k)=2 mu sigma2 k(k-1)+4 v k.

These formulas follow from conditional variance of the increment
2 sum_(i=1)^J X_(k-1,i). They separate shared-branch correlation from randomness
in which root pairs survive. Differentiating the generating functions gives
all ordered-port counts, since r-port Burnside counts equal E X^r or E Y^r.
These count recurrences alone do not supply the transport kernels of P2.

## P4. Paley-7 and Paley-11 have the same critical variance, different ports

Let p be a prime congruent to3 modulo4, h=(p-1)/2, and orient x->y when
y-x is a nonzero square modulo p. Use the affine subgroup

    G_p={x->a x+b : a a nonzero square, b arbitrary}.

It has order ph and preserves the tournament. Reflection x->-x reverses it.
For p=3,7,11, `paley.py` checks by unrestricted isomorphism search that this
is the full automorphism group and by exhaustive subset checks that the base
is indecomposable. The intrinsic modular-decomposition mechanism then gives
the full wreath actions at every depth for those three bases. For other p,
the theorem is stated for the specified affine wreath action.

More explicitly, for an indecomposable strong quotient the substituted copies
are its maximal proper strong modules. Every automorphism permutes those
intrinsic modules through Aut(Q), and its restrictions are independent child
automorphisms. Conversely all such maps preserve substitution. This gives the
wreath equality inductively, using the inherited modular-decomposition theorem,
rather than extrapolating a large-group census from the small bases.

Its root fixed-count law is

    P(N=0)=2/p,
    P(N=1)=(p-3)/(p-1),
    P(N=p)=2/[p(p-1)].

These are respectively nonidentity translations, nonunit multipliers, and
identity. Consequently E N=1 and Var(N)=2 for every such p. The first three
ordered-port orbit counts of A_(p,k) are

    a_(k,1)=1,
    a_(k,2)=2k+1,
    a_(k,3)=6k^2+2(p-2)k+1.

The third moment detects the quotient that the variance misses. In the root
centralizer, c(x)=a x, so J=h with probability1/h and J=0 otherwise. Hence

    G_k(s)=G_(k-1)(s)[(h-1)/h + F_(k-1)(s^2)^h/h],
    c_(k,1)=2k+1,
    c_(k,2)=8k^2+4k(h-1)+1.

At p=3 this recovers the inherited 8k^2+1. At p=7 it is 8k^2+8k+1;
at p=11 it is 8k^2+16k+1. The receipt computes exact moments through order6
and depth8, without enumerating the giant groups.

## P5. Exact holonomy for these larger atoms

Keep the inherited definitions theta(g)=tgt, K={g:theta(g)=g^-1}, C=Fix(theta),
L=<K>, and marking holonomy Hol=C intersect L. The equality is proved in the
inherited FULL_HOLONOMY, not assumed from a single nonflat triangle.

First, L is normal in G. The elements delta(g)=g^-1 theta(g) lie in K, and
every element of K is of this form by taking its inverse square root.
The identity delta(xy)=y^-1 delta(x)y delta(y) proves normality. Thus G/L is
the largest quotient on which theta acts trivially.

On the odd abelianization G_ab, the plus and minus eigensubgroups of theta
intersect trivially. The image of C is plus, and that of L is minus. Therefore

    Hol is contained in [G,G].

This necessary obstruction applies to any odd group with an involution.

For p>3 in P4, G_p has commutator subgroup equal to its translations: a
nonidentity square multiplier has commutator (a-1)b with translation b, which
ranges over the field. Thus (G_p)_ab=H, the multiplier group of order h.
Transitivity of the root action gives inductively

    (A_(p,k))_ab = H^k.

Explicit coordinates chi_j are the products of the local affine multipliers
at all tree nodes of depth j, 0<=j<k. Composition only permutes the factors,
so each chi_j is a homomorphism. Digit negation leaves all chi_j unchanged.
It follows that L is contained in ker(chi)=[A_(p,k),A_(p,k)].

For the reverse inclusion, a nontrivial root translation belongs to K, and
L is normal. In A_(p,k)/L this translation is identity. Conjugation by it
cycles the p child factors, so all their images agree. Distinct child factors
commute, forcing their common image to be abelian. The remaining root
multiplier group is abelian and merely permutes those equal child images.
Hence the quotient is abelian, so the entire commutator subgroup lies in L.
We have proved

    L=[A_(p,k),A_(p,k)],
    Hol=ker(chi restricted to C_(p,k)).

Finally chi maps C_(p,k) onto H^k. A local scalar at the all-zero-prefix node
of depth j commutes with digit negation and changes precisely chi_j. Therefore

    [C_(p,k):Hol] = h^k,                p>3.

The lower bound is thus an **equality**. With

    E_(p,k)=((p^k-1)/(p-1)-k)/2,
    |C_(p,k)|=h^k (ph)^E_(p,k),
    |Hol|=(ph)^E_(p,k).

At depth2: Paley7 has centralizer order83349 and holonomy order9261; Paley11
has centralizer order12582109375 and holonomy order503284375. At depth1 the
connection is flat although the centralizer has order3 or5. The literal base
triangle computation is included in `paley.py`; the all-depth conclusions
come from the proof, not exhaustive large-group enumeration.

For p=3, h=1 and the inherited theorem instead gives L=A_k and Hol=C_k.
Its constructive bound 3[1+(k-1)3^k] remains valid in that family. It cannot
be exported as a repair bound for arbitrary centralizer elements of the
larger Paley families: h^k classes of corrections are unreachable by loops.

**Incoming improvement credited.** Session94's
[portrait protocol](../session94_recursive_ports_crt_clocks_20261029/loop_protocol/LINEAR_LOOP_PROTOCOL.md)
and the concurrent
[linear repair](../crt_weighted_ports_20261028/LINEAR_REPAIR.md) lower the
ternary bound to 3(3^k-2k-1)/4 legal re-marking edges: at27 vertices,15 rather
than165. They construct endpoints of conjugated individual triangles while
retaining full port permutations. This improves the reachable C3 case; it
does not remove the Paley multiplier obstruction proved above.

## P6. A mirrored spine, with the probability law kept explicit

The classical size-biased critical Galton-Watson construction chooses a spine
and gives its other children independent ordinary descendant trees; see
[Lyons, Pemantle and Peres, Section2](https://arxiv.org/html/math/0404083v1).
For offspring0/3, the size-biased spine has two independent side branches at
each level. Its generation count has generating function

    s product_(j=0)^(k-1) F_j(s)^2 = s F_k'(s).

The ternary reversal centralizer instead has **identical paired side
branches**, with generating function

    s product_(j=0)^(k-1) F_j(s^2).

Both means are 2k+1. Their variances are respectively 2k(k-1) and 4k(k-1),
and their second moments are 6k^2+2k+1 and 8k^2+1. This is an elementary
consequence of independence versus duplicate use of the same X_j, whose
variance is2j. The extra variance has a concrete source-sharing mechanism.
It does not give an arithmetic orbit distribution for Collatz or a density
law for prime-valued Lucas numbers.

**Next targets.** Classify the fixed-pair polynomial h_C alongside the full
root transporter action for other indecomposable self-converse atoms. Determine
which theta-trivial quotients can be nonabelian, where the abelian obstruction
would be incomplete. Minimize weighted kernels only after fixing which native
complement contexts can query their off-diagonal entries.
