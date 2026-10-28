# Proofs, observer languages, and exact boundaries

**Status: PROVED elementary / CITED inherited / FINITE-EXACT.**
The proofs received same-agent alternate-path review; no independent audit
or Lean formalization is claimed. A strong atom means a strongly connected
tournament in the ordered-join decomposition, not necessarily a modular prime.

## K1. Which of the 125 charts a named attachment can see

Let T5 have vertices Z/5 and arcs t -> t+1,t+2. Let
R=C3[T5,T5,T5], with vertices (i,t), i in Z/3. Its internal rotation
subgroup is K=(Z/5)^3, acting by (i,t) -> (i,t+a_i). The inherited
attachment diag(J5,J5,J5) is fixed by all K on one side.

For a named binary attachment matrix M with these 15 columns, call module i
active when some row is nonconstant on its five columns. Then

    Stab_K(M) = {a : a_i=0 on every active module}.

Indeed a nonzero rotation of a five-cycle is transitive, so a binary word
fixed by it must be constant. The conditions apply separately to every row
and module, proving necessity and sufficiency. With k active modules the
kernel has size 5^(3-k), and the orbit has size 5^k. A single binary row
has binom(3,k)*30^k*2^(3-k) possibilities of each kind: 8,360,5400,27000.
This classifies all 32768 named one-vertex attachments, not just random probes.
Allowing every named row separates all125. Even one nonconstant probe,
combined with cyclic module re-marking, eventually sees every coordinate.

Vertices and orientation are intrinsic tournament vertices and arcs;
attachments prescribe one arc direction per cross pair and have no ties.
The quotient loses the inactive phase coordinates. Its required sidecar
is exactly the phases later probes activate. No arithmetic tie is forced
into a tournament.

## K2. A genuine 25-state quotient, realized by strong attachments

Attach a new T5 hub with vertices u to R, setting

    (i,t) -> u  iff  u-t=a_i mod5.

All other cross arcs go from the hub to R. Write H_a for this tournament.
R and the hub are strong, and there are cross arcs in both directions, so
their union is strong. Each a gives an actual twenty-vertex strong atom.
The observation keeps the R labels and the cyclic hub type but forgets only
the hub's absolute rotation. Its cross-block response is therefore the orbit
of these fifteen singleton rows under a common cyclic column shift.

Two responses agree iff a-b is constant in its three coordinates: column
rotation adds the same amount to every a_i, and equality of the singleton
rows forces that same amount in every module. Consequently the exact
observational quotient is

    K/D,  D=<(1,1,1)>,  |K/D|=25.

It is minimal for this response: 25 different responses already occur.
Coordinates (a1-a0,a2-a0) suffice. This is a quotient of the fixed,
typed attachment response, not a claim that the H_a are 25 unlabelled
tournament classes.

There is a useful field description. Take F25=F5[omega] with
omega^2+omega+1=0; the polynomial has no root in F5. Define

    S=a0+a1+a2,
    Z=a0+omega*a1+omega^2*a2.

The map a -> (S,Z) is bijective: ker Z=D and S(t,t,t)=3t. Thus K is
F5 x F25 as an additive group, and the relative attachment retains Z.
Cyclic module rotation rho(a)=(a2,a0,a1) acts by Z -> omega*Z.
The chart transformation c(a)=-(a0,a2,a1) acts by

    (S,Z) -> (-S,-Z^5).

Here omega^5=omega^2. All translations of K, rho, and c descend to the
25-state phase quotient. To realize **global converse** on the actual H_a,
one must also retain the cross-edge polarity bit: canonical relabelling
(i,t) -> (-i,-t), u -> -u changes matching cross arcs to their complement.
The code checks this actual converse identity on all125 attachments.
The phase quotient has25 states per polarity, rather than a claim that
the original matching-only family is itself converse closed.

Among linear quotients of K invariant under rho and c, there are only
four possible kernel sizes:1,5,25,125. The rho-module splits into the
constant line and an irreducible two-dimensional module (the two coprime
factors of x^3-1 over F5); every invariant submodule is a sum of these.
Both summands are c-invariant. The exact enumeration of all64 subspaces
agrees. The possible quotient sizes are125,25,5,1. The five-state
sum-only quotient is algebraically invariant but loses the actual hub
responses, so it is not a valid further compression here.

## K3. Native re-marking destroys the discarded common phase

For b in F5 let P_b be the named vertex sequence

    0,1,...,14, then the hub cycle starting at 4+b.

The first fifteen vertices follow a Hamilton path in R; the hub segment
also follows its directed cycle. The only uncertain arc is the bridge
(2,4) -> hub(4+b), so

    P_b is legal in H_a  iff  a2=b.

For a=(j,j,j), the25-state response is constant, whereas the five guards
have the identity 5-by-5 truth table. The relative coordinates together
with these guards recover a uniquely. Thus a carrier retaining these
observations must distinguish all125, including before the subsequent
complement operation is evaluated.

Native complement tau_P reverses all arcs except the consecutive arcs
of the chosen directed Hamilton path. It is not global converse.
Already tau_(0,1,2,3,4)(T5) has outdegree multiset (1,2,2,2,3),
where T5 has (2,2,2,2,2). Therefore the cyclic-five family leaves the
restricted chart model under native complement.

The strongest survivor is the restricted relative quotient in K2, with
its polarity sidecar and declared chart changes. The first failed
implication in a proposed full reduction would be “a hub chart is
unobservable now, hence unobservable to later named guards.” P_b supplies
the missing absolute coordinate explicitly.
The all-order faithful named carrier and the aligned source/partner
representation are inherited from
[strong-atom C1](../strong_atom_continuation_vitali_20261026/README.md)
and [THM-4807 / named-language route](../legality_language_goldbach_free_perfect_companions_20261026.md);
they are not new open conjectures.

## K4. Frobenius and LCA on actual Paley7 attachments

Let P7 have vertices Z/7 and arcs u -> u+1,u+2,u+4. The21 maps

    u -> 2^a*u+b,  (b,a) in F7 x Z/3

are automorphisms. Exhausting all5040 permutations proves they are all
the automorphisms (FINITE-EXACT). Under u -> alpha^u in F8*, where
alpha^3+alpha+1=0, the same maps become

    x -> alpha^b * x^(2^a).

Pinning the owner vertex0, or field element1, forces b=0 and leaves
Frobenius3. A literal matching attachment between two P7 cells has
cross arcs i -> j iff j=g(i), with the other directions reversed.
All21 such fourteen-vertex tournaments are strong; they are concrete
attachments realizing the chart group rather than abstract field examples.

For g=(b,a), h=(c,d), composition is

    g*h=(b+2^a*c,a+d).

On a connected tree of charted cells define directed gain
g_uv=g_v*g_u^-1. Choose the root chart as identity. Each child chart is
g_uv times its parent chart, so every chart and every off-tree gain is
recovered. Common right multiplication of all charts is the forgotten
gauge. There are exactly21^(r-1) relative configurations on r cells,
so the tree is also information-minimal if every matching response may
be queried. In the pinned language replace21 by3.

For arbitrary endpoints, multiply along the path through their LCA;
reverse traversal uses inverse gains. All9261 four-cell configurations
and their55566 pair matching maps are checked. Global converse normalizes
the group by (b,a) -> (-b,a), with a cross-polarity sidecar as in K2.
Native complement along the identity path gives degree multiset
(2,3,3,3,3,3,4) and exits the Paley family.

There is a useful correction to an overly strong normality requirement.
The actual pinned orbit {alpha,alpha^2,alpha^4} is not a normal basis:
its sum is zero. Nevertheless its three trace profiles

    (Tr(alpha^(2^j)*y))_(j=0,1,2)

are011,101,110, and distinguish all three permitted y. On the whole
field the map is not injective: the polynomial masks2 and3 collide.
Thus normality is sufficient for a full-field linear reconstruction,
but not necessary for reconstruction on a specified orbit. The promise
is essential and must survive every operation. This is the actual-atom
application of the earlier F64/F256 work, not a claim that those fields
are automorphism charts of every strong atom.

## K5. Generate classes using a minimum-root certificate

Assume one representative and its actual automorphism action are known
for every order-(n-1) tournament. A pointed n-tournament is encoded by
its deleted-root parent class and the orbit of the root's outneighbour
mask under that parent's automorphism group. This encoding is complete:
an isomorphism of parents carrying the mask extends to a root-preserving
isomorphism, and restriction proves the converse.

Order pointed types lexicographically by

    (root outdegree, parent class index, canonical mask index).

Generate one extension for each parent-mask orbit, and accept it exactly
when its new root has minimum pointed type among all its vertices.
Every ordinary isomorphism class has a minimum pointed type; all roots
with that type are in one automorphism orbit. Each pointed type is
generated once. Hence exactly one extension per ordinary class is
accepted. This proves completeness and uniqueness for every n,
conditional only on the exact previous bank.

A root of higher outdegree can be skipped without classifying its
deletion. If the new root is not of minimum outdegree, reject immediately.
Otherwise stop as soon as a smaller pointed key is found. For accepted
children, if t vertices have the same pointed key as the new root,

    |Aut(child)| = t * |Stab_Aut(parent)(mask)|.

This is orbit-stabilizer, and the actual child representative is retained
for the next level. An automorphism order alone is not a sufficient
attachment carrier.

This is an application of classical canonical construction, not a
priority claim; see [McKay, Isomorph-free exhaustive generation,
1998, with errata](https://users.cecs.anu.edu.au/~bdm/papers/orderly.pdf).
Here the gain is avoiding the full root-switch matrix and exploiting a
cheap invariant before exact parent isomorphism. At order9, the prior
15449472 root queries become278461, a deterministic factor of about55.5.
The generation phase took19.877 seconds in this run; timings exclude
the subsequent full reference-bank bijection and are not controlled
cross-session speed ratios.

## K6. One exact structural DAG, three counting channels

For an odd partition lambda of m, let k be its number of parts,
z_lambda its centralizer denominator, and e(lambda) the number of
unordered-pair orbits. Put

    A_m(q,z) = sum_(lambda odd) q^e(lambda)*z^k/z_lambda.

The ordinary tournament count is A_m(2,1). Doubling every odd cycle
gives the allowed anti-automorphism cycles from THM-283:

    SC(2m)=A_m(4,1),   SC(2m+1)=A_m(4,2).

Indeed doubling sends the pair exponent to2e+k and the centralizer
denominator to2^k*z_lambda; one additional fixed vertex adds k more
pair orbits. These are classical self-converse formulas, also recorded
in [OEIS A002785](https://oeis.org/A002785).

Reuse the inherited defect DAG state (ell,S,D), where S is the unfilled
size, r is the number of previously chosen cycles, and
D_d=d*r-sum_past gcd(d,p). Its gauge bound is
b=max_(odd d>=3) ceil(D_d/(d-1)), or0 if empty.
For a history weight W, the ordinary integral coordinate is

    Y=oddpart(u!)*2^(S*(r-b))*W,  u=m-S.

For q=2^beta,z=2^gamma instead use exponent beta*S*(r-b) in this
coordinate. An edge adding k cycles of length ell has the same odd
factor T as before. If its original binary shift is J, its new shift is

    J_(beta,gamma) = beta*J + (beta-1)*v2(k!) + gamma*k.

To derive this, the combinatorial exponent and gauge changes both
scale by beta, while the denominator k! occurs only once. The original
J already subtracts v2(k!), so the displayed correction restores the
single subtraction. The z weight adds gamma*k. For beta>=1,gamma>=0,
all these shifts are nonnegative integers. Initial and terminal gauges
are1, and final division by oddpart(m!) recovers A_m.

Thus graph construction, reachable defects and odd transfer factors
are shared by (beta,gamma)=(1,0),(2,0),(2,1). No active defect states
are identified beyond the previously proved forced-skip contraction.
The three channels are scalar action counts, not a reconstruction of
individual tournament attachments or of all native metagraph edges.

## K7. A precise general rule for removing invisible coordinates

Let V be a finite abelian group. Suppose all translations are allowed,
along with a monoid M of endomorphisms, and outputs f:V -> O_f. Define
the period subgroup

    Per(f)={h : f(x+h)=f(x) for every x},
    H=intersection_(f,A in M) A^-1 Per(f).

Then two states have identical outputs after every allowed word iff
their difference lies in H. For sufficiency every affine word has linear
part A, so the output difference is A*h. For necessity, use an arbitrary
final translation to place the tested state at any x; failure of the
period condition supplies a distinguishing output. H is stable under M,
so V/H is the smallest such closed observation carrier.

For linear outputs L this becomes the familiar intersection of ker(LA).
For the K2 relative response, H=D. Adding the named K3 guards makes H=0.
For binary square-class syndrome tests this computes the unobserved
relation subspace. K4 uses group-valued transport instead of a linear
abelian carrier; it must not be folded into the same theorem without
changing its hypotheses.

K5 also identifies a different optimization: full future-state equality
can require all information even though a one-time accept/reject decision
has a short certificate. A smaller sufficient certificate is not a
smaller closed carrier. This distinction lets faithful data coexist with
early stopping and exact enumeration.
