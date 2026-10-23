# Nonfunctional does not mean arbitrary: the group heap and a prefix bank

23 October 2026. **PROVED:** the reduction and the group-chain statements
below. **CITED:** polynomial feasibility of explicitly supplied Mal'tsev
constraints. **FINITE-EXACT:** the native recursive-nine-cell control in
`nonabelian.json`. The classical constraint algorithm and stabilizer-chain
idea are recovered tools, not new algorithms claimed by this packet.

## 1. The stronger existing route

For a group G define the ternary operation

    mu(x,y,z)=x y^-1 z.

It satisfies mu(x,y,y)=mu(y,y,x)=x. If R=rH is a left coset of a subgroup
of a product of groups, then coordinatewise

    mu(rh1,rh2,rh3)=r h1 h2^-1 h3 in rH.             (1)

Thus every actual nonempty cross-block comparison relation from the main
proof is preserved by the SAME operation. This includes nonabelian groups
and nonfunctional projection kernels. Empty relations are preserved too.

For different cell groups, use their tagged disjoint union as a domain.
On three arguments of one tag use that group's mu. For mixed tags set
mu(x,y,z)=z if x=y, and x otherwise. The two identities still hold.
Every coordinate of a well-typed relation has a fixed tag, so (1) proves
preservation there. Add a unary type constraint at each variable.
Pointwise actual-owner requirements are empty sets or cosets of point
stabilizers; prescribed whole frames are singleton cosets. They are also
preserved. An arbitrary permitted subset of frames need not be preserved.

**CITED consequence.** The finite constraint system is a Mal'tsev CSP.
Dalmau, *Mal'tsev Constraints Made Simple*, Theorem 1, gives tractability
and its proof supplies a solution. For our binary constraints the algorithm
is polynomial also in the explicit domain and operation tables: retain the
displayed domain-size factors in its loops, with projection arity at most
three. See [the primary preprint](https://eccc.weizmann.ac.il/report/2004/097/download/),
sections 3–4, especially `Nonempty`, `Fix-values`, and `Next-beta`.

This complexity statement uses listed groups and listed relations (and a
supplied operation, directly computable here). Listing a cell's whole
automorphism group can be exponentially expensive in its vertex count.
No polynomial algorithm for unstructured tournament isomorphism follows.
The current Python packet does NOT implement the cited CSP algorithm.

The reusable connection is precise: actual block equality -> cosets ->
closure under a ternary operation -> compact feasibility. Replacing a
whole-cell map by separate edge choices destroys the shared variable before
the theorem can be applied. Adding arbitrary palette restrictions can
destroy closure. Those are different failure mechanisms.

For a cheap hostile, {0,1} in the additive group C3 is not closed:
mu(1,0,1)=2. Even the three-colour inequality relation fails: the allowed
pairs (0,1),(0,2),(1,2) combine to the forbidden pair (1,1). Thus the
coset theorem does not by itself solve arbitrary graph-colouring contexts.

## 2. A complete prefix bank for all compatible frames

Suppose S is the nonempty set of full source-to-target frame tuples. Pick
r in S. Every edge equation, pulled back to the whole product of cell
groups, is a coset of a subgroup. Their nonempty intersection is therefore

    S=rH,       H <= G1 x ... x Gt.                  (2)

The same is true with the allowed coset mark constraints. In actual block
comparison H is exactly the group of source automorphisms respecting the
named cells and marks, expressed in the supplied charts.

Set H^(i)={h in H: h_1=...=h_i=1}, with H^(0)=H. For each i let

    K_i=pr_i H^(i-1) <= G_i.

For every a in K_i choose one whole tuple L_(i,a) in H^(i-1) with
coordinate i equal to a; choose the identity for a=1. These are compatible
FULL tuples, not freely chosen local edge lifts.

**Theorem.** Every h in H has a unique ordered decomposition

    h=L_(1,a1) L_(2,a2) ... L_(t,at),     a_i in K_i. (3)

In particular |H|=product_i |K_i|. The bank needs at most
sum_i(|G_i|-1) nonidentity tuples, each of length t. Together with r and
the group operations, it represents every compatible frame without storing
the Cartesian product of all frames.

**Proof.** Starting with h, choose a1=h_1 and replace h by
L_(1,a1)^-1 h. The first coordinate is now identity. At stage i the
remainder is in H^(i-1); choose its coordinate a_i and remove L_(i,a_i)
on the left. Earlier coordinates stay identity. The last remainder is
identity. This proves existence and forces every choice successively, so
it also proves uniqueness. Conversely every displayed product lies in H.
Equivalently H^(i-1)/H^(i) is isomorphic to K_i. No commutativity is used.

This is a standard stabilizer-chain/transversal mechanism in the present
source-and-port coordinates. Its size bound is in the SUM of the explicit
group sizes; it is not logarithmic in arbitrarily large implicit groups.

To extract a complete bank with a feasibility oracle, normalize by r,
then for each i and a in G_i ask for a homogeneous solution whose first
i-1 coordinates are identity and whose ith coordinate is a. If feasible,
retain a full witness. These are again constraints of the same type. A
decision oracle can recover a witness by successively fixing coordinates,
using only polynomially many additional queries in the explicit domains.
This is an oracle-relative extraction proof, not a claim that the finite
checker implements the Mal'tsev algorithm.

A submitted bank alone only proves a subgroup of the solution group unless
its coordinate images K_i have been certified complete. This is the exact
coverage obligation; a few successful witnesses cannot certify it.

## 3. A native nonabelian control from the incoming recursive nine-cell

Use R9^lex=C3[C3,C3,C3], already proved in the incoming
[odd-cell completion](../odd_cell_completion_20261023/NETWORKS.md).
Its automorphisms act on (a,x) in F3^2 as

    g(a,x)=(a+s,x+t_a),  g=(s,t0,t1,t2),  |G|=81.

For composition g after h, the outer coordinate is s_g+s_h and the
inner coordinate at a is t_h(a)+t_g(a+s_h). This is nonabelian. The
cross-block I3 tensor J3 is preserved by two frames exactly when their
outer s coordinates agree. Its one-sided kernel has size 27.

For t cells joined by a connected support with these blocks, incoming
work proves |H|=3^(3t+1). In the prefix bank,

    |K_1|=81,       |K_i|=27 for i>1.                (4)

For the first lift, propagate its outer rotation through all cells and
use its inner triple only in the first cell. Each later lift has outer
rotation zero everywhere and its chosen inner triple only in that cell.
This describes the complete bank for every t. There are 80+26(t-1)
nonidentity whole tuples, compared with 3^(3t+1) compatible frames.

The finite control exhausts all 81^3=531441 candidate three-cell tuples,
checks the actual binary-block relation, obtains 59049 solutions, and
verifies unique decomposition/reconstruction for every solution. It also
checks that the displayed 81 maps are exactly the automorphisms among
all 9! permutations, and compares the group multiplication against
actual permutation composition. A twisted triangle has no solution
despite every edge separately having 2187 permitted pairs. The flat
case cannot substitute for this hostile.

This control uses the inherited family as an authenticated input to the
newly recovered general route. It does not relabel the incoming group
order or the classical transversal theorem as discoveries.

## 4. Where the new boundary is

For supplied finite cell groups and authentic comparison cosets, generic
compatibility is already covered by the cited theory; cyclicity is needed
only for this packet's much more explicit two-congruence implementation.
The useful open targets are compact native authentication, compression
relative to implicit group descriptions, and exact observations after
compression. Arbitrary integer cost optimization or descent does not
follow from feasibility. The character test in `CHARGE.md` handles the
specified abelian observations; nonabelian observation tables require
their own retained representation data.

The zero-slot analogy is then an operation rather than a numerical match:
remove an earlier coordinate while retaining the rest of the source.
The prefix bank removes a compatible group lift; Fibonacci expansion
removes a geometric block while retaining its zero-length charge. In
both cases the remainder is part of the certificate, and neither rule
licenses replacing a shared remainder with an independently chosen copy.
