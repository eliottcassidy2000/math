# Folding self-converse tournaments into smaller paired tournaments

**PROVED explicit inverse constructions; FINITE-EXACT controls.** No literature-priority claim. This is a carrier-level strengthening of the counting identity in the companion packet, not a claim that a scalar count determines metagraph edges.

## 1. Statement with the correct morphisms

For m>=0, let P_m be the groupoid of ordered pairs(A,B) of tournaments on the **same** m-element vertex set, with simultaneous vertex relabelings as isomorphisms. Let P_m^col additionally carry one binary color on each vertex, preserved by isomorphisms.

Let T_n^anti be the groupoid of tournaments T equipped with a specified involutory anti-automorphism rho. Its arrows are tournament isomorphisms intertwining rho. There are explicit equivalences

`T_(2m)^anti ≃ P_m`,
`T_(2m+1)^anti ≃ P_m^col`.

Forgetting rho induces a bijection on isomorphism classes onto the self-converse tournament classes. It generally does **not** preserve automorphism groups. Thus the equivalence is for marked objects; the unmarked assertion is a class bijection, not a groupoid equivalence.

## 2. Internal winners provide a canonical transversal

An anti-involution fixes at most one vertex: two fixed vertices would force their connecting arc to reverse. Therefore rho has m two-element orbits, with one additional fixed center z at odd order.

Inside each pair choose its unique internal winner u_i and loser v_i, so u_i->v_i and rho(u_i)=v_i. This choice is intrinsic to the marked tournament. Let A be the tournament induced on the winners. For distinct i,j define

`g_ij = [u_i -> v_j]`.

The anti-automorphism relation gives

`u_i -> v_j iff u_j -> v_i`,

so g is a symmetric graph on the set of pairs. Define the second tournament B by flipping A precisely on the edges of g. Equivalently `B_ij=A_ij XOR g_ij`; symmetry of g ensures B_ji=1-B_ij. At odd order add color `c_i=[z -> u_i]`.

The construction lives on the set of rho-orbits and needs no arbitrary ordering of them. Ordering is used only to serialize a finite matrix. Any isomorphism intertwining rho sends internal winners to internal winners, so it induces a simultaneous isomorphism of(A,B), with colors when present.

## 3. Exact inverse inflation

Given(A,B), put g_ij=A_ij XOR B_ij and create vertices u_i,v_i. Direct

- each internal arc u_i->v_i;
- winner arcs according to A;
- loser arcs opposite to A;
- both u_i->v_j and u_j->v_i precisely when g_ij=1.

The last rule determines both cross arcs and their reverses. Swap u_i,v_i to obtain rho. For an odd-order input add center z, with z->u_i and v_i->z iff c_i=1, and the reverse two arcs otherwise. This is a tournament and rho is an involutory anti-automorphism.

Folding recovers A,B,c literally. Conversely inflation of a folded marked tournament is isomorphic to the original by the internal-winner/loser identification. An isomorphism of the pair uniquely lifts by u_i->u_(sigma i), v_i->v_(sigma i), and z->z. This proves full faithfulness and the stated groupoid equivalences.

## 4. Why forgetting the involution gives a class bijection

Every finite tournament automorphism group has odd order. Indeed an automorphism of order two would swap some pair of vertices and reverse its arc; Cauchy's theorem rules out an even group order. For a self-converse tournament of order at least two, the group of preserving and reversing isomorphisms has twice this odd order. Its Sylow two-subgroups therefore have order two, and their nonidentity elements are involutory anti-automorphisms.

Sylow conjugacy says all such involutions are conjugate in the signed group. The conjugator can be chosen preserving: if g is reversing and g rho g^-1=sigma, replace it by g rho, which preserves orientation and has the same conjugation effect. Thus every self-converse class admits an anti-involution, and all choices give isomorphic marked objects. Empty and singleton tournaments satisfy the statement directly.

Any isomorphism between two unmarked inflations can likewise be adjusted by an automorphism of the target to intertwine the chosen involutions. Consequently two pair objects inflate to the same unmarked tournament class iff they were already isomorphic as pair objects. This proves the unmarked class bijection.

**Automorphism sidecar.** The pair automorphism group is

`Aut(A,B[,c]) ≅ C_Aut(T)(rho)`,

not all of Aut(T). The number of anti-involutions is

`|Aut(T)| / |Aut(A,B[,c])|`.

For a directed triangle, the folded object is one colored vertex with trivial automorphism group, while Aut(C3) has order3. This is the smallest obstruction to asserting an unmarked categorical equivalence. Full Aut(T) remains recoverable by reconstructing T; it is not equal to the smaller pair stabilizer by default.

## 5. Recursive counting consequence

For an odd-part partition lambda with multiplicities m_d, let

`z_lambda=product d^m_d m_d!`,
`E(lambda)=sum m_d(d-1)/2 + sum d*binom(m_d,2)`
`           + sum_(d<e) gcd(d,e)m_d m_e`.

Define the Burnside specialization

`F_m(q,x)=sum_(lambda odd, |lambda|=m) q^E(lambda) x^(sum m_d)/z_lambda`.

The [counting packet](../counting/THEORY.md) proves

`A_m=F_m(2,1)`,
`S_(2m)=F_m(4,1)`, `S_(2m+1)=F_m(4,2)`,
`M_n=(A_n+S_n)/2`.

Here A counts ordinary unlabelled tournaments, S self-converse classes, and M converse-merged nodes. The fold explains the q=4 factor structurally: each unordered pair has two tournament bits. A coloring contributes2 per permutation cycle. A permutation with an even cycle cannot fix either tournament; the remaining odd-type orbit count is E.

The two tournaments share their vertex set and are quotiented by **simultaneous** relabeling. Their count is not A_m squared: already S4=2 whereas A2²=1; S6=12 whereas A3²=4. Relative alignment is the missing information. This is a useful recursive carrier of half the order, not a closed recurrence on bare scalar counts.

Strong-component recurrences and exact pointed growth in the companion packet provide further recursive structure. None of these identities assigns the retained literal blue/black line multiplicities from a count alone. Reconstructing the inflated graph and marked-path fibres restores those data.

## 6. Exact checks and scope

Run `python3 check_folding.py` and `python3 -O check_folding.py`. No external package or producer import is used.

The checker independently partitions all labelled tournaments through order6 into permutation orbits, finds every anti-involution in each self-converse class, verifies independence of the folded class, and checks the centralizer/automorphism relation. It also enumerates ordered-pair and colored-pair orbits through m=4, giving

`S_(0..9)=1,1,1,2,2,8,12,88,176,2752`.

Each representative inflates and folds back exactly. The C3 automorphism loss is an explicit negative control. Failure gates remain active under optimization; `--test-failing-gate` must exit1.

The all-order statements follow from the inverse constructions and group argument, not from the finite census. This packet supplies no canonical labeled choice of rho, no equality between pair and full automorphism groups, and no autonomous merged-node successor map. It makes no priority assertion beyond the displayed proof and verification.
