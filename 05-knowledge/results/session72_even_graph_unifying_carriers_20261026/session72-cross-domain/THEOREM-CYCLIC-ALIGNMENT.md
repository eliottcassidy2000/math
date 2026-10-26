# Synchronized cyclic alignment: exact all-arity descent and its sparse tests

**Status: PROVED elementary synthesis, with FINITE-EXACT controls.** The pair case, normal-basis dictionary and carousel source family are recovered from Session68. The contribution here is an explicit all-arity compatibility theorem, its sharp prime-layer criterion for selecting fewer pair tests, and an all-arity native tournament realization. The arithmetic ingredients are classical; no literature-priority claim is made. Independent audit is recorded separately by the parent session.

## 1. The complete phase theorem

Let `r>=1`, `d_1,...,d_r>=1`, `L=lcm(d_1,...,d_r)`, and let `m` be a positive multiple of `L`. On

    X = product_i Z/d_i Z

let `C_m` act by simultaneous translation: `t.x=(x_i+t)_i`. This is the exact operation, not independent translations in the factors. Put `g_ij=gcd(d_i,d_j)` and retain the **named relative phases**

    r_ij(x)=x_j-x_i mod g_ij, i<j.                         (1)

Extend by `r_ji=-r_ij` and `r_ii=0`.

**Theorem CA1.**

1. Every orbit has size `L` and stabilizer order `m/L`; there are

       N=(product_i d_i)/L                               (2)

   orbits.
2. Two states have the same complete phase record `(1)` if and only if they lie in one simultaneous-translation orbit.
3. A proposed phase record is realizable if and only if, for every triple of distinct indices,

       r_ij+r_jk=r_ik mod gcd(d_i,d_j,d_k).                (3)

   Thus the complete orbit quotient is exactly a finite set of compatible, modulus-labelled edge cocycles. Neither independent local orbit types nor arbitrary unconstrained edge labels give that quotient.

**Proof.** A shift fixes a tuple exactly when every `d_i` divides it, equivalently `L` divides it; this proves (2). If two records agree, their coordinate differences `delta_i` satisfy `delta_i=delta_j mod gcd(d_i,d_j)`. The generalized Chinese remainder theorem supplies a common `t` with `t=delta_i mod d_i` for every i. Conversely a common shift preserves all differences.

For completeness, the generalized CRT assertion follows prime by prime: for each prime p choose an index a whose `v_p(d_a)` is maximal; the desired t modulo the maximal p-power is `delta_a`, and all smaller-p-power equations follow from pairwise compatibility. Combine these coprime prime-power congruences by ordinary CRT.

Necessity of (3) is subtraction. To prove sufficiency, fix p and write `v_i=v_p(d_i)`. Choose an anchor a with maximal v_i. Set

    x_i = r_ai mod p^v_i.

For any i,j, (3) applied to a,i,j gives `x_j-x_i=r_ij mod p^min(v_i,v_j)`. Cases with an index equal to a are immediate from antisymmetry. For each i, combine its residues over distinct primes by CRT. This constructs x satisfying every equation (1), including the empty-prime case d_i=1. QED.

**A real correction of the tempting negative claim.** In this cyclic model, all labelled pairwise relative phases **are sufficient** for every arity. It would be false to claim that an inherently three-way deterministic coordinate is missing after these complete pair phases have been retained. Independent necklace types omit (1); pairwise probability marginals are a different observer again.

## 2. The economical version: exactly which pair tests are necessary?

Let E be a simple graph on the r named factors. Keep only the phase records r_ij on edges of E. For each prime p and exponent a>=1 define

    I_(p,a)={i:p^a divides d_i}.

Ignore empty sets; a singleton is connected.

**Theorem CA2 (prime-layer connectivity criterion).** The edge-restricted record separates simultaneous-translation orbits if and only if every induced graph `E[I_(p,a)]` is connected.

**Proof, sufficiency.** If two states have equal E-records, their differences delta satisfy `delta_i=delta_j mod p^a` on every edge in `E[I_(p,a)]`. Connectivity propagates equality through that layer. For each i,j and every p-power dividing gcd(d_i,d_j), they therefore agree modulo that p-power. CA1 now gives a common shift.

**Proof, necessity.** Suppose `E[I_(p,a)]` has at least two components, and choose one component C. Write `L=p^v M`, with gcd(p,M)=1. Define one state to be zero and the other by

    delta_i=p^(a-1) M mod d_i if i belongs to C;
    delta_i=0 otherwise.                                (4)

For any retained edge crossing C, the other endpoint is outside `I_(p,a)`, since no edge connects distinct induced components. Its p-valuation is at most a-1. Thus (4) vanishes modulo the p-part of that edge's gcd; M kills every other prime part. Every retained edge sees zero. But a vertex in C and a vertex in another layer component differ nontrivially modulo p^a, so no common shift exists. This proves necessity. QED.

This criterion concerns completeness on **actual phase states**. It does not say arbitrary sparse edge labels automatically admit a phase state; CA1's full-record compatibility or a separate lifting test must still be enforced.

### Sharp hostiles and positive boundary

- A star through a short-period factor can fail: periods `(2,4,4)`, states `(0,0,0)` and `(0,0,2)`, retain equal two anchor phases but differ on the omitted mod-4 edge. All three factors being connected in the ordinary graph is insufficient.
- The odd-period version is `(3,9,9)` with states `(0,0,0)` and `(0,0,3)`.
- Periods `(6,10,15)` require all three pair edges: the prime-2, prime-3 and prime-5 layers are the three different pairs. No spanning tree suffices.
- Within the **odd carousel scope**, periods `(15,21,35)` inside m=105 require all three edges, now for primes 3,5,7. There are `15*21*35/105=105` joint alignment classes.
- When all periods are equal to d, CA2 reduces to ordinary connectedness: a spanning tree's full mod-d differences suffice. Cyclic transport has an economical frame exactly when the arithmetic layers permit it.

## 3. The field interpretation is an actual intertwiner

For a prime power q and m>=1, choose a normal basis

    beta_i=alpha^(q^i), 0<=i<m, of F_(q^m)/F_q.

The coordinate map `(epsilon_i)->sum epsilon_i beta_i` conjugates cyclic coordinate rotation to q-Frobenius. The normal-basis theorem is the only imported existence statement: [Keith Conrad, Linear Independence of Characters, Theorem 3.7](https://kconrad.math.uconn.edu/blurbs/galoistheory/linearchar.pdf). The choice of normal basis is retained; an arbitrary polynomial basis need not intertwine these actions.

Fix r Frobenius orbits of respective exact degrees d_i dividing m, and choose one reference element in each. The choices of conjugates are exactly X above. A **joint** Frobenius orbit is therefore classified by CA1, with sparse sufficient observations classified by CA2. Replacing each element separately by its minimal polynomial loses these N possible relative placements.

The compatible algebra statement is

    tensor_(i=1)^r F_(q^d_i) over F_q
        is isomorphic to product_(j=1)^N F_(q^L).        (5)

**Proof of (5).** For two finite fields of degrees a,b, write the first as `F_q[t]/f(t)` for a degree-a irreducible. Its a roots lie in one q-Frobenius orbit. Over `F_(q^b)`, q^b-Frobenius partitions those roots into gcd(a,b) orbits of length a/gcd(a,b). Separability makes the corresponding factors distinct. Polynomial CRT gives gcd(a,b) fields, each of total degree lcm(a,b) over F_q. Iterating this two-field result proves (5), since a current factor of degree L_j yields gcd(L_j,d_(j+1)) new factors; the total count telescopes to product(d_i)/L. Component embeddings are indexed by the same simultaneous-Frobenius phase orbits. The isomorphism is not a canonical numbering of components without references. QED.

This is neither a claim that tournament composition is field multiplication nor a transfer of an arithmetic survival theorem. It is an exact correspondence of a specified finite action, its joint orbit quotient, and its stabilizers.

## 4. Native all-arity realization in tournaments

This extends Session68's two-component construction without claiming a classifier for arbitrary tournaments.

Let m>=3 be odd, write m=2k+1, and let R_m be the carousel tournament on Z/m with `i->j` iff `1<=j-i mod m<=k`. Let `A=1->C3` on four labelled vertices and `A^1=A^op`, `A^0=A`. For an r-bit colour e define

    B_e=A^(e_1) -> A^(e_2) -> ... -> A^(e_r).             (6)

Here arrows denote ordered joins, not tensor products. For each of r binary words epsilon_j of length m, put e_i=(epsilon_1(i),...,epsilon_r(i)), and form

    T=R_m[B_(e_0),...,B_(e_(m-1))].                      (7)

**Theorem CA3.** The 2^r colours B_e are pairwise nonisomorphic. T has `4rm` vertices, and two such T are isomorphic exactly when their r words differ by one common rotation. If their least periods are d_1,...,d_r, then

    Aut(T) = C3^(rm) semidirect C_(m/L),
    |Aut(T)|=3^(rm) m/L.                                 (8)

Fixing the r individual necklace classes therefore leaves exactly N joint tournament types, distinguished by CA1's relative phases; CA2 gives exactly which pair phases suffice to recover the joint type.

**Proof.** The strong-component sizes of A are `(1,3)`, those of A^op are `(3,1)`. Strong components of an ordered join occur in their displayed linear order. Thus the 2r-term strong-component-size sequence of B_e decodes e, by grouping consecutive positions in pairs. Its automorphism group is C3^r, with no permutation of strong components.

For self-containment, R_m is prime: if a proper module M has at least two vertices, an outside vertex confines M to one of its in/out cyclic intervals of length k. Lift that interval consecutively to integers, take a<b in M, and let z=a+k+1 modulo m. This z is outside the interval, and z->a but b->z, contradicting that M is a module. Also an arc i->i+d has k-d common out-neighbours; maximizers d=1 are the directed successor cycle. Hence every automorphism is a rotation, and every rotation is allowed.

In (7), a proper module meeting more than one outer block has a met-index set that is a module of R_m, and hence would meet all blocks. If it is proper, an omitted vertex in one met block would have to see vertices of every other block uniformly, forcing its carousel vertex to dominate or lose to all other carousel vertices. Regularity forbids that. Thus every proper module is contained in one block. The displayed blocks are consequently canonical maximal proper modules. Any isomorphism induces one carousel rotation and must preserve every colour; conversely every colour-preserving rotation lifts. The rotation stabilizer has order m/L; independent local C3 rotations give its kernel and the fixed local labels provide a splitting. QED.

**Scopes.** Even m is excluded from the native construction: an involution cannot be a tournament automorphism because it reverses the arc of a swapped pair. CA1, CA2 and (5) remain valid for even periods as arithmetic statements. Vertex number, labelled colour positions and permitted common rotations are retained. No arbitrary wreath-action classification follows.

## 5. Why this does not erase the other obstructions

1. **Multidimensional group boundary.** Let G=F3^2 act on X=F3^3 by `(u,v).(x,y,z)=(x+u,y+v,z+u+v)`. Every individual orbit and every pair orbit is full, but `z-x-y` distinguishes three joint orbits. This is an abstract finite-action hostile, **not** asserted to be the full automorphism action of a tournament. The cyclic hypothesis matters, not just abelianness. As a separate nonabelian control, A4 on four points is transitive on ordered distinct pairs, while `(0,1,2)` and `(0,2,1)` have distinct triple orbits. A4 is not a tournament automorphism group.
2. **Law boundary.** Deterministic recovery from coherent pair differences does not recover an arbitrary distribution from pair marginals. Uniform parity laws on a Boolean k-cube have equal every-proper-coordinate marginals and different full conjunction probability. THM-2237 proves this; Sessions70/71 realize the conjunction in actual tournament/Penrose guard protocols. These are different state and observer categories.
3. **Realization boundary.** Flat seam transport produces globally consistent relative placements but not disjoint interiors. The five-square E,N,W,S path from Session71 is flat and overlaps its first and last pieces. Even a complete phase/frame classification does not discharge a geometric acceptance predicate.
4. **Counting boundary.** Equal numbers of orbits do not construct an equivariant bijection or identify stabilizers. In (7) there is an actual source map and canonical block recovery; in the normal-basis dictionary there is an actual intertwiner. These are what make proof transport valid, rather than coincident cardinalities.

## 6. Finite-exact checks

`check.py` independently constructs phase orbits, enumerates all possible complete-edge labels in its declared small universes, reconstructs compatible phases prime by prime, and tests every edge subgraph against the layer criterion. It also constructs tournament adjacency matrices and recovers outer blocks by minimal-module closure; selected automorphism groups are enumerated directly from adjacency rather than inferred from (8).

The tensor control is a literal multiplication-preserving F2-linear isomorphism `F8 tensor F8 tensor F8 -> F8^9`: a 27-by-27 matrix has rank27, and every pair of 27 tensor basis elements has matching products. All arithmetic is exact, using no repository implementation and no optimization-sensitive assertions. The JSON receipts give the complete finite counts; these tests corroborate, not replace, the all-parameter proofs above.
