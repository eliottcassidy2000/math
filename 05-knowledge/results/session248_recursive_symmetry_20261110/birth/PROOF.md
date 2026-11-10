# Session248: symmetry exchange, compatible rebirth, and reversible growth

**Status: PROVED elementary propositions below; inherited barrier explicitly separated; OPEN directions at end.** No literature-priority claim. This is a mechanism-building lane, not a re-audit of the inherited results.

## Inheritance and concept board

Closest inherited mechanism: `triangle_return_170807_20261103/TRIANGLE.md`, §6 (TR5), charges the first birth that breaks an actual old orbit. Hostile: equal Hamilton counts or isomorphic vertex-deleted cards do not specify a compatible rooted extension. Corrected near miss: old automorphisms preserving the incidence are only the newborn stabilizer, not all child automorphisms. Least-used sidecars recovered here:

* THM-2716, `01-canon/theorems/THM-2716-c4-arm-transporter-groupoid-and-relative-degree-holotopy-boundary.md`: transporter fibres are bitorsors, not canonically chosen arrows.
* THM-3345, `01-canon/theorems/THM-3345-prime-xor-ancestry-path-groupoid-and-source-dependent-berggren-cost.md`: source-dependent transport labels can be essential even when scalar or ambient holonomy is flat.
* THM-4987, `01-canon/theorems/THM-4987-subgroup-path-orbit-gaps-and-small-symmetry-cost.md`: subgroup index is a real multiplicative cost in path-orbit quotients; carousel Q7 is an important hostile outside the low quotient window.

Live concepts: incidence fracture; compatible rebirth; groupoid mass; reversible growth; transport holonomy; symmetry rupture cost. Anchor is the exact birth law, niche is groupoid detailed balance, wildcard is distinguishing scalar flatness from transporter holonomy.

## 1. Exact exchange law

Let T be an n-vertex tournament. Add x with S=N_U^+(x) contained in V(T), obtaining U. Write

    G=Aut(T), K=Stab_G(S), A=Aut(U),
    r=|G.S|=[G:K], o=|A.x|=[A:A_x].

**Proposition B1 (PROVED).** Restriction identifies A_x with K, and

    |Aut(U)| / |Aut(T)| = o/r.                         (B1)

Thus old incidence fracture costs r, compatible replacement of the newborn supplies o, and their quotient is the net symmetry change. Both integers are odd. In particular growth occurs iff o>r, equality iff o=r, and decay iff o<r. If a birth genuinely breaks old symmetry (r>1), its growth factor is at most the largest odd integer at most n+1, divided by 3.

Proof: an automorphism fixing x restricts to an automorphism preserving its out-neighborhood; every such automorphism extends uniquely by fixing x. Apply orbit-stabilizer twice. Tournament automorphism groups have odd order: an involution would swap some vertex pair, reversing its directed edge. Cauchy's theorem then excludes even order. Hence both indices are odd. This also proves the bound.

Interpretation: a large child automorphism group need not be inherited. If r=1, old symmetries survive but may become a proper point stabilizer. If o=1, every child symmetry is inherited and symmetry can only decrease. The integers o,r are NOT determined by the two group orders alone.

## 2. What compatible replacement actually means

For y in V(U), let C(x,y) consist of isomorphisms

    f: U-x -> U-y

such that for every v != x,

    x -> v iff y -> f(v).                             (B2)

Equivalently, f transports the entire deleted-vertex incidence word, including the coordinate for y when y != x.

**Proposition B2 (PROVED).** C(x,y) is nonempty exactly when y belongs to A.x. Its arrows extend uniquely by x->y to automorphisms of U. Every nonempty C(x,y) is an A_y-A_x bitorsor, so |C(x,y)|=|K|. Consequently

    |A| = sum_y |C(x,y)| = o |K|.                     (B3)

Proof: (B2) checks precisely the only arcs not checked by the deleted-card isomorphism. Extension therefore preserves every arc. The inverse implication is restriction. Postcomposition and precomposition by stabilizers act freely and transitively on a nonempty transporter fibre.

This gives a rigorous replacement for the invalid shortcut “isomorphic deletion cards imply movable newborn.” Bare deletion isomorphisms are candidates; incidence compatibility is the acceptance condition. No reconstruction conjecture is being assumed.

### Holonomy and loss ledger

Objects are vertices of U; arrows x->y are these compatible extended permutations. This is exactly the action groupoid A acting on V(U). Composition is composition of FULL extensions followed, if desired, by deletion restriction, not naive composition of differently punctured maps.

A loop based at x has holonomy in A_x=K. Picking one transporter s_y:x->y for each reachable y gives flat selected arrows s_z s_y^{-1}, but does not make all arrows flat or produce a canonical choice. For example take U=C3[C3], three directed-triangle blocks in cyclic orientation. A rotation in a block not containing x is a nontrivial h fixing x. For any a moving x to y, the two compatible arrows a and h a^{-1} form a loop with product h. This is holonomy of chosen transport, not an intrinsic obstruction to the existence of a flat section.

The coarse deletion-support graph forgets compatibility, multiplicities, and stabilizer action. The orbit partition forgets the bitorsors. Restoring complete extension permutations is the sidecar. This is the precise analogue of THM-2716's warning about identifying distinct arms, without importing its C4 parity mechanism into odd-order tournament groups.

## 3. Sharp, arbitrarily large gains and losses

Let R_q, q=2m+1>=3, be the carousel on Z/q with i->j iff j-i belongs to {1,...,m}.

**Proposition B3 (PROVED).** R_q-0 is rigid and Aut(R_q) is the rotation group C_q. Hence the birth reconstructing R_q from R_q-0 has r=1, o=q, and symmetry gain q.

Proof: in R_q-0, vertices 1,...,m have outdegree m and vertices m+1,...,2m have outdegree m-1. Each degree class induces a transitive tournament, which has only the identity automorphism. Thus the deletion is rigid. Rotations act transitively on R_q, and a point stabilizer embeds into the rigid deletion, so the whole group has order q.

This attains the absolute one-birth bound for every odd child size. For even child size N>=4, take a source singleton followed by R_{N-1}, and insert the newborn in the carousel. The parent is rigid and the child has group C_{N-1}, attaining the maximal possible odd orbit N-1. Thus the universal bound is sharp at every child size >=3.

**Proposition B4 (PROVED).** For q=2m+1>=5, a birth at R_q whose newborn beats exactly one old vertex produces a rigid child. Its symmetry ratio is 1/q.

Proof: the newborn has outdegree 1; old vertices have outdegree m or m+1, so it is uniquely recognized and fixed. The old rotation group stabilizing a singleton is trivial. Thus o=1, r=q, K=1.

The same cyclic symmetry can therefore be created from a rigid puncture and destroyed by one deliberately asymmetric incidence. Neither group order nor vertex count is monotone under induced growth.

## 4. Groupoid mass gives exact birth/death detailed balance

For unlabelled isomorphism classes T at size n and U at size n+1, fix a labelled representative of T and define

    a(T,U) = number of subsets S of V(T) whose birth has type U,
    d(U,T) = number of vertices y of U with U-y isomorphic to T.

These are integer multiplicities, not Boolean support indicators.

**Proposition B5 (PROVED).**

    a(T,U)/|Aut(T)| = d(U,T)/|Aut(U)|.                (B4)

Proof: split the vertices y counted by d into Aut(U)-orbits O_j. Each orbit corresponds to exactly one Aut(T)-orbit of birth incidence subsets: choose an isomorphism T->U-y, and transport y's incidence; different choices differ by Aut(T). Conversely incidence subsets in one old orbit define the same pointed child. If K_j is the pointed-child stabilizer, its incidence orbit has size |Aut(T)|/|K_j|, and |O_j|=|Aut(U)|/|K_j|. Sum these equalities. Notice that all deletion cards of type T contribute, potentially in several incompatible rooted orbits; B4 does not collapse them.

This is the groupoid cardinality of the span of pointed births, expressed from its two ends. It recovers the exact old/new multiplicity factor omitted by a simple unweighted metagraph.

Define kernels

    B_n(T,U)=a(T,U)/2^n,
    D_{n+1}(U,T)=d(U,T)/(n+1),
    pi_n(T)=n! / (2^(n(n-1)/2) |Aut(T)|).

Here pi_n is the isomorphism-class distribution of a uniformly random labelled tournament, hence sums to 1.

**Corollary B6 (PROVED).**

    pi_n(T) B_n(T,U) = pi_{n+1}(U) D_{n+1}(U,T).       (B5)

Therefore birth followed by uniform deletion gives a Markov chain
P_n=B_n D_{n+1} on size-n classes, reversible with stationary measure pi_n. It is irreducible: by successive replacement one may replace all vertices and choose all edges of any target. It has self-loops by deleting the newborn. These elementary facts do not estimate a mixing time.

Important distinction: uniform incidence births do NOT produce a uniform distribution on unlabelled classes. They produce the reciprocal-automorphism weighting pi. Also B5 compares adjacent levels; pi_n alone is not a stationary distribution for a chain whose size always increases.

## 5. Scalar flatness versus transport memory

For a literal birth chain T_0 subset T_1 subset ... subset T_k, compute o_i,r_i from its actual pointed steps. B1 telescopes:

    log |Aut(T_k)| - log |Aut(T_0)|
       = sum_i (log o_i - log r_i).                  (B6)

Thus two birth paths with the same endpoint isomorphism types have equal net scalar symmetry gain, even if their separate fracture and replacement budgets differ. Assign inverse ratios to death edges: their product is 1 around every closed metagraph walk. This scalar 1-cocycle is exact. It does NOT force compatible transport arrows around a rooted loop to multiply to the identity: those can retain the nontrivial stabilizer element of §2.

This separation connects to THM-3345: flat scalar endpoint data do not determine the actual source-dependent word. Here the scalar potential is log automorphism order, while the retained information is a permutation of the marked parent. The analogy preserves compositional dependence; it does not identify the Berggren tree with a tournament graph.

## 6. Interaction with inherited symmetry barriers

**INHERITED:** TR5 says that for a nontrivial subgroup G<=Aut(T), g=|G|, and mu_G=min H(T-{u,v}) over distinct vertices in a common G-orbit, every finite induced extension V obeys

    H(V)=H(T), or H(V)-H(T)>=2 min(g,mu_G).

We do not re-prove or re-test that theorem here. B1 supplies its missing group-theoretic companion. Preserving the original action is an incidence condition, not a condition on total child group order. A birth may destroy the old subgroup and still acquire new symmetries through compatible replacements. Therefore an observed large |Aut(U)| cannot certify that TR5's original group survived. Conversely survival only guarantees K contains that group, not that the newborn is fixed by all of Aut(U).

Write h(T)=H(T)/|Aut(T)|, an integer by the inherited free action on Hamilton words. B1 gives the exact identity

    h(U)/h(T) = [H(U)/H(T)] r/o.                      (B7)

Thus increased Hamilton count and increased raw symmetry can pull the path-orbit count in opposite directions. THM-4987's subgroup-index cost is r's natural companion, but its low-quotient classifications cannot be applied solely from B7: their exact hypotheses and inherited path-count bounds remain necessary.

## 7. Next mechanisms, explicitly OPEN

1. Quantify the distribution of the pair (r,o), rather than merely the ratio, under uniform incidence births. Can old group orbit structure predict symmetry destruction while deletion transporters predict exceptional large gains?
2. Develop a marked metagraph retaining transporter bitorsors. A scalar reversible chain is already supplied by B5; a representation-valued refinement should retain monodromy, but needs a specified marking space and transition rule before any spectral claim.
3. Separate irreducibility from efficient mixing for P_n. B5 gives a canonical measure, not a polynomial mixing theorem.
4. **RESOLVED later in this same lane, §9:** families where one step both breaks old symmetry (r>1) and gains net symmetry (o>r). The carousel examples isolate the currencies; cyclic block substitution below combines them.

## Reproduction and scope

`python3 05-knowledge/results/session248_recursive_symmetry_20261110/birth/probe.py`
checks the newly defined ratios and weighted balance on four small named parents (TT2, C3, a punctured R5, R5), using exact permutations and all their birth incidences. This is a cheap diagnostic, not a census or a re-verification of inherited claims. The proofs above are independent of that finite check.

## 8. Cross-lane identity: retained cores are conditional birth histories

**PROVED.** This section connects the one-birth chain to the bridge lane's retain-k/resample-exterior chain. All orientations in the resampled exterior, including edges between two exterior vertices, are independently fair. Retained-core edges alone are held fixed.

Fix 0<=k<=n and put M=binom(n,2)-binom(k,2). Let

* d_k(U,T) count k-element subsets W of V(U) with U[W] isomorphic to T;
* e_{k,n}(T,U) count completions of a fixed labelled T on the first k vertices by n-k distinct new labelled vertices, resulting in unlabelled type U.

Write a_T=|Aut(T)| and a_U=|Aut(U)|, with the empty tournament having automorphism order 1. The multistep multiplicity identity is

    e_{k,n}(T,U) a_U = (n-k)! d_k(U,T) a_T.           (B8)

Proof: count pairs consisting of a labelled completion of T and an isomorphism to the fixed representative U. The left side counts by completion. On the right choose the image W of the first k vertices (d_k choices), an isomorphism from T to U[W] (a_T choices), and a bijection from the remaining n-k labels to the exterior ((n-k)! choices). For n=k+1 this is B4.

Define

    D_{n:k}(U,T)=d_k(U,T)/binom(n,k),
    B_{k:n}(T,U)=e_{k,n}(T,U)/2^M.                    (B9)

Deleting one uniformly chosen vertex at each stage leaves a uniformly chosen k-subset: every survivor set has (n-k)! deletion orders. Adding uniformly random incidence subsets successively exposes every new edge once. Thus, as row-stochastic matrix products,

    D_{n:k}=D_n D_{n-1} ... D_{k+1},
    B_{k:n}=B_k B_{k+1} ... B_{n-1}.                 (B10)

Products are identities when k=n. By B8,

    pi_n(U) D_{n:k}(U,T)=pi_k(T) B_{k:n}(T,U).        (B11)

This is also the exact reverse conditional. Sample U from the uniform-labelled class law pi_n and then a uniform k-subset. Conditional on its induced isomorphism type being T, the distribution of U is B_{k:n}(T,·), not a uniform distribution among extending unlabelled classes. In fact sum_U pi_n(U)D_{n:k}(U,T)=pi_k(T), so Bayes' rule gives B11 directly.

### The retain-and-resample factorization

The bridge kernel on size-n classes is exactly

    R_{n,k}=D_{n:k} B_{k:n}.                          (B12)

Indeed, choose the retained k-subset of the current labelled tournament, keep its actual internal orientations, and resample every other edge. Conditional on its induced class T, relabelling the core changes no resulting unlabelled probabilities, so the extension law is B_{k:n}(T,·). Summing over T proves B12. It follows that

    pi_n(U) R_{n,k}(U,V)
      = sum_T pi_k(T) B_{k:n}(T,U) B_{k:n}(T,V),       (B13)

which is symmetric in U,V. This is a two-step conditional resampling identity, not merely a common stationary measure.

The distinction between the two orders is useful. Our original one-birth/delete chain is P_n=B_n D_{n+1}, whereas retaining n-1 vertices at size n gives R_{n,n-1}=D_n B_{n-1}. More generally define the excursion chain Q_{k,n}=B_{k:n}D_{n:k}, which temporarily grows a k-vertex tournament to n vertices and returns to k. The nonzero spectra of Q_{k,n} and R_{n,k} coincide, including algebraic multiplicities, since they are the products of the same two rectangular matrices in opposite orders. This relates the birth-history and exterior-resampling viewpoints without claiming they are the same chain or state space.

Both chains have nonnegative spectrum. Explicitly for any real f on size-n classes,

    <f,R_{n,k}f>_{pi_n}
      = sum_T pi_k(T) [sum_U B_{k:n}(T,U)f(U)]^2 >=0. (B14)

Thus R is a positive semidefinite reversible contraction. The rank bound rank(R_{n,k})<=number of k-vertex isomorphism classes is immediate. B14 does not provide sharp eigenvalues or a mixing-time estimate; the bridge lane's spectral description adds information beyond this factorization.

### Marking changes mass: the exact orbit correction

An unlabelled subset-marked type (U,W), |W|=k, has automorphism group A_W, the SETWISE stabilizer. Under a uniform labelled U and a uniform retained k-subset, its probability is

    pi_n(U) |Aut(U).W|/binom(n,k)
       = n!/[2^binom(n,2) binom(n,k) |A_W|].         (B15)

In particular rooted types (U,x) have mass

    (n-1)!/[2^binom(n,2) |A_x|].                    (B16)

Uniform root choice is uniform on vertices, not on root-orbit types. Formula B16 is precisely why the stabilizer, rather than the full automorphism group, weights pointed births. Summing B15 over subset orbits recovers pi_n(U).

If instead the k retained vertices are individually ordered, use the POINTWISE stabilizer A_(w1,...,wk) and denominator (n)_k=n!/(n-k)!; the marked-type mass becomes

    (n-k)!/[2^binom(n,2) |A_(w1,...,wk)|].           (B17)

Setwise and pointwise marking are not interchangeable. If the retained ordered tournament is fixed as a particular labelled T, conditioning divides B17 by its prior probability 2^-binom(k,2), giving (n-k)!/[2^M |A_(w1,...,wk)|] for each compatible ordered-marked extension type. This recovers the extension law at the fully marked level and explicitly identifies the sidecar forgotten when one keeps only the unrooted tournament class.

### Exact spectral bridge and graded commutation

The [spectral proof](../bridges/PROOF.md), using inherited signed-even-graph coordinates, strengthens B14: on unlabelled classes R_{n,k} has eigenvalues binom(k,v)/binom(n,v), with multiplicity t_v-t_{v-1}, where t_v is the number of v-vertex tournament classes and t_{-1}=0. Thus its rank is exactly t_k, not just at most t_k. These diagonalization statements belong to that lane; B8–B17 supply their reversible probabilistic interpretation.

A direct cross-lane identity also relates the adjacent up-down and down-up chains:

    B_n D_{n+1} = [I + n D_n B_{n-1}]/(n+1).         (B18)

Proof: following the birth, uniform deletion removes the newborn with probability 1/(n+1), returning the starting tournament. Otherwise it removes a uniform old vertex; the newborn's edges to all n-1 survivors are independently fair, exactly the one-vertex refresh D_n B_{n-1}. Hence B18. The common eigenbasis in the bridge proof then gives eigenvalues 1-v/(n+1) for B_n D_{n+1}, versus 1-v/n for D_n B_{n-1}, at the same order n and with the same support multiplicities. This is an exact graded birth/death commutation law, not merely a comparison of stationary laws.

## 9. Degree-recognizable births and simultaneous fracture/gain

### Regular children cannot fracture the parent's group

**Proposition B7 (PROVED).** If U is regular of degree m and T=U-x, then every automorphism of T extends fixing x. Consequently r=1, and the symmetry ratio is exactly o.

Indeed, for v in N_U^+(x) one has d_T^+(v)=m, while for v->x one has d_T^+(v)=m-1. Thus the deleted parent's own degree function recovers the newborn incidence. Every parent automorphism preserves it. More generally the same conclusion holds whenever the two sets of parent outdegrees attained on N_U^+(x) and its complement are disjoint. This sufficient condition need not be necessary.

Thus mixed fracture and net gain (r>1 and o>r) requires an irregular child. Restricting a search to regular symmetric banks would hide this regime entirely.

### A concrete infinite mixed regime

**Proposition B8 (PROVED).** For every odd q>=5 there is a one-birth pair with r=3, o=q and symmetry ratio q/3>1.

Take B=R_q-0, rigid by Proposition B3, and form

    T=C3[B,B,B],
    U=C3[R_q,B,B],

where brackets mean homogeneous substitution into the three cyclically oriented blocks. The new vertex completes the first B to R_q. The child is irregular: its blocks have sizes q,q-1,q-1.

Here Aut(T)=C3 and Aut(U)=C_q. To prove this without assuming a modular-decomposition theorem, observe that any proper module meeting two cyclic blocks must meet no such pair: if it contains a in A and b in B, every c in C distinguishes them (b->c->a), so all C lies in the module. Any vertex of A outside the module then distinguishes b from c; similarly every remaining vertex of B does. Thus the module would be the whole tournament. The same argument holds for each cyclic pair. Therefore the three substituted blocks are exactly the maximal proper modules and must be permuted by every automorphism.

For T, the equal rigid blocks admit precisely the three cyclic block rotations. For U, the unique block of size q forces every block fixed, and the only automorphisms are the C_q rotations inside R_q. The newborn orbit has size q, and its stabilizer is trivial; the parent's C3 rotates the exceptional incidence pattern into three distinct patterns. Hence o=q and r=3.

At q=5 this gives a completely explicit 12-to-13-vertex step with |Aut(T)|=3 and |Aut(U)|=5. The two actual automorphism groups even have coprime orders: every nontrivial old symmetry is destroyed while a larger, entirely new group appears. The step attains the general fracture-growth bound o/r<=q/3 relative to its actual newborn orbit q, although it does not attain the size-only bound at child size 3q-2.

This resolves OPEN direction 4 of §7 positively. Directions 1–3 remain open. It also sharpens the TR5 interpretation: total symmetry gain is not evidence of persistence of any old nontrivial subgroup.
