# The Walsh bijection preserves the whole Sylow p action

**Status: PROVED.** This strengthens `PRIME_VALUATION_THEOREM.md`. For each fixed odd prime p, the constructive p-unit Walsh matching preserves the S_n-conjugacy class of a Sylow p-subgroup of the automorphism group, not only its order. It follows that HYP-9340(a) holds for **every odd prime-power-order subgroup H<=S_n**. Mixed-prime odd-order H remain outside this proof.

## 1. The signed orbit-localization lemma

Fix an even graph X on [n]. For each labelled graph Y in its orbit choose

    epsilon(Y)=(-1)^r(g,X), where gX=Y,

and r(g,X) counts, modulo two, the increasing-label edges of X whose endpoints g reverses. Because X is sign-even, epsilon(Y) is independent of the chosen transporter g. Set

    w_T(Y)=epsilon(Y) phi_Y(T) in {+1,-1}.

Then Phi_X(T)=sum_Y w_T(Y).

**Lemma.** If P is any subgroup of Aut(T), the weight w_T(Y) is constant on every P-orbit of graph copies. If P is a p-group and p does not divide Phi_X(T), at least one graph copy Y is P-invariant.

**Proof.** For h in P the cocycle identity gives

    epsilon(hY)=epsilon(Y)(-1)^r(h,Y).

Transport of the oriented edge product gives

    phi_hY(hT)=(-1)^r(h,Y) phi_Y(T).

Since hT=T, multiplying these two formulas proves w_T(hY)=w_T(Y). Therefore each P-orbit contributes its size times a single sign. A P-orbit has p-power size; every non-singleton orbit contributes zero modulo p. Consequently

    Phi_X(T) = sum_{Y: P Y=Y} w_T(Y) mod p.             (1)

A nonzero left side forces at least one fixed copy Y. Thus P<=Aut(Y). □

This is the ordinary p-group orbit-counting mechanism applied to the **signed** terms. The sign transport must be retained; unsigned edge monomials need not be constant on these orbits.

## 2. Action-preserving bijection

Let beta_p be any perfect matching in the support of the Walsh matrix modulo p, for instance the deterministic cofactor construction in the main note. For every tournament T choose a Sylow p-subgroup P_T of Aut(T). The lemma supplies a copy Y of beta_p(T) with

    P_T <= Aut(Y).

The main theorem proves

    v_p(|Aut T|)=v_p(|Aut beta_p(T)|).

Thus P_T already has the full p-power order of Aut(Y), so it is a Sylow p-subgroup of Aut(Y). Transporting Y to the chosen canonical representative of beta_p(T) identifies the two Sylow subgroups up to conjugacy in S_n.

**Theorem.** For every n>=0 and every odd prime p there is an explicitly computable bijection on isomorphism classes

    beta_p: tournament types -> sign-even graph types

such that all of the following hold:

1. p does not divide Phi_beta_p(T)(T);
2. a Sylow p-subgroup of Aut(T) is conjugate in S_n to a Sylow p-subgroup of Aut(beta_p(T));
3. after an explicitly searchable relabelling of the output graph, the same selected P_T acts on both T and its output.

The canonical subgroup-action datum is well defined: all Sylow p-subgroups of one automorphism group are conjugate within that group, hence within S_n. The transporter in item 3 can be found by enumerating the graph copies and testing P_T-invariance, whose existence (1) guarantees. This gives a finite constructive algorithm, not an efficient-in-n claim or a canonical equivariant choice of a Sylow subgroup.

## 3. Every p-group symmetry predicate, simultaneously

For an embedded p-group H<=S_n define

    tau_T(H) = number of tournament types whose automorphism group
               contains an S_n-conjugate of H,

and tau_E(H) analogously for sign-even graphs.

Every p-subgroup of a finite group lies in a Sylow p-subgroup. It follows that the predicate “Aut(object) contains a conjugate of H” is determined by the S_n-conjugacy class of the object's Sylow p-subgroup. The theorem therefore gives

    tau_T(H)=tau_E(H)                                    (2)

for **all p-groups H<=S_n**, for each odd prime p. A single beta_p preserves all these predicates simultaneously, including arbitrary finite Boolean combinations of them. In particular, it preserves existence of several specified p-subgroups even when they have different orbit structures on [n].

This closes the prime-power portion of HYP-9340(a), extending the inherited THM-4721 Theorem S from a Sylow subgroup of ambient S_n to arbitrary p-subgroups of ambient S_n. It includes the previously highlighted H=C3 acting as a single 3-cycle on six vertices, as well as diagonal C3 actions and C3-by-C3 actions. It does not identify conjugacy classes of H *inside the full automorphism group*, counts of embeddings, or fusion of such subgroups.

## 4. Refined block structure modulo p

Assign each row and column its Sylow p-subgroup's S_n-conjugacy class. The lemma implies that a p-unit matrix entry can occur only if the row subgroup is conjugate to a subgroup of the column subgroup. Ordering by subgroup order makes the matrix block upper triangular. The equality of p-valuation multiplicities from the main theorem aligns each order layer; inside a fixed order layer, the lemma forces zero entries between distinct S_n-conjugacy classes. Invertibility then forces equality of the row and column multiplicities of every Sylow-action class and invertibility of each corresponding diagonal block.

For clarity, the last assertion does not assume equal class multiplicities in advance. At a fixed order layer, the rectangular subblocks for different conjugacy classes have no off-diagonal connections. A square invertible layer can have no row or column rank defect; each nonzero class block must therefore be square and invertible. This yields a direct algorithm matching inside the Sylow-action blocks.

The initial cofactor algorithm also necessarily remains inside those blocks: any complete p-unit matching has equal p-orders at every pair, and then the containment from the lemma becomes conjugacy.

## 5. Exact failure boundaries

- This is one prime at a time. At n=5 the empty graph has Sylow 3- and Sylow 5-symmetries, while no five-vertex tournament automorphism group has order divisible by 15. No single bijection can preserve both prime-group predicates simultaneously.
- Preserving a Sylow p action does not preserve the whole automorphism group. At n=3 and p=3 the directed triangle is paired with the empty graph, whose automorphism groups are C3 and S3.
- It does not produce a labelled S_n-equivariant map. Canonical representatives, a selected Sylow subgroup, and a searched output transporter are additional choices. The known full-automorphism equivariance obstructions survive.
- The p-group hypothesis is essential to this proof. For a mixed-prime group, nontrivial orbit sizes need not all be divisible by one common prime, so equation (1) no longer forces a globally fixed copy. The remainder of HYP-9340(a) for mixed-prime odd groups stays OPEN here.

## 6. Source and novelty boundary

The linear Walsh basis and weighted orthogonality are inherited from THM-4721, with their elementary derivation in `PRIME_VALUATION_THEOREM.md`. The prime-local matching and the orbit-localization implication above are the additional steps. The argument uses only the standard Sylow containment/conjugacy theorems and p-group orbit sizes; no external research theorem about p-permutation modules is being claimed or imported. No global literature-priority claim is made.

## 7. Direct finite action controls

`walsh_sylow.py` checks the localization congruence by explicitly enumerating P-fixed signed graph copies for **every** p-unit matrix entry whose tournament has nontrivial p-Sylow subgroup, for p=3,5,7 and every n<=6. There are respectively 1,4,13,97 such entries for n=3,4,5,6. For all 23 matched nontrivial cases across these orders, it saves an actual vertex permutation carrying the canonical graph to a P-fixed copy and verifies that P has the full output p-power order.

At n=6 the matched nontrivial actions are: ten C3 groups with vertex-orbit sizes (3,1,1,1), two C3 groups with orbit sizes (3,3), one C3-by-C3 group with orbit sizes (3,3), and two C5 groups with orbit sizes (5,1). The order distinguishes the two (3,3) action cases. This bounded verifier uses the checked fact that tournament automorphism orders through six are 1,3,5,9, so each nontrivial Sylow subgroup here is the whole automorphism group; the all-order proof does not use that shortcut.

Normal and optimized runs agree. The complete witnesses are in `walsh_sylow.json` and the short run transcript in `walsh_sylow.out`.
