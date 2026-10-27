# Individual tournament symmetry: the decorated modular stabilizer

Status: **PROVED deduction from cited canonical decomposition; FINITE-EXACT checks**.
No literature-novelty claim. This packet is not a prime census and does not
recover individual groups from the aggregate prime cycle index.

## Inheritance and boundary

The structural input is Gallai's canonical decomposition, as stated in
Theorem 16 of Boussaïri–Chergui–Ille–Zaidi, *3-uniform hypergraphs: modular
decomposition and realization by tournaments* (2018),
https://arxiv.org/html/1805.04607v1 (primary text checked 2026-10-27).
Maximal proper strong modules give either a transitive quotient or a prime
quotient. The transitive root is the intrinsic ordered SCC decomposition.
The [prime-cycle-index recovery](../session79_tournament_enumeration_20261027/PRIME_CYCLE_INDEX.md) supplies the aggregate species side; the
[corrected composition work](../tournament_composition_classification_20261026/README.md) supplies the C3 hostile.
These results concern ordinary directed tournaments, not converse classes.

Concept board: canonical blocks; child isomorphism colors; quotient permutation
action; coherent transport; actual vertex cycles. The anchor is the exact group,
the niche is its cycle index, and the wildcard is symmetry breaking by colors.
The missing coordinate in the C3 near miss is ordinary child isotype, not child
group order, abstract group, or even child vertex-action cycle index.

## 1. Exact split sequence

Let T=Q[A_i : i in I] be a canonical root decomposition with nonempty children.
Put G_i=Aut(A_i), and

    H = {h in Aut(Q) : A_i is isomorphic to A_{h(i)} for every i}.

There is a canonical exact sequence

    1 -> product_i G_i -> Aut(T) -> H -> 1.                 (1)

It splits, but the following splitting depends on choices. For each child
isotype c choose one representative R_c and isomorphisms
phi_i : R_c -> A_i. Lift h by

    s(h)(phi_i(x)) = phi_{h(i)}(x).                        (2)

Then s(hk)=s(h)s(k), so

    Aut(T) is isomorphic to (product_i G_i) semidirect H.  (3)

In representative coordinates H simply permutes factors belonging to the same
isotype. On actual blocks it also transports by phi_j phi_i^{-1}. Changing
all phi_i conjugates this chosen complement by an element of the kernel.
We claim no intrinsically preferred splitting or preferred complement.

**Proof.** Strong modules are defined intrinsically and preserved by every
automorphism. Hence an automorphism permutes canonical root blocks, induces an
automorphism of Q, and restricts to isomorphisms between its source and target
children. Its induced permutation is therefore in H. Conversely every such
block permutation and every choice of child isomorphisms preserve internal
arcs and cross-block arcs, hence give an automorphism. The kernel consists
exactly of independent child automorphisms. Formula (2) gives surjectivity and
a multiplicative section. This proves all statements.

At an ordered SCC node H=1, even if several components are isomorphic, hence
Aut(T)=product_i Aut(A_i). Repeated equal children never justify arbitrary
permutations. At a prime node with all children equal to A the answer is the
imprimitive wreath product Aut(A) wr Aut(Q), with **the actual quotient action**.
For mixed child types, (3), not a product of full symmetric wreath products,
is the safe general statement: a single quotient element may couple motions
on several color classes or orbits.

## 2. Order and actual action

    |Aut(T)| = |H| product_i |Aut(A_i)|.                   (4)

Recursing down the canonical tree yields

    |Aut(T)| = product_{internal nodes v} |H_v|.          (5)

Linear-node factors are one. The leaves are singletons. Therefore T is rigid
iff every local decorated stabilizer H_v is trivial. This is an iff, not merely
a sufficient condition. All factors and hence the answer are odd: a tournament
automorphism cannot have an even cycle, since a suitable power exchanges an
antipodal pair and reverses its arc; an even-order finite group has an involution.

To construct generators, recursively embed generators of every G_i, and lift
generators of H by (2). This yields the full action on the original vertices.
To construct just the abstract group, quotient action on factors and child
abstract groups suffice *once H and compatible identifications are already
supplied*. Child group isomorphism does not supply the child isotype colors
needed to calculate H.

## 3. Cycle transport and normalized cycle-index recursion

For a specific lift g above h, write f_i:A_i -> A_{h(i)} for its restriction.
For a cycle C=(i_0 i_1 ... i_{ell-1}) of h, the return map is

    r_C = f_{i_{ell-1}} ... f_{i_1} f_{i_0} in Aut(A_{i_0}).

Every d-cycle of r_C gives exactly one ell*d-cycle of g on the union of the
blocks in C. Proof: one must make a multiple of ell steps to return to the
starting block, and exactly d such rounds to return to the starting vertex.
This describes individual permutations, not only averages.

Define Z_A=(1/|Aut A|) sum_g product_d p_d^{c_d(g)} and
psi_ell Z_A = Z_A(p_ell,p_{2ell},...). Then

    Z_T = (1/|H|) sum_{h in H}
                  product_{C cycle of h} psi_{|C|} Z_{A_{i_C}}.   (6)

Each h-cycle stays within one child isotype, so any representative i_C works.
For fixed h and one length-ell cycle, all |G|^ell tuples of internal maps are
equiprobable; each possible return map occurs |G|^{ell-1} times. Different
cycles are independent. Averaging the preceding cycle transport proves (6).
For linear roots this reduces to product_i Z_{A_i}.

One more action statistic follows without computing all permutations. If
r(A) is the number of vertex orbits of Aut(A), then

    r(T) = sum_{O an H-orbit on I} r(A_{i_O}).              (7)

Indeed the kernel acts independently within blocks; H identifies the matching
child orbits across precisely its block orbits. This uses coherent transport;
an automorphism of a child preserves each of its own group orbits.

## 4. Sufficient carriers, and what each forgets

* Exact tournament: canonical rooted modular tree, quotient tournaments with
  addressed child slots, and leaf attachments. This also recovers every group.
* Exact permutation group: recursively the decorated quotient actions H_v,
  child vertex actions, and coherent equal-type transports. Up to relabelling,
  representative coordinates suffice. Quotient arc sets are dispensable once
  their *correct decorated stabilizers* are already supplied.
* Individual normalized cycle index: H at each node with the type decorating
  every slot, plus each child cycle index. Equation (6) suffices; child groups
  themselves are unnecessary after H is known.
* Group order alone: the local integers |H_v| suffice via (5).
* Aggregate unrefined family counts: the summed prime cycle index of Session79
  suffices by species substitution. It does **not** supply any particular Q,
  its decorated stabilizer, the individual Aut(T), or converse information.

The maps group action -> cycle index -> order (identity coefficient 1/|G|)
lose information, and none is inverted here. Equal child actions or cycle
indices cannot be used as interchangeable colors for calculating H.

## 5. Hostile strengthened: same C3 and child actions, different parent exponent

Let A=1 -> C3 and A*=C3 -> 1. Their scores are (1,1,1,3) and (0,2,2,2),
so they are nonisomorphic. Their automorphism vertex actions are nevertheless
isomorphic: a fixed point and one cyclic orbit of length three. Put

    T=C3[A,A,A],  U=C3[A*,A,A].

The ordinary C3 hypergraphs agree on the same twelve labels. Internal cyclic
triples survive reversal, triples using exactly two blocks are transitive, and
transversal triples are unchanged. There are 3+4^3=67 such triples.
For T the root H=C3; for U its distinguished color forces H=1. Thus

    Aut(T) = C3 wr C3 (order 81),   Aut(U) = C3^3 (order 27).

Write z=(p1^4+2*p1*p3)/3. Formula (6) gives

    Z_T = (z^3 + 2*(p3^4+2*p3*p9)/3)/3,
    Z_U = z^3.

In particular T has 36 automorphisms of cycle type (3,9), and U has none.
Aut(T) has exponent nine while Aut(U) has exponent three: the upper rotation
combines with a nontrivial return rotation to multiply cycle lengths. This is
an explicit symmetry-creation mechanism beyond multiplying group orders.

More generally replacing child colors while leaving all child groups unchanged
can destroy or create a parent's quotient stabilizer. If a new coloring refines
an old coloring on fixed Q, H_new <= H_old and, provided child group orders stay
unchanged, the group-order ratio is [H_old:H_new]. Arbitrary recolorings need not
be refinements and their stabilizers need not be nested. The transport kernel
never compensates for an inadmissible child isomorphism.

## 6. Exact validation and stopping boundary

Run:

    python3 check_groups.py
    python3 -O check_groups.py

The checker independently extracts all strong modules by subset tests; constructs
the root quotient; recursively enumerates all valid child maps; and compares the
resulting full permutation set against independent global arc-preserving,
degree-pruned backtracking. It tests **all 1,099 labelled tournaments of orders
1 through 5**, without isotype filtering. For every case, the cycle-index formula
is also compared with the cycles of the actual reconstructed group. It separately
checks both order-12 hostiles against global backtracking, including identical
labelled C3 hypergraphs. Normal and optimized outputs agree byte for byte.
No assertion statements control validity. No external binary or package is used.
This finite check supports, but does not replace, the general proof above.


## 7. Independent audit and the canonical-partition guard

A separately written audit checked all 1,098 labelled tournaments of orders
two through five, both complete twelve-vertex groups and cycle indices, and
the exact image/kernel sizes 3/27 and 1/27. Normal and optimized outputs agreed.

Canonicality is essential for identifying the whole group. In
C3[C3,C3,C3], split just the first child into three singleton blocks and keep
the other two children whole. This is a modular partition but not an invariant
one: a top-level rotation sends its singleton blocks into an unsplit child.
The corresponding block-preserving formula gives order 27, whereas the whole
automorphism group has order 81. For an arbitrary modular partition the formula
computes its setwise stabilizer, not automatically the whole tournament group.
The independent checker verifies this counterexample explicitly.
