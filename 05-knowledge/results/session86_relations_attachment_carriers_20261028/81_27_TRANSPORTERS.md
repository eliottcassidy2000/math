# Session86: identical local isotropy, missing port transporters

**Status:** PROVED elementary statements for the explicit family below; FINITE-EXACT independent verification. No new minimal-order, literature-priority, arbitrary-Frobenius, or universal reconstruction claim.

## 1. Pinned inheritance and loss ledger

Read-only public repository source pin: `5ed1835235b36a57b8cd11eb705aa259b6bc0f08`. Paths below are relative to the repository root.

- `05-knowledge/results/session81_counterexample_families_20261028/LIFTING.md`, Sections 1, 3, 4: source/sink-converse leaves, module run reversal, complete subset-H/path-cover equality, and the inherited 81/27 seed.
- Same packet, `OBSERVER_GRAMMAR.md`: named restriction and uniform substitution preserve the exact cyclic-triple observer quotient; arbitrary orientation-sensitive attachments are outside that grammar.
- `01-canon/MISTAKES.md`, “Cyclic triples lose child colours before they lose path counts” and “Session67 incoming audit”: descendant reversal changes ordinary child type and the ancestor's color stabilizer. Neither abstract child groups nor cyclic-triple data determine this stabilizer.
- `05-knowledge/results/session85_arithmetic_codes_free_bits_20261028/LINEAR_GLUING_INFORMATION.md`, Sections 2–5: marginal joins retain supported equations, not hidden joint constraints; affine rank needs actual coefficient rows and right-hand sides.

Closest mechanism: canonical decorated substitution. Canonical hostile: 81/27. Corrected near miss: retaining the abstract quotient and all child groups without retaining child identity and transport. Useful sidecar: **transporter sets**, not only stabilizer groups. Board: cyclic triples / endpoint readouts / module reversal / role-colored ports / transporters / invariant versus fixed relations.

## 2. The exact pair and an elementary full-group proof

Use vertices `(i,q)` for `i=0,1,2`, `q=0,1,2,3`, encoded `4i+q`. Let `A0` have source `q=0` and cycle `1→2→3→1`; let `A1` be its **literal converse on the same labels**, including reversal of the three-cycle. Set

    T  = C3[A0,A0,A0],
    T' = C3[A0,A0,A1],

with block arcs `0→1→2→0`. In both tournaments the cyclic triples are exactly the three within-leaf cycles plus all `4^3` triples selecting one vertex from each block: **67 identical labelled triples**. Triples with two vertices in one block and one outside cannot be cyclic.

Write `c_i=(4i+1,4i+2,4i+3)` and let `rho` advance the block index while retaining `q`. These are actual arc-preserving permutations wherever used.

For T, the three source vertices have degree 7 and the remaining nine vertices degree 5. Any automorphism preserves the source set and acts there by a cyclic rotation. Each cycle vertex's block is recovered from its relations to the three sources (there is exactly one source it beats), so no other block mixing is possible. Within each recovered block the automorphism group is exactly its cycle rotations. Hence

    G=Aut(T)=<c_0,c_1,c_2> ⋊ <rho> ≅ C3^3 ⋊ C3,
    |G|=81.

For T', vertex 8 has the unique degree 4; vertices 9,10,11 have degree 6; vertices 0,4 have degree 7; the other six vertices have degree 5. The degree-7 pair is individually fixed because its induced arc is `0→4`. Their incidences separate the two degree-5 cycles; the degree-6 cycle is already a separate class. Thus

    H=Aut(T')=<c_0,c_1,c_2> ≅ C3^3,
    |H|=27.

This proof establishes the complete groups without importing modular-decomposition exhaustiveness. Unrestricted arc-preserving backtracking independently verifies it. G has element-order counts `1:1, 3:44, 9:36`; H has `1:1, 3:26`.

## 3. All nonempty pointwise stabilizers agree; transporters do not

The two displayed groups are literal permutation groups on the same twelve labels, and `H<G`. Every element of `G\H` moves all three blocks cyclically and therefore fixes no vertex. Consequently, for **every nonempty named subset S**,

    {g∈G : g(s)=s for all s∈S}
      = {h∈H : h(s)=s for all s∈S}.

In particular all actual named vertex-stabilizer subgroups, all nonempty ordered-tuple stabilizer subgroups (including repeated entries), and all their orders agree. This is stronger than an equality of abstract local groups.

But the transporter from vertex 0 to vertex 4 has size 27 in G and is empty in H. Local isotropy does not retain arrows between different named ports. The orbit partitions are

    G: {0,4,8}, {1,2,3,5,6,7,9,10,11};
    H: {0}, {1,2,3}, {4}, {5,6,7}, {8}, {9,10,11}.

For any nonempty ordered tuple x, orbit–stabilizer gives `|Gx|=3|Hx|`. Equivalently, the residual `G/H=C3` acts freely on H-orbits of nonempty tuples. Thus for any G-stable collection of such tuples, its number of H-orbits is three times its number of G-orbits. The finite distinct-port checks are:

| Ordered distinct ports | T orbit count | T' orbit count | T orbit-size histogram | T' orbit-size histogram |
|---:|---:|---:|---|---|
| 1 | 2 | 6 | `3:1, 9:1` | `1:3, 3:3` |
| 2 | 12 | 36 | `3:2, 9:8, 27:2` | `1:6, 3:24, 9:6` |
| 3 | 72 | 216 | `3:2, 9:38, 27:30, 81:2` | `1:6, 3:114, 9:90, 27:6` |

A histogram entry `s:k` means k orbits of size s. **Do not replace “pointwise” by “setwise”:** the nonempty set `{0,4,8}` is stabilized setwise by the outer rotation in G, but H has no such rotation. Do not include the empty tuple, whose stabilizers are G and H themselves.

## 4. Explicit role-colored attachments distinguish the sources

Keep all twelve original vertices blue. For an ordered tuple of distinct ports `(p_1,...,p_r)`, add uniquely colored markers `z_1,...,z_r`. Marker `z_i` beats exactly `p_i` among the blue vertices; every other blue vertex beats `z_i`. Orient marker-marker arcs by their fixed role order.

A color-preserving attachment isomorphism must preserve the blue source, fix each marker role, and carry the designated tuple to the target tuple. Conversely every such source transporter extends uniquely. Therefore attachment isomorphism types correspond exactly to the tuple orbits above. One red role gives **2 versus 6** types; red/green roles give **12 versus 36**. At every fixed named tuple, however, the attachment automorphism groups agree literally under restriction to the blue labels, by Section 3. Counts of self-symmetries do not reconstruct the attachment groupoid.

There is also a direct scalar path readout, not merely a distinction of isomorphism classes. Attach one red marker at **named port 8**. In T's attachment `(8,9,z)` is cyclic; in T' it is transitive. Independent Hamiltonian endpoint DP gives

    H(T with red marker at 8)  = 1,395,279,
    H(T' with red marker at 8) = 1,061,775.

This does not contradict the inherited restriction–substitution theorem: the marker treats vertices within the third leaf nonuniformly, so it breaks that leaf's module property. It is an explicitly richer operation.

## 5. What unmarked path data still preserve

Reversing the third leaf module maps every directed path to a path by reversing each maximal run inside that module. It preserves the exact path vertex set. Applied componentwise it preserves any unordered path cover and the component vertex sets. Hence every named induced subset has the same H and every path-cover count `pc(S,d)` in T and T'.

The independent check compares **all 4,096 subset H values** and the entire subset-path-cover table for all eight binary leaf words. It obtains `H(T)=H(T')=462,915`; the full-set path-cover counts for `d=1,...,12` are

    462915, 4505949, 13435794, 18254079, 13393350, 5803110,
    1558116, 264627, 28305, 1839, 66, 1.

Endpoint information is not preserved. For T the endpoint-count vector is `(31185,41040,41040,41040)` repeated three times. In T' the first eight entries agree, while the last four are `(47142,35721,35721,35721)`. The common canonical tree, root quotient, abstract child groups and child group cycle indices therefore do not repair the lost endpoint or ordinary-type data.

## 6. Invariant distribution is not pointwise fixed alignment

On `X=F2^3`, let `sigma(x0,x1,x2)=(x2,x0,x1)`. The uniform law on all eight points is cyclically invariant. Its orbit sizes are `1,1,3,3`; it certainly does **not** force each sampled alignment to be fixed. The pointwise fixed constraint is

    sigma x=x  iff  x0=x1=x2,

and only this additional constraint reduces eight possibilities to two. It has rank two over F2. Invariance of a family means `sigma R=R`; invariance of a law means `mu(sigma^-1 E)=mu(E)`. Neither says `R⊆Fix(sigma)` or `mu(Fix(sigma))=1`. The distinction applies to an actual Frobenius action as well: a Frobenius-stable set need not consist of Frobenius-fixed points. This three-bit model is a cyclic-action control, not an unproved identification of arbitrary port coordinates with a field's normal basis.

A stronger same-form/marginal hostile uses the two disjoint relations

    R_even={x : x0+x1+x2=0},
    R_odd ={x : x0+x1+x2=1}.

Both are sigma-stable, have size four, use the **same ambient form** `B(x,y)=Σ x_i y_i`, and sigma preserves that form. Their uniform laws have identical full pairwise marginals on every coordinate pair. Joining even just the bags `{0,1}` and `{1,2}`—a valid running-intersection tree—returns all eight assignments. Retaining all three pair projections still returns eight. The missing one-bit joint parity, including its right-hand side, is not supplied by the ambient form, cyclic invariance, or bag/LCA structure.

The repaired positive statement is conditional: if actual owned constraints and shared variable identities are retained, compatible joining is exact for those constraints; if an actual pointwise symmetry equation is imposed, its fixed space can be computed. No automatic inference “Frobenius + LCA recovers joint ports” follows from these structures' mere presence. For orbit/stabilizer targets, retain the actual action and transporters on witnesses as well as the relation, not just marginal tables or isotropy summaries.

## 7. Reproduction and finite scope

From this directory:

    python3 -B check_hostiles.py > check_hostiles.json
    python3 -B -O check_hostiles.py > check_hostiles.optimized.json
    cmp check_hostiles.json check_hostiles.optimized.json

The standard-library checker uses explicit exception guards, not `assert`. It builds adjacency matrices independently of the pinned implementation, performs unrestricted automorphism backtracking for all eight order-12 sources, computes every subset H, derives unordered covers by anchored set partitions, verifies all nonempty named pointwise stabilizer sets, enumerates ordered-port orbits through three ports, and separately recomputes colored attachment automorphisms. The source emits deterministic receipts to standard output. All finite universes and exclusions are explicit; the all-subset-stabilizer and invariant-versus-fixed proofs do not rely on extrapolating a finite search.
