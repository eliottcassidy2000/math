# Every prime-eight quotient converse is visible through arbitrary children and ambient

**Status: finite-exact theorem with unbounded child size, pending independent packet review.** The finite quotient classification and the three-anchor certificate are exact. Neither an all-order quotient theorem nor unrestricted named-response injectivity is claimed.

## 1. Statement and physical convention

Let Q be any prime tournament on eight quotient names, and let A_0,...,A_7 be arbitrary nonempty tournaments on fixed disjoint physical child sets. Compare

    T=Q[A_i],  U=Q^op[B_i],

where B_i uses exactly the same physical names as A_i and has the same complete ordered-path-cover profile. Then a child-uniform two-port query distinguishes the complete responses in every fixed ambient tournament. One fixed query suffices for all but two quotient isomorphism types; a fixed three-query battery suffices for each of those two types, independently of all child orders and orientations.

Ports x,y have the fixed edge x->y. Bits a_i,c_i mean all physical names in child i point to x,y respectively. To lift the query through a fixed ambient root, both ports inherit that root's incidence to every exterior name. Source and target names and port bits remain aligned. A relabeling used to choose a quotient representative transports the entire child sets and query masks together, never each response independently.

This is an arbitrary-child and arbitrary-fixed-ambient quotient-converse theorem. Its original-module family consequence is stated separately below. An arbitrary sequence using newly created modules is outside scope.

## 2. General inherited endpoint and singleton-sensitive anchor mechanism

The pinned S188 dependencies, at commit `8f14b96ab020f7d5ceeb88e7bda413bb6f8bc9b7`, are under `05-knowledge/results/session188_structural_bridges_20261104/`:

- `parity_profile_transfer/PROOF.md`: nonzero quotient endpoint parity, universal cover residues, and the complete modulo-four ambient lift;
- `uniform_block_starts/PROOF.md`: the exact uniform-block start identity;
- `mod8_transfer/PROOF.md` and `quotient_cycle_parity/PROOF.md`: the complete higher-cover congruences and odd-order cycle correction;
- `module_decomposition/PROOF.md` and `monotone_dual_queries/PROOF.md`: hereditary original-module monotonicity and full-profile converse-complement duality.

The singleton-sensitive combination below is developed in the separate Session189 prime-seven profile packet; we reproduce its complete argument so this proof does not assume that quotient-order-seven computation.

Write S_i(Q),E_i(Q) for the integer quotient Hamilton-start and Hamilton-end counts, p_i=S_i+E_i modulo two. If p_i=1, the actual query (a,c)=(e_i,0) has contrast two modulo four for arbitrary nonempty children, and complete augmented-profile contrast (2,0,0,...) modulo four. It therefore survives every fixed ambient.

It remains to treat quotients with p=0. Let v_j be one exactly when child j is nonsingleton, and V its bitmask. Let R_ij(Q) count directed words starting at i, visiting every quotient name once and j one additional time. Repeated copies are indistinguishable. Let Q_i^+ add one port y with j->y exactly when j=i, and let K_ij(Q_i^+) count full words visiting every old name and y once, with j one additional time. Set

    t_i(V)=S_i-E_i+2 sum_j v_j[R_ij(Q)-R_ij(Q^op)] modulo4,

    d_i(V)=H(Q_i^+)-H((Q^op)_i^+)
             +2 sum_j v_j[K_ij(Q_i^+)-K_ij((Q^op)_i^+)] modulo4.

Then the exact physical-query residues are

    Delta(0,e_i)=d_i(V) modulo4,
    Delta(e_i,0)=2t_i(V)+2d_i(V) modulo8.                 (1)

**Why arbitrary children are allowed.** Divide each complete child profile by its odd Hamilton count. Modulo four, its first coordinate is one; its second is2v_j; every coordinate from four onward vanishes; and its third is even. Two repeated children vanish modulo four. A three-run child has two nonempty disjoint return segments. Their interchange preserves all arcs, ports and endpoint constraints and is fixed-point-free, so its even profile factor kills that term. Thus the normalized aggregate start skew is t_i, and the one-port contrast is d_i. Query(0,e_i) has x as a universal source in its augmented core; deleting x identifies its count with that one-port extension.

The common normalization factor L=product_j H(A_j) is odd. Each d_i is even, and p_i=0 implies t_i is even. Consequently L fixes d_i modulo four and fixes2t_i modulo eight. Apply the exact inherited identity

    Delta(e_i,0)-2Delta(0,e_i)=2[S_(A_i)(T)-S_(A_i)(U)]

to obtain(1). This argument concerns complete physical responses, not a quadratic truncation or independent endpoint marginals.

**Why every fixed ambient is allowed.** A d_i=2 query has full profile gap(2,0,...) modulo four because all higher cover residues depend only on the common augmented order. For the second query, its Hamilton residue is0 or4 modulo eight. If the augmented order is even, the inherited higher-cover formulas make every higher contrast zero modulo eight. If the order is odd, the only vertex entering y is x, so every Hamilton cycle contains x->y. Deleting both ports is a bijection to old Hamilton paths ending in child i. The source-target cycle parity difference is p_i=0, independently of child sizes. Thus the odd-order correction also vanishes. Every h>=4 ordered-cover count is separately zero modulo eight by h! divisibility.

Therefore the second anchor always has full profile contrast(4b,0,...) modulo eight for its bit b. Under fixed ambient substitution these pure first-coordinate contrasts multiply by the odd Hamilton count of the contracted ambient, preserving the stated bits. Nonzero local scalar gaps are not lifted without this profile argument.

## 3. Exactly two endpoint-zero prime-eight types

The complete degree-partition canonical enumeration gives6,880 eight-name isomorphism classes, of which4,008 are prime. Of those,4,006 have nonzero endpoint parity and are covered by the first criterion of Section2. The only endpoint-zero representatives are

    Q0=(4,17,42,131,45,75,31,119),
    Q1=(16,5,41,67,46,139,55,95).

Their exact Hamilton count is249. Their second-endpoint residues in the all-nontrivial case are

    sigma(Q0)=(2,0,0,0,2,2,2,0),
    sigma(Q1)=(0,2,2,2,0,0,0,2).

The finite classification retains every prime representative and every exact integer start/end count. No zero test is treated as a response collision.

These are not newly invented obstructions: they are exactly the S188 endpoint obstruction `(238,4,240,118,195,146,34,74)` and its converse, up to explicit checked physical relabeling. The target-index-to-original-index permutations are respectively

    known obstruction -> Q1: (1,6,5,7,2,4,3,0),
    its converse -> Q0: (0,3,4,2,7,5,6,1).

Thus the finite structural residual is one known converse pair, rather than a new unclassified large family.

## 4. Three affine anchor bits cannot vanish together

For Q0, the exact word counts in Section2 give, over F_2,

    Delta(0,32)/2 = v_2+v_4                 modulo2,
    Delta(0,128)/2 = 1+v_3+v_5+v_6         modulo2,
    Delta(1,0)/4 = v_2+v_3+v_4+v_5+v_6     modulo2.       (2)

The first two equations specify contrasts modulo four; the third specifies a contrast modulo eight, known to be divisible by four. Their three right sides sum to one. Hence at least one of the three fixed actual queries `(0,32),(0,128),(1,0)` separates, for every singleton-position vector and every allowed child profile.

For Q1, the corresponding certificate is

    Delta(0,1)/2 = 1+v_1+v_2+v_4           modulo2,
    Delta(0,4)/2 = v_3+v_5                 modulo2,
    Delta(128,0)/4 = v_1+v_2+v_3+v_4+v_5   modulo2.       (3)

Again the three right sides sum to one. Thus `(0,1),(0,4),(128,0)` is a fixed battery. All equations remain true for the lifted ambient bits by Section2. These identities are stronger than merely retaining a different witness for each of512 type patterns.

The coefficients in(2)-(3) are finite exact certificates on length-eight through length-ten quotient words. Their validity is independently checked by multiset-prefix dynamic programming and recursive completion in the audit. The affine argument, not a child-order cutoff, extends them to arbitrary children. All512 patterns and their chosen literal physical queries are retained as controls.

## 5. Completeness without a 2^28 labeled brute-force scan

`cover.json` proves coverage inductively from the empty tournament. At level n, take every retained(n-1)-name representative and every one of its2^(n-1) possible new-vertex incidence masks. For every such rooted extension, the file retains its target n-name representative and an explicit permutation carrying the extension to that target.

The numbers of representatives at orders0 through8 are

    1,1,1,2,4,12,56,456,6880.

The number of extension records at levels1 through8 is

    1,2,4,16,64,384,3584,58368.

Every record is checked edge by edge. In particular the456 seven-name representatives are proved to cover their universe by the earlier levels, rather than assumed from a table. Any n-name tournament can delete one vertex, map its remainder to a retained parent, and then appears among that parent's complete rooted extensions. This proves coverage at the next level.

For canonicalization, sort vertices by outdegree and enumerate all permutations within equal-degree groups. The lexicographically least resulting row tuple is an isomorphism invariant and is identical for two isomorphic tournaments: every isomorphism preserves degrees and transports the entire allowed permutation set. Distinct resulting canonical tuples are therefore distinct classes. The independent review additionally uses a different minimum integer edge-code convention and automorphism/orbit totals to check irredundancy.

Every retained eight-name representative is independently tested against all nontrivial proper vertex subsets for modules. The4,008 prime records have every integer endpoint count checked again by all8! permutations. This checks64,128 exact marginal coordinates, not only the zero-parity histogram. The audit additionally retains full joint start/end matrices.

## 6. Original-module consequences and precise boundary

For a fixed prime-eight source core and a fixed ambient, take the one- or three-query baseline battery above and add the core converse-complements

    J(a,c)=(255-c,255-a),

while retaining the same exterior port incidences. Complete augmented-profile duality makes the two baseline ambient contrasts opposite over the integers. At least one baseline query is nonzero, so one member of its pair has negative source-minus-baseline contrast.

Consider an XOR family of modules of the original source, contained inside this core, that selects the core root oddly. Every proper original module lies inside one prime child. Hereditary original-module monotonicity makes each target child profile coordinatewise no smaller; the whole-core reversal merely converses each interior and does not change its cover profile. Because the chosen query is child-uniform, its run-substitution coefficients are nonnegative. The strict negative contrast cannot be canceled by these profile increases. Hence a battery of at most six actual queries detects every such root-containing family from its source. The negative member may depend on the fixed ambient; the battery does not depend on the chosen family.

This is source-relative detection, not all-pairs injectivity among arbitrary family outputs. Operations on exterior names and newly created modules are not included in this root-family statement. The separate highest-changed-node restoration argument is needed to integrate this theorem with arbitrary original-module families at different tree nodes. This packet alone does not replace that argument or claim unrestricted named reconstruction.

## 7. Full physical banks and reproducibility

The verifier retains the complete1,369-coordinate named Boolean basis for each eight-name residual and its converse, including every zero, cubic and quartic coefficient. `physical_records.json` also stores all262,144 exact full response values across the four named banks, indexed by a+256c. These are full integer Hamilton counts on actual ten-name extensions. The full cubes are recovered losslessly from the complete named polynomial and independently checked at additional nonbasis queries. They are not full banks for arbitrary large child substitutions.

Literal vertex-subset DP verifies a selected actual query for all512 singleton/pair position patterns, plus eight cyclic-ambient controls, on6,548 total physical extensions of orders at most20. Counts fit unsigned64-bit integers since20!<2^64. All original child names, physical port masks, exact response values and complete smaller named banks remain retained. No independent relabeling or endpoint-product approximation is used.

Run `python3 verify.py > normal.json` and `python3 -O verify.py > optimized.json`, then compare. The programs use the Python standard library and the already installed C++17 compiler; no installation is performed. Checks use explicit exceptions, not optimization-removable assertions. The final verification record pins all source and data bytes.
