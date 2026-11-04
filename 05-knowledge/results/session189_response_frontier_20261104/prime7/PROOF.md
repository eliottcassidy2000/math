# Two all-cover queries resolve the even-child carousel-seven residual

**Status: PROVED by a structural profile reduction and a finite exact certificate, independently checked by literal physical interpolation.** The assertion covers four explicitly listed seven-name quotients, arbitrary even child tournaments of orders at least four, profile-preserving target replacements, and every fixed ambient. Completeness of any prime-seven classification is a separate dependency and is not asserted here.

## 1. Physical setting

The fixed outgoing quotient rows, in the indicated order of physical child sets, are

- Q0 = (8,5,17,70,43,15,55),
- Q1 = (14,28,56,112,97,67,7),
- Q2 = (14,52,88,98,41,69,19),
- Q3 = (66,76,25,49,35,7,60).

Here Q1 is the cyclic carousel: i points to i+1,i+2,i+3 modulo seven. For each i take an arbitrary tournament A_i of even order at least four. Let T=Q[A_i] and U=Q^op[B_i], where B_i uses exactly the same physical names as A_i and has the same complete ordered-path-cover profile. The special case B_i=A_i is allowed.

Actual ports x,y have the fixed edge x->y. Query masks a,c mean that every physical name in child i points to x or y exactly when bit i of a or c is one. Physical child names are never identified or independently relabelled between T and U. Let q_h count ordered spanning covers by h nonempty directed paths; q_1 is the Hamilton count.

For every fixed ambient root substitution, there is a child-uniform physical query that distinguishes T from U. Both ports receive precisely the original core root's incidences to every exterior name. At most two fixed core queries are needed for each listed quotient. These assertions concern complete augmented-tournament responses, including all cubic and quartic terms.

## 2. A general affine reduction for every ordered-cover coordinate

The complete normalized modulo-sixteen upper domain for an even child of order at least four is inherited from [S188 even cyclic profiles](https://github.com/eliottcassidy2000/math/blob/8f14b96ab020f7d5ceeb88e7bda413bb6f8bc9b7/05-knowledge/results/session188_structural_bridges_20261104/even_cyclic_mod16/PROOF.md):

    q(A_i)/H(A_i) = (0,1,14+8e_i,4+8f_i,8,8g_i,0,...),
    e_i,f_i,g_i in {0,1}, modulo sixteen.

Every H(A_i) is odd. The domain includes every actual profile; its formal vectors need not all be realized. We retain the common odd normalization scalar Lambda=product_i H(A_i).

For any quotient Q, query z, visit vector k, and h>=1, let K_h(Q,z;k) count cut-aware run words: each child i occurs k_i times, each port occurs once, and h-1 cuts split the word into h nonempty paths. An uncut adjacent pair must follow an actual arc, and cannot be the same child twice. A cut allows any adjacent symbols, including equal child symbols. Child runs across all path slots are filled by one ordered child cover. Therefore the normalized augmented q_h is exactly

    Phi_h(Q,z;r) = sum_k K_h(Q,z;k) product_i r_(i,k_i) modulo sixteen.

Use the actual TT4 baseline r=(0,1,14,4,8,0,...). Expanding around it gives

    Phi_h = B_h + 8 sum_i [e_i D_(h,i,2)+f_i D_(h,i,3)+g_i D_(h,i,5)],
    D_(h,i,k) = K_h(Q,z;1,...,1,k,1,...,1) modulo two.       (1)

This holds simultaneously for all ordered-cover coordinates. A product of two perturbations has factor64 and vanishes modulo sixteen. In a first-order term, every other child must occur once, because all its higher baseline coordinates are even. No return-segment involution across path cuts is assumed. This is why all three kinds of derivatives are retained, even when some later contrast entries vanish.

The baseline only needs visit vectors with sum_i(k_i-1)<=3, with k_i in1..4. Indeed baseline valuations at k=1,2,3,4 are0,1,2,3. There are120 such vectors at seven children. The seven single-child k=5 vectors supply the remaining derivative terms. Consequently127 small-model shapes determine (1), for arbitrary child orders.

This affine all-cover reduction is valid for arbitrary quotient order. Nonvanishing below is the finite certificate specific to the four listed quotients.

## 3. Explicit full-profile identities

Write Delta q_h(z)=q_h(T+x+y)-q_h(U+x+y), with the same physical query. The complete exact residue bank and the independent physical interpolation give the following identities modulo sixteen:

- Q0, z=(1,0): Delta q(z)=(8,0,0,0,0,...).
- Q3, z=(3,1): Delta q(z)=(8,0,0,0,0,...).
- Q1, z0=(0,3): Delta q(z0)=(8(e_3+e_5),0,0,0,0,...).
- Q1, z1=(4,68): Delta q(z1)=(8(1+e_3+e_5),0,0,0,0,...).
- Q2, z0=(0,3): Delta q(z0)=(8(e_0+e_1),0,0,0,0,...).
- Q2, z1=(4,68): Delta q(z1)=(8(1+e_0+e_1),0,0,0,0,...).

For h=1,...,5 these are full affine identities, including all21 independent formal profile bits. The common odd Lambda fixes zero and eight, so the displayed formulas already describe the actual physical contrasts. For h>=6, sixteen divides h!, hence each ordered-cover coordinate separately vanishes modulo sixteen. This establishes the complete infinite profile identities, rather than only their first coordinate.

In particular exactly one of the two displayed queries for Q1 or Q2 has Hamilton contrast eight. All higher profile contrasts vanish at both queries. The pair works without knowing which of the two child-profile bits is present. The carousel-seven all-E branch is therefore resolved by two fixed actual queries, uniformly over all child sizes and orientations.

## 4. Every fixed ambient and root-containing original-module families

Substitute the core for a fixed root of any tournament R, and lift both ports with that root's exterior incidences. The augmented core remains a common module. The usual run substitution writes the ambient Hamilton contrast as sum_h c_h Delta q_h, with c_1=H(R) odd. At a selected residue-eight query this is8H(R)=8 modulo sixteen. Exterior size and orientation are unrestricted.

This lift uses the entire profile identities in Section3. A general nonzero modulo-sixteen Hamilton difference alone would not imply it.

A source-relative root-family consequence follows from the separately audited S188 hereditary original-module monotonicity and converse-complement duality. For each query z in the fixed one- or two-query baseline battery, include Jz=(127-c,127-a) inside the core, leaving exterior port incidences fixed. Complete baseline profile duality makes the two ambient baseline contrasts opposite over the integers. At least one original query has nonzero contrast, so one member of its dual pair has negative source-minus-baseline contrast. For an original-module XOR family containing this prime core root oddly, every proper selected original module lies inside a child, and each target child profile is coordinatewise no smaller. The chosen uniform query has nonnegative run-substitution coefficients, so the strict negative contrast persists. The resulting battery has at most four queries.

This is source-relative detection for families contained in the core, not all-pairs injectivity among arbitrary orbit outputs. Exterior operations require a separate highest-node restoration argument. No claim of unrestricted named-response reconstruction or an exact collision is made.

## 5. Full data, independent verification, and reproduction

`even_bank.cpp` counts the first five ordered-cover coordinates with a cut-aware multiset automaton. It retains both orientations of all four quotients at every one of the841 Boolean query basis points with at most two a bits and two c bits. Each point includes five TT4 baseline residues and all105 derivative parities (three child multiplicities, seven labels, five cover coordinates). `even_bank.txt` is the complete bank; `coefficients.json` retains its complete Boolean coefficient bank, including zero, cubic and quartic coordinates. Every channel has bidegree at most(2,2), because each actual port occurs exactly once and contributes at most its incoming and outgoing incidence. This supplies a lossless representation of every full block-query vector, not an endpoint-only approximation.

For independent verification, `physical.cpp` uses literal named vertex-subset dynamic programming. It counts complete ordered covers on actual tournament extensions, without any quotient-run compression. `independent.py` evaluates1,524 physical instances, on at most13 vertices, and uses exact triangular inversion of the TT_n profiles q_k=k!S(n,k) to recover all127 run shapes at each of the six selected queries and both orientations. Every recovered integer is nonnegative and every inversion division is checked. It reproduces every first-five baseline coordinate and all21 profile derivatives at each selected query. Complete literal values and recovered coefficients are retained in `independent.json`.

No actual28-old-name Hamilton enumeration is claimed: exact small-model interpolation is the independently justified route to the arbitrary-size coefficient functional. The graph-order bound13 also keeps all integer physical counts below13! binomial(12,4)<2^64.

Run `python3 verify.py > normal.json` and `python3 -O verify.py > optimized.json`; their outputs must agree. The standard-library programs use an already installed C++17 compiler and perform no installation. Checks remain active under optimization. Frozen allowlist: PROOF.md, rows.json, verify.py, even_bank.cpp, physical.cpp, independent.py, even_queries.json, even_bank.txt, coefficients.json, independent.json, normal.json, optimized.json. Exploratory endpoint files are not dependencies.
