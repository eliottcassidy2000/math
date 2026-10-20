# Complete mixed-prime symmetry census through ten vertices

**Status:** PROVED: the elementary completeness reduction and least-degree consequence have passed fresh proof audit. Counts are FINITE-EXACT and have passed a newly written independent full-isomorphism and full-automorphism-generator implementation, in addition to the original Burnside replay. The prime-degree classification also has the elementary proof below. No unproved hypothesis or computer-algebra database is used. This staged report does not itself assign a new canon theorem identifier.

**Universe:** every subgroup H of S_n of odd order divisible by at least two distinct primes, for every 1 ≤ n ≤ 10, up to S_n-conjugacy. Fixed vertices are part of the permutation representation. Different embeddings of the same abstract group are not identified unless conjugate in S_n.

This is a complete bounded census, not the earlier selected-group list. It yields ten degree/class pairs: none through six vertices; one at seven; two at eight; two at nine; five at ten.

## 1. What is counted

A graph is **sign-even** when every automorphism has edge-reversal sign +1, with sign computed as the parity of edges whose endpoints are reversed relative to an ordering. This is not the Eulerian-degree condition.

For N=N_(S_n)(H), write:

- tau_T(H): full S_n-isomorphism types of tournaments admitting a conjugate of H;
- tau_E(H): full S_n-isomorphism types of sign-even graphs admitting a conjugate of H;
- R(H): N-orbits of H-fixed tournaments, equivalently N-orbits of relatively sign-even H-fixed graphs;
- m_T and m_E: numbers of N-orbits inside each full tournament or full sign-even graph type;
- Z: number of relatively-even N-orbits whose graphs are globally odd.

The histogram notation `1:31` means 31 full types, each with one relative orbit. The unrestricted multiplicity m_x(H) is the number of Aut(x)-conjugacy classes of embedded S_n-conjugates of H in Aut(x). For globally odd graphs only the **relatively-even** embedding classes contribute to Z; the table does not call that restricted number the unrestricted m_x(H).

Set F_T=sum_T(m_T-1) and F_E=sum_E(m_E-1). The exact defect balance is

    tau_T - tau_E = F_E + Z - F_T.

This is the inherited relative Walsh identity with the two independent losses made explicit: forgetting H can fuse presentations, and forgetting N(H) can expose a new odd automorphism.

## 2. Complete table

Let F21=C7⋊C3, with x↦x+1 and x↦2x on F7. “Fixed” means its usual seven-point action plus n−7 fixed vertices. The diagonal F21 acts on the last three vertices through its quotient F21/C7≅C3. C15 is the independent disjoint 5-cycle and 3-cycle action.

| n | H and vertex orbits | order H | order N(H) | pair orbits d | fixed objects each 2^d | R | tau_T | tau_E | m_T histogram | m_E histogram | Z |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---|---|---:|
| 7 | F21, 7 | 21 | 42 | 1 | 2 | 1 | 1 | 1 | 1:1 | 1:1 | 0 |
| 8 | F21, 7+1 | 21 | 42 | 2 | 4 | 2 | 2 | 2 | 1:2 | 1:2 | 0 |
| 8 | C15, 5+3 | 15 | 120 | 4 | 16 | 2 | 2 | 2 | 1:2 | 1:2 | 0 |
| 9 | F21, 7+1+1 | 21 | 84 | 4 | 16 | 4 | 4 | 4 | 1:4 | 1:4 | 0 |
| 9 | C15, 5+3+1 | 15 | 120 | 6 | 64 | 8 | 8 | 8 | 1:8 | 1:8 | 0 |
| 10 | F21, 7+1+1+1 | 21 | 252 | 7 | 128 | 12 | 12 | 12 | 1:12 | 1:12 | 0 |
| 10 | C15, 5+3+1+1 | 15 | 240 | 9 | 512 | 32 | 32 | 31 | 1:32 | 1:31 | 1 |
| 10 | C21=C7×C3, 7+3 | 21 | 252 | 5 | 32 | 4 | 4 | 4 | 1:4 | 1:4 | 0 |
| 10 | diagonal F21, 7+3 | 21 | 126 | 3 | 8 | 4 | 2 | 2 | 2:2 | 1:2 | 2 |
| 10 | F21×C3, 7+3 | 63 | 252 | 3 | 8 | 2 | 2 | 2 | 1:2 | 1:2 | 0 |

Only the displayed C15 action fails individual equality. It is also the only class not covered by the inherited odd-pronormal-base/one-quotient-prime theorem.

## 3. Elementary completeness proof

### 3.1 Odd transitive groups and long prime cycles

Every odd-order permutation group admits an invariant tournament: it cannot interchange the endpoints of an unordered pair, since that would give an element of even order. Choose one orientation on each pair orbit and propagate it.

If the group is transitive of odd degree d, every vertex in such an invariant tournament has the same outdegree, necessarily (d−1)/2.

**Long-cycle lemma.** A transitive odd-order group of degree d cannot contain a p-cycle with d/2<p<d.

Indeed, a vertex outside that cycle is fixed by it. Its out-neighborhood either contains all p cycle vertices or none. Its outdegree is therefore at least p or at most d−p−1. Neither range contains (d−1)/2. This is an elementary group obstruction extracted from the invariant tournament.

For degree nine, Cauchy's theorem would turn any factor 5 or 7 in the group order into a single 5- or 7-cycle. The lemma excludes both. The only remaining odd prime divisor of 9! is 3, so every odd transitive subgroup of S9 is a 3-group. We do not need to enumerate the various transitive 3-groups.

### 3.2 Prime degrees three, five, seven

An odd transitive subgroup at prime degree p contains a p-cycle. After conjugating that cycle to translation, an invariant tournament is a circulant: its connection set S chooses one element of each pair {a,−a} in F_p^×.

At p=3 the image is C3.

At p=5 all four connection sets are multiplier images of {1,2}. The out-neighbors and in-neighbors of 0 each induce a directed edge, so an automorphism fixing 0 fixes all vertices. Translation is transitive, hence the full tournament automorphism group is C5. Every odd transitive image of degree five is therefore C5.

At p=7 the eight connection sets split under multiplication into the six images of {1,2,3} and the two images of {1,2,4}. In the first case the three out-neighbors and the three in-neighbors of 0 induce transitive tournaments, so the vertex stabilizer is trivial and the full group is C7. In the second case the out-neighbors form a directed triangle, so their action has order at most three. If those out-neighbors are fixed pointwise, their adjacency signatures distinguish all three in-neighbors, making the action faithful. Multipliers 1,2,4 realize the three rotations. Thus the full automorphism group is exactly F21. Any transitive subgroup of it containing C7 is C7 or F21.

The verifier separately exhausts all four and eight cyclic tournaments. It obtains automorphism-order histograms {5:4} and {7:6,21:2}, verifies every automorphism is affine, and checks the refined backtracking result against all p! vertex permutations. This is a positive control for the elementary argument, not an assumption that the subgroup list is complete.

### 3.3 Assemble the vertex orbits

All orbits of an odd group have odd size. Through degree ten the possible nontrivial sizes are 3,5,7,9, with transitive images as just proved. The group injects into the product of its orbit images.

- If a nine-orbit occurs, only one vertex remains. The whole group is a 3-group.
- With no seven-orbit, a mixed-prime group needs a five-orbit and a three-orbit. Their images are C5 and C3. A subdirect subgroup of this coprime product is the whole C5×C3=C15. Up to two fixed vertices may be added. Two five-orbits instead give a 5-group; any collection of three-orbits gives a 3-group.
- A seven-orbit leaves at most three vertices. If they are all fixed, the sole mixed-prime image is F21.
- With orbit sizes 7+3, the image on seven points is C7 or F21 and that on three is C3. A subdirect subgroup of C7×C3 is the whole C21. For F21×C3, the kernel of projection onto F21 is either the whole second C3, giving F21×C3, or trivial. In the latter case the subgroup is the graph of a surjection F21→C3. Since F21/[F21,F21]≅C3, there are two such surjections. A transposition on the three-point orbit interchanges them, so they give one ambient-conjugacy class, the diagonal F21.

This proves exhaustiveness. Distinct rows at a fixed n are distinguished by order, vertex-orbit sizes, abelianness, or projection onto the seven-point orbit. Thus there are no duplicates.

## 4. Least possible failure degree

**Corollary.** Ten is the least n for which there exists an odd-order subgroup H≤S_n with tau_T(H)≠tau_E(H). At n=10 exactly one ambient-conjugacy class of odd-order subgroups fails: C15 with vertex-orbit sizes 5+3+1+1.

For p-groups, use the inherited all-p-group theorem in

`05-knowledge/results/prime_local_walsh_bijection_20261020/SYLOW_ACTION_ADDENDUM.md`.

The trivial group is covered by the same statement (or the original tournament/sign-even-graph equinumerosity). Every remaining odd subgroup through degree ten belongs to the complete table. The unique failure is already proved independently in

`05-knowledge/results/c15_symmetry_counterexample_20261020/C15_ORDER_TEN_COUNTEREXAMPLE.md`.

The new content is the complete reduction, which removes the earlier qualification that degree ten was only the first failure found after a selected-group census. No least-degree statement beyond this scope is inferred from a selected list.

## 5. Pronormal-base coverage, including the nonpronormal survivor

The inherited theorem is

`05-knowledge/results/pronormal_quotient_bijections_20261020/PRONORMAL_QUOTIENT_THEOREM.md`, Theorems 1 and 3 and Lemma 4.

It applies when K is odd and pronormal in S_n, K≤H≤N(K), and H/K is a q-group for one odd prime q.

| class | base K | quotient prime | H itself pronormal? | explanation |
|---|---|---:|---|---|
| F21 with 0–3 fixed vertices | C7 | 3 | yes | C7 is ambient Sylow; multiplier C3 is normal in C6×S_(n−7) |
| C15 at n=8,9 | C5 | 3 | yes | C5 is ambient Sylow; the C3 image is Sylow in S3 or S4 |
| C21 on 7+3 | C7 | 3 | yes | root C3 is normal in C6×S3 |
| F21×C3 on 7+3 | C7 | 3 | yes | quotient C3² is normal in C6×S3 |
| diagonal F21 on 7+3 | C7 | 3 | no | diagonal quotient line is a 3-group but is not pronormal |
| C15 at n=10 | none | none | no | all possible nontrivial proper bases fail ambient pronormality |

For the last row, the only nontrivial proper subgroups are the displayed C3 and C5. Each is nonpronormal: conjugate its support to a disjoint support of the same size, obtaining different factors of an abelian direct product. That product cannot conjugate the factors. K=1 cannot give a prime-power quotient, and K=H is impossible because the sign-escape witness contradicts the pronormal descent lemma. Thus this particular theorem has no alternative base hidden inside H.

Eight of the ten rows have pronormal H. A ninth has a pronormal base but nonpronormal H. Only the tenth is outside the theorem, and it is the unique failure. This is a bounded observation, not a necessity theorem for larger degrees.

## 6. Two different mechanisms at degree ten

### 6.1 Uncompensated escape: C15

The unique extra relatively-even globally-odd graph is K_(5,5). Under the verifier's labeling, the five-cycle acts on {0,1,2,3,4}; the three-cycle acts on {5,6,7}; vertices 8 and 9 join the latter orbit's independent twin class. Swapping the two five-vertex parts is an odd automorphism and is outside N(H).

Every tournament and every globally even graph has relative multiplicity one, so F_T=F_E=0 and Z=1. The balance gives tau_T−tau_E=1. This reproduces the inherited structural obstruction, with a different convenient labeling.

### 6.2 Compensated escape: diagonal F21

Write A=F7 and B={7,8,9}. The H-fixed pair orbits are the pairs within A, the pairs within B, and the cross pairs. The H-fixed cube is therefore exactly the F21×C3-fixed cube, although the two groups are not conjugate and their normalizers have different orders.

The two tournament types are the two directed joins of a Paley seven-tournament and a directed triangle. Each splits into two H-normalizer orbits, so m_T=2 on both types and F_T=2. The two globally even graph types are the empty graph and K_(7,3), each with m_E=1, so F_E=0.

The remaining two relatively-even graph types are

- K3 disjoint union 7 isolated vertices;
- K3 joined to an independent seven-vertex set.

A transposition on the three-vertex set is an odd automorphism of each: it reverses the single edge between its two moved vertices, while all contributions involving the seven-set pair off. This transposition does not normalize the diagonal subgroup. Thus Z=2 and the defect balance gives 0=0+2−2.

This explicitly separates three claims:

1. H itself is pronormal: false here;
2. N(H)-evenness implies global evenness: false here;
3. tau_T(H)=tau_E(H): true here, with the inherited pronormal C7 base controlling the relevant quotient instead.

The strongest survivor is therefore not “no sign escape.” It is the actual weighted balance, or the separately sufficient descent theorem at a suitable base.

## 7. Reproduction and verification

Requirements: Python 3 standard library only. No GAP, nauty, external graph database, installation, or network is needed.

    python3 mixed_symmetry_census.py > census.json
    python3 -O mixed_symmetry_census.py > census.optimized.json
    cmp census.json census.optimized.json
    python3 independent_replay.py > independent_replay.json
    python3 independent_automorphism_audit.py > independent_automorphism_audit.json
    python3 -O independent_automorphism_audit.py > independent_automorphism_audit.optimized.json
    cmp independent_automorphism_audit.json independent_automorphism_audit.optimized.json

The main path enumerates every H-fixed graph and tournament from pair orbits, exactly enumerates N(H) from the intrinsic vertex-orbit partition, and partitions the fixed cube into normalizer orbits. It then uses full directed/undirected isomorphism backtracking to recover S_n-types, retaining the number of relative orbits inside each full type. Graph signs are tested by the intrinsic twin quotient: true twins give an explicit odd transposition; otherwise all independent-twin quotient automorphisms are tested and lifted. Permutations inside independent twin classes are sign-even, so this does not omit a sign source.

The independent path classifies **all** fixed labelled objects directly before normalizer reduction, including graphs that the main path excluded as relatively odd. It then computes ordinary and signed Burnside sums independently of orbit partitioning, compares the results, and repeats the entire normalizer/fusion/sign census after reversing all vertex labels. It is algorithmically independent in enumeration order and Burnside/orbit aggregation, but shares the full-isomorphism and twin-sign primitives; it is not claimed to be a wholly independent software implementation.

The main script also compares the small prime-degree automorphism search to brute-force S5 and S7 enumeration. All explicit checks use runtime exceptions, so optimized Python does not remove them. Exact masks, generators, full-type representatives, multiplicities, and escaping odd automorphisms are in the JSON files. `SHA256SUMS.txt` pins the staged artifacts.

The fresh `independent_automorphism_audit.py` imports neither earlier verifier.
It reconstructs every fixed object directly from the certificate generators,
uses its own constrained adjacency search for full isomorphism, and computes
full graph automorphism generators by a point-stabilizer chain. At each level
it finds one transporter for every attainable image of the next vertex while
fixing the prefix. These are coset representatives of the next point stabilizer;
induction proves their union generates the entire automorphism group. Testing
the edge-sign character on these generators is therefore complete, including
graphs with automorphism group S10. It does not use twin reduction or presume
that an automorphism normalizes H. Its 14,122 runtime checks pass in normal and
optimized Python, and all ten rows agree exactly. The original main and
Burnside paths retain respectively 306 and 621 checks.

See `AUDIT_AND_NEXT_TARGET.md` for the precise copy-transitivity descent
criterion, the signed copy-orbit formula, and the remaining mixed-prime target.

## 8. Inheritance, concept board, and boundaries

Inheritance: the closest proved mechanism is relative fixed-cube Walsh duality plus the pronormal-base/quotient theorem; the canonical hostile is C15 at n=10; the corrected near miss is mistaking the older selected subgroup list for a complete census; the least-used sidecar is the multiplicity of an H-presentation inside a full isomorphism type.

Concept board used here:

1. Orbit-size decomposition: turns subgroup classification into a bounded list of transitive images.
2. Regular invariant tournament: excludes long prime cycles in transitive odd actions.
3. Pronormal base: controls fusion and sign only at the chosen base, not automatically at H.
4. Relative multiplicity: measures the quotient loss from N(H) to S_n.
5. Sign escape: can be uncompensated (C15) or exactly compensated (diagonal F21).

The map is from H-fixed labelled objects to normalizer orbits and then full types. It preserves the ambient symmetry predicate, loses the embedded copy of H, and requires relative multiplicity and sign escape as sidecars. The decisive hostile is the diagonal row: equality survives while both individual descent properties fail.

The census does not classify odd subgroups beyond n=10, characterize every equality case at larger degree, or produce a single bijection preserving all odd-prime predicates simultaneously. The inherited five-vertex joint-prime obstruction remains unchanged. No literature-priority claim is made; the completeness proof and finite program are supplied explicitly.
