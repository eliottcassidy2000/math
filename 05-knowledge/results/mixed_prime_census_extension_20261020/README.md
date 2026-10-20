# Odd mixed-prime symmetry classes at eleven and twelve vertices

**Status:** Structural classification and compensated-escape claims **PROVED**, with independent proof review. The 18-row numerical census is **FINITE-EXACT**; all rows passed an uncompressed pair-action and original unbucketed-isomorphism replay. See [verification scope](VERIFICATION.md). This is not an extrapolation of the earlier selected-subgroup experiments.

**Universe:** odd-order subgroups H≤S_n whose order has at least two distinct prime divisors, up to ambient S_n-conjugacy, for n=11 and n=12. Fixed vertices and the actual permutation embedding are retained.

The complete classification has **8 classes at n=11 and 10 at n=12**. Exactly one class at each degree fails tournament/sign-even-graph equinumerosity: the cyclic C15 action with one five-orbit, one three-orbit, and all remaining vertices fixed. Its counts are 172 versus 168 at n=11 and 1376 versus 1360 at n=12. These finite counts are separate from the structural classification proof.

The principal new mechanism is an infinite family with **exactly compensated sign escape**: a diagonal Cp⋊C3 and the independent product (Cp⋊C3)×C3 have the same fixed cube, but different normalizers. This forces tournament presentation multiplicity two and one escaping relative graph orbit per surviving type.

## 1. Dependencies, definitions, and boundaries

The predecessor is [Complete census through ten](../odd_mixed_prime_census_20261020/README.md), including its elementary degree-3,5,7,9 transitive-image classification and the defect bookkeeping. Its proof does not use a subgroup database, Burnside's prime-degree theorem, or the odd-order theorem.

The proved descent mechanism is [Pronormal bases and quotient refinement](../pronormal_quotient_bijections_20261020/PRONORMAL_QUOTIENT_THEOREM.md), especially Lemmas 2 and 4 and Theorem 3. The earlier C15 failure is [C15 at ten vertices](../c15_symmetry_counterexample_20261020/C15_ORDER_TEN_COUNTEREXAMPLE.md).

For H≤S_n and N=N_(S_n)(H), use:

- tau_T, tau_E: full S_n-isomorphism types of tournaments and sign-even graphs admitting an ambient conjugate of H;
- R: N-orbits of H-fixed tournaments, equal by the relative Walsh theorem to N-orbits of relatively sign-even H-fixed graphs;
- F_T, F_E: sum of m−1 over full tournament or full globally sign-even graph types, where m is the number of their N-orbits;
- Z: number of relatively sign-even N-orbits whose graphs are globally odd.

Then tau_T−tau_E=F_E+Z−F_T. “Sign-even” means trivial edge-reversal character on every automorphism, not the Eulerian-degree condition. Z counts relative orbits, not necessarily distinct global graph types.

For globally odd graphs, the number entering Z is a restricted relative-even embedding multiplicity. It is not silently identified with the unrestricted number of all embedded H-classes.

## 2. An elementary degree-eleven transitive-image theorem

**Claim.** An odd transitive subgroup of S11 is conjugate to C11 or F55=C11⋊C5 in its usual affine action.

Every odd permutation group has an invariant tournament: no group element can exchange the two points of a pair, so orient one pair in each orbit and transport the orientation. A transitive degree-11 subgroup contains an 11-cycle by Cauchy's theorem. Label this cycle as translation on F11. An invariant tournament is therefore T(S), with u→v iff v−u∈S, where S contains one member of each pair {a,−a}. There are 2^5=32 such S.

Multiplication partitions these sets into the following four orbits:

| representative S | orbit size | score sequence inside the five out-neighbors of 0 | full automorphism order |
|---|---:|---|---:|
| {1,2,3,4,5} | 10 | 0,1,2,3,4 | 11 |
| {1,2,3,4,6} | 10 | 1,1,2,3,3 | 11 |
| {1,2,3,5,7} | 10 | 1,2,2,2,3 | 11 |
| {1,3,4,5,9} | 2 | 2,2,2,2,2 | 55 |

Here is a proof of the orbit assertion, not an appeal to an enumeration. A nonidentity multiplier stabilizing a skew S has odd order, hence order five. Its two nonzero orbits are the quadratic residues and nonresidues. Thus the first three sets have multiplier stabilizer one, while the final set has stabilizer five. Their displayed neighborhood score sequences distinguish the four multiplier orbits. Their sizes sum to 32, exhausting all skew S.

It remains to bound the full vertex stabilizer, without presuming all automorphisms are affine. An automorphism fixing 0 preserves its out-neighbor set S. In the first row that subtournament is transitive. In the second row its two score-3 vertices and two score-1 vertices each induce a directed edge, so each is fixed, as is the unique score-2 vertex. In the third row vertices 2 and 3 have unique scores three and one; the remaining vertices 1,5,7 are distinguished by their adjacency to (2,3), respectively (1,1),(0,0),(0,1). Thus the first three out-neighbor subtournaments are rigid.

For the fourth row the out-neighbor tournament is regular on five vertices. An automorphism of any regular five-tournament fixing one vertex fixes its two out-neighbors and its two in-neighbors, since each pair induces a directed edge. Its automorphism group therefore has order at most five. The five quadratic-residue multipliers realize all five possibilities here.

In every row, an automorphism fixing 0 and all out-neighbors fixes all in-neighbors, because their adjacency signatures to the ordered out-neighbors are distinct:

| S in increasing order | signatures of the increasing in-neighbors |
|---|---|
| 1,2,3,4,5 | 00000, 10000, 11000, 11100, 11110 |
| 1,2,3,4,6 | 00001, 01000, 10100, 11010, 11100 |
| 1,2,3,5,7 | 00011, 01001, 01000, 10110, 11000 |
| 1,3,4,5,9 | 01010, 00101, 10010, 10001, 01100 |

Translation now gives the full group orders 11,11,11,55. The displayed affine groups realize those orders, so these are exactly the automorphism groups. A transitive subgroup of C11 or F55 has order divisible by 11 and hence is C11 or F55. All order-11 subgroups of F55 are its unique Sylow-11 subgroup. This proves the claim.

`prime11_certificate.py` checks all 32 connection sets independently. For each one it enumerates only the 5! permutations of the out-neighbor set, extends them uniquely through these signatures, and verifies the full adjacency matrix. It obtains 30 automorphism groups of order 11 and two of order 55. No full S11 search is needed or claimed.

## 3. Complete subgroup assembly

All vertex-orbit sizes of an odd group are odd. Through degree twelve the nontrivial possibilities are 3,5,7,9,11. The predecessor proves:

- degree three has image C3;
- degree five has image C5;
- degree seven has image C7 or F21=C7⋊C3;
- every odd transitive degree-nine image is a 3-group.

The last fact follows from the invariant-tournament long-cycle lemma: in an odd transitive degree d action, a single p-cycle with d/2<p<d is impossible, since a vertex outside its support has outdegree either at least p or at most d−p−1, never (d−1)/2. At degree nine this excludes the prime divisors 5 and 7.

The group embeds as a subdirect subgroup of the product of its orbit images. The following cases exhaust the possibilities.

1. **An eleven-orbit.** At most one fixed vertex remains. The only mixed-prime group is F55.
2. **A nine-orbit.** At most one three-orbit remains, and every image is a 3-group. There is no mixed-prime group.
3. **A seven-orbit.** The remaining nontrivial orbit is absent, one three-orbit, or, only at n=12, one five-orbit. With no other nontrivial orbit one gets F21. With a three-orbit one gets C21, diagonal F21, or F21×C3, exactly as in the predecessor. With a five-orbit, coprime image orders force C35 or F21×C5.
4. **No orbit larger than five.** Two five-orbits cannot be accompanied by a three-orbit before n=13, so such groups are 5-groups. Mixed primes therefore require exactly one five-orbit and either one or two three-orbits. One three-orbit gives C15. With two three-orbits the 3-image is a subdirect subgroup of C3×C3: it is either a diagonal line or the full product. All nonzero-slope diagonal lines are conjugate by independently reversing a triple. Coprime order with C5 then gives C5×C3 acting on 5+3+3, or C5×C3² on the same orbit partition.

The elementary coprime-product step can be seen from orders: a subgroup surjecting onto groups of coprime orders a and b has order divisible by both, hence by ab, so must equal their product. For a subgroup inside C5×(C3×C3), first project onto its actual 3-image, then apply this argument.

The diagonal F21 is the graph of a surjection F21→C3; its two choices are conjugate by a triple reversal. Distinct rows below are distinguished by group order, vertex-orbit partition, or the order of the seven-orbit image. Therefore the list has no duplicates as well as no omissions.

## 4. Complete finite table

All unspecified vertices are fixed. “Diagonal C15” means a five-cycle independent of one synchronized pair of three-cycles; it is not the single-three-orbit C15 action. The order-45 group is C5×C3², **not cyclic C45**.

| n | H | order H | order N(H) | pair orbits d | R | tau_T | tau_E | F_T | F_E | Z |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 11 | F21 on 7 | 21 | 1008 | 11 | 48 | 48 | 48 | 0 | 0 | 0 |
| 11 | C15 on 5+3 | 15 | 720 | 13 | 176 | 172 | 168 | 4 | 4 | 4 |
| 11 | C21 on 7+3 | 21 | 252 | 7 | 16 | 16 | 16 | 0 | 0 | 0 |
| 11 | diagonal F21 on 7+3 | 21 | 126 | 5 | 16 | 8 | 8 | 8 | 0 | 8 |
| 11 | F21×C3 on 7+3 | 63 | 252 | 5 | 8 | 8 | 8 | 0 | 0 | 0 |
| 11 | F55 on 11 | 55 | 110 | 1 | 1 | 1 | 1 | 0 | 0 | 0 |
| 11 | diagonal C15 on 5+3+3 | 15 | 720 | 9 | 16 | 12 | 12 | 4 | 0 | 4 |
| 11 | C5×C3² on 5+3+3 | 45 | 1440 | 7 | 4 | 4 | 4 | 0 | 0 | 0 |
| 12 | F21 on 7 | 21 | 5040 | 16 | 296 | 296 | 296 | 0 | 0 | 0 |
| 12 | C15 on 5+3 | 15 | 2880 | 18 | 1408 | 1376 | 1360 | 32 | 34 | 14 |
| 12 | C21 on 7+3 | 21 | 504 | 10 | 64 | 64 | 64 | 0 | 0 | 0 |
| 12 | diagonal F21 on 7+3 | 21 | 252 | 8 | 64 | 32 | 32 | 32 | 0 | 32 |
| 12 | F21×C3 on 7+3 | 63 | 504 | 8 | 32 | 32 | 32 | 0 | 0 | 0 |
| 12 | F55 on 11 | 55 | 110 | 2 | 2 | 2 | 2 | 0 | 0 | 0 |
| 12 | diagonal C15 on 5+3+3 | 15 | 720 | 12 | 128 | 96 | 96 | 32 | 0 | 32 |
| 12 | C5×C3² on 5+3+3 | 45 | 1440 | 10 | 32 | 32 | 32 | 0 | 0 | 0 |
| 12 | C35 on 7+5 | 35 | 840 | 6 | 4 | 4 | 4 | 0 | 0 | 0 |
| 12 | F21×C5 on 7+5 | 105 | 840 | 4 | 2 | 2 | 2 | 0 | 0 | 0 |

Each fixed graph/tournament cube has 2^d labelled objects. The largest is dimension 18, not all labelled graphs on twelve vertices.

At n=11 the single-three-orbit C15 row has multiplicity histograms m_T=(1:168,2:4), m_E=(1:164,2:4), with four escaping types, all of multiplicity one. At n=12 these are m_T=(1:1344,2:32), m_E=(1:1326,2:34), and fourteen escaping types, all of multiplicity one. Thus the degree-twelve deficit 16 is **not** merely the number 14 of escaping types: target fusion contributes two more than source fusion. The exact balance is 16=34+14−32.

Combining the predecessor, the classification above, and the inherited all-p-group theorem, the only odd-subgroup failure classes through degree twelve are C15 on 5+3 with 2,3,4 fixed vertices, respectively. The least failure degree remains ten. No n≥13 classification is inferred.

## 5. Infinite compensated-escape family

Let p≡1 mod3 be prime, let a∈F_p^× have order three, and let f∈{0,1,2}. On A=F_p, B=F3, and f individually fixed vertices, define

- t: x↦x+1 on A;
- a: x↦ax on A;
- b: y↦y+1 on B;
- D=〈t,ab〉, the diagonal Cp⋊C3;
- P=〈t,a,b〉=(Cp⋊C3)×C3.

Write a_p for the number of affine-multiplier orbits of Cp⋊C3-invariant tournaments on A. This is a precisely defined finite number, not an assumption of uniqueness. In particular a_7=1.

**Family claim.** The groups D and P have the same full symmetry predicates on graphs and tournaments. Put c_f=2,8,32 for f=0,1,2. Then

- tau_T(D)=tau_E(D)=tau_T(P)=tau_E(P)=c_f a_p;
- every tournament type for D has m_T=2;
- every globally sign-even graph type for D has m_E=1;
- F_T(D)=Z(D)=c_f a_p and F_E(D)=0;
- P has no fusion or sign escape: F_T(P)=F_E(P)=Z(P)=0.

### 5.1 Same cube, different normalizers

The pair orbits for D and P coincide. On A, both project to Cp⋊C3; on B, both project to C3. The cross pairs A×B form one orbit because t moves freely on A while the quotient moves transitively on B. Pairs involving fixed vertices also have the same orbits. Odd-order groups never reverse a pair, so equality of pair orbits gives equality of the invariant tournament cube as well as the graph cube.

The local normalizer of Cp⋊C3 on A is AGL(1,p): Cp is its unique Sylow-p subgroup, hence characteristic, and every affine map normalizes Cp⋊C3. The local normalizer of C3 on B is S3. Conjugation by AGL(1,p) induces the identity on the quotient (Cp⋊C3)/Cp≅C3, since the multiplier group is abelian. Therefore a normalizer of D must centralize the action of C3 on B. Consequently

    N(D)=AGL(1,p)×C3×S_f,
    N(P)=AGL(1,p)×S3×S_f.

Their index is two, represented by any transposition in B. Fixed and nontrivial vertex orbits cannot be exchanged, so no omitted block permutation enlarges these groups.

### 5.2 Why the product has full descent

The subgroup Cp is an ambient Sylow-p subgroup because n=p+3+f<2p. Hence it is pronormal. In its normalizer quotient C_(p−1)×S_(3+f), P/Cp is the product of a normal multiplier C3 and a Sylow-3 subgroup of S_(3+f). It is pronormal: the normal factor stays fixed under conjugation, and Sylow conjugacy supplies a conjugator in the generated group for the other factor. The inherited pronormal-tower lemma makes P pronormal in S_n. Its normalizer orbits therefore coincide with full types, and relative evenness coincides with global evenness.

### 5.3 Exact compensation

Every tournament automorphism has odd order: an even-order permutation has an involution power, and a nonidentity involution exchanges the endpoints of at least one pair. Thus a tournament stabilizer in N(P) has trivial image in N(P)/N(D)≅C2. Every N(P)-tournament orbit splits into exactly two N(D)-orbits. Pronormal descent for P gives m_T(D)=2.

A P-invariant graph induces either an independent triple or K3 on B. All three vertices have identical external neighbors. When B is independent, a transposition in B is an even automorphism, and the N(P)-orbit remains one N(D)-orbit. Relative evenness for N(D) and N(P) is the same in this case, because that even transposition represents the missing coset. By P-descent, such a relatively-even graph is globally even.

Changing only the three internal edges on B to K3 preserves every N(D)-sign, since N(D) induces only C3 on B. It creates an odd automorphism, namely a transposition in that true-twin triple. Conversely every N(D)-relatively-even graph with K3 on B arises from this toggle. The toggle commutes with the normalizer action and gives exactly one extra relatively-even N(D)-orbit for each globally even type. It follows that m_E(D)=1 and Z(D)=tau_E(P). There is no assertion that all these extra relative orbits necessarily give different full odd graph types in the entire infinite family.

### 5.4 Counting the surviving types

The affine action on A supplies a_p internal tournament types, and S3 makes the internal triangle unique. Outside these blocks, a tournament is a binary orientation for A versus B, two binary incidences for each fixed vertex, and a tournament on the f fixed vertices. For f=0 and f=1 this gives 2 and 8. For f=2 there are 2^(1+4+1)=64 choices before exchanging the two fixed vertices. That exchange fixes no tournament because it reverses their mutual arc, hence gives 32 types. These choices are independent of the internal types. Pronormality prevents any additional full-isomorphism fusion.

For p=7 the formula gives 2,8,32 at n=10,11,12 and exactly explains every entry of the diagonal F21 rows. It is stronger than observing numerical equality: it identifies both losses and their exact compensation.

### 5.5 The exact relative law survives without the small-f restriction

For every f≥0, the same cube and normalizer arguments remain valid, even when P is no longer pronormal. Consequently the full symmetry predicates still agree and, for each admitted full type,

    m_T(D)=2 m_T(P),       m_E(D)=m_E(P) on globally even graphs.

The independent-triple/K3 toggle still doubles the relatively-even orbit set. Therefore

    R(D)=2R(P),
    F_T(D)=2F_T(P)+tau_T(P),
    F_E(D)=F_E(P),
    Z(D)=R(P)+Z(P).

These identities imply the same final defect for D and P, as they must from their identical symmetry predicates. They do not force that common defect to vanish. The restriction f≤2 in the family claim supplies the additional pronormal descent that sets all three P-defects to zero; it is a sufficient scope, not a claimed sharp equality boundary.

The separate [family scope check](FAMILY_SCOPE.md) proves equality for 0≤f≤p−4 by the inherited quotient theorem and records four explicitly selected probes beyond the complete census. It keeps those probes separate from the 18-row classification.

## 6. Other structural count checks

F21 on seven vertices with f fixed vertices, for 0≤f≤5, has exactly b_f types on both sides, where b_f counts two-colored tournaments on f vertices. The two colors record the orientation between a fixed vertex and the seven-block. Cp=C7 is ambient Sylow and the multiplier C3 is normal in the normalizer quotient, so F21 is pronormal and there is no hidden fusion/sign escape.

If λ=(l_1,...,l_k) is an odd partition of f and z_λ=∏_j j^(m_j)m_j!, Burnside gives

    b_f = sum over odd λ of 2^[k + sum_i (l_i−1)/2 + sum_(i<j) gcd(l_i,l_j)] / z_λ.

An even cycle in a permutation would prohibit a fixed tournament; for all-odd cycles the exponent counts the k color orbits and the unordered-pair orbits. This gives b_0,...,b_5=1,2,4,12,48,296. Thus the two new F21-fixed counts are independently explained.

The same formula with 2^k replaced by c^k counts c-colored tournaments; call it b_f(c). The relative orbit count for C15 on 5+3 with f fixed vertices is exactly 2 b_f(4): each fixed vertex has four possible incidence colors to the two nontrivial blocks, their mutual orientation supplies a factor two, and each internal block has one normalizer-orbit of tournaments. This proves R=176 and R=1408 independently of cube enumeration, since b_3(4)=88 and b_4(4)=704. Likewise R(P)=2 a_p b_f(4) and R(D)=4 a_p b_f(4) for every f≥0; no full-isomorphism descent is implicit in these relative formulas. `relative_burnside_control.py` checks twelve rows by partition arithmetic alone, without importing the cube or isomorphism engine.

The seven-vertex C7 sector has two internal tournament types, while F21 has one and C5 has one. Hence the C21-on-7+3 counts are twice the product-family counts: 16 and 64. C35 on 7+5 has 2×1×2=4 types, while F21×C5 has 1×1×2=2. Both seven-plus-five groups are pronormal: start with ambient Sylow C7 and use respectively a Sylow C5, or normal C3 times Sylow C5, in C6×S5. F55 is likewise pronormal using ambient Sylow C11 and the normal multiplier C5, yielding one or two types with zero or one fixed vertex.

The equality of the two new 5+3+3 embeddings is presently recorded as FINITE-EXACT here, not falsely attributed to the ambient-Sylow sufficient theorem. A structural uniform transfer theorem for those embeddings remains a useful next target.

## 7. Computational certificate and replay boundaries

Requirements: Python 3 standard library. No GAP, nauty, subgroup database, package installation, or web source is used.

    python3 prime11_certificate.py > prime11_certificate.json
    python3 -O prime11_certificate.py > prime11_certificate.optimized.json
    cmp prime11_certificate.json prime11_certificate.optimized.json
    python3 extension_census.py > extension_counts.json
    python3 -O extension_census.py > extension_counts.optimized.json
    cmp extension_counts.json extension_counts.optimized.json
    python3 pair_action_replay.py > pair_action_replay.json
    python3 relative_burnside_control.py > relative_burnside_control.json

`fixed_cube_engine.py` is an unchanged, checksum-pinned copy of the predecessor's engine. `extension_census.py` replaces three parts:

1. It generates every normalizer using the proved local affine normalizers, all permutations of fixed vertices, and both exchanges of repeated three-orbits, then tests conjugation of every group generator. Projection to an intrinsic vertex orbit must normalize its orbit image, so these candidates exhaust the normalizer.
2. It represents each normalizer action by a permutation of pair-orbit coordinates, the affine tournament flip, and the graph sign functional. It discards only actions identical in all three coordinates, retaining their multiplicities for full stabilizer orders. This is exact signed-action compression, not a quotient that forgets signs.
3. It buckets potential full isomorphisms by a canonical equitable-refinement invariant, then applies the predecessor's unrestricted full-isomorphism search inside each bucket. Bucket equality is only a necessary test, never the final equivalence relation.

`pair_action_replay.py` recomputes all 18 rows with the predecessor's uncompressed vertex-pair actions and original unbucketed full-isomorphism classifier. It compares every count, all three multiplicity histograms, all fusion/escape defects, normalizer orders, and pair-orbit counts. This replay is independent of the compressed action, its deduplication, the coordinate enumeration order, and the invariant buckets; it shares normalizer generation, pair-orbit construction, full-isomorphism search, and twin-quotient global signs. It is not advertised as a wholly independent implementation.

The JSON retains generators, full-type representative masks, every multiplicity, and an explicit odd automorphism for each escaped odd type. Every check is a runtime exception, so Python -O cannot disable verification.

## 8. Research inheritance, map, and next stopping boundary

- Closest mechanism: relative Walsh balance plus pronormal-base descent.
- Canonical hostile: C15 on 5+3+1+1, with the K_(5,5) odd part swap.
- Corrected near miss: identifying Z with the entire deficit; n=12 gives 16≠14.
- Least-used sidecar: the index-two normalizer extension and the internal-triple toggle.

Concept board: orbit-image classification; elementary neighborhood rigidity; signed coordinate compression; normalizer index; source/target fusion; twin-sign escape. The subgroup list is assembled structurally before any cube is enumerated. The diagonal family then isolates a mechanism that the raw count alone would hide.

The map is H-fixed labelled object → N(H)-orbit → full S_n-type. It preserves the ambient symmetry predicate and loses the chosen embedded H. The restored sidecars are presentation multiplicity and escaped sign. In the infinite diagonal family, the normalizer index controls the first loss and a commuting K3 toggle controls the second.

The next concrete boundary is n=13, where a new degree-13 prime action and the repeated-five-orbit partition 5+5+3 appear together. The present file does not classify them. A second target is an intrinsic proof for the synchronized 5+3+3 family, rather than further blind enumeration. No literature-priority or universal mixed-prime equality claim is made.
