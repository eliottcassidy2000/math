# Creator minima versus substitution primes

Status: elementary PROVED statements, inherited CITED inputs, and separately marked FINITE-EXACT census. Repository baseline: f870c2abf. Proofs, finite evidence, and inherited classification are distinguished below.

## Useful conclusions

1. **PROVED, with explicit finite witness:** even the set of *all* first creators is not a substitution basis for all tournaments. The earliest obstruction is order seven, not the two order-eight creators of 63. The critical prime W7 has H=31, which is already realized at order six. A prime cannot be assembled by any genuinely nontrivial substitution.
2. **FINITE-EXACT:** the two creator orphans of 63 are not new substitution primes. Both are substitutions into W5, the prime first creator of 9, with three directed-pair children and two singleton children. Their size vectors in the explicit W5 convention below are (1,2,2,2,1) and (2,2,2,1,1).
3. **PROVED, all orders:** no finite list of tournament seeds generates all tournaments using substitution (including ordinal sums). The W_(2k+1) family supplies arbitrarily large substitution-prime obstructions. Allowing H-neutral vertex additions does not repair this claim.
4. **PROVED, all orders:** H(W_(2k+1)) = (4^(k+1)+6k+5)/9, so its two-vertex critical lift obeys H(W_(n+2)) = 4 H(W_n) - n. This gives a concrete operation and scalar calculus for the critical obstruction, rather than unrestricted insertion.
5. **CITED bridge, not claimed new:** Schmerl–Trotter/THM-4169 says the only prime tournaments with no prime one-vertex deletion are the three odd-order critical families. Thus a prime-focused grammar must distinguish ordinary prime-preserving attachments from critical-family lifts. An all-order prime recursion can use the critical families as parameterized constructors and prime-preserving one-vertex extension, with substitution for decomposable objects. This is structurally restricted, but it is not the stronger conjecture that all prime residuals have bounded cut defect, and it does not prove the existing five-operation grammar at all orders.

## 1. Substitution irreducibility theorem (PROVED)

A module M is a vertex set such that each outside vertex either dominates every member of M or is dominated by every member of M. A tournament is prime here when it has order at least three and no proper module of size at least two. Keep the directed pair as a separate elementary constructor.

In Q[A_1,...,A_m], every child A_i is a module. If the output T is prime and m>=2, every A_i has size one; hence Q is isomorphic to T. If m=1, the operation returns its sole child. Therefore, in any finite expression built from seeds by substitution, every prime output must already occur as a seed (after deleting identity wrappers).

Conversely, the directed pair, singleton, and all prime tournaments generate every tournament. Induct on order. If T is prime it is supplied. Otherwise choose any proper module M with 2<=|M|<|T|. Contract M to one vertex, obtaining the strictly smaller quotient Q. Then T=Q[... ,T[M], ...], and both factors are smaller. Order two is the separate base. No uniqueness assertion about the chosen decomposition is needed for this proof.

Thus substitution irreducibles, not value-minimal first creators, form the necessary and sufficient primitive set for substitution-only generation. This is a statement about this operation set, not about every imaginable grammar.

Ordinal sum is substitution into a transitive quotient, and therefore adds no escape. H-neutral insertion, by inherited THM-1900(B), inserts a singleton strong component, hence is also an ordinal-sum operation. It adds no escape either.

## 2. Elementary critical family and count (PROVED)

Let W_(2k+1), k>=1, have chain vertices 0,...,2k-1 with i->j iff i<j, and a hub z with z->i for even i and i->z for odd i.

### Primality

A module M not containing z must be an interval of the chain: any omitted vertex between its minimum and maximum distinguishes those extremes. If |M|>=2 this interval contains both parities, so z distinguishes its vertices. Hence no such nontrivial module exists.

If M contains z and at least one chain vertex, then 0 must belong to M: otherwise z->0 but every positive chain member of M loses to 0, violating the module property. Similarly 2k-1 must belong to M. An outside even vertex is beaten by z but beats 2k-1; an outside odd vertex beats z but loses to 0. Neither is possible, so M is the whole tournament.

Thus W_(2k+1) is prime at every odd order >=3. This alone proves the no-finite-substitution-basis result.

### Criticality

Deleting z leaves a transitive tournament. Deleting an interior chain vertex i leaves the adjacent chain vertices i-1,i+1 as a two-element module: they have equal parity and every remaining chain vertex sees them equally. Deleting 0 leaves 1 as a source, so the complement of 1 is a nontrivial module when k>=2. Deleting 2k-1 similarly leaves 2k-2 as a sink. Hence W_(2k+1) is critical for k>=2. Order three is the seed exception, since a two-vertex deletion card is not called prime here.

### Hamiltonian paths

The paths with z first or last contribute two. A path with interior z is an increasing chain on a subset, followed by z, followed by the increasing complementary chain. Let a be the predecessor and b the successor of z. Necessarily a is odd and b even.

If a<b, the partition can cover all chain vertices only when b=a+1, giving k-1 paths. If a>b, every vertex below b must be in the first chain, every vertex above a in the second chain, and the a-b-1 intervening vertices may be assigned freely. For a-b=2d+1 there are k-d choices of (a,b), each contributing 4^d paths. Therefore

H(W_(2k+1)) = k+1 + sum_(d=0)^(k-1) (k-d)4^d
             = (4^(k+1)+6k+5)/9.

The polynomial-geometric sum follows by induction or differentiation of the finite geometric series. Subtraction gives H(W_(n+2))=4H(W_n)-n. Deleting the last two chain vertices of W_(2k+3) leaves W_(2k+1), so this is a genuine nested two-vertex construction, and every step after the seed lands in a prime with no prime one-vertex deletion.

No claim is made that W_n is a non-first-creator for every n>=7. That would require additional realization bounds or constructions.

## 3. Explicit smallest obstruction (FINITE-EXACT plus the preceding proof)

W7 has H=31 by the formula. There is a six-vertex tournament with outgoing-neighbor bitmasks

    [2,12,33,21,7,27]

and H=31. The standalone DP in probe.py verifies this integer independently of the generator's Hamilton counts. In upper-triangle encoding (successive pairs (0,1),(0,2),...; one means i->j), its code is 6241. Hence W7 is not a first creator. Its primality is proved in section 2, so section 1 excludes it from the substitution hull of *all* first creators, even allowing all representatives in every creator fiber.

All prime classes through order six are first creators. Through order eight the census is:

| order | classes | prime classes | noncreator prime classes | critical noncreator primes |
|---:|---:|---:|---:|---:|
| 3 | 2 | 1 | 0 | 0 |
| 4 | 4 | 0 | 0 | 0 |
| 5 | 12 | 3 | 0 | 0 |
| 6 | 56 | 15 | 0 | 0 |
| 7 | 456 | 197 | 8 | 1 |
| 8 | 6880 | 4008 | 214 | 0 |

The eight order-seven obstructions have H-multiset {27,31,31,41,41,41,43,43}; only W7 has no prime deletion. This separates the critical obstruction from the seven ordinary prime attachments. This finite witness is not a counterexample to THM-4942's richer finite grammar, which already allows first-creator vertex extensions.

## 4. Resolving the 63 genealogy orphans (FINITE-EXACT)

Exhaustive class generation gives exactly two order-eight H=63 classes, with upper-triangle codes 135528708 and 36175876. Independently checking *every* proper vertex subset gives their complete nontrivial module lists:

- 135528708: {0,1}, {2,4}, {5,7};
- 36175876: {1,2}, {4,6}, {5,7}.

Contract these disjoint pairs. Each quotient is prime, has order five and H=9, and explicit isomorphism testing identifies it with the W5 convention of section 2. In W5 vertex order (0,1,2,3,z), the respective child blocks are

- ({6}, {5,7}, {2,4}, {0,1}, {3}), sizes (1,2,2,2,1);
- ({5,7}, {4,6}, {1,2}, {0}, {3}), sizes (2,2,2,1,1).

Each pair induces the unique directed pair. W5 is a first creator of 9: H(W5)=9, and the complete census through order four has maximum H=5. Both orphans are therefore compositions from the first creator W5 and elementary directed pairs; a directed pair itself is an ordinal sum of singleton seeds.

Their deletion H multisets are both {1,11,25,25,33,33,37,37}. Their failure to have a first-creator deletion is a failure of the one-vertex *value-birth* genealogy, not a new modular primitive. This is exactly why the two notions of fundamental structure should not be conflated.

## 5. Reproduction and evidence boundaries

Files:
- classgen.c is copied verbatim from snapshot f870c2abf, path 05-knowledge/results/first_creators_genealogy_20261107/classgen.c. Its class enumeration is inherited, not independently invented here.
- probe.py independently implements subset DP and checks all proper subsets for modules. It checks every generated class through order eight, not a sample. It also checks W_n through n=13 against the proved formula and criticality.
- `probe.out` retains results; `classes_1.txt` through `classes_8.txt` are deterministic regenerated inputs, built in isolation by the packet's `verify.py`.

Commands, from this directory:

    gcc -O3 classgen.c -o classgen
    ./classgen 8 classes
    python3 probe.py > probe.out

The class counts 1,1,2,4,12,56,456,6880 agree with the established census. Independent DP reproduces H on every class. The module check is exhaustive over all subsets, independent of the inherited classifier. The formulas above are elementary all-order proofs, not extrapolations from these checks.

Sources inherited and inspected:
- THM-4942, 01-canon/theorems/THM-4942-first-creators-at-order-eleven-attachment-form-generation-of-all-tournaments-connectivity-floors.md and 05-knowledge/results/first_creator_generation_20261108/PROOF.md.
- THM-4922, 05-knowledge/results/first_creators_genealogy_20261107/PROOF.md (first-creator definitions, orphan genealogy, all first creators need not be known strong).
- THM-4169, 01-canon/theorems/THM-4169-prime-parent-one-vertex-augmentation-and-quartic-johnson-transfer.md, especially section 4; Schmerl–Trotter DOI https://doi.org/10.1016/0012-365X(93)90516-V.
- 05-knowledge/results/session188_structural_bridges_20261104/critical_quotient_synthesis/PROOF.md records the independently inspected critical-family source definitions.

OPEN: Whether all prime noncreators admit a bounded-cut or bounded-defect derivation from creator parents at arbitrary order; whether the finite five-operation grammar persists; whether the three critical families have predictable first-creator status at every order. None is resolved by substitution irreducibility alone.

## 6. All Fibonacci fans have exactly two fixed prime kernels (PROVED)

This is a stronger all-order interpretation of the 63 calculation. Credit THM-4748, `01-canon/theorems/THM-4748-two-reversed-arcs-exact-path-counts-the-seven-mod-eight-lines-and-the-mersenne-side-three-arcs-from-transitivity.md`, for the Fibonacci fan construction F_k and H(F_k)=2^k-1, k>=4. The construction is not new here. In its convention the chain is 1,...,k+1 and hub 0 loses exactly to a=floor(k/2), b=k+1-a, and k+1.

Let TT_s be a transitive s-tournament; write singleton children as 1. In the W-family convention of section 2 (chain first and hub last), for every r>=2,

    F_(2r)   = W5[TT_(r-1), TT_2, TT_(r-1), 1, 1],
    F_(2r+1) = W7[TT_(r-1), 1, 1, 1, TT_(r-1), 1, 1].

Equalities mean explicit isomorphism, not equality of H alone.

Proof: in the even case, a=r and b=r+1 are consecutive. The chain's maximal constant-hub-orientation runs are

    [1,r-1], [r,r+1], [r+2,2r], {2r+1},

of lengths (r-1,2,r-1,1). The hub alternately beats, loses to, beats, loses to these four runs. Each run is a module because it is a chain interval and the hub treats all its vertices equally. Contracting yields W5. For odd k=2r+1, the losses are r,r+2,2r+2, so the six runs have lengths

    (r-1,1,1,1,r-1,1),

and alternate in the same manner. Contracting gives W7. This proves both formulas at every r>=2.

Consequences:
- The Mersenne fan construction of all 2^k-1, k>=4, does not require a growing set of modular primes: W5 and W7 suffice with transitive inflation.
- Its even exponent branch uses the first creator W5 (H9); its odd exponent branch uses the structurally indispensable noncreator W7 (H31).
- F6 is W5[TT2,TT2,TT2,1,1], exactly one of the two 63 orphans; its converse gives the other size vector recorded above.
- The known first-appearance exception 63 is thus on the simpler first-creator-kernel branch, whereas the 31 witness is the first genuine new modular kernel on the odd branch.
- This does not show every first creator or every tournament uses only these two primes. Section 1 proves any finite substitution basis fails globally.

Prior-art check: THM-4922 names the orphan pair as F6 and its converse, crediting THM-4748; THM-344 establishes their complete odd-cycle conflict graph. THM-4748 explicitly credits earlier one-apex signatures 1001100 and 1100110. Those signatures already encode the runs used here. The inspected S16/S17, session224, THM-344 and THM-4748 sources did not state this W5/W7 modular-kernel reinterpretation. This is a recovered consequence of their explicit construction, not a claim that the fan or its count is new.

The final portion of probe.py constructs both sides directly, verifies identical adjacency matrices and H=2^k-1 for every 4<=k<=12, independently of canon's symbolic count.

## 7. Creator-rooted structural grammar (CITED CONSEQUENCE)

This is an all-order consequence of the cited critical-prime classification, with elementary explicit pair deletions. It is not a new classification or a proof of THM-4942's narrower five-operation conjecture.

**Grammar.** Allow a singleton and directed pair, together with the full first-creator library as optional stopping leaves. The operations are:

1. Proper substitution: contract a proper module, and reconstruct its strictly smaller quotient and child.
2. A one-vertex extension whose parent and output are both prime.
3. A two-vertex extension along a single critical family T_(2r-1)->T_(2r+1), U_(2r-1)->U_(2r+1), or W_(2r-1)->W_(2r+1), with the explicit embeddings below. There are no arbitrary two-vertex attachments in this clause.

**Theorem (CITED CONSEQUENCE).** Every finite tournament has a derivation in this grammar, with strictly decreasing order in the reverse derivation.

Proof: proceed by induction on order. Stop at any first creator, and at the singleton or directed pair. A nonprime tournament of order at least three has a proper nontrivial module, so section 1 gives operation 1 with smaller inputs. For a prime tournament, either some vertex deletion is prime, giving operation 2, or it is critical. Schmerl–Trotter classifies every critical prime of order at least five into T,U,W. Each admits the same-family prime two-vertex deletion explicitly established below, giving operation 3. The only order-three prime is C3, a first creator. There is no order-four prime. The three order-five primes are W5,U5,T5, with H=9,13,15, all first creators since max H(4)=5. These base facts are also checked by the exhaustive small census. Thus the induction terminates.

**Exact citation boundary.** The classification is imported from Schmerl–Trotter, Discrete Mathematics 113 (1993), 191–205, DOI https://doi.org/10.1016/0012-365X(93)90516-V, as already imported in THM-4169 section 4. The source definitions and their correspondence were checked previously in `05-knowledge/results/session188_structural_bridges_20261104/critical_quotient_synthesis/PROOF.md`. The pair-deletion properties need not be left as an additional unverified citation: they follow directly from those explicit definitions, as follows.

- **T_(2r+1):** vertices 0,...,2r, and i->j iff j-i modulo 2r+1 is in {1,...,r}. This is exactly the cyclic family defined in `session188_structural_bridges_20261104/ordered_module_skew/PROOF.md`, section 6. Delete vertices r and 2r, retaining the circular order. Every retained vertex's interval of r original successors contains exactly one deleted vertex: the sole r-consecutive interval missing both is {0,...,r-1}, whose predecessor 2r was deleted; no r-consecutive interval contains both deleted vertices. Thus the retained outneighbors are exactly the first r-1 successors in the new circular order. The result is T_(2r-1).
- **U_(2r+1):** use A_0,...,A_(r-1), B_0,...,B_r, with A_j->A_i for i<j, B_i->B_j for i<j, and A_i->B_j iff j<=i. This is exactly `session188_structural_bridges_20261104/staircase_prime_quotient/PROOF.md`, section 1. Delete A_(r-1) and B_r; the identical defining inequalities on the remaining names give U_(2r-1).
- **W_(2r+1):** in section 2's convention, delete chain vertices 2r-2 and 2r-1 and retain the hub. The remaining chain and its alternating hub incidences give W_(2r-1).

All three statements hold for r>=2, with the order-three family members isomorphic to C3. The checker now verifies literal adjacency equality for these embeddings at every r=2,...,8; the proofs supply the unbounded scope.

**What this achieves and does not achieve.** Insertions occur only inside the prime part of a substitution decomposition; a two-vertex jump is allowed only in one of three explicit critical families. This is a genuine structural restriction compared with unrestricted vertex insertion. Optional creator leaves let a derivation stop as soon as a value-minimal realization is reached. Nevertheless, the *full* creator library is redundant for this grammar: singleton/directed pair/C3 already suffice once the cited critical lifts and prime-preserving extensions are allowed. Neither this derivation nor its stopping point is asserted unique. It proves no minimality of the creator library, no all-order completeness of the H spectrum, and no bound on distance to a creator leaf or on cut defect. Those arithmetic-strengthened assertions remain separate open questions.

Final extended-check receipt: after adding the section 7 checks, normal and optimized runs both completed successfully and their entire outputs are byte-identical. The current SHA-256 of `probe.out`, `probe.normal.out`, and `probe.optimized.out` is `65b5417a90e21987ab56b477657b8745be80ca9df42829d000368f0906ded19c`. The order-five check corrected a draft base-value typo: U5 has H=13, not 11; the first-creator conclusion is unchanged.
