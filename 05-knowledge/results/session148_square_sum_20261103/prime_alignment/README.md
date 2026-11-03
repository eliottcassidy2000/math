# Prime shells and square-sum graphs: explicit alignment, not identification

## Provenance and scope

This is a finite comparison, not a moonshine explanation or a natural isomorphism. Session132 source `05-knowledge/results/session132_fractran_tournament_recursion_20261031/prime_graphs/PROOFS.md` establishes the first two maps and original five groups. Session131 `.../moonshine/README.md` provides the sporadic-group support, exact PRIMEGAME coordinate roles, and relabeling control. The new third map and finite witnesses below are supplied as a separately verifiable comparison. `alignment.csv` preserves all 25 rows and fields; `paths.json` contains all 25 Hamilton paths, `verify.py` checks them and exhausts all 2≤n≤24.

The 25 actual primes less than 100 are 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97.

The input groups recovered verbatim are:

- A={2,3,5,7,11,13}; B={17,19,23,31,47}; C={29,41,59,71}; D={37,43,53,61,71}; E={73,79,83,89,97}.
- There are 25 entries but only 24 distinct primes: 71 belongs to both C and D, while 67 is missing. Do not silently repair. If the purpose was to have disjoint groups covering primes <100, replacing D's 71 by 67 is one possible and especially natural conditional repair, because then D∪E is exactly the non-Ogg primes; it is not evidence of intent. Swapping C's 71 instead is a different set-theoretic repair but destroys the clean C interpretation.
- From the published order factorizations, A=π(Fi22), A∪B=π(Baby Monster), A∪B∪C=π(Monster). Also A={p:p² divides |Monster|}. The Ogg primes are all 15 in A∪B∪C. Absolute/Ogg supersingular means every geometric supersingular j-invariant in characteristic p lies in F_p. This is equivalent to genus zero of X_0^+(p), NOT X_0(p), and is distinct from supersingular reduction of a particular elliptic curve.
- Sources: [Fi22 ATLAS](https://brauer.maths.qmul.ac.uk/Atlas/v3/group/Fi22), [Baby Monster ATLAS](https://brauer.maths.qmul.ac.uk/Atlas/v3/spor/B), [Monster ATLAS](https://brauer.maths.qmul.ac.uk/Atlas/v3/spor/M), [Duncan–Ono](https://arxiv.org/html/1411.5354v2).

## The three indexing maps

G_n has vertex set {1,...,n}, and distinct x,y are adjacent iff x+y is a perfect square. A square-sum permutation corresponds to a Hamilton path, not a graph itself. Let p_i be the i-th prime and q_i the i-th nontrivial admissible order in OEIS A090461. Then compare:

1. Prime order: p_i ↦ G_(p_i).
2. Raw rank: p_i ↦ G_i.
3. Admissible rank: p_i ↦ G_(q_i), where the OEIS-style list begins q_i=15,16,17,23,25,26,...,45 for i≤25.

For map 3 the graph-theoretic trivial G_1 also has a Hamilton path: under a convention that includes it, q_1=1 and every subsequent pairing shifts. OEIS A090461 prints 15 as its first term, so its quoted initial list is adopted explicitly; this is a convention, not a claim that the vacuous singleton is mathematically impossible. See [OEIS A090461](https://oeis.org/A090461). Gerbicz's Jan. 21, 2018 OEIS comment reports a proof that all n≥25 have paths and all n≥32 have cycles, linking the original Mersenne Forum post. That linked post returned HTTP 403 in this check, so its algorithm was not independently audited here. Repository THM-4851 likewise credits these global claims as cited classical results. Our finite table needs only its own exhaustive checks and witnesses.

The map 3 paired rows p_i:q_i are 2:15,3:16,5:17,7:23,11:25,13:26,17:27,19:28,23:29,29:30,31:31,37:32,41:33,43:34,47:35,53:36,59:37,61:38,67:39,71:40,73:41,79:42,83:43,89:44,97:45. This rank matching follows from a convention, not arithmetic or group theory.

## Finite checks, theorems, and negative controls

- G_n is connected iff n=1 or n≥14. A is exactly those first-25 prime orders p for which G_p is disconnected. Proof: G_14 is a tree, all n≥15 attach to an earlier vertex by a square root, and finite checks handle n=2..13. The source gives the explicit tree and covering intervals.
- G_n is a forest for n≤14, has girth 5 for 15≤n≤29, and girth 3 for n≥30: the 5-cycle 1-3-6-10-15-1 exists in G_15 and triangle 6-19-30-6 exists in G_30. Session132 certifies the absence of triangles and 4-cycles through G_29 and uses induced inclusions. Thus map 3 sends the first nine primes through 23 to girth 5 and the remaining 16 to girth 3. Six Ogg primes 29,31,41,47,59,71 map to girth 3 alongside all 10 non-Ogg primes, so girth does not classify Ogg status.
- For map 2, within 1≤n≤25 the nontrivial Hamilton-path orders are 15,16,17,23,25, so the corresponding prime labels are 47,53,59,83,97; n=1 is the trivial exception. There are no Hamilton cycles in this range. For map 1, prime orders 2,3,5,7,11,13 are disconnected; 19 is connected but has no Hamilton path; 17,23,29,31 have a Hamilton path but no cycle; 37,41,43,47,53,59,61,67,71,73,79,83,89,97 have checked Hamilton cycles. These are finite Session132 certificates, not an all-prime theorem.
- In map 1 the same arrival-root set R(41)=R(43)=R(47)={7,8,9} crosses C, D, and B. In map 3, p=31,37 map to G_31 and G_32 yet one is Ogg and one not; their pairwise group/class labels cannot be inferred merely from graph order. More strongly, the finite 5/3 girth test exhibits the six explicit Ogg false negatives listed above. These negate specified proposed classifiers, not all conceivable relationships.
- Independent checker constructs adjacencies by integer square testing, exhaustively searches Hamilton paths on n=2..24 by bitmask recursion with exact memoization, validates explicit paths on the 25 q_i, and tests a missing vertex and a broken edge as negative controls. Run `python verify.py` and `python -O verify.py` from any directory. It uses no `assert`; both outputs passed. The checker also independently constructs graph adjacency lists, component counts, edge totals, and girths, and compares the regenerated table to the stored CSV without writing files.

## PRIMEGAME and why equal labels do not force a bridge

The ordered fractions are 17/91,78/85,19/51,23/38,29/33,77/29,95/23,77/19,1/17,11/13,13/11,15/2,1/7,55/1. Multiply the current integer by the first fraction yielding an integer, and repeat. Conway's initial seed 2 produces successive pure powers 2^p for *all primes* p; the initial seed itself is 2^1 and is not a prime output. [Conway's original chapter](https://gwern.net/doc/cs/computable/1987-conway.pdf).

The prime coordinates appearing in the printed program are exactly first ten primes, all Ogg: data 2,3,5,7 and control 11,13,17,19,23,29. This role split does not equal A/B/C; control 11 and 13 are in A, 17,19,23 in B and 29 in C. Under an injective prime renaming applied consistently to fractions and seed, divisibility and first-enabled priority are preserved, giving a trajectory conjugacy. In particular 29→67 retains the machine's structure but replaces an Ogg coordinate by a non-Ogg coordinate. This is not an automorphism of the unchanged printed program and does not rule out every richer connection, but it disproves an inference from its pure register/guard structure alone.
