# Determination of the deficit-two center polynomial

Status: complete producer proof with finite-exact residual certificates; independent residual replay pending. The preceding degree-28 reduction has already been independently approved. 2026-10-20. No remote writes.

## Result

The only strong self-converse tournament with Delta=H/|Aut|-B=2 is the unique strong tournament S4 on four vertices. Consequently

    D(z)=z^4,
    deg D=4,
    d_4=1 and d_j=0 for every j!=4.

This conclusion combines the audited all-order cutoff with a complete residual enumeration, not an extrapolation from a census through a fixed small order.

## 1. Already approved structural reduction

The companion PREFIX_CONNECTIVITY_CUTOFF.md (frozen SHA256 f4b36f9ce5a7549311a1d2a1d54cf140f59395e6f32177c8fb274d0b5980cdf4) proves:

1. A half-supported clique of t Hamilton words with even-order pairwise relative permutations forces Delta>=2(t-1).
2. In particular, if an induced strong four-set S has three vertices beating v outside S and T-S is strong, then Delta>=4, provided n>=8. Rotate a directed four-cycle on S and append a Hamilton path of T-S starting at v.
3. Every strong self-converse deficit-two tournament of n>=8 has vertex connectivity at most4. The local configuration needed here is checked exhaustively on six vertices, and Moon's primary theorem extends it to all larger strong tournaments.
4. For any Hamilton cycle successor R, write m=|Aut(T) intersect <R>| and q=n/m. The cycle-cut lemma gives q<=4, and m is odd. The q rotation orbits are regular tournaments of order m and linked by directed perfect matchings.
5. For q>=2, connectivity is at least ceil(m/2). Thus m<=7. The q=1 case is removed at every unknown order using its weaker ceil(n/3) bound, the odd-order cycle obstruction and the already complete census through9.

The only previously unknown orders are therefore10,12,14,15,20,21,28. The inherited exact census through9 gives the unique positive S4 example and no others.

## 2. The m=7 endpoint reduction

For q>=2, m=7, and connectivity<=4, both minimum outdegree and minimum indegree of a self-converse T are exactly4.

Proof. Connectivity is at least4, so a minimum cut has size s=4. In the cut split X->Y, a rotation orbit meeting both sides costs at least ceil(7/3)=3 removed vertices. There can be only one such orbit, O_0; no split orbit would contradict the surviving perfect matchings around the orbit cycle. The other orbit labels cannot have a Y-to-X transition. If both X and Y labels occur, the matching inequalities in section5.1 of the cutoff proof imply7<=s, impossible. Suppose all other labels are X. Write s_i for the deleted vertices in O_i and a_0,b_0 for the two surviving sides in O_0. We have s_0>=3, b_0<=s_1 and a_0<=s_0. The choice s_0=4 would force b_0=0, so s_0=3,s_1=1, a_0=3,b_0=1 and all other deletions vanish. Thus Y is globally a singleton, whose outdegree is at most4. Its internal regular outdegree is3 and the Hamilton cycle supplies an external outgoing edge, so its outdegree is exactly4. The all-Y case is dual. Self-converseness exchanges the in/outdegree inventories, giving both minima4. QED.

The independent cutoff auditor separately approved this refinement in section6 of session15-deficit-cutoff-audit/AUDIT.md.

## 3. Complete periodic presentation space

Label the chosen Hamilton cycle 0->1->...->n-1->0. Then R^q:i->i+q modulo n is an automorphism.

For odd q=3, the cycle-cut lemma gives exactly one blue cut orbit. Rotate the labels until that cut is the standard path; its reversal r:i->n-1-i is an anti-automorphism.

For q=4, there are exactly two blue cut orbits, so the same normalization supplies r. These cases therefore have both R^q automorphic and r anti-automorphic.

For q=2, a Hamilton cycle may have no blue cut; we make NO reflection assumption and enumerate the whole R^2-invariant cycle space. This avoids the known dangerous even-order gap.

Generation uses one Boolean variable per orbit of unordered vertex pairs under R^q and, when available, r. Carry an orientation parity with each pair: R^q preserves arcs, r reverses arcs, and sorting an ordered pair flips orientation. The generator propagates these equations and rejects inconsistent orbits. The n prescribed cycle arcs force their variables. Every remaining variable is enumerated independently.

At m=7 and q=3 or4, the degree lemma supplies an orbit containing outdegree4 vertices. Each vertex already has internal outdegree3, and its successor along the Hamilton cycle is external. Therefore it has exactly that ONE external outneighbor. For each possible orbit residue low=0,...,q-1, force vertex low to beat its cycle successor and lose to every other external vertex. Rotation imposes the same condition throughout that orbit. Inconsistent residue cases are empty and discarded. A candidate with low indegree also has low outdegree because it is self-converse, so this restriction loses none.

The complete residual presentation counts are:

    n10,q2:        7 free variables,        128 presentations
    n12,q4:       11 free variables,      2,048 presentations
    n14,q2:       11 free variables,      2,048 presentations
    n15,q3:       12 free variables,      4,096 presentations
    n20,q4:       21 free variables,  2,097,152 presentations
    n21,q3,low0:   6 free variables,         64 presentations
    n21,q3,low2:   6 free variables,         64 presentations
    n28,q4,low0:  12 free variables,      4,096 presentations
    n28,q4,low1:  12 free variables,      4,096 presentations
    n28,q4,low2:  12 free variables,      4,096 presentations
    n28,q4,low3:  12 free variables,      4,096 presentations

Total: 2,121,984 presentations. These are normalized presentations, NOT pairwise nonisomorphic tournament classes. Duplicates are harmless for exhaustive exclusion.

## 4. Every presentation has a direct prefix obstruction

For every presentation, the verifier finds a four-set S and terminal v outside S such that:

- T[S] is strong (equivalently, every internal outdegree is1 or2);
- at least three vertices of S beat v;
- T-S is strongly connected, verified by independent forward and backward reachability.

All 2,121,984 presentations pass this obstruction test. There are zero survivors. Since n>=10 in every residual family, the half-support hypothesis n>=8 holds. Section1 therefore gives Delta>=4 in each case, contradicting the putative Delta2.

No Hamilton-path counts, heuristic automorphism estimates, floating-point comparisons, isomorphism pruning or random samples are used in the exclusion. The only group data imposed are the necessary literal successor symmetries. The search records an actual small witness for each presentation.

Together with the exact orders through9, this proves D=z^4.

## 5. Reproduction and compact certificates

All files are local staging files; the frozen preceding session is unchanged.

    python3 generate_residual_frames.py
    g++ -O3 -std=c++20 scan_residual_frames.cpp -o scan_residual_frames
    for f in frames_n*.txt; do ./scan_residual_frames "$f"; done

The descriptor files encode fixed/free orientation variables as two possible outneighbor-mask vectors per variable. The scan enumerates the free variables in binary-reflected Gray order: at step i the assignment mask is i xor(i>>1). It records one five-byte witness per presentation: the four increasing vertex indices of S followed by v. The raw binary streams are also gzip-compressed with reproducible timestamps. witness_certificate_index.json gives record counts and both raw/compressed SHA256s. A checker can reconstruct the tournament for each step and verify the five-byte witness without repeating the search.

The complete run is residual_all_scan.jsonl. The exhaustive smaller-family checks and a separately selected514-frame order20 test additionally verify tournament completeness, cycle arcs, rotation symmetry, optional anti-symmetry, lowdegree conditions, and an independent cut/reachability implementation of the prefix predicate; these are residual_sample_checks.json. Those checks are validation only and are NOT substituted for the complete order20 enumeration.

The small last-subset search positions explain the short runtime: the n20 scan never examines beyond the41st four-subset in lexicographic order, and every n28 low0/low1 frame uses the first four vertices. Runtime is not mathematical evidence; the full frame coverage and checkable certificates are.

## 6. Exact leaf-series consequence

Insert D=z^4 into the inherited leaf reduction:

    MixedLeaf(z)=z^4/(1-z^2-z^6),

    Leaf(z)=1/[2(1-z-z^3)]
            +(2z^4-1-z-z^3)/[2(1-z^2-z^6)].

Equivalently, every mixed paired-black leaf is a unique rigid singleton/C3 word A, followed by S4, followed by A^op. This classifies mixed leaves, not all bridges or the wiring of the nonleaf core.

The mixed-leaf coefficients are exactly m_n=0 for odd n and m_(2r)=a_(r-2), where a_0=1, a_j=0 for j<0, and a_j=a_(j-1)+a_(j-3). In particular the mixed recurrence m_n=m_(n-2)+m_(n-6) holds for every n>4, with the single forcing term at n=4. No unknown eventual onset remains.
