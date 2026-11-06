# Independent audit of the paired carriers and full-state obstruction

6 November 2026. **PASS for both stated all-order results.** Four physical O/Z carrier orbits supply simultaneous staggered residue certificates. A separate start-two family proves that no fixed finite O/Z/J seed bank covers every full nine-coordinate state. Neither result settles completeness of the pair (B,s), F(8,137), or an eventual odd-spectrum threshold.

## Frozen inputs and conventions

The reviewed proof has SHA-256 777fec9e14b9219855268c157b3579b1bc5b76ba8ae1f194a8e569dc8d1f57a2. The producer verifier has hash d643529d29b0548317ba8ac77185c46a8b07e9d02c9ee9c1462a6a59d04b5a6f. Both frozen producer receipts have hash 5cc9a1e61d8ca0b4efe2f260391fd6dd1527c742f438e030f7e5a2fddc7fb55a.

The inherited nine-coordinate transfer is retained at dependencies/ROOTED_COVER_CLOSURE.md, hash 8e78951f9dd15decff0d23cfb3ebe64f358b3c9d4561b269b527a1e9aa902168. The conditional staggered criterion is retained at dependencies/STAGGERED_PROFILE_CRITERION.md, hash cf97dd951c7b641f236f3e931b9b3e2bce23125fe31f746769ec1c4034cfd196. No discovery-bank completeness is used.

The coordinates h,s,e,L,R,d,C,b,T use labeled cover bins and allow exactly the stated ordinary bins to be empty. In particular C counts two ordinary bins on U, b counts them on U-root, and T has one nonempty root-starting bin plus two optional empty ordinary bins. Re-marking under J is essential; J is the same underlying addition as Z but has a different root.

## The four actual coupled carriers

All eight displayed adjacency lists and roots are valid tournaments. Independent exact counting reproduces every stated seven-coordinate seed profile and the four common adjusted rows:

    (X,e,V,d)=(677,27,1061,49), (759,29,927,29),
                (1547,101,1691,109), (821,33,1079,33).

The fixed j pairs are (3,3), (4,3), (4,3), (5,5). For each pair, the first source has s=1+3j and the second has s=2+3j. Their unadjusted B values need not agree. What agrees is X=B+2j together with e,V=R+b,d. This is the coupled invariant needed for the two staggered residue tracks.

O and Z retain the marked root and preserve s on these nontrivial inputs. Therefore j stays fixed and remains in [0,8]. The inherited physical transfers give

    O: (X,e,V,d) -> (X,0,2V,d),
    Z: (X,e,V,d) -> (X+V-e,d,2V-d,d).

Applying any identical finite word preserves the common adjusted row. Each resulting source is an actual tournament, not a point admitted merely by coordinate inequalities. The fixed four-core readout has H=3B+2s=3X+2r. The independent checker constructs that graph literally and confirms both residues on 120 cases.

Every adjusted base is odd because B is the Hamilton count of a physical edge extension and X differs by 2j. O preserves X. Z increases it: R>=d>=e and b>0, so V-e>0. Thus the lower bound X>=137 is preserved from each listed seed. This establishes certificates on the four word orbits, not on every odd base.

## Pure-Z boundaries and the J hostile

For a pure-Z word of length k>=1, the exact solution is

    V_k=2^k(V_0-d)+d,  e_k=d,
    X_k=X_0+(2^k-1)(V_0-d)+d-e_0.

Substitution gives the four sequences 1012*2^k-313, 898*2^k-139, 1582*2^k-27, and 1046*2^k-225. The second and fourth extend to k=0 because their initial e equals d; the first and third do not. These endpoint qualifications are correct and independently checked.

The O/Z restriction is real. For the 759 pair, independent counting gives (C,T)=(4384,2713) and (4552,2797). Applying J to the actual graphs gives adjusted bases 4379 and 4547. Their formerly common row separates. Equality of (X,e,V,d) is therefore not enough for arbitrary J continuation.

## The three physical faces

The quantities e, d-e and f=R-h+e-d are nonnegative for every physical rooted tournament, including the singleton. Root deletion injects ending paths into paths on U-root, proving d>=e.

For f, count R-covers with a designated root-starting path Q and an ordinary path P. Cutting any Hamilton path immediately before the root gives h distinct such covers, recoverable by concatenating P then Q. The covers with Q equal to the singleton root give d further covers. Their intersection consists exactly of Hamilton paths ending at the root, so has size e. The remaining covers have Q of length at least two, P nonempty, and root->last(P). Their number is R-h-d+e=f. This proves the stated interpretation, not merely an inequality on a relaxed state cone.

Every O child has e=0. Every Z child has d-e=0. Every J child satisfies f=0 by the exact all-input formulas h'=h+R-e, e'=h-e, R'=h+R and d'=h. These face identities include singleton inputs. The separate singleton exception for the start coordinate of Z does not affect them.

## The strict-interior family and grammar conclusion

Starting from a rooted edge, repeated J gives U_n. The newly marked root has only the previous root as an out-neighbor. The previous root retains its unique start under the underlying nontrivial Z addition. Thus both designated vertices have outdegree and start count one. U_3 is a directed triangle; adding a vertex which beats the old root and loses to the other vertices preserves strong connectivity, proving U_n strong for n>=3.

W_n adds a new root v beating precisely those two unit-start vertices. Every v-starting Hamilton path next enters one of them, so its start count is exactly two. The old strong U_n reaches a vertex outside those two, which beats v, while v reaches U_n through a unit-start vertex. Hence W_n is strong. Its ending count is positive, and

    d(W_n)-e(W_n)=end_0(U_n)+end_1(U_n)>0.

For n>=4, removing the latest root from U_n leaves the strong U_(n-1). A Hamilton path of that deletion ending at the preceding root, together with Q=(v,latest root), gives the required remaining cover counted by f. Therefore f>0. The elementary fact that every vertex of a strong tournament starts and ends a Hamilton path can be proved by root-preserving insertion from a singleton, or its converse; no Hamilton-cycle theorem is required.

Every full state of W_n lies outside the three possible last-letter faces. It cannot equal the state of any nonempty O/Z/J word, even if that generated graph is not isomorphic to W_n. Moreover d(W_n)=H(U_n)>=n, since every vertex of U_n starts at least one Hamilton path. These deletion counts are unbounded, so the family contains infinitely many distinct full states. A finite seed bank cannot contain them all as empty-word states. The finite-bank conclusion follows.

This argument does not exclude a generated graph having the same (B,s) while different e,R,d or other coordinates. It is a full-state obstruction, not profile-pair incompleteness and not a refutation of F. The proof's restriction is explicit and correct. No assumption that the selected finite carrier seeds represent all physical sources enters either theorem.

## Independent computation

check.cpp uses append-at-end subset DP with a separate root-start table. This differs from the producer's remove-a-vertex recurrence. It preserves the seeds' original root labels instead of first relabeling every root to zero. It computes optional-empty bin covers directly and independently counts the covers constituting f by literal permutations and cuts on the small controls.

Normal and debug builds agree on:

- All eight complete seed profiles and 504 actual carrier graphs from all O/Z words of length at most five.
- 120 literal four-core readout graphs, including both residue outputs.
- The four pure-Z formulas through k=20, with their exact k=0 boundaries, and the two J-hostile graphs.
- All 285 marked inputs through order four and their 855 actual O/Z/J children, checking the complete inherited nine-coordinate transfers and all three output faces.
- 1,140 direct remaining-cover interpretations of f on those inputs and children.
- All eleven stated W_n cases for 4<=n<=14, with strongness, the two unit starts/outdegrees, s=2, both endpoint identities and strict interior values. An extra n=3 control has f=2; the proof's n>=4 range is sufficient, not asserted sharp.

There are 1,810 exact profile computations. The largest graph has 15 vertices. Even the three-bin quantities are bounded by binom(17,2)*15!, safely below the signed 64-bit limit. All checks use explicit failures in every build. The producer's optimized and debug programs are replayed and reproduce both frozen receipts.

## Reproduction and scope of this packet

    g++ -std=c++17 -O3 -Wall -Wextra check.cpp -o check
    ./check > NORMAL.json
    g++ -std=c++17 -O0 -g -D_GLIBCXX_ASSERTIONS -Wall -Wextra check.cpp -o check_debug
    ./check_debug > DEBUG.json
    cmp NORMAL.json DEBUG.json

Producer replays use PRODUCER_VERIFIER.cpp. The optional density addendum and any later cited order cap are separate results and are not part of this main proof approval. The included dependency provenance uses immutable published commits; all local SHA-256 and Git blob identities match. Portable publication should preserve these reviewed bytes. Status promotion and publication remain separate checks.
