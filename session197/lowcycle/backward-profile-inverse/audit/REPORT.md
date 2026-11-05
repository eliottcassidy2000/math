# Independent audit of the backward start-one profile inverse

**PASS, no correction.** This binds the final frozen snapshot in `session197-lowcycle/`:

- BACKWARD_PROFILE_INVERSE.md:8a07f1bb59f7ee3af9c172a1cfec8551c87125747eeb95dddef7335df33df1fb.
- verify_backward_profiles.py:3f11b79a5470736e01ce1596a5d7b5fe9b75be9a78489cc285a9cfc2e22cf4b4.
- backward_profile_certificate.json:740db5c2caf7f0c1ea218195ba175d1d502622d76dc775d07483b1c4d2c113c9.

The result is exact all-order membership and a physical witness recipe for the complete strong rooted start-one profile class, plus the two explicit count-one bases. It gives an exact scalar C1 criterion after existentially checking its bounded response domain. It does not give a cofinite scalar-coverage theorem or an inverse for arbitrary roots, C2/C3 profiles, or all integer cone points.

## 1. All-order proof audit

The inherited one-seam grammar is correctly used as a physical completeness theorem, not merely as a forward semigroup of formal matrices. Every gadget is determined by a finite binary word, and its transitive root-deletion and unique ending path guarantee start-one preservation and strongness except at the two named small bases.

For a final L1 step at h>3, the predecessor is nontrivial strong; singleton and two-vertex transient predecessors could only produce count1 or3. The predecessor cone gives h<=3e. At equality the output slack R-h is exactly its positive predecessor ending count, so the boundary R=h is excluded. For L>=2, the predecessor h is the new e. Even L gives h>=3e, with equality only for L2 and predecessor R=e; then the output R equals h. Odd L>=3 gives h>=5e. These are exhaustive possibilities and prove the exact numerical last-step partition, including the tie case. No inference from arbitrary C1 roots is used.

The singleton inverse is the exact solution of the L1 equations. Its necessary growth filter is valid for every nontrivial strong predecessor. It is deliberately followed by reachability testing; it is not treated as an inverse-closed cone map.

In the nonsingleton branch, R_old=(R-Qe)/B must be integral and at least e. Since Q>=1 and B>=L, this gives precisely L<=floor(R/e)-1, including the correct upper-end off-by-one. Odd L forces the old endpoint count. Even L genuinely erases that count, so testing all integers0 through e is necessary; in particular0 must remain available for the two-vertex transient. Recursive acceptance then rejects impossible endpoints. There is no division by an unchecked zero, since this branch is entered only with e>=1.

Every recursive predecessor has smaller h: at most2h/3 for the singleton branch and at most h/3 for a nonsingleton. The only nongrowing count-one transition is absorbed into the explicit base state. At count3 both possible presentations give exactly the triangle triple and have predecessor count1. Thus accepting only the three stated base triples is complete, and returning just one triangle recipe loses no profile or witness.

Soundness follows by physical reconstruction from an accepted smaller predecessor. Completeness follows by taking an actual last gadget in the full rooted grammar; the numerical partition chooses its type and the finite candidates include its index and old endpoint. The induction therefore establishes an all-order decision procedure, not a bounded census result.

For the deletion sidecar, deleting the newest gadget root leaves the predecessor followed by a transitive tail, so the actual deletion count is exactly predecessor h. The last-step partition and the singleton inverse determine that h from the three output responses. The triangle ambiguity has equal predecessor counts. This proves the stated nonlinear piecewise readout and genuine deletion uniqueness for reachable triples. It is not a linear quotient of the unrestricted four-response interface.

The numerical complexity statement is appropriately scoped. There are polynomially many valid triples below a numerical H bound. At each nonsingleton state there are O(H/e) gadget indices and at most e+1 endpoint alternatives, hence O(H) candidate calls. Invalid candidates may be rejected immediately; memoizing them too does not invalidate polynomial numerical complexity. The argument does not claim polynomial complexity in log H.

## 2. Independent complete physical comparison through H101

Rather than regenerating the source's forward grammar, the audit uses the independently regenerated complete C<=3 graph bank from the earlier C3 census review:

- Full bank `session197-bridges/c3-independent-records.csv`, SHAa9df9d85412f37e65a1059db21a70a981ae30f98e553ee58229786bad8e0038c.
- Its completeness review `session197-bridges/C3_COMPLETE_CENSUS_REVIEW.md`, SHA0b6d87063bd9a655925a31189962fcd00931ea03d72de31849e9b1223b7d73d3.

The supplied `c1-inputs.csv` is the exact C1/order<=9 extraction, containing232 graphs; the checker adds the unique directed triangle explicitly. Completeness below H101 follows also from the already reviewed strong Busch bound:101^3=1030301<1953125=5^9, so no strong graph of order10 or larger can intervene. This finite comparison attributes bank completeness to the earlier independent review; it is not a new unlabelled generation claim.

The standalone C++ checker enumerates all full vertex permutations and all cuts of each word. It directly counts Hamilton paths, starts, ends, ordered two-covers and root-deletion paths. Across233 actual C1 graphs it finds357 start-one roots with H<=101 and324 distinct triples after adding the two count-one bases. It checks that every repeated physical triple has the same deletion value.

The independent Python inverse deliberately tries both possible last-step types, rather than using the newly claimed numerical partition. It compares every358904 integer triple with odd h<=101,0<=e<=h,1<=R<=2h against the complete physical set. They agree exactly. It then verifies the last-step partition and nonlinear deletion formula on all324 accepted profiles.

Its gadget constants do not use the source binary recurrence. For an m-bit gadget word w, Q is1 plus the sum of2^(m-i-1) over base positions i with w_i=1. The Hamilton count is obtained by placing the increasing base before/after its apex: empty sides give the two boundary terms; with last left index i and first right index j, there is one assignment when j=i+1, or2^(i-j-1) assignments when i>j, subject to w_i=0,w_j=1. The m0 case is handled separately. These are literal slot formulas and remain efficient for the four larger recipes.

Every independently returned small recipe is constructed as an actual tournament by disjoint gadget blocks and the sole reverse seam. The C++ permutation counter verifies all324 reconstructed graphs, their three responses, start count1 and deletion count. The four large returned source recipes are also replayed with the independent gadget constants; no direct adjacency count is claimed on their larger graphs.

## 3. Replay and files

The source verifier and its three named proof dependencies were copied into `source/` before execution, because the verifier writes its receipt. Both normal and Python -O source replays reproduce the frozen receipt exactly. Original source files remain unchanged.

Independent reproduction from this audit directory:

    g++ -O3 -std=c++17 physical.cpp -o physical
    ./physical c1-inputs.csv physical-profiles.csv
    python3 check.py
    python3 -O check.py
    ./physical c1-inputs.csv physical-profiles-replay.csv recipes.csv

The Python checker reads source-normal.json only to obtain the four large already-returned recipes; its358904 small membership comparisons and physical reconstruction are independent of that source output. An assertion-enabled unoptimized C++ replay and the optimized replay agree, including the entire physical profile CSV and all324 returned graphs.

Frozen audit hashes:

- physical.cpp:e0a8d43345e1cf187d3291b1ac8d8db3097a785667aff76e7fc632b1cb89a116.
- check.py:df7bbef9e8f9cab402b7f3070f0dc604af05a35be40c1840a832a6a0d7107ea5.
- c1-inputs.csv:c26cec3efde1cde2acc99b3f5d86c074246cacaf33d3f823aba89dfd24a596e8.
- physical-profiles.csv:d3f2a5ca30fa54b8960b0d903869ffbf045f1b5a0e7ca099c80ac0b68b36fc3c.
- independent-normal.json and independent-optimized.json:fb8530ee4a751646fc352b7003720642a5d64e678df8bd216eb565748e4b5132.
- recipes.csv:66bdc6ab9cde4dd07d58735d5a6e7744d083ed361d4b097f3845347d6f5a771f.
- physical-recipes.txt and physical-debug.txt:de14669cdcbbcfb810cbc17a0a429427615bf82d556388353cdf47288a6649b4.
- source-normal.json and source-optimized.json:740db5c2caf7f0c1ea218195ba175d1d502622d76dc775d07483b1c4d2c113c9.

All-order validity comes from the inverse induction and the inherited physical grammar. These independent physical and arithmetic controls test boundaries and implementation without supplying an extrapolation argument.
