# Independent audit: Session203 rooted profile order bound

**Mathematical correctness: PASS. No required corrections.** No theorem status promotion or publication was performed.

## Frozen object

Before analysis, copied PROOF.md, verify.py, NORMAL.json, OPTIMIZED.json and SHA256SUMS into this directory. The proof hash is `f6366dd37201d6c4e655ed1bf9a25f1a80d440270bcf3bd8fee59d768ca51810`; verifier hash is `edca19d57370e8937f913221a155d9c30eaa1dd1d41375f1424e38ce38187e8d`. Both original receipt hashes are `a5adac75bcde063212467d9a446c89b7eaeb22cf5d053519e3f67c0bc5a8fee6`. These agree with the supplied manifest. Audit applies exactly to that proof and definition of L.

## Proof audit

1. Rooted Hamilton existence is valid. If the visited set is proper, strong connectivity gives an outgoing arc. When the last visited vertex loses to the new vertex but some visited vertex beats it, a consecutive beat/lose transition supplies an insertion position without changing the starting vertex. Reversal gives specified ending vertices. Distinct starting vertices imply H(D)>=|D|, including singleton D.
2. Formula (3) holds when U-v is nonempty: component order forces the first post-root vertex to be in D1; the rest of each component path is independent. The proof uses it only for n>=3, so no empty-deletion exception is overlooked.
3. Strong deletion case: every root out-neighbor in U-v contributes at least one start, giving outdegree(v)<=s. Subsets of in-neighbors define different designated vertex sets; arbitrary Hamilton paths within each selected set and its complement form valid covers. Empty selected set and empty complement have the specified conventions. Thus L>=2^indegree(v).
4. Split deletion case: a terminal component D must contain w->v; otherwise D has no outgoing arc in U. All earlier components dominate D. Concatenating any selected earlier vertices, a fixed Hamilton path of D ending at w, then v is legal, and the complementary earlier vertices form Q. This is genuinely injective by the P vertex set. Because there are at least two components, H(D) is a factor of q, and a positive integer leading sum gives s>=q>=H(D)>=|D|. The exponent direction is correct.
5. Weak-tail identity is a bijection, not merely an inequality. P ending at v cannot visit a later component and return; Q contains all later components in forced order. Removing their suffix leaves a unique internal cover, including the case of empty internal Q. Conversely every internal cover extends in exactly q ways. Hence h,s,L,B all multiply by q, so deleting singleton components preserves BOTH coordinates.
6. Every retained non-singleton strong component has size and Hamilton count at least three. Sum of their sizes is bounded by the product q of their Hamilton counts. The root singleton case is separate. In the non-singleton root case, s0+q<=qs0+1 and floor(log2 B0)<=floor(log2(qB0))-1 for q>=3 together yield precisely the claimed order cap; no off-by-one defect occurs. The integer logarithm is legitimate because B>=L>=1.
7. Singleton profile is (h,s,L,B)=(1,1,1,1); exponent is -1 and poses no problem. The positive-start rooted edge is the source-rooted edge, normalizing to the singleton. The sink-rooted edge has s=0 and is excluded. There is no strong two-vertex tournament.
8. The target consequence substitutes s=r+3j<=26 and positive B=X-2j<=X, giving cap 27+floor(log2 X). It is a per-base finite-universe guarantee only. It does not establish occupancy, F(8,137), a uniform finite check for all X, or an efficient algorithm.

## Reproduction and independent checks

Commands (run in this directory):

    python3 verify.py > RERUN.json
    python3 -O verify.py > RERUN_OPTIMIZED.json
    cmp NORMAL.json RERUN.json
    cmp NORMAL.json RERUN_OPTIMIZED.json
    python3 independent.py > INDEPENDENT.json

Both original verifier reruns matched the supplied NORMAL.json byte-for-byte. This includes all 1,099 root-zero labeled tournaments through order 5, 2,952 literal injected covers, and four explicit weak-tail examples of orders 9–12.

The independent script does not import the candidate verifier. It enumerates full permutations rather than subset dynamic programming: a cover has a unique full permutation P followed by Q, split immediately after the root, and requires every consecutive edge except possibly the split edge. Its SCC routine uses individual DFS reachability, rather than bitset transitive closure. It checks every root of all labeled tournaments through order 5 and 128 seeded distinct graphs at each of orders 6 and 7 (seed 203137). All 7,069 rooted graph checks passed, including 5,230 positive-root cases and 4,083 strong-root cases. Checked full profile preservation, component product, normalized cap and strong exponential inequality. Finite tests support, but do not replace, the proof.

## Scope and checkout preservation

Read repository AGENTS.md and research protocol. The shared checkout is dirty and substantially behind its cached origin/main. It was not edited, synchronized, staged, committed, or published. This audit deliberately concerns the frozen candidate rather than claiming repository freshness or theorem-number availability. All intentional audit artifacts are confined to this scratch directory. No mathematical blocker remains.
