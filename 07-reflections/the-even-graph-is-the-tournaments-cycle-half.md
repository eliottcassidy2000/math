# The Even Graph Is the Tournament's Cycle Half — One Cut⊕Cycle Object Through Every Lens

*mac-mini-2026-06-22-S38. A synthesis of the even-graph ↔ tournament lenses the project has
accumulated (cycle-space bijection, H=I(Ω,2), the E_n/G_n metagraph duality, nowhere-zero flows,
spectral, perfect-graph), grounded in one frame: the **Cut ⊕ Cycle decomposition** of K_n's edge
space. Owner's prompt: "really understand how their structure shines through different frames."*

## The one frame everything else is a view of

The edge space of K_n over GF(2) splits canonically (graph homology):

> **E(K_n) = Cut(n−1) ⊕ Cycle(C(n−1,2)).**

Verified dims (S38): |E(K_n)| = (n−1) + C(n−1,2) for all n. Fixing the base Hamiltonian path P_0
n→n−1→…→1 **chooses the cut summand**: the n−1 base-path arcs are a spanning tree = the cut/score
space; the C(n−1,2) non-base arcs (**the tiles**) are the cycle space. So:

| | dimension | what it carries | object |
|---|---|---|---|
| **Cut** | n−1 | score sequence / hierarchy / Rédei transitive backbone | the base path P_0 |
| **Cycle** | C(n−1,2)=m | odd cycles / rotation / H | **the even graph** |

The 2^m tilings = 2^m even graphs (cycle-space bijection: even graph = XOR of the fundamental
cycles of the chosen tiles, cycle_space_bijection_s20ge). **The even graph is not "related to" the
tournament — it is precisely the tournament's cycle half, read off after the score (cut half) is
quotiented out.** This is why CLAUDE.md's "base-path = cut, wiggly = cycle" is the whole story.

## Each lens = the same Cut⊕Cycle object from one angle

**1. Score / Rédei (the CUT side).** The score sequence is a cut-space functional (out-degrees =
coboundary of a potential). Rédei's "H is odd" is a cycle-side statement (H counts cycle structure);
the transitive tournament is pure cut (cycle part = 0, H = 1). The hierarchy lives entirely in Cut.

**2. OCF / conflict graph — H = I(Ω,2) (the CYCLE side, VERIFIED).** H(T) = I(Ω(T),2), Ω = the
odd-cycle conflict graph (independent set = vertex-disjoint odd-cycle collection). Verified on **all
1096 tournaments n=3,4,5**. H reads ONLY the cycle half: it is the independence polynomial, at the
"two-graph" value 2, of the odd cycles. The even graph (cycle half) and Ω are the same information.

**3. The metagraph duality E_n vs G_n.** Both are Q_m/S_n, but E_n quotients the **cycle half**
(even-graph iso classes: V=2,3,7,16,54 for n=3..7) and G_n quotients the **whole** tournament. The
80/10/8/2 tile-flip correlation (even-graphs-as-first-class) measures how much a cycle-side flip also
moves the cut-side class. E_n is the DUAL of G_n exactly because it is G_n with the cut forgotten.

**4. Nowhere-zero flows / Tutte tension (the duality made literal).** Over GF(2): an even graph IS a
flow (degree-even = cycle space = flow space); the score potential IS a tension (cut = coboundary).
Cut⊕Cycle = Tension⊕Flow. So lrc-as-nowhere-zero-flows (oracle S537o) is saying the LRC obstruction
lives on the **flow = cycle = even-graph** side — the same side as H. The apex-7 is a cycle-side
phenomenon, not a cut-side one.

**5. Spectral.** The Laplacian eigenbasis respects Cut⊕Cycle; "additive energy = spectral 4th moment"
(my S29) and "eigenvalues of the merged metagraph" are cycle-side spectra. The score/cut side is the
trivial/Perron part.

**6. Perfect-graph — the apex-7, where the two halves touch (C_5 ↔ K_3).** This is the sharpest
cross-frame image, and it is now two-sided:
- **Cycle side:** C_5 = XOR of 3 pairwise-vertex-conflicting triangles (S37, all 5 decompositions in
  K_5). E_7 loses chordality exactly by gaining C_5 odd holes at the apex prime.
- **Conflict side:** H=7 ⟺ Ω=K_3 (3 pairwise-conflicting triangles), and kps's **THM-200** shows
  Ω=K_3 is FORBIDDEN because three pairwise-sharing triangles force a common vertex → a fifth vertex
  carrying a directed **C_5** → a fourth odd cycle. **Forbidding the K_3 (conflict side) literally
  manufactures a C_5 (cycle side).** The pentagon and the triangle are the two faces of the single
  apex-7 obstruction, glued along the Cut⊕Cycle seam. H = I(Ω,2): the clique value I(K_3,2)=7 is the
  unique forbidden one (kps S31g), and its forbiddenness is the C_5.

## What this buys

- **Where does a problem live?** Ask: cut or cycle. Score, hierarchy, transitivity → cut. H, odd
  cycles, even graphs, LRC apex-7, the {7,21} gaps → cycle. The two are orthogonal; don't look for H
  in the score or for the apex-7 in the hierarchy.
- **The even graph carries H — but only with the cut fixed.** In the tiling model (base path = the
  standard hierarchy, i.e. cut held fixed), H is a function of the even graph: the cycle half
  determines the tournament, so H = f(even graph). CAVEAT (discipline): this is the *labeled,
  fixed-base-path* level. The iso quotient E_n uses a *different* S_n action than G_n: the class
  counts already differ (V(G_n)=2,4,12,56,456 vs V(E_n)=2,3,7,16,54 for n=3..7), and the bridge matrix
  B[tourn,even] is full-rank but not square — so the cycle-only quotient does not separate tournament
  classes the way G_n does. H is therefore a G_n invariant; the cycle half determines H once you fix
  the cut, but it does not survive the cycle-only quotient to become an E_n invariant. VERIFIED (S38):
  the 64 tilings at n=5 fall into exactly 7 even-graph iso classes (= V(E_5)), and **5 of the 7 carry
  more than one H value** — one |E|=4 class spans H ∈ {5,9,11,13,15}. So forgetting the cut genuinely
  scrambles H; the even graph is the cycle half, but H needs the cut fixed to be read off it.
- **The apex-7 is one object.** {7,21} forbidden-H, E_7 non-chordality (C_5 holes), H=7's K_3, and
  (via flow=cycle) the LRC(14) apex are not four coincidences at 7 — they are the Cut⊕Cycle seam
  showing the same odd obstruction in cut-quotient (E_7), conflict (K_3), value (7=I(K_3,2)), and
  flow (LRC) coordinates.

Verified core (S38): the Cut⊕Cycle dims; H=I(Ω,2) on 1096 tournaments. Synthesis links: §3–§6 read
existing results (even-graphs-as-first-class, kps S31g/THM-200, oracle S537o, my S37/S29) through the
one frame. Related: [[seven-is-the-unique-forbidden-clique-value]], [[even-graphs-as-first-class]],
[[even-graphs-through-the-metagraph]], HYP-2880 (C_5=K_3), HYP-2852 (additive energy=spectral).
