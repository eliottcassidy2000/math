# Chain-run lifting strengthens the component theorem

**Later closure in this session:** [EXTENSION4.md](EXTENSION4.md) proves the universal quantitative and single-cycle inequalities. OPEN labels for those two targets below record the earlier research stage; Chan–Pak remains open.

**PROVED all orders:** the quantitative defect inequality for every disjoint union of chain graphs, with arbitrary tournaments inside both shores. **OPEN:** the single-cycle domination conjecture introduced below. This is a continuation of `PATH_COVER_INTERFACES.md`, not a replacement of its biclique proof.

## 1. Even-path inventory

For a tournament U with a fixed partition A⊔B, let P(U;A,B) count all directed simple paths with an even, positive number of vertices, beginning in A and ending in B. These paths need not span U. Closing any such path through a new vertex z with B→z→A gives a directed odd cycle containing z. Inherited OCF therefore gives

    Δ=H(U+z)−H(U) ≥ 2P(U;A,B).                           (1)

This lower bound uses only single-cycle packings and leaves every unused vertex outside all cycles. It does not need any complementary H weights.

## 2. Chain-run lifting lemma

Suppose the bipartite incidence graph G of A→B edges is a chain graph. Choose an interleaved linear order of A⊔B in which an A vertex precedes a B vertex exactly when they are adjacent in G. Let E be the number of even-size nonempty vertex subsets whose first selected vertex is A and last selected vertex is B in this order.

Then, regardless of U's internal orientations in A and B,

    P(U;A,B) ≥ E = (d(G)−1)/2.                            (2)

**Proof of the inequality.** For each counted subset, read its selected vertices in the interleaved order and split them into maximal monochromatic runs. Each run lies entirely in A or entirely in B; choose a Hamilton path of that induced tournament canonically. Concatenate those paths in run order. Every edge joining consecutive runs points forward: both endpoints have different colours, and every vertex of the earlier run precedes every vertex of the next run in the chosen chain representation. The resulting directed path starts in A, ends in B, and uses an even number of vertices. Distinct subsets give distinct paths by their vertex sets.

**Proof of the identity.** Replace U by the completely transitive tournament in the chosen chain order. Every even A-to-B directed path is exactly one of these selected subsets. Every odd directed cycle of U+z contains z, so no two can coexist in a packing; its OCF count is exactly 1+2E. The inherited chain adapter says this count equals d(G). Therefore E=(d−1)/2. This uses only the already established exact transitive adapter, not the conjectural quantitative inequality.

Combining (1) and (2) proves the conjectured Δ≥d−1 for all chain G, including arbitrary nontransitive choices inside either shore. The earlier zero-defect equality theorem alone did not imply this strengthened assertion.

## 3. Disjoint chain-component theorem

Suppose G is a disjoint union of chain graphs G₁,…,G_k, and isolated vertices. Set E_i=(d(G_i)−1)/2. Then

    P(U;A,B) ≥ Σ_{∅≠J⊆[k]} |J|! ∏_{i∈J} E_i,
    Δ ≥ 2 Σ_{∅≠J⊆[k]} |J|! ∏_{i∈J} E_i
      ≥ d(G)−1.                                         (3)

For each component choose one of the canonically lifted paths in (2). Each path has even order, starts in A_i and ends in B_i. No graph edges join distinct components, so all B_i→A_j edges exist for i≠j. Any ordering of selected component paths concatenates to a directed even A-to-B path. The resulting path recovers the component order and each component path, so this is an injection. Product expansion and j!≥2^{j−1} finish the proof exactly as in the biclique case.

The same strict factorial surplus from three nontrivial components holds. The argument even proves the stronger single-cycle inequality (1)→(3). Complete bipartite graphs are only the one-run A^pB^q special case of this lifting.

If this class is merely induced on W⊂V(U), the OCF construction can be combined with arbitrary packings on U−W. Therefore the general induced witness bound is

    Δ ≥ 2 H(U−W) Σ_{∅≠J⊆[k]} |J|! ∏ E_i.               (4)

This is useful for any tournament whose marked graph contains a large induced chain-cluster witness. It does not assert d(G[W])≥d(G), and so does not prove the full inequality from a small witness alone.

## 4. Stronger remaining conjecture and its exact missing step

A cleaner possible target than reserving old Hamilton paths is

    d(G) ≤ 1+2P(U;A,B)                                  (SC)

for every bipartite partition of every tournament U. If true, (SC) implies the Session260 quantitative defect conjecture directly, without any complementary cycle packings. Sections 2–3 prove (SC) for disjoint unions of chain graphs.

The missing universal map is now precise: inject half the nontrivial dominating configurations, with an explicit pairing, into simple even directed A-to-B paths. Arbitrary bipartite graphs need not admit the chain order used above. Taking a biclique cover instead of actual components also fails: a cross-component graph edge can reverse the seam B_i→A_j. The existence of a Hamilton path on a chosen subset does not guarantee both required endpoint colours. These are concrete lost coordinates, not a proof that the universal map is impossible.

A 500-case seeded scout across rooted orders 4..8 found no violation. This is **FINITE-EXACT evidence only**, deliberately a cheap hostility test rather than a large range extension. No universal claim is promoted from it. The path-count DP is in `chain_component_probe.py`, and `single_cycle_scout.json` retains the scope and seed.

## 5. Third random pull and transfer ledger

After the component theorem and the initial idea of chain-run lifting, another uncontrolled pull used seed **2611706**, the same 4,130-path canonical-theorem pool as Wave2 at commit f92e070445e43c59f1252f2979db564df88d7b60, and sample size two. Actual paths:

- `01-canon/theorems/THM-3147-length-singleton-endpoint-jet-facet-observer.md`
- `01-canon/theorems/THM-2613-canonical-root-diagonal-opposite-shift-section.md`

The first source distinguishes a compressed observer from a positive transport. Map tested: count a path by its number of runs, versus retain its run colours and endpoint constraints. Preserved by the full run record: legal seams. Destroyed by the run count alone: which colours meet at a seam. Cheapest test: a two-run AB block works, BA does not close through z. Outcome: retain the literal run-level lift, not just a length statistic. No symmetric-function positivity is imported.

The second source constructs a canonical section only after retaining physical origins. Map tested: choose one Hamilton path inside every monochromatic run using fixed vertex labels. Preserved: vertex subset and boundary colours. Lost: naturality under arbitrary relabelling, which is not needed for an injection. Cheapest test: changing the chosen path inside a run leaves every cross-run edge legal. Outcome: a deterministic local section suffices; it is not claimed canonical in an isomorphism-invariant sense. No arithmetic phase or LRC bound is transferred.

These pulls did not cause the earlier chain insight; they sharpened its proof's loss ledger and clarified which quotient is safe. The record distinguishes chronological inspiration from subsequent comparison.

## 6. Reproduction and consequence boundary

Run `python3 05-knowledge/results/session261_creative_chanpak_20261111/bridges/chain_component_probe.py`. It checks 84 actual constructions with seven specified word-component shapes, arbitrary seeded internal-shore orientations, their component path counts, and the factorial lower bound. The separate scout is labelled conjectural. The proved result is structural and unbounded, not the finite range.

This extension controls positive defect on a larger graph class; it does not cancel the defect or provide all desired Hamilton counts. The two-apex family in the companion note is a different operation and retains its own path-cover memory. Chan–Pak remains OPEN.
