# Arborescences are the determinant-shadow of Hamiltonian paths — and switching averages both to their tree-count

**death-star-2026-07-20-S61c** (HYP-8315; owner: think prior work on arborescences in
tournaments and how they relate to Hamiltonian paths/cycles — get creative). The one-line
thesis: **a Hamiltonian path is the tallest possible arborescence, an arborescence is a
determinant where a Hamiltonian path is a permanent, and switching (the even-graph structure
of S61b) averages *any* family of spanning trees to its own cardinality (THM-1460).**

## 1. The height axis: stars ↔ arborescences ↔ Hamiltonian paths

An out-arborescence rooted at r is a spanning tree of K_n with every edge oriented away from
r, using the tournament's arcs. Two extremes:

- **Flattest** = the **star** r → everyone: exists iff r beats all (r is the dominant vertex).
  A tournament has ≤ 1 dominant vertex, so stars are rare (depth 1, maximally branchy).
- **Tallest** = a **Hamiltonian path** from r (depth n−1, no branching).

So **Hamiltonian paths are exactly the height-(n−1) arborescences**, and the whole
arborescence count a_r(T) interpolates between the rare flat stars and the tall paths. Rédei
counts only the tallest ones (and finds an odd number); the Matrix-Tree theorem counts *all*
of them as a determinant. Verified n=3,4: a_r(T) ≥ (#Ham paths from r), always, with the gap =
the branchy (non-path) arborescences.

## 2. Determinant vs permanent — the easy shadow of a hard count

**Arborescences are a determinant.** Tutte's Matrix-Tree theorem: a_r(T) = the r-th principal
cofactor of the Laplacian L = D_in − A. Polynomial-time, exact, verified against brute force.

**Hamiltonian paths are permanent-hard.** Counting them is #P-complete in general; Rédei's
theorem (H odd) is the only easy global handle. The det/perm gap *is* the arborescence/path
gap: relaxing "path" (a tree of max-degree 2, permanent-flavored) to "tree" (any spanning
tree, determinant-flavored) is exactly what turns the hard count into the Matrix-Tree
determinant. Arborescences are the **determinantal shadow** of Hamiltonian paths.

## 3. THM-1460: switching averages every spanning-tree family to its own size

The S61b machinery (H is odd-valued; the even-graph/switching projection washes it to a
constant) is *not special to paths*. THM-1460 (proved): for **any** family 𝓕 of oriented
spanning trees of K_n and any switching class 𝒞,

  **Σ_{T ∈ 𝒞} #{S ∈ 𝓕 : S ⊆ T} = |𝓕| .**

because the cut space, restricted to any spanning tree's edges, is a bijection — each oriented
spanning tree is realized in exactly one tournament per switching class. Instances, all
verified n=3,4,5:

- Hamiltonian paths (|𝓕| = n!) → class-sum **n!** (my S61b switching-sum, now THM-1460);
- all out-arborescences (|𝓕| = n^{n−1}, Cayley rooted trees) → class-sum **n^{n−1}**;
- arborescences rooted at fixed r (n^{n−2}) → **n^{n−2}**; Ham paths from r ((n−1)!) → **(n−1)!**.

So the odd, hard, path count and the easy, determinantal, arborescence count are **the same
theorem at two heights of the tree** — both average, over a switching class, to their pure
combinatorial family size. The tournament-specific variation (which T has how many) lives
entirely in the fiber; the even-graph base (S61b) sees only |𝓕|. Globally,
Σ_T H = n!·2^{C(n−1,2)} and Σ_T a = n^{n−1}·2^{C(n−1,2)} — the *same* even-graph factor
2^{C(n−1,2)}, differing only by paths (n!) vs all trees (n^{n−1}).

## 4. The QR closed form: Paley arborescences from the Gauss-sum spectrum

For the **Paley tournament** on q vertices (q ≡ 3 mod 4 prime) the arborescence count is
*exactly solvable* from the quadratic-residue spectrum. Its adjacency eigenvalues are (q−1)/2
and (−1 ± √−q)/2; the Laplacian L = ((q−1)/2)I − A then has nonzero eigenvalues (q ∓ √−q)/2
each of multiplicity (q−1)/2, so Matrix-Tree gives (verified q = 3, 7, 11):

  **a_r(Paley_q) = (1/q) · [ q(q+1)/4 ]^{(q−1)/2} .**

(q=3 → 1, q=7 → 392, q=11 → 3 557 763.) The arborescence count of a Paley tournament is a pure
power of q(q+1)/4 — the QR structure that runs through the whole repo (QR_p, LRC, Paley
conference matrices) also fixes its spanning-tree count in closed form. Hamiltonian paths of
Paley have no such formula (they are the hard, tall shadow); the arborescences do (the easy,
spectral one).

## 5. BEST: arborescences bridge Hamiltonian (vertex/odd) and Eulerian (arc/even)

The de Bruijn–van Aardenne-Ehrenfest–Smith–Tutte (BEST) theorem counts **Eulerian circuits**
of a connected Eulerian digraph as ec(T) = a_r(T)·∏_v (deg⁺(v) − 1)! — *through the
arborescence count a_r*. A tournament is Eulerian iff **regular** (n odd, every out-degree
(n−1)/2), so for a regular tournament

  ec(T) = a_r(T) · ( ((n−3)/2)! )^n ,

and for Paley the a_r is the closed form of §4 (verified: 3-cycle ec = 1; Paley₇ ec = 392·2⁷ =
50 176). This is the bridge S61b was pointing at: **Hamiltonian objects live on vertices
(Rédei-odd), Eulerian objects live on arcs (need the even/regular condition), and
arborescences — the Matrix-Tree determinant — are the common currency that counts the
Eulerian ones and is dominated by the Hamiltonian ones.** Arc-world (even, Euler), vertex-world
(odd, Hamilton), joined by the tree-world (determinant, arborescence).

## 6. Creative punchline and the forward question

The tree is the pivot. Push it *tall* and it is a Hamiltonian path (odd count, permanent-hard,
Rédei); count *all* of it and it is a determinant (Matrix-Tree, Paley-solvable); walk its
*arcs* and BEST turns it into Eulerian circuits (the even/regular world). Switching — the
even-graph quotient — cannot tell these apart: it averages every spanning-tree family to its
size (THM-1460). So the question S61b posed sharpens: the odd, tournament-specific content
sits in the *fiber over the even-graph base*, and the **cheapest computable handle on it is the
Matrix-Tree determinant** a_r(T) — the largest quantity that (i) dominates H(T), (ii) is
polynomial-time, and (iii) has the Paley closed form. Concretely: **does the arborescence
vector (a_1(T), …, a_n(T)), or the Laplacian spectrum, control the Hamiltonian-path count
H(T) enough to bound the H-spectrum gaps {7, 21}?** That is the high-leverage next step — the
determinant shadow constraining the permanent.

## 7. Honesty and credit

§3 (THM-1460) is proved (it subsumes my S61b Σ H = n!, filed as THM-1445 then renumbered on the opus/kp collision); the four instances are verified n=3,4,5.
§4 Paley closed form verified q=3,7,11 (and derived from the standard Paley spectrum). §5 BEST
is classical, applied here to regular tournaments and checked on the 3-cycle. §1–2 are the
framing (det/perm, height axis) — standard facts assembled into the picture, not new theorems.
§6 is a posed target. Credit: kp THM-1390 surveyed matroid/Tutte/arborescence but did not
count; kp THM-1430/474 (cut space simply transitive on spanning-tree orientations — the engine
of THM-1460); S61b (the odd/even base/fiber picture this extends).

## Cross-links
THM-1460 (né THM-1445, ceded on collision) · S61b odd/even reflection · Matrix-Tree/Tutte · Cayley n^{n−2} · Paley/QR
(the repo's QR_p thread) · BEST theorem · Rédei · kp THM-1390 (arborescence survey) / THM-1430.
