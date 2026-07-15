---
id: THM-793
title: The Mode-B line tower — the three objects (tilings, nodes, lines) and their five maps commute down the strip-legs+apex projection; the axis current lives on the fiber; SC = A000570 so merged nodes = (A000568+A000570)/2; the transitive node is a fixed point of the descent kernel; self-line census with the parity-forced blue-odd vanishing
status: PROVED (the tower: fiber count, both commutations, induced-subtournament, leg-measurability — verified exactly at n = 5→3, 6→4, 7→5; the one-line proofs below) + CENSUS-EXACT n = 3..7 (SC ≡ A000570; self-lines; merged self-loops; descent kernels) + CONJECTURE (black self-lines = SC/2 at odd n ≥ 5; two further data points needed)
source: opus-2026-07-14-S306 (owner directive: track nodes, tilings, blue/black edges and all their relations as recursively as possible)
depends_on:
  - THM-790   # the leg law (the fiber-measurability input)
  - THM-787   # the flow tables
related: [HYP-6860, everything-is-the-triangle (Mode B), merged-metagraph-invariants]
verification: 04-computation/mode_b_line_tower_census_opus_S306.py
  (+ 05-knowledge/results/mode_b_line_tower_census_opus_S306.out)
---

# THM-793 — the Mode-B line tower

**The three objects and five maps.** Tilings T_n (2^m of them); nodes = iso
classes (A000568(n)), merged nodes = classes mod reversal; lines = the d=m
edges (t ↔ t̄), coloured blue/black by grid symmetry. The maps: π (tiling →
node), κ (the flip, tiling → tiling), σ (the grid reflection), the line-to-
node-pair map, and now the DESCENT

> **p: T_n → T_{n−2}, forget the two legs and the apex** (the 2n−5 tiles
> (x,1), (n,y), (n,1)); the interior tiles are exactly the staircase of
> {2, …, n−1}.

## (1) The tower (PROVED; each part one line, verified exactly at n = 5,6,7)

- **(a) Fiber:** p is exactly 2^{2n−5}-to-1 (the forgotten bits are free).
- **(b) κ-equivariance:** p∘κ_n = κ_{n−2}∘p (flipping all tiles restricts).
- **(c) σ-equivariance:** p∘σ_n = σ_{n−2}∘p (σ maps interior to interior:
  σ(x,y) has coordinates n+1−y, n+1−x, interior iff (x,y) is). Hence blue
  descends to blue, and on blue fibers p is 2^{n−2}-to-1
  ((m+f)/2 − (m′+f′)/2 = n−2).
- **(d) Tournament meaning:** T(p(t)) = the induced subtournament of T(t) on
  the interior vertices {2, …, n−1} (path arcs and tiles both restrict).
- **(e) THE CURRENT LIVES ON THE FIBER:** Δx_n(t) = 8(e₁ − e_n) (the leg law)
  depends ONLY on the forgotten coordinates. Equivalently the layer drop GF
  factors as **GF_n(z) = 2^{m(n−2)} · (1+z)^{n−3}(1+z⁻¹)^{n−3}(z+z⁻¹)** — the
  2-power IS the child tiling count — and the blue GF as
  **#blue(n−2) · (z+z⁻¹)^{n−2}**.

So the whole line-layer of level n is a bundle over the line-layer of level
n−2: base = the child tilings with their own κ/σ/line structure intact,
fiber = the leg+apex data carrying ALL the axis current. Iterating p strips
the staircase two vertices at a time — the Mode-B recursion of
everything-is-the-triangle, now carrying the metagraph's edge structure with
it. (Mode A — forget one leg row — is also κ-equivariant with fiber 2^{n−2},
but does NOT commute with σ; the blue/black structure only descends the
Mode-B tower.)

## (2) The census (EXACT, n = 3..7)

| n | tilings | nodes | SC | merged | lines | blue | black | self-blue | self-black | merged-self-blue | merged-self-black |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 3 | 2 | 2 | 2 | 2 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |
| 4 | 8 | 4 | 2 | 3 | 4 | 2 | 2 | 1 | 0 | 1 | 0 |
| 5 | 64 | 12 | 8 | 10 | 32 | 8 | 24 | 0 | 4 | 0 | 4 |
| 6 | 1024 | 56 | 12 | 34 | 512 | 32 | 480 | 2 | 6 | 2 | 24 |
| 7 | 32768 | 456 | 88 | 272 | 16384 | 256 | 16128 | 0 | 44 | 0 | 114 |

- **SC(n) = A000570(n)** (self-converse tournaments) exactly — so
  **merged nodes = (A000568 + A000570)/2**: 2, 3, 10, 34, 272, and the n=8
  prediction 3,528 (with SC(8) = 176).
- Closed counting laws: tilings 2^m; lines 2^{m−1}; blue lines 2^{(m+f)/2−1}
  (the half-tiling model); fiber(node)·|Aut| = H; Σ fibers = 2^m.
- **Blue self-lines vanish at odd n — PROVED** (a self-line has Δx = 0; blue
  lines are never level at odd n by the parity law). At even n they exist
  (1 at n=4, 2 at n=6).
- Self-lines are level (Δx = 0) by the leg law — self-lines ⊆ {e₁ = e_n} —
  a proved necessary condition on the fiber coordinate.
- **Observed coincidence (CONJECTURE):** black self-lines = SC(n)/2 at odd
  n = 5, 7 (4 = 8/2, 44 = 88/2), and self-blue + self-black = SC/2 fails at
  even n (8 ≠ 6 at n=6) — the odd-n law needs n = 9 or a proof; the natural
  route: a self-line's class contains t and t̄, and the flip's score action
  (reversal + leg defect) pairs such classes with self-converse structure.

## (3) The descent kernel (the node-level shadow of p)

p induces a Markov kernel K_n: node(n) → distribution over node(n−2)
(project the fiber). Exact facts (n = 5→3, 6→4, 7→5):

- **The transitive node is a FIXED POINT: K(transitive) = point mass on the
  child transitive** (the all-forward tiling restricts to the all-forward
  tiling — its fiber is a single tiling at every level). The one node the
  owner's program started from persists identically down the whole tower.
- Kernel support sizes: max 2 (of 2), 4 (of 4), 11 (of 12) — at n=7 NO class
  sees all twelve children (support ≤ 11); full-support rows: 8/12, 31/56,
  0/456. The kernel is genuinely sparse at the top and localizes with n —
  the node-level recursion is a proper tree-like refinement, not mixing.

## (4) The relations, assembled (the owner's requested dictionary)

```text
            kappa (flip; carries the current 8(e1-en))
   T_n  <-------------------------------------------->  T_n
    |  \___ sigma (blue = fixed sector)                   |
    | pi                                                  | pi
    v                                                     v
  nodes  <---- lines (blue/black; self-lines level) ---> nodes
    |                                                     |
    | reversal quotient                                   |
    v                                                     v
 merged nodes  <---- merged lines (self-loops absorb reversal pairs) ----
        AND every layer of this diagram maps onto the (n-2) diagram via p,
        equivariantly for kappa and sigma, with the axis current carried
        entirely by the p-fiber (legs + apex).
```

Next checks named: n=9 for the SC/2 self-line law (needs the invariant-
certified classifier one level up, 2^28 tilings — vectorizable); the kernel's
axis-compatibility statistics; the quarter model as the σ-quotient of the
fiber itself.
