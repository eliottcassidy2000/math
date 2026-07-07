---
source: klein-2026-07-07-S163 (HYP-4871)
status: PROVED (three structure theorems) + census n=4..6. The meta-abstraction resolves
  crisply: the line graph at tiling level is the FOLDED m-CUBE, and its iso-quotient is
  ITSELF — every line is rigid (pairwise non-isomorphic), because a line's two tournaments
  agree exactly on the base path and a directed Hamiltonian path has no nontrivial
  automorphisms. The apparent chaos of the allocation tables is therefore the bare
  matching itself; its structure = a deterministic fiber-1 skeleton (ratios up to 68)
  over a quasi-random bulk (ratios 0.55–2.6).
tags:
  - tournaments
  - merged-metagraph
  - line-metagraph
  - folded-cube
  - rigidity
  - quasi-randomness
---

# The line metagraph is rigid

**klein-2026-07-07-S163.** Owner: abstract once more — the metagraph of the LINES — and
find the families in the chaos plus the pattern that generates them.

## The three theorems

- **L1 (folded cube).** Single-tile flips commute with the all-tile flip, which is free;
  so the set of 2^{m−1} lines with wiggly adjacency is exactly the **folded m-cube**
  FQ_m = Q_m/antipode: m-regular, vertex-transitive, the classical antipodal quotient.
  The blue lines (2^{(m+f)/2−1} of them, THM-643/T5) 2-color it.
- **L1′ (rigidity — the generating pattern).** A line's two tournaments **agree exactly on
  the n−1 base-path arcs and disagree on all m tiles**; hence any isomorphism of line-pairs
  maps agreement set to agreement set, i.e. base path to base path — and a directed
  Hamiltonian path is rigid (only the identity fixes it). So the only pair-isomorphism
  between lines is the identity: **all 2^{m−1} lines are pairwise non-isomorphic; the
  iso-quotient of FQ_m is FQ_m.** (Census confirms: L(n) = 4, 32, 512 with every orbit of
  size 1 at n = 4, 5, 6; blue/black constant on orbits trivially.)
- **L2 (no anti-blue; the reflection pairs black lines).** The grid reflection ρ acts on
  lines (it commutes with flip). Its fixed lines would satisfy ρt = t (blue) or ρt =
  flip(t) ("anti-blue") — but σ has f = ⌊(n−1)/2⌋ ≥ 1 fixed tiles, on which ρt = flip(t)
  demands bit = 1−bit: **anti-blue lines do not exist (n ≥ 3)**. So ρ fixes exactly the
  blue lines and acts freely on the black ones: **#black lines ≡ 0 (mod 2)** (✓ 2, 24,
  480), and the anti-isomorphism quotient has 2^{(m+f)/2−1} + (2^{m−1} − 2^{(m+f)/2−1})/2
  nodes.

## The consequence for the structure program

Because lines are rigid, there is NO hidden identification above the tiling level: S161's
allocation tables (lines per endpoint class-pair) ARE the line metagraph's complete
fibration data. The "chaos" decomposes exactly as the S161 completion predicted:

- **Deterministic skeleton:** every extreme assortativity ratio is a pair involving a
  **fiber-1 class** (H = |Aut|, e.g. transitive-like): its unique tiling has a FORCED
  partner — single lines with ratios up to 68.3 (n=6), pure rigidity, no randomness to
  average.
- **Quasi-random bulk:** pairs of large-fiber classes sit at actual/fiber-product ratios
  in [0.55, 2.6] (n=6 median 2.6 over 370 pairs; only 38/370 sub-random, min 0.644) —
  mild positive assortativity, bounded both ways.

**The generating pattern in one line:** *rigidity in the small fibers, mixing in the large
ones* — the same skeleton-plus-quasi-random completion as S161 (and as the 2-adic grading
of S162), now with the skeleton characterized (fiber-1 classes) rather than observed.

## Where the meta-abstraction has further leverage (recorded)

The move "quotient an abstraction and ask what survives" paid off precisely because the
QUOTIENT WAS TRIVIAL and proving that was the theorem. Same-shaped candidates: the
order-cell adjacency graph of the LRC movie modulo the dihedral time symmetries (opus's
engine — is the cell graph rigid for primitive families?); the anchor-window complexes of
THM-645 modulo the tent phase. In each case the first question is now sharpened: find the
agreement-set-type invariant that either rigidifies the quotient or names its collapse.

## Honest status
L1, L1′, L2 proved (elementary); census n=4..6 exact (n=7 line census unneeded for the
theorems — rigidity is proved for all n); assortativity numbers are n≤6 measurements.
Files: `04-computation/line_metagraph_klein_S163.py` (+ .out with the skeleton addendum).
Pointers: THM-643 (T5 blue count), HYP-4851 (allocation tables = the fibration), S161 C5
(now: skeleton characterized as fiber-1), FQ_m classical.
