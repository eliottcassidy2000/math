---
source: klein-2026-07-07-S161 (HYP-4851)
status: PROVED (five parity/structure theorems, one-line proofs from Rédei + the grid
  involution) + COMPLETE CENSUS to n=7 (0 violations) + new conjectures (even-n blue
  self-loops; pure-blue never adjacent to pure-blue for n>=4; transitive always pure-blue).
  Headline micro-theorem: every self-complementary tournament has an ODD number of
  Aut-orbits of anti-reversible Hamiltonian paths — a Rédei-type refinement.
tags:
  - tournaments
  - merged-metagraph
  - blue-black-lines
  - redei
  - parity
  - tiling-model
---

# The blue/black line structure IS Rédei parity

**klein-2026-07-07-S161.** Owner: even/odd duality of the blue/black lines — line counts,
node-type counts, allocation vs n, toward complete structural knowledge of the tiling
fibration. Extends klein-S75's line accounting (n ≤ 6) with proofs and the n=7 layer.
Concurrent: mac-mini-S46 claimed THM-643 for the canon statement minutes after my claim —
this note supplies proofs + data; coordinate there.

## The five theorems (proofs complete)

Let ρ = the grid reflection on tilings ((x,y) ↦ (n+1−y, n+1−x), bits carried), and note
`ρ(t)'s tournament = (relabel v ↦ n+1−v applied to T(t))^op`, so ρ maps the fiber of class
[T] bijectively onto the fiber of [T^op]. Let flip = complement every tile (the line
partner), which commutes with ρ and is fixed-point-free.

- **T1 (odd fibers).** Every unmerged class fiber is odd: fiber × |Aut| = H (orbit-
  stabilizer, LEM-003), H is odd (Rédei), |Aut| is odd (classical). ∎
- **T2 (SC ⟺ gridsym present).** If [T] is SC, its fiber is ρ-stable, and an involution
  on an odd set has a fixed point: some tiling is gridsym — so SC nodes are pure-blue or
  MIXED, never pure black. If [T] is not SC, a gridsym tiling would force [T] = [T^op]:
  zero gridsym — NS-merged nodes are ALWAYS pure black. (CLAUDE.md's "verified n=3..7"
  claim, now proved for all n.) ∎
- **T3 (the owner's parity slogan).** At an SC node: #gridsym = bx + 2·bs ≡ fiber ≡ 1
  (mod 2) ⟹ blue-cross degree bx is ODD; #non-gridsym = kx + 2·ks ≡ 0 ⟹ black-cross kx
  EVEN. At an NS node: blue absent, fiber = 2|C| even. "Blues contribute odd amounts,
  blacks even" — at every self-complementary node. ∎
- **T4 (tripartite laws).** Blue endpoints are gridsym ⟹ their classes contain gridsym ⟹
  blue never touches pure-black. Black endpoints are non-gridsym ⟹ black never touches
  pure-blue. ∎
- **T5 (blue count).** #gridsym tilings = 2^{(m+f)/2} (f = ⌊(n−1)/2⌋ σ-fixed tiles);
  flip preserves gridsym and is free ⟹ **#blue lines = 2^{(m+f)/2 − 1}**, #black =
  2^{m−1} − 2^{(m+f)/2−1}. ∎

**Headline corollary (a Rédei-type refinement).** A gridsym tiling of class [T] is an
Aut(T)-orbit of Hamiltonian paths reversed by an anti-automorphism. T2+T3 say: *every
self-complementary tournament has an ODD number of Aut-orbits of anti-reversible
Hamiltonian paths* (in particular ≥ 1); non-self-complementary tournaments have none.
Rédei counts all paths odd; this refines the count on the anti-symmetric locus.

## The census (n = 3..7; 0 violations of T1–T5)

| n | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|
| blue lines | 1 | 2 | 8 | 32 | 256 |
| black lines | 0 | 2 | 24 | 480 | 16128 |
| pure-blue | 2 | 1 | 3 | 2 | 4 |
| mixed | 0 | 1 | 5 | 10 | 84 |
| pure-black | 0 | 1 | 2 | 22 | 184 |
| blue SELF (all at mixed) | 0 | 1 | 0 | 2 | 0 |
| black SELF (MX, PB) | 0 | (0,0) | (4,0) | (2,24) | (18,96) |

n=7 allocation: blue-cross = 250 (MX,MX) + 6 (MX,PBlue); black-cross = 858 (MX,MX) +
5044 (MX,PB) + 10112 (PB,PB). Fibers: pure-blue tiny (1..3), mixed 5..159, pure-black
(merged) 2..306; totals 6 + 7302 + 25460 = 2^15 ✓.

## New conjectures (the niche insights, each checkable at n=8)

- **C1 (even-n blue self-loops).** Blue self-loops exist ONLY at even n (0,1,0,2,0 for
  n=3..7): a gridsym tiling with flip(t) ≅ t requires even n. Candidate mechanism: a
  score-multiset/parity obstruction under total tile reversal at odd n. Count at even n:
  1, 2 → 4 at n=8? (2^{n/2−2}?)
- **C2 (pure-blue anti-adjacency, n ≥ 4).** No line connects two pure-blue classes: every
  pure-blue tiling's flip-partner lies in a MIXED class (checked all n=4..7; n=3 is the
  exception — its unique blue line joins the two pure-blues). Suggests pure-blue classes
  are "extreme points" whose complements are deep in the bulk.
- **C3 (transitive is always pure-blue).** The transitive class has fiber 1 (H=1, |Aut|=1)
  and its unique tiling (all tiles forward) is gridsym — PROVED, giving pure-blue(n) ≥ 1
  and anchoring the principal line's H=1 end in the blue skeleton.
- **C4 (blue-cross concentration).** Blue cross lines concentrate on (mixed, mixed) as n
  grows (3/8 → 28/32 → 250/256): the pure-blue attachment count (5, 2, 6) stays O(1)-ish.
  Equivalently: almost every anti-reversible path-pair connects two mixed SC classes.
- **C5 (black bulk scaling).** Black-cross mass ratios (MX,MX):(MX,PB):(PB,PB) =
  858 : 5044 : 10112 at n=7 — compare fiber-mass products (7302², 2·7302·25460, 25460²)
  normalized: the line matching is close to fiber-proportional ACROSS types (a
  quasi-randomness statement for flip at the class level) but with measured positive
  assortativity — quantify at n=8.

## The unifying picture

The tiling space 2^m fibers over classes with odd fibers (T1); flip is a free involution
whose quotient (the lines) carries a 2-coloring by the ρ-fixed locus (blue = the
anti-symmetric stratum); SC-ness of a class ⟺ its fiber meets that stratum (T2), with the
parities (T3) and adjacency prohibitions (T4) forced purely by the two commuting
involutions (ρ, flip) and Rédei's oddness. What full structure still needs: the class-level
flip-matching multiplicities (the allocation tables) beyond parity — C5 says they are
near-fiber-proportional with measurable assortativity, so the completion is a
quasi-randomness theorem for the complement-tiling map plus the O(1) exceptional blue
skeleton (C1–C4). Then "which tilings map where" is determined: fibers by H/|Aut| (known),
the gridsym sublocus by T2/T3, and partners by the (near-random + skeleton) matching.

## Files
`04-computation/merged_metagraph_lines_n7_klein_S161.py` (+ .out; refinement
canonicalizer validated against brute force at n ≤ 5 via identical class counts 4/12/56 —
and NS-merged 184 at n=7 matches canon). Extends `merged_metagraph_line_accounting_klein
.py` (S75). Pointers: CLAUDE.md blue/black section (three claims now proved), LEM-003,
MISTAKE-033 (cross-class lines), mac-mini-S46 (THM-643 claim — proofs and data here feed
it), A000568.
