---
id: THM-584
title: The complement involution R is the ANTIPODAL MAP of the tournament arc-hypercube Q_{C(n,2)}; the iso-class metagraph adjacency is the S_n-quotient of Q_{C(n,2)}, so its R-even/R-odd (eps=+-1) eigenspace split equals the EVEN/ODD hypercube-level split, with eigenvalues =C(n,2) (mod 4) on the R-even block and =C(n,2)-2 (mod 4) on the R-odd block; the Perron mode is always R-even
status: PROVED (elementary representation theory of Q_d under coordinate-permutation S_n and the antipodal involution); VERIFIED exactly n=3..6 (block spectra, mod-4 residue law, commutation A o R = R o A all hold; A000568 class counts and SC counts 2,2,8,12 reproduced)
source: klein-2026-06-29-S1
depends_on:
  - HYP-3538   # the R +-1 eigenspace organizing principle (this gives its spectral skeleton on the metagraph)
related:
  - THM-583    # the witness-side half-system (eps=-1 stored in phi); this is the metagraph-side analogue
  - HYP-2685   # half-tiling address quotient: the discarded "side" bit IS the R-odd coordinate
  - HYP-2080   # self-complementary / Burnside counts (SC = R-even fixed points)
  - HYP-3539   # the n=6 face-compression overhead is NOT the R-odd coordinate (refuted) -- bulk packing
related_objects:
  - merged metagraph G_n/Z_2 = the R-even (+1) projection; V_merged = (A000568+SC)/2 = dim of R-even block
results:
  - 04-computation/r-eigenspace-merged-compression.py
  - 05-knowledge/results/r-eigenspace-merged-compression.out
  - 04-computation/r-block-spectra-antipodal.py
  - 05-knowledge/results/r-block-spectra-antipodal.out
---

# THM-584 — Complement is the antipodal map; the metagraph spectrum splits by hypercube level parity

## Setup

Let `d = C(n,2)`. Identify a labeled tournament on `n` vertices with a vertex of the hypercube
`Q_d = {0,1}^d`: one coordinate per unordered pair, the bit recording the arc orientation. Then:

- **One single-arc reversal = one hypercube edge.** The *arc-flip graph* on the `2^d` labeled
  tournaments is exactly `Q_d`, a `d`-regular graph.
- **`S_n` acts by permuting the `d` coordinates** (relabeling vertices permutes the pairs). This
  action **preserves Hamming weight** (it permutes coordinates), i.e. it preserves hypercube *levels*.
- **The complement / reversal involution `R` (`T -> T^op`) reverses every arc = flips every bit =
  the ANTIPODAL MAP `x -> x ⊕ 1` of `Q_d`.**

The iso-class metagraph adjacency `A` is the `S_n`-quotient of `Q_d`: in the basis of iso-class
indicator functions (= the `S_n`-invariant functions on `Q_d`),
`A[i][j] = #{ arcs e : flip_e(rep_i) ∈ class j }` (well-defined: `S_n` is transitive on each class
and commutes with arc-flips). Row sums are `d`. `A` is the standard tournament dominance-reversal
metagraph at the iso-class level.

## Statement

1. **Spectrum lives on hypercube levels.** Every eigenvalue of `A` has the form `d - 2k` for some
   level `0 <= k <= d` (the `S_n`-invariant subspace decomposes into level eigenspaces of `Q_d`).
2. **Complement = antipodal acts as `(-1)^k` on level `k`.** Hence `A` commutes with `R`
   (`A∘R = R∘A`, verified exactly n=3..6), and the iso-class space splits
   `V = V_+ ⊕ V_-` into the `R`-even (`eps=+1`) and `R`-odd (`eps=-1`) eigenspaces with
   - `V_+` = sum of **even-level** invariant eigenspaces, `dim V_+ = V_merged = (A000568(n)+SC(n))/2`;
   - `V_-` = sum of **odd-level** invariant eigenspaces, `dim V_- = (A000568(n)-SC(n))/2` = #NS complement-pairs.
3. **mod-4 residue law.** Because the level-`k` eigenvalue is `d - 2k`:
   - on `V_+` (k even): every eigenvalue `≡ d (mod 4)`;
   - on `V_-` (k odd): every eigenvalue `≡ d - 2 (mod 4)`.
4. **The Perron / bulk mode `d` is always `R`-even** (it is level `k=0`, the all-ones invariant).
   So the dominant "bulk" of the metagraph dynamics is the `eps=+1` part — exactly the
   SOS/Brouwer-provable side of the HYP-3538 organizing principle.

## Proof

`Q_d` has the orthonormal Fourier basis `{χ_S : S ⊆ [d]}`, `χ_S` an eigenvector of the adjacency at
eigenvalue `d - 2|S|`. `S_n` permutes coordinates, hence permutes the `χ_S` within each level
`|S| = k`; the level-`k` eigenspace is `S_n`-invariant. The metagraph `A` is the restriction of the
`Q_d` adjacency to the `S_n`-invariant subspace (class-indicator functions are exactly the
`S_n`-invariant functions), so its eigenvalues are a sub-multiset of `{d-2k}` with multiplicity =
`dim(invariants at level k)`. The antipodal map sends `χ_S -> (-1)^{|S|} χ_S` (since `χ_S(x⊕1) =
(-1)^{|S|} χ_S(x)`), so it acts as `(-1)^k` on level `k`; restricted to invariants this is exactly
the complement action `R`. Therefore `V_+` = even-level invariants, `V_-` = odd-level invariants,
and a level-`k` eigenvalue `d-2k` is `≡ d (mod 4)` iff `k` even. The block dimensions follow from
`R` being a fixed-point-free-on-NS involution on the `A000568(n)` classes (SC = fixed points), giving
`dim V_+ = SC + (A-SC)/2 = (A+SC)/2` and `dim V_- = (A-SC)/2`. ∎

## Verification (n=3..6, exact)

| n | d=C(n,2) | A000568 | SC | dim V_+ (=V_merged) | dim V_- | R-even residue d%4 | R-odd residue (d-2)%4 |
|---|----|----|----|----|----|----|----|
| 3 | 3 | 2 | 2 | 2 | 0 | 3 | — |
| 4 | 6 | 4 | 2 | 3 | 1 | 2 | 0 |
| 5 | 10 | 12 | 8 | 10 | 2 | 2 | 0 |
| 6 | 15 | 56 | 12 | 34 | 22 | 3 | 1 |

Observed block spectra (klein-S1):
- n=4: `V_+ = {-2, 2, 6}` (all `≡2 mod4`), `V_- = {0}` (`≡0 mod4`).
- n=5: `V_+ = {-6, -2×3, 2×4, 6, 10}` (all `≡2`), `V_- = {0, 4}` (all `≡0`).
- n=6: `V_+` eigenvalues all `≡3 (mod4)` (`{-9,-5,-1,3,7,11,15}`), `V_- ` all `≡1 (mod4)`
  (`{-7,-3,1,5,9}`). Perron `15 ∈ V_+`.

`commute A∘R == R∘A` returns True at every n=3..6.

## The Klein-`V_4` base case (n=4)

The four classes `T=(0,1,2,3)`, `+=(0,2,2,2)`, `−=(1,1,1,3)`, `S=(1,1,2,2)` are the Klein four-group
`V_4=(Z_2)^2` under XOR of two bits: `x` = "source destroyed", `y` = "sink destroyed". `R` is the
coordinate swap `x<->y` (reversing arcs swaps sources and sinks). The split of THM-584 here is concrete:

- **`R`-even coordinate** `u = #boundary defects ∈ {0,1,2}` (T:0, ±:1, S:2) — `R`-invariant; this is
  the merged-class label `G_4/Z_2 = {T, [±], S}` and the level-parity grading.
- **`R`-odd coordinate** `w = (source-killed) − (sink-killed) ∈ {−1,0,+1}` — zero on SC `{T,S}`,
  `±1` on the NS pair `{+,−}`, and `w -> −w` under complement. This is the signed `eps=-1` content.

So `V_4` is the smallest nontrivial instance of THM-584: `Q_6 / S_4`, antipodal = swap, `u`/`w` =
even/odd level coordinates.

## Significance

THM-584 gives the **metagraph-side realization of the HYP-3538 organizing principle**, complementary
to THM-583's witness-side (`f` on the half-system). The same `eps=±1` of the same `R`:
- the **merged metagraph `G_n/Z_2` is the `R`-even projection** (its dimension `V_merged` is the
  `eps=+1` block dimension), so quotienting by complement = discarding the odd-level coordinate;
- **complement = the antipodal map**, which is precisely the map Borsuk–Ulam is about — placing the
  project's topological certificate on the same footing as the hypercube level grading;
- the **mod-4 residue law** is a hard fingerprint: any claimed metagraph eigenvalue must satisfy it,
  a free correctness check on spectral computations.

Caveat preserved (MISTAKE-033): `R` (reverse ALL arcs = antipodal of `Q_d`) is NOT the *tiling*
complement (flip all *tiles*, fixing the base path — the d=m blue/black line). THM-584 is about `R`
on the full `Q_d` arc-cube, not the fixed-base-path tiling cube `Q_{C(n-1,2)}`.
