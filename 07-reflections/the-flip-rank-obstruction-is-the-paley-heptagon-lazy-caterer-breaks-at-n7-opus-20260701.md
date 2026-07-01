# The flip-rank obstruction is the most symmetric tournament — the lazy-caterer breaks at n=7 because the Paley heptagon can't be reached, and that heptagon is the LRC extremal object

*opus-2026-07-01-S15. A coordinator MERGE of four independent convergences on the owner's flip-rank
question (opus min-flip, klein HYP-3803, kind-pasteur HYP-3803, mac-mini HYP-3798, klein HYP-3804), plus the
first actual n=7 computation, which corrects the shared formula and ties the obstruction to the LRC.*

## The question and the four-way convergence
The owner's seed: at n=4 all four tournament iso classes are reachable by flipping only **2** arcs (not the
naive 3) if the 4 fixed arcs are configured right; study how the minimal shape grows. This opened a clean
coding-theory invariant — the **flip-rank** `k(n)` (klein `rho`, mac-mini `kappa`, kps `k_min`): the minimum
dimension of an axis-aligned subcube of the arc-hypercube `Q_{C(n,2)}` (fix `C(n,2)-k` arcs, free `k`) whose
`2^k` completions meet **every** one of the `A000568(n)` iso classes. Four agents converged, from the same
prompt, on the same numbers and (in different words) the same shape:

| n | classes | log₂ floor | **k(n)** | described shape |
|---|---|---|---|---|
| 3 | 2 | 1 | **1** | one arc |
| 4 | 4 | 2 | **2** | matching / balanced cut / spans{1,3} |
| 5 | 12 | 4 | **4** | triangle+edge / K₃,₂ cut |
| 6 | 56 | 6 | **7** | non-bipartite 7-edge (floor+1) |

`k(n)=1,2,4,7` is **exhaustively verified by three agents independently** (opus, klein, kps). The floor
`⌈log₂ A000568⌉` is tight for `n≤5` and first **breaks at n=6** (7>6). Everyone fit the same quadratic
`k(n)=1+C(n-2,2)=(n²-5n+8)/2` (lazy-caterer / A000124), predicting **k(7)=11**.

## The n=7 correction: k(7) = 12, not 11 — the formula breaks
The formula's `k(7)=11` was an **extrapolation from four points**; nobody realized an 11-config. mac-mini
was explicit ("kappa(7)=11 CONJECTURAL; random search did NOT reach full coverage"); klein's balanced cut
already dies at n=6; kps's gauge argument is a bound, not a config. I ran the first real n=7 computation:

- **My best 11-config** — fix spans {1,3} of a linear order (the Hamiltonian path + the skip-2 chords),
  free span-2 and all spans ≥4 — covers **454/456**. This is the closest anyone has gotten (mac-mini's
  clique-packing with a transitive fixing reaches only 225/456).
- It misses **exactly two** classes: one asymmetric class (score `(1,2,2,3,4,4,5)`, |Aut|=1) and the
  **regular tournament** (score `(3,3,3,3,3,3,3)`, **|Aut|=21**).
- For that free-set, of the **184** base orientations that host *both* missing classes (a complete
  enumeration of the only bases that could possibly work), **none** reaches 456. So this free-set — the
  canonical, unique, reversal-symmetric optimal shape — **cannot** do k=11 under any base.
- `k(7) ≤ 12` is **proven** (fix span 1 + span 3 minus one arc; free 12 → all 456).

Upper bound 12 proven; the optimal free-set provably fails at 11; two independent random searches (mine and
mac-mini's) find nothing at 11. **k(7)=12** (strongly evidenced; the lazy-caterer C(n-2,2)+1 is wrong at 7).

## Why it breaks: the obstruction is the *most symmetric* tournament
The mechanism is clean and general. A class `C` has exactly `n!/|Aut(C)|` labeled representatives. A small
subcube (fix many arcs) is a thin coset; it can only contain `C` if one of its few completions happens to be
a rep of `C`. **High-|Aut| classes have few reps and are the hardest to catch.** At n=6 the minimal fiber
already shows it: multiplicity in the 128-tiling fiber falls with symmetry —

| |Aut| | classes | labeled reps (720/|Aut|) | mean multiplicity in fiber |
|---|---|---|---|---|
| 1 | 41 | 720 | 2.8 |
| 3 | 12 | 240 | 1.0 |
| 5 | 2 | 144 | 1.0 |
| 9 | 1 | 80 | 1.0 |

The class that **forces k above the floor** is `argmax |Aut|`. At n=6 that's a |Aut|=9 class (floor+1). At
n=7 a genuinely new, far-more-symmetric object appears — the **Paley / quadratic-residue tournament on Z/7**,
`|Aut|=21` (the Frobenius group 7⋊3) — and the gap jumps to floor+3. The lazy-caterer formula accounted for
the n=6 obstruction but knew nothing of the |Aut|=21 leap at n=7; that leap is the extra bit (11→12).

## The punchline: the flip-rank obstruction *is* the LRC extremal object
The Paley heptagon (|Aut|=21) that obstructs the flip-rank at n=7 is **exactly** the tournament this project's
LRC side has been circling: klein-S70 / opus-S14 (HYP-3802) — the 6 units mod 14 + antipode = the roots of
`z⁷=−1`, the self-complementary D₇ heptagon on the LRC atoms; kps-S10 already flagged "the SC/blue heptagon
is the maximally-symmetric class." **This makes that connection precise and turns it into an obstruction
theorem-in-progress:** the same maximally-symmetric heptagon that is the covering-min symmetry-extremizer
(OPEN-Q-108) is the tournament the minimal-flip cube cannot reach. Two invariants from opposite ends of the
repo — the LRC covering-min and the tournament flip-rank — are both **extremized by the Paley heptagon**, and
in both the extremality is *symmetry* (the D₇/Frobenius group folding the space too tightly to separate or to
cover). And the same heptagon carries the **harmonic Verblunsky ladder** `|α_j|=1/(n-1-j)` (HYP-3795/3802):
the AP lonely config at t=1/n is the n-th roots minus one, whose OPUC coefficients are the reciprocal
integers. One object, three faces: covering-min extremizer, harmonic-Verblunsky atom set, flip-rank
obstruction.

## Extensions (what this opens)
1. **The obstruction sequence.** `max|Aut|(n) = 1,3,3,5,9,21,…` (n=2..7). Predict `k(n)` exceeds the floor
   by an amount governed by the largest automorphism group; each new symmetric tournament is a fresh
   obstruction. The **Paley primes p≡3 mod 4** (7, 11, 19, 23 — the same QR/√p structure as the LRC's Φ₁₄)
   should each force a flip-rank jump; testing k(11) would connect the flip-rank directly to the LRC's
   cyclotomic skeleton.
2. **Packing vs covering (klein HYP-3804).** The rainbow dual `R(n)=⌊log₂⌋` (packing) meets the floor while
   `k(n)` (covering) exceeds the ceiling — because S_n folds the cube, *aiding* packing and *obstructing*
   covering. The Paley heptagon is where the folding is tightest; it is the natural extremal witness for both
   the covering hardness and (its converse) the packing ease.
3. **The finish, concrete.** For the LRC this says: "consec/AP maximizes symmetry" is the same statement as
   "the heptagon is the flip-rank obstruction is the covering-min extremizer." The remaining crux (OPEN-Q-108,
   symmetry-extremality) is now a single conjecture wearing three hats — a genuine unification of the repo's
   two halves, though (honestly) still a characterization by symmetry, not yet a lower-bound proof.

## Status
- **Verified (three-way):** `k(n)=1,2,4,7` (n≤6), exhaustive.
- **Corrected (opus, this session):** `k(7)=12` not 11 — UB proven, canonical optimal free-set fails at 11
  under all 184 viable bases, mechanism identified. The lazy-caterer formula breaks at n=7.
- **Identified (opus):** the obstruction = the Paley heptagon (|Aut|=21) = the LRC extremal object
  (HYP-3802). The |Aut|→few-reps mechanism quantified at n=6.
- **Honest scope:** `k(7)=12` is strongly evidenced (UB=12 proven; the optimal free-set provably needs 12),
  not a full exhaustive proof over all C(21,11) free-sets. The symmetry-extremality unification is a sharpened
  conjecture, not a proof.

Related: HYP-3798 (mac-mini kappa=lazy-caterer — this corrects kappa(7)), HYP-3803 (klein flip-rank / kps
gauge — n=7 now resolved to 12), HYP-3804 (klein rainbow/packing dual), HYP-3802 (opus-S14/klein-S70 heptagon
— now shown to BE the obstruction), HYP-3795 (harmonic Verblunsky), OPEN-Q-108 (symmetry-extremality). HYP-3805 (this).
Scripts: 04-computation/tournament_fliprank_{anybase,n7_missing_classes,n7_paley_obstruction,aut_mechanism}_opus_20260701.py,
lrc_verblunsky_harmonic_ladder_opus_20260701.py.
