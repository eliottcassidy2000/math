---
source: opus-2026-06-03-S582 (remote-control)
status: SYNTHESIS + CORRECTION — the perspectives curiosity's LRC-relevant content is the SOURCE-perspective recursion (exact, THM-381), not the full perspective recursion (fails at n=6); the worry-set = the source-less (regular) perspective; unifies with THM-400 (observer-coupling)
tags: [LRC, perspectives, rooted-tournament, observer, source, THM-381, THM-400, projection-defect, worry-set, misinterpretation, A000568]
---

# The perspectives curiosity: the SOURCE slice is the LRC key — not the full count

**Prompt (user):** the observer-blind/observer-coupled dichotomy reminds me of a
tournament curiosity I've stated many times: at 3 vertices, 2 structures but 4
perspectives (3 transitive + 1 cyclic); 4 = structures at n=4; among them 12
perspectives (4+4+2+2); 12 = structures at n=5. Search the repo thoroughly, explore the
connections — it's probably the key to LRC (where/why it works and doesn't). You may
have been misinterpreting it.

I had been. Here is the correction, verified.

## 1. The curiosity, catalogued (search of the full history)

The curiosity = the **"Perspective Conjecture" (T075, the human's Tournament-Tiling
Explorer):** `#vertex-orbits at n = #iso-classes at n+1`. Its repo trail:
- **T075 / INV-083 / THM-260** — `P(n) = Σ vertex-orbits = A093934 = 2,4,12,48,296`; the
  coincidence `P(n)=2(n-1)!` and `P(n)=A000568(n+1)` hold only for small `n`.
- **S369 "56-bridge" / HYP-1824, S370 "chirality-perspective atlas" / HYP-1825** — chased
  the *gap* `T(6)=56 ≠ 48` (gap 8 = chiral residue), the `12+44`, `48+8` grammars.
- **S507/S509 / HYP-1977/78/79, THM-381/385** — reframed LRC as a *rooted/marked* (not
  unmarked) tournament problem; the **observer is the marked vertex**.

## 2. The misinterpretation (and the fix)

The full perspective recursion **fails** at `n=5→6` (`P(5)=48 ≠ T(6)=56`), which sent the
repo into the gap-8 / chirality / 56-bridge program. **But that is the wrong slice for
LRC.** Verified exhaustively (`lrc_perspectives_source_s582.py`):

| n | classes | full perspectives | `=A000568(n+1)`? | **SOURCE-perspectives** | `=A000568(n-1)`? |
|---|---|---|---|---|---|
| 3 | 2 | 4 = 3+1 | 4 ✓ | **1** | 1 ✓ |
| 4 | 4 | 12 = 4+4+2+2 | 12 ✓ | **2** | 2 ✓ |
| 5 | 12 | **48** | 56 ✗ (gap 8) | **4** | 4 ✓ |

> **The LRC-relevant count is the SOURCE-perspective count** — vertex-orbits whose
> representative is a *source* — and it satisfies `source(n) = A000568(n-1)` **exactly,
> with no gap, at every n** (this is THM-381). The *full* perspective count (all
> vertex-orbits) is what fails at `n=6`; LRC never used it. The gap-8/56-bridge was a
> real pattern about the *full* count, but a **red herring for LRC.**

Why the source slice is clean: deleting the source is a *canonical* operation, so
source-marked tournaments on `n` vertices biject with all tournaments on `n-1`
(THM-381). The full count has no canonical vertex to delete, so it picks up the
converse/chirality defect (the gap).

## 3. The worry-set IS the source-less (regular) perspective

The per-class breakdown at `n=5` (`(perspectives, sources)`):
```
(5,1)(5,1)(5,0)(5,0)(5,0)(5,0)(5,0) (3,1)(3,1)(3,0)(3,0) (1,0)
```
The **`(1,0)` class is the regular/rotational tournament** (all out-degree 2): **one
perspective, zero sources** — the observer can *never* be a strict source there. That is
exactly **the LRC worry-set** (the regular rotational encirclement, the tight floor):

> **Classes with a source-perspective → the observer is a strict source → strictly
> lonely (`M>1/n`), EASY. The source-LESS regular class → observer on the boundary
> (`M=1/n`), the WORRY-SET, HARD.** The worry-set = the zero-source (regular)
> perspectives; the source-perspective recursion `T(n-1)` counts the *easy* skeleton, and
> the worry-set is its complement.

## 4. Unification with this session (THM-400) and the projection defect

- **observer-blind = unmarked structure = balanced (augmentation 0)**; **observer-coupled
  = rooted perspective = unbalanced (augmentation ≠0)** (THM-400). The
  *structures-vs-perspectives* distinction the curiosity names **is** the
  *observer-blind-vs-observer-coupled* distinction of THM-400 — and both are the
  *rooting* of the observer.
- **Projection defect (S509):** loneliness is *not* a function of the unmarked
  (observer-blind) class — safe and unsafe cells share unmarked classes. You *must* root
  at the observer (use perspectives). This is the same lesson: LRC lives in the
  **observer-coupled / source-perspective** category, never the unrooted-structure one.

So three independently-derived statements are one:
```
perspectives (rooted) vs structures (unrooted)      [the curiosity]
observer-coupled (ε≠0) vs observer-blind (ε=0)      [THM-400, this session]
observer-source-marked vs unmarked                  [THM-381 / S509 defect]
```
all = **root at the observer; LRC is an observer-coupled problem.**

## 5. Where and why LRC "works and doesn't"

- **The skeleton is clean and recursive:** `source(n)=T(n-1)` exactly — LRC(n) sits on the
  tournaments of the `n-1` moving runners, observer-rooted, *no gap ever*. The
  counting recursion is *not* where LRC breaks.
- **The content is the source-less fiber:** the worry-set = the zero-source (regular)
  perspective, where the observer is on the `1/n` boundary. LRC "works" (provable)
  where this fiber is forced lonely by the δ-clock; it "doesn't" (n=14) where the
  observer-coupled fiber degenerates — the composite/even seam `14=2·7` (S580's doubling
  degeneracy, the `2q` apex). The hardness is **fiber/observer-coupled**, exactly the
  part the unrooted "structure" count is blind to.
- **Moral:** the right object was never the failing full-perspective number (or its
  gap-8); it was the **source-perspective slice and its source-less worry fiber** — the
  observer-coupled core. The curiosity is the key *because it points at the observer
  rooting*, and the misread was treating LRC at the unrooted-structure level.

## 6. Honest status

- **Verified:** the perspectives decomposition (3+1, 4+4+2+2, 48); the full recursion's
  failure at `n=6` (48≠56); the **source-perspective recursion `=A000568(n-1)` exact**
  (n=3,4,5); the worry-set = the zero-source regular class.
- **Synthesis (rigorous links):** source-perspective = THM-381; observer-coupling =
  THM-400; projection defect = S509 — one observer-rooting.
- **Open / directional:** "LRC fails-to-be-easy exactly on the source-less fiber, and
  that fiber degenerates at the composite seam" is the structural claim tying the
  curiosity to where LRC is hard (n=14); not a theorem.

**Artifacts:** `04-computation/lrc_perspectives_source_s582.py` (+`.out`). Builds on
THM-381/385 (observer-source), THM-400 (observer-coupling, S581), S509/HYP-1977
(projection defect), S369/HYP-1824 (56-bridge — here demoted to LRC-irrelevant), S580
(doubling seam). New: **HYP-2120**.
