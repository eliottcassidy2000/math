---
source: opus-2026-06-02-S558 (remote-control)
status: SURVEY/SYNTHESIS — literature-grounded map of LRC proof methods by n, why each works at its range, and how the repo's tools map on
tags: [LRC, survey, methodology, proof-history, finite-checking, sieve, chromatic-number, n14, corrections]
---

# LRC proof methodologies by n: the complete map, and why each lives where it does

**Prompt (user):** go through the proof methodologies available at each size up
to 14; understand exactly which are available at what ranges of n and *why*.

This is a literature survey (verified June 2026) cross-mapped to the repo's own
tools. **Two corrections to our prior framing fall out and matter a lot.**

## 0. Convention (now pinned)

Standard modern statement uses **`k` = number of non-zero integer speeds, gap
`1/(k+1)`, total runners `k+1`**. The repo's **`n` = total runners = `k+1`**, gap
`1/n` — they agree. So repo `n=14` ⟺ `k=13` non-zero speeds. WLOG integer,
gcd-1 speeds (scaling invariance).

## 1. The complete map (n = total runners)

| n (=k+1) | year | authors | method | regime |
|---|---|---|---|---|
| 2,3 | — | folklore | elementary (three-gap / continued fractions) | structural |
| 4 | 1972 | Betke–Wills; Cusick | geometry of numbers / view-obstruction | structural |
| 5 | 1984 | Cusick–Pomerance | computer-assisted Diophantine (later elementary; flows: Bienia et al.) | structural |
| 6 | 2001 | Bohman–Holzman–Kleitman (Renault 2004 simpler) | combinatorial averaging over chosen times | structural |
| 7 | 2008 | Barajas–Serra | **regular/circular chromatic number of distance (Cayley) graphs** | structural |
| 8 | 2025 | Rosenfeld | **finite-checking + divisibility sieve** | computational |
| 9 | 2025 | Rosenfeld; Trakulthongchai | refined divisibility sieve | computational |
| 10 | 2025 | Trakulthongchai | refined sieve (small-multiplier lifts) | computational |
| 11,12,13 | Apr 2026 | Sungkawichai–Trakulthongchai (arXiv 2604.23906) | intermediate sieves + **polynomial method** | computational |
| **14** | — | **OPEN (immediate frontier)** | — | — |

**Correction #1 to our prior sessions:** LRC is proved through **n = 13**, not
n = 7. Our "even-fold uses the proven LRC(7)" (S554) was leaving enormous power on
the table — see §4.

## 2. The structural era (n ≤ 7): bespoke, each a new idea

Each case is settled by exploiting *low-dimensional* special structure, and each
needs a genuinely different technique:

- **n=3:** one effective rotation; the three-gap theorem / a single
  equally-spaced grid meeting the far-band. (The repo proves this two ways:
  Fourier/Legendre-mod-3, S526; geometric center-grid, S522o.)
- **n=4:** view-obstruction — LRC ⟺ a line on the torus misses a box; in low
  dimension the lattice/geometry-of-numbers argument closes (with a finite check).
- **n=5:** Diophantine case analysis, originally computer-assisted; the residual
  resonant families are finite and checkable.
- **n=6:** Bohman–Holzman–Kleitman *averaging*: average a "far" indicator over a
  cleverly chosen finite set of times and bound the overlaps so the runners can't
  all be near at once. Combinatorial, not arithmetic.
- **n=7:** Barajas–Serra: recast as the **circular chromatic number of the
  distance graph `G(ℤ,D)`** (a Cayley graph); LRC ⟺ a multiplier circular
  `(k+1)`-coloring exists. A deep structural classification of these colorings.

**Why the structural era caps at 7.** These are not one method but a ladder of
escalating bespoke arguments; the case analysis / coloring classification grows
combinatorially with the number of incommensurate rotations, and *the geometric
codimension of the obstruction rises by one per runner*. The repo's own geometric
attempt makes the cap mechanism explicit (S522o): the "center trick" reduces `n`
runners to "`n-2` others far on a 1-D grid," a **codimension-(n-3) hit** — automatic
at n=3, a 1-D-grid-vs-2-D-box gap at n=4, and only worse after. No single witness
suffices; you must iterate, and the iteration is what explodes. By n=8 no one has
found a closed structural argument — the field switches tools.

## 3. The finite-checking era (n ≥ 8): one uniform machine

A single pipeline replaces bespoke cleverness with a (huge) finite computation:

1. **Reduce to bounded integer speeds.** Tao (2018): suffices to check integer
   speeds of size `n^{O(n²)}`. Malikiosis–Santos–Schymura (2025): sharpened — a
   counterexample has `∏ uᵢ < B_k`, `B_k = (\binom{k+1}{2}^{k-1}/k)^k` (a
   *linearly-exponential* bound), via the **LR-zonotope** convex-geometric
   reformulation (Henze–Malikiosis).
2. **Divisibility sieve (Rosenfeld 2025).** For a prime `p`: if *every* residue
   tuple mod `pℓ` that avoids multiples of `p` admits a witness `t ∈ (1/ℓp)ℤ`,
   then **any counterexample must have all speeds divisible by `p`**.
3. **Prime-product contradiction.** Run the sieve over a set `P` of primes for
   which it succeeds ⟹ a counterexample has `∏uᵢ ≥ ∏_{p∈P} p`. If
   `∏_{p∈P} p ≥ B_k`, this contradicts step 1 ⟹ **no counterexample**.
4. **Accelerators (Trakulthongchai; Sungkawichai–Trakulthongchai 2026):**
   intermediate sieves (lift by small `c=2,3`, project back to `ℤ_p`), and a
   **polynomial method** that proves the single hardest tuple — the tight AP
   `(1,…,k)` — is "proper" *analytically, with no computation*, **whenever `k+1`
   is an odd prime** (Combinatorial-Nullstellensatz-style: two indicator
   polynomials over `ℤ_{k+1}` forced equal).

**Why this era starts at 8 and not earlier:** below 8 the structural proofs
already exist and are cleaner; the finite bound `B_k`, while "small" in theory, is
astronomically large, so the sieve is only worth it once no closed argument is in
reach.

## 4. Why **n = 14** is exactly the wall — and it is the 2·7 reason

From the Apr-2026 paper's own account, the bottleneck at `k=13` (n=14) is the
`c=(k+1)` lift, and **the reason is that `k+1 = 14 = 2·7` is composite**:

- The **polynomial-method shortcut requires `k+1` to be an odd prime.** For
  `k=12`, `k+1 = 13` is prime → the tight tuple `(1,…,12)` is handled *for free*,
  no computation. For `k=10`, `k+1 = 11` prime → same.
- For `k=13`, `k+1 = 14` is **not** an odd prime → the shortcut fails → one must
  run the full, exponentially expensive `c=14` computational lift. (k=12 already
  cost ~40 days on 10 cores; k=13 needs far more.)

So the literature's obstruction at n=14 is **precisely the compositeness of
`14 = 2·7`** — the same structural fact the repo's threads keep circling:
"the 7 impossibility" (S554/S555), the even-fold `n=14 → n=7` (S554), the mod-7
CRT singleton (oracle-S552o), the doubled-prime `n=2p` analyses. We were not
chasing a metaphor: **`2·7` is the actual reason the modern machine stalls at 14.**

## 5. The repo's tools, mapped onto the literature

| repo tool | literature equivalent |
|---|---|
| covering reformulation (S525) | the standard gap/forbidden-interval cover |
| Fourier character sum, proves n=3 (S526) | Fourier/measure lower-bound method |
| geometric center-grid, proves n=3 (S522o) | view-obstruction / grid witness |
| **sieve THM-369** (`q∤v ⇒ ‖v/q‖≥1/q`, witness `a/q`) | **= Rosenfeld's divisibility sieve** — the modern engine, independently rediscovered |
| "counterexample needs a multiple of each `q≤n`" (S551) | the divisibility necessary condition feeding the prime-product contradiction |
| sieve incompleteness at bounded modulus (S551, HYP-2052) | why the sieve must be paired with the MSS *finite bound* (a counterexample's speeds are bounded) — the two halves of the real proof |
| even-fold `M(S) ≤ M(fold)` (S554, HYP-2056) | a genuine reduction; **should use proven LRC(13)**, not LRC(7) |
| pinch lemma, radius `r/s` (S557, HYP-2059) | extremizer/critical-time structure (tight tuple analysis) |
| spectral gap `2/(2n-1)` (oracle-S552) | margin of the tight extremiser |

**Correction #2 (a real strengthening).** The even-fold sends the even speeds
`{2u}` to `{u}` (fold) with `M(S) ≤ M(fold)`. A primitive 13-set has at most 12
even speeds, so `|fold| ≤ 12` and **`M(fold) ≥ 1/13` by the now-proven LRC(13)** —
for *every* config, not just `e≤6`. So the even half of LRC@14 is *entirely*
protected by an established theorem; the whole residual is the odd-runner coupling
(the antipodal split, S554), now sitting on top of a `1/13`-safe even base. Our
`e ≤ 6` restriction was an artifact of using LRC(7).

## 6. What this means for attacking n=14

The honest landscape:
- **The published route to 14 is computational** (extend the `c=14` sieve lift),
  blocked by `14=2·7` composite. A repo contribution there = a better algorithm
  for the `c=14` improper-tuple computation, or an *algebraic* substitute for the
  polynomial method when `k+1` is twice an odd prime (`k+1=2q`). This is the
  single most leveraged target, and it is exactly our `2·7` wheelhouse.
- **A structural (non-computational) proof at 14** would be a major result; our
  pinch/even-fold/spectral-gap tools are aimed here but face the codimension wall
  (§2) that has held since n=8.
- Either way, the precise object to crack is the **tight tuple `(1,…,13)` at
  `k+1=14`**, where the prime-based polynomial method fails — i.e. find the
  `2q`-analogue of "`(1,…,k)` is proper."

**Sources:** Wikipedia "Lonely runner conjecture"; Perarnau–Serra, *The Lonely
Runner Conjecture turns 60* (arXiv:2409.20160); Tao (2018); Malikiosis–Santos–
Schymura, *Linearly-exponential checking…* (arXiv:2411.06903); Rosenfeld, *…eight
runners* (arXiv:2509.14111); *Nine and ten lonely runners* (arXiv:2511.22427);
*…nine runners* (arXiv:2512.01912); Sungkawichai–Trakulthongchai, *Eleven, twelve,
and thirteen lonely runners* (arXiv:2604.23906); Quanta, *New Strides…* (2026-03).
Repo: S505 (route landscape), S430 (lens atlas), S359 (distance-graph coloring),
S522o/S526 (small-n proofs), S554 (even-fold), S557 (pinch), oracle-S552/S552o.
