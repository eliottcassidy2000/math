---
source: opus-2026-06-07-S703 (user seed: derived-series monodromy of roots↔coefficients; "no quintic
  = two 3-subsets sharing one element")
status: SYNTHESIS — the user's solvability/derived-series picture, verified (THM-436), pushed onto the
  repo's famous problems. (1) The monodromy of the roots↔coeffs cover (S699p) is graded by the DERIVED
  SERIES of S_n; threshold n=5 = A_5 perfect = two cyclic triangles sharing ONE vertex (3+3−1=5) =
  the round tournament C_5 = the LRC n=5 cyclotomic witness. (2) Abel–Ruffini (derived series never
  reaches 1) MIRRORS THM-406 (no finite Bonferroni; all-orders cancellation) — both "no finite-depth
  tower." (3) CONJECTURAL FRAMEWORK (HYP-2303): the LRC witness hierarchy (clock⊃shell⊃pair-sum,
  THM-420/430) is a RADICAL TOWER; worry-set (cyclotomic/abelian) = solvable bottom = TIGHT; residual
  R(n) = perfect/unsolvable core. (4) THE INVERSION: in Galois solvable=easy; in LRC solvable=
  cyclotomic=TIGHT=hard. (5) THE ICOSAHEDRAL BRIDGE (classical + S699h): A_5 = icosahedral rotation
  group; Klein solves the quintic via the icosahedron; the repo's A_5 unit-distance Cayley graph
  (spherical HN, S699h) sits literally on the quintic's group. (6) HN "5" resonance: honest/speculative.
tags: [galois, solvability, derived-series, abel-ruffini, quintic, A5, icosahedron, klein, monodromy,
  round-tournament, C5, cyclic-triangle, lonely-runner, worry-set, witness-tower, vitali-wall,
  hadwiger-nelson, unit-distance, tournament-analysis, FTA, inversion]
---

# The solvability tower: Galois ↔ LRC ↔ the icosahedron

**Prompt (user):** the monodromy of the roots↔coefficients cover, graded by nested commutators; degree
3 kills double commutators, degree 4 triple, but at degree 5 "two subsets of 3 solutions each sharing
one element" allow arbitrarily deep nesting — no quintic formula. "Take this further to work on famous
problems."

The user's picture is exactly the **derived series of `S_n`** acting on the root-fiber of the FTA
cover (my S699p monodromy = Galois). It is correct and now verified (THM-436): largest `k` with
`S_n^{(k)}≠1` is `n−2` for `n≤4`, **∞** for `n≥5`. This reflection takes it onto the repo's problems.

## 1. The threshold IS a tournament: C₅

The "two 3-subsets sharing one element" is the overlap `3+3−1=5`. A 3-subset with a cyclic order is a
**cyclic triangle** (3-cycle) — the atom of both `A_n` (3-cycles generate it) and tournament theory
(the rock-paper-scissors triangle). Two cyclic triangles sharing exactly one vertex first fit on **5
vertices**, and the smallest tournament realising it is the **round tournament `C₅`** (verified: 5 of
its 10 triples are cyclic). `A₅` is *perfect* (`[A₅,A₅]=A₅`) precisely because `[(123),(345)]` is a
3-cycle — the overlapping pair regenerates all of `A₅`. So:

> **The quintic's unsolvability is the appearance of two overlapping cyclic triangles, first at `C₅`
> — and `C₅` is the LRC `n=5` cyclotomic worry-set witness (THM-403).** The repo's smallest
> "interesting" tournament and the smallest unsolvable Galois group are the same object on 5 points.

This is the Tournament-Analysis reading (directive): the OCF/Rédei world is built from odd cycles, and
the 3-cycle is the smallest; the derived series of the vertex-symmetry group becomes perfect exactly
when two 3-cycles can interlock, at `n=5`.

## 2. Abel–Ruffini mirrors the Vitali wall (THM-406)

The deepest structural echo. Abel–Ruffini says the derived series of `S_n` (`n≥5`) **never reaches
`1`** — no finite tower of radicals (each radical = one abelian/commutator step) expresses the roots.
THM-406 says the LRC covering-depth inclusion–exclusion **cancels to all orders** — no finite
Bonferroni certificate (each Bonferroni order = one inclusion–exclusion step) certifies loneliness.

> **Both are "a finite-depth tower fails because of a perfect (depth-∞) subobject."** `A_n` is the
> perfect group; the all-orders overlap is the perfect measure-cancellation (the Vitali wall). The
> "no quintic formula" and the "no finite LRC certificate" are the same shape of obstruction. This
> sharpens the repo's Vitali-wall theme: it is the *solvability wall* in disguise.

## 3. The conjectural payoff: the witness hierarchy as a radical tower (HYP-2303)

The LRC witness hierarchy I/we built — **clock `t=1/k` ⊃ shell `t=m/(2n−1)` ⊃ pair-sum `t=m/(a+b)`**
(THM-420/430) — looks like a **radical tower**: each level adjoins a coarser symmetry (cyclotomic
roots of unity at the clock; the antipodal `σ`/quadratic at the shell, S702; pair-sums above). The
proposal:

> **A config is certifiable by the witness tower iff its local witness-monodromy is *solvable*; the
> worry-set (cyclotomic, abelian monodromy, THM-403/HYP-2282) is the *solvable bottom*; the residual
> `R(n)` (S700) is the *perfect / unsolvable core* where no finite tower certifies.** The depth of the
> tower needed = the derived length of the local monodromy.

This is conjectural (HYP-2303) — the verified content is the group theory and the C₅ realisation, not
an LRC equivalence. But it predicts where to look: `R(14)`'s difficulty should be a *non-solvable*
local monodromy, tied to the `3³` shell tower (THM-428) — the prime-power depth = the commutator depth.

## 4. The inversion (solvable = tight = hard)

A genuine subtlety worth stating. In Galois theory **solvable = good** (a radical formula exists). In
LRC the **solvable/cyclotomic worry-set is the *bad* case**: abelian monodromy = roots of unity =
rigidity = `M=1/n` *exactly* (the tight obstruction). The "unsolvable"/generic configs are **loose**
(lots of room; loneliness is generic in measure). So:

> **Rigidity (= solvability = cyclotomic) is what pins LRC to the floor.** The hard LRC cases are the
> *algebraically simplest* (most symmetric) ones — the opposite of Galois, where the simplest
> (solvable) are the *solvable* ones. The LRC worry-set is "too symmetric to be loose."

The same inversion appears in HN: the cyclotomic/Eisenstein lattice is the χ=3 *floor* (most
symmetric, easiest to colour), and escaping it (Heegner/non-cyclotomic, S699m) is what *raises* χ —
symmetry = low complexity there too, but in the colouring it lowers the invariant. The two problems
read the symmetry↔complexity axis with opposite signs.

## 5. The icosahedral bridge — quintic ↔ A₅ ↔ unit distance (classical + S699h)

The classical fact (Klein, *Lectures on the Icosahedron*, 1884): the rotation group of the
**icosahedron is `A₅`** (order 60), and the general quintic is solved by reducing to the *icosahedral
equation* (a non-radical "hypergeometric/theta" solution). So `A₅` — the first unsolvable group — is
*the symmetry of a Platonic solid in `ℝ³`*.

The repo already has this group in a metric problem: **S699h's unit-distance Cayley graph on `A₅` is a
spherical Hadwiger–Nelson problem**, and `A₅ ≅` icosahedral. So:

> **The quintic's group `A₅` is literally the symmetry group of an `ℝ³` unit-distance / HN object
> (the icosahedron / the `A₅` Cayley graph, S699h).** The "no quintic formula" obstruction and the
> spherical-HN chromatic problem live on the *same* `A₅`. Klein's icosahedral solution suggests the
> probe: is the `A₅` spherical-HN colouring controlled by the *icosahedral equation* / the `A₅`
> invariant theory (degrees `2,6,10,15`)? — a concrete, classical handle the repo has not used.

## 6. The Hadwiger–Nelson "5" (honest/speculative)

Tempting: `χ(ℝ²) ≥ 5` (de Grey) and the quintic threshold both single out **5**. Is the first
*combinatorially forced* chromatic jump the first *unsolvable* degree? I have **no proof or mechanism**
— de Grey's 5 is a finite-gadget count, not obviously an `A₅`/derived-series fact. What is *defensible*
is the resonance with the repo's field tower (S699m): `χ=2,3,4,…` adjoins Heegner rotations; the
*Galois groups* of those Heegner fields are abelian (class number 1) — **solvable**. If a chromatic
jump ever needed a *non-solvable* rotation field, that would be a genuine `A₅`-type wall. Recorded only
as a question; the verified content of this session is §§1–2 and the group theory.

## 7. Honest status

- **Proved/verified (THM-436):** the derived-series thresholds; `A_n` perfect for `n≥5` via overlapping
  3-cycles (the 5-point cause); the `C₅` tournament realisation; the FTA/monodromy dictionary.
- **Established analogy:** Abel–Ruffini ↔ THM-406 Vitali wall (both "no finite-depth tower");
  the `A₅` ↔ icosahedron ↔ S699h unit-distance bridge (classical + repo).
- **Conjecture (HYP-2303):** the LRC witness hierarchy = a radical tower; worry-set solvable ⟺ tight;
  residual = unsolvable core. **Not proved; resolves no open case.**
- **Speculative:** the HN "5" ↔ quintic "5" — a question, not a claim.

**Artifacts:** `04-computation/galois_solvability_tower_s703.py` (+`.out`). Theorem **THM-436**. New
**HYP-2303**. Builds on S699p/HYP-2282 (monodromy=Galois), S699l (FTA n↔n+1), THM-403 (cyclotomic
worry-set = round tournament), THM-406 (Vitali wall), THM-420/430 (witness hierarchy), S700 (residual),
S699h (A₅ icosahedral unit-distance), S699m (Heegner tower).
