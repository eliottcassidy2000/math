---
source: opus-2026-06-08-S708 (user: Explore covering systems / Owens min-modulus 42; Work On the
  min-degree-3 even-cycle question)
status: RESULT + EXPLORATION + bridge. The even-cycle question answers YES (THM-443, classical;
  min degree ≥3 ⟹ even cycle ≥4, via longest path + pigeonhole; even-cycle-free = odd cactus = min
  degree ≤2; exhaustively verified n≤7). Covering systems: Owens' min-modulus 42 is the largest KNOWN,
  Hough (2015) proved min-modulus is BOUNDED (≤10^16, later 616 BBMST) — so Erdős's "arbitrarily large
  min modulus?" is NO. The two prompts are one theme: COVERING / PARITY. LRC failure = the danger
  residues COVER ℤ/q for every q = a persistent (interval) covering system; the S704 unbounded depth q*
  = the covering never breaks; the worry-set = the tight (boundary) cover. The directed even-cycle
  problem ties to Pfaffian orientations / Pólya permanent (per=±det), bridging to S707; a strong
  tournament on ≥4 vertices has an even dicycle (Moon). HYP-2313.
tags: [even-cycle, min-degree, longest-path, cactus, odd-cycle, covering-system, owens, hough,
  erdos-min-modulus, lonely-runner, danger-covering, pfaffian-orientation, polya-permanent,
  even-dicycle, parity, cut-cycle, tournaments, moon-pancyclic]
---

# Even cycles, covering systems, and the one lens under both: parity-covering

The user paired two things — **explore** covering systems (Owens' minimum modulus 42) and **work on**
whether minimum degree ≥3 forces an even cycle. They are the same theme read in two domains:
*covering* (residues covering ℤ; danger arcs covering the circle; cycles covering the cycle space),
graded by *parity*.

## 1. The even-cycle question: YES, cleanly (THM-443)

Every finite graph with `δ≥3` has an even cycle of length `≥4`. The proof is the longest-path
argument: the endpoint of a longest path has all `≥3` of its neighbours on the path, three positions
fall into two parities, and the same-parity pair bounds an even cycle `(p_j−p_i)+2 ≥ 4`. The sharp
converse: **even-cycle-free ⟺ every block is an edge or an odd cycle (a cactus of odd cycles) ⟺ min
degree ≤2** — because any 2-connected non-cycle has a theta, and every theta has an even cycle (the
three path-sums `a+b,b+c,a+c` sum to even, so one is even). Exhaustively verified over all graphs on
`≤7` vertices: all 173 with `δ≥3` have an even cycle; even-cycle-free ones cap at `δ=2`.

## 2. Covering systems: Owens' 42 and Hough's wall

A **covering system** is a finite set of congruences `a_i (mod m_i)` (distinct `m_i ≥ 2`) covering
every integer. **Erdős's minimum-modulus problem:** can `min m_i` be arbitrarily large?
- **Owens (2014)** constructed a covering system with **minimum modulus 42** — the largest *known*
  minimum modulus (the best lower bound on the extremal value).
- **Hough (2015)** proved the answer is **NO**: every covering system has `min m_i ≤ 10^{16}`. The
  **distortion method** of Balister–Bollobás–Morris–Sahasrabudhe–Tiba lowered the bound to **616** and
  reproved Hough. So the extremal minimum modulus lies in `[42, 616]` — open.

The mechanism (Hough): the density a covering can place while keeping small moduli scarce is bounded;
the residues "run out of room." This is a **covering-saturation** phenomenon.

## 3. The bridge: LRC is a danger-covering system (the S704 wall in covering language)

For integer speeds `S` and rational witnesses (S704: the optimal tick is `t=m/q`), the **danger set**
of runner `v` at resolution `q`, threshold `1/n`, is the residue set
`D_v = { m ∈ ℤ/q : ‖v m/q‖ < 1/n }` — an *arc* of residues. A **lonely tick at resolution `q`** exists
iff `⋃_v D_v ≠ ℤ/q`. So:

> **LRC failure (`M(S) < 1/n`) ⟺ the danger arcs `D_v` COVER `ℤ/q` for *every* `q`** — a persistent
> (interval) covering system that never breaks. The S704 "**unbounded cyclotomic depth `q*`**" is
> exactly "the covering survives to arbitrarily fine `q`," and the **worry-set** is the *tight*
> covering (covers everything but a measure-zero boundary tick, where `M=1/n` is hit with equality —
> verified: AP₅ does not strictly cover `ℤ/5`).

This is a reframing, not a new theorem (the `D_v` are arcs, an *interval* covering system, not Erdős's
single-residue kind). But it puts LRC and covering systems in one language: **LRC asks whether a
danger-covering system can persist at all resolutions; Hough-type saturation is the analogue of "the
covering must eventually break."** The repo's Vitali/depth wall (THM-406/S704) is a
covering-saturation statement.

## 4. The directed twin: even dicycles = Pfaffian orientations (bridge to S707)

The **directed** even-cycle problem — does a digraph have a directed cycle of even length — is the
core of **Pólya's permanent problem**: a digraph has **no even directed cycle ⟺ it admits a Pfaffian
orientation ⟺ `per(A)=±det(A)`** (poly-time, Robertson–Seymour–Thomas / McCuaig). So the undirected
`δ≥3` theorem and last session's **Pfaffian** thread (S707: `Pf²=det(I+2A)`, max-H ⟺ `|Pf|=1`) are two
faces of the same even-cycle/parity object. In **tournaments**: a strong component of size `≥4` has an
even directed cycle (Moon — strong tournaments are vertex-pancyclic, lengths `3..n`); a transitive
tournament has none (verified n=3..6). So "even dicycle in a tournament" = "a strong component of size
≥4" = the cyclic (non-transitive) part — the same strong-component decomposition that makes `H`
multiplicative (S531).

## 5. The unifying parity-covering lens

> Three covering problems, one parity grading:
> - **integers:** residues covering ℤ (covering systems; Hough's wall);
> - **the circle:** danger arcs covering `ℝ/ℤ` (LRC; the S704 depth wall);
> - **the cycle space:** cycles covering / realising parities (even-cycle theorem; Pfaffian orientations).
>
> Each asks "can a covering persist / which parities are forced," and each has a saturation answer:
> covering systems must break (Hough); `δ≥3` forces an even cycle (THM-443); a Pfaffian orientation
> exists iff no even dicycle. The repo's **OCF / Rédei (odd structures)** and **even graphs `E_n`** are
> the tournament face of this parity grading.

## 6. Honest status

- **THM-443** is classical (clean proof recorded, characterization, exhaustive `n≤7`).
- **Covering systems:** Owens 42 / Hough ≤10^16 / BBMST 616 are cited facts (not my results).
- **New/repo-relevant:** the **LRC = persistent danger-covering** reframing of the S704 depth wall
  (a lens, not a theorem); the **even-dicycle ↔ Pfaffian-orientation** bridge tying THM-443 to S707;
  the parity-covering unification. **No open case resolved** (LRC, Erdős min-modulus both stay open).
- Recorded the lens as **HYP-2313** (the covering-saturation analogy between Hough's wall and the LRC
  depth wall), flagged as a framing to be made precise — does a Hough-type density/saturation argument
  bound the LRC cyclotomic depth `q*`?

**Artifacts:** `04-computation/even_cycle_mindeg3_and_covering_s708.py` (+`.out`). Theorem **THM-443**.
New **HYP-2313**. Builds on THM-443, S704 (depth wall), THM-406 (Vitali wall), S707/THM-174 (Pfaffian),
S531 (strong components / H), Owens/Hough/BBMST, Pólya/RST/McCuaig, Moon.
