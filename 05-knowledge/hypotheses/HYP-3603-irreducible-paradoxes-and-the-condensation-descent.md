---
id: HYP-3603
title: The IRREDUCIBLE PARADOXES are the STRONGLY CONNECTED tournaments (iso atoms 1,0,1,1,6,35 for n=1..6; the 3-cycle is the unique minimal one, n=2 has ZERO -- 2 things can't paradox), and the CONDENSATION (every tournament = a UNIQUE TOTAL ORDER of strongly-connected blocks = ORDER + one irreducible paradox per block) is the TOURNAMENT-SIDE DESCENT -- the exact structural analog of the LRC 2-adic descent (THM-580): both PEEL the orderable skeleton (the transitive condensation / the even 2-part) and EXPOSE the irreducible cores (the SC blocks / the odd Z_7-cores), finitizing a large object into a FINITE family of irreducible atoms where the true minimum lives (the SC tournament types / the doublet gap 4cos^2(3pi/7), klein-S17 HYP-3598)
status: COMPUTED + verified (SC counts & condensation block-compositions, n=1..6, exhaustive). The condensation<->descent parallel is a structural correspondence (both = peel-orderable/expose-irreducible). Complements klein-S17 (the LRC-apex finite family = 127 Z_7-cores); this is the tournament-side finite family (the SC atoms).
source: mac-mini-2026-06-30-S35
related:
  - HYP-3599  # mac-mini S34: tournaments = intransitivity among n things (the paradoxes ARE the intransitivity)
  - HYP-3598  # klein-S17: the descent's finite families = all 127 Z_7-cores (the LRC-apex side of this)
  - THM-580   # the LRC 2-adic descent (the move the condensation mirrors)
  - THM-590   # the apex odd-cycle gap (the minimum on the LRC finite family)
  - THM-588   # cyclicity = the intransitivity count; SC blocks carry it
results:
  - 04-computation/irreducible_paradoxes_and_descent_macmini_20260630.py
  - 05-knowledge/results/irreducible_paradoxes_and_descent_macmini_20260630.out
---

# HYP-3603 -- irreducible paradoxes and the condensation descent

Working the owner's three-part ask (more irreducible paradoxes; the descent; the finite families) on the
intransitivity reframe (HYP-3599). klein-S17 (HYP-3598) gave the LRC-apex finite family (the 127 Z_7-cores);
this is the tournament-side, and the structural link between them.

## The irreducible paradoxes = the strongly connected tournaments
A paradox (intransitivity) is **irreducible** if the things cannot be split into "these all beat those" --
i.e. the tournament is **strongly connected** (every thing reaches every other through the dominance
relation). Iso-type counts (exhaustive n=1..6):
> `#irreducible paradoxes = 1, 0, 1, 1, 6, 35` (n=1..6); continues `353, 6008, ...`
- `n=2`: **zero** -- two things cannot paradox (one beats the other, always orderable).
- `n=3`: **one** -- the 3-cycle, the unique minimal paradox (the Condorcet atom).
- The atoms are FEW and fully classifiable. (Lead: verify the iso sequence in OEIS.)
These are exactly the tournaments where the OCF's `H` is genuinely irreducible -- the transitive (orderable)
tournament has `H=1`; an SC tournament has `H>1` (Moon: it has a Hamiltonian cycle).

## The condensation IS the tournament-side descent
Every tournament has a **condensation**: a unique TOTAL ORDER of its strongly-connected blocks (the
condensation tournament on the blocks is transitive). So
> **tournament = (a ranking of blocks) + (an irreducible paradox inside each block) = ORDER + IRREDUCIBLE
> INTRANSITIVITY.**
Verified (n=4,5,6) by block-composition: the transitive tournament is the UNIQUE all-singletons composition
`(1,1,...,1)` (zero paradox); the fully-irreducible single block `(n)` accounts for the `6` (n=5) / `35`
(n=6) SC types; the rest mix. This is a **descent**: peel the orderable condensation (the block ranking),
expose the irreducible SC cores.

## The two descents are the same move
| | TOURNAMENT condensation | LRC 2-adic descent (THM-580, klein-S17) |
|---|---|---|
| object | a tournament | meas(lonely S), S an infinite covering family |
| PEEL (orderable) | the transitive condensation (block ranking) | the even 2-part `E/2` (the doubling skeleton) |
| EXPOSE (irreducible) | the strongly-connected blocks | the odd `Z_7`-cores `O_j` |
| finite family of atoms | SC tournament types `0,0,1,1,6,35,..` | the 127 nonempty `Z_7`-cores |
| the true minimum | the 3-cycle (minimal atom) | the doublet gap `4cos^2(3pi/7)>0` (THM-590) |
Both take a large/infinite object and **finitize it to a finite family of irreducible paradoxes**, where a
genuine minimum is attained (klein-S16: a provable infimum needs a finite family; the descent supplies it).
"Knowing the finite families" = knowing the irreducible cores: the SC atoms (tournament side), the 127
`Z_7`-cores (LRC apex side, klein-S17).

## Why this matters
- It gives **"the descent"** a clean conceptual identity: it is a CONDENSATION -- peel the orderable,
  expose the irreducible. The 2-adic `E/2` peeling is the LRC's version of taking the condensation.
- It gives **"the finite families"** a clean identity: the irreducible-paradox atoms. Finite and
  classifiable (SC tournaments; `Z_7`-cores).
- It explains why the proof lives there: the orderable part carries no obstruction (it is a coboundary,
  HYP-3599); ALL the content is in the irreducible cores, and there are finitely many, so a minimum exists
  and is attained at the smallest atom (the 3-cycle / the doublet `C_7`).
- Forward: classify the SC atoms by their cyclicity / `H` / spectral gap, and ask which SC atom is the
  tournament-side image of the binding LRC doublet (the apex `C_7`); the descent should land the binding
  covering on the minimal-gap irreducible core.

## What it buys
The irreducible paradoxes are named and counted (SC tournaments, atoms `0,0,1,1,6,35`); the descent is
identified as a condensation (peel orderable / expose irreducible); the finite families are the irreducible
cores on both sides (SC atoms; 127 `Z_7`-cores). The proof lives on the finite family of irreducible
paradoxes -- exactly where klein-S16/S17 and THM-590 put it.
