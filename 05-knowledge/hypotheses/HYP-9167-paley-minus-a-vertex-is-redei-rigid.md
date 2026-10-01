---
id: HYP-9167
title: "Paley minus a vertex is Redei-rigid: for every prime power q = 3 (mod 4), every arc of QR_q minus a vertex lies on an odd number of Hamiltonian paths. More broadly, a tournament on N >= 3 vertices with every arc on an odd number of HPs exists iff N = 2 (mod 4)"
status: >
  OPEN. Proved parts (THM-4524):
  - the antipodal arcs {u, -u} of QR_q minus 0 are odd (an F_q^* orbit
    of odd size (q-1)/2, with sum_e c(e) = (q-2) H odd);
  - all-odd forces N = 1, 2 (mod 4).
  FINITE-EXACT support:
  - QR_q minus a vertex is all-odd for q = 7, 11, 19, 23, 27 (N = 6, 10,
    18, 22, 26);
  - all-odd classes: exactly one at N = 6, none at N = 5 and N = 9,
    exactly two at N = 10 (exhaustive censuses).
  First open cases of (b):
  - N = 13 (non-existence; N = 1 mod 4 is allowed by parity);
  - N = 14 (existence; 15 is not a prime power; none of the 128
    circulants on Z_15 minus a vertex is all-odd).
source: collatz-procgen-20260922 session, selfie lane (2026-10-01), Conjectures C2 and C5; promoted with THM-4524
related:
  - 01-canon/theorems/THM-4524-selfie-tournaments-loop-gauge-arc-parity-and-the-odd-mallows-sloane-partner.md
  - 05-knowledge/results/procgen_selfie_20261001_selfie_tournaments.md
---

# HYP-9167 — Paley minus a vertex is Rédei-rigid

**(a)** For every prime power `q = 3 (mod 4)`, the tournament `QR_q - 0` (`N = q - 1 = 2 mod 4`) has `c(e)` odd for every arc `e`. Here `c(e)` is the number of Hamiltonian paths through `e`.

**(b)** For `N >= 3`, some tournament on `N` vertices has every arc on an odd number of HPs iff `N = 2 (mod 4)`.

## Why it matters

Rédei's theorem says `H` is odd. Deleting an arc `e` changes `H` by `c(e)`, so "all arcs odd" means that no single arc can be shaved off while keeping Rédei's parity. These are the parity-rigid tournaments.

The opposite extreme is proved: every Cayley tournament on an odd abelian group, including `QR_q` itself, has every `c(e)` even (THM-4524, C1). So deleting one vertex from Paley turns all-even into all-odd. The problem is to understand that flip.

## What is known

- **The antipodal orbit is proved odd.** The multiplicative group `F_q^*` acts on the HPs of `QR_q - 0`: residues act as automorphisms, and non-residues act as anti-automorphisms followed by reversal. Every arc orbit has size `q - 1` except the antipodal orbit, whose size `(q-1)/2` is odd. Since `sum_e c(e) = (q-2) H` is odd, the antipodal arcs are odd.
- **End counts.** `end(u)` is odd exactly for the non-residues `u`. This is proved by arc-transitivity of `QR_q`; the orchestrator also checked it for `q = 7, 11, 19`.
- **The open part.** It remains to show the non-antipodal orbits are odd: `(q-3)/2` orbits of size `q - 1`. The evidence is the five values `q = 7, 11, 19, 23, 27`.

## Limits of the data

Exhaustive data reach only `N = 10`:
- `N = 9` has none (although parity allows it);
- `N = 10` has exactly two classes: `QR_11 - v` and one rigid class with `H = 3929`.

So (b) has thin evidence on the non-existence side at `N = 1 (mod 4)`.
