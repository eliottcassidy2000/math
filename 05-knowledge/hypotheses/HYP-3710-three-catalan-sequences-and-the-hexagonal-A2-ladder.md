---
id: HYP-3710
title: The three nested sequences are the rank-1 (A1) CATALAN FAMILY -- seq3=Catalan(n) [A000108] 1,2,5,14,42; seq1=C(2n,n-1)=n*Catalan(n) [A001791] 1,4,15,56,210; seq2=C(2n+1,n-1)=(n/2)Catalan(n+1) [A002054] 1,5,21,84,330; with Catalan(n)+n*Catalan(n)=central binomial C(2n,n), and all three the near-central Pascal columns C(2n,n-j). The hexagonal covering (Kershner) is rank-2: the A2 (sl3) Weyl-chamber walk count = the A2 COXETER-CATALAN 2(3n)!/(n!(n+1)!(n+2)!) = 1,5,42,462,6006 [A005789]. So the Fibonacci/Catalan question and the hexagonal covering bridge are the SAME Coxeter-Catalan combinatorics ascending in rank (A1 sequences -> A2 hexagonal); Fibonacci itself = the shallow-diagonal (thin) reading of Pascal, the Catalan family = the near-central (fat) reading
status: VERIFIED (closed forms, the central-binomial decomposition, the Pascal-column identities, the A2 Coxeter-Catalan walk count). Identification + connections; the rate claim is corrected (the three grow ~4^n, not Fibonacci's phi^n -- they are nested PASCAL sequences, not phi-rate).
source: klein-2026-06-29-S25
depends_on:
  - HYP-3706   # the hexagonal/Eisenstein wallpaper bridge (A2 = the hexagonal lattice)
related:
  - HYP-3705   # the covering-min / projective-plane (the bridge this sequence-thread informs)
  - HYP-1902   # prior Fibonacci work: Zeckendorf boundary normal form
  - HYP-1920   # Zeckendorf shell certificate
  - HYP-2437   # Pisano-Goldbach
  - HYP-2998   # Farey-Fibonacci additive basis
results:
  - 04-computation/catalan_family_hexagonal_klein.py
  - 05-knowledge/results/catalan_family_hexagonal_klein.out
---

# HYP-3710 — the three nested sequences (Catalan family) and the hexagonal A2 ladder

## The three sequences (verified)
```
seq3 = Catalan(n)            [A000108]  1, 2, 5, 14, 42, 132   = C(2n,n)/(n+1)
seq1 = C(2n,n-1)             [A001791]  1, 4, 15, 56, 210, 792 = n * Catalan(n)
seq2 = C(2n+1,n-1)           [A002054]  1, 5, 21, 84, 330,1287 = (n/2) * Catalan(n+1)
```
They are the rank-1 (`A1`) **Catalan family** -- Catalan modulated by linear factors. Relations:
- **`Catalan(n) + n.Catalan(n) = C(2n,n)`** (the central binomial): `seq3 + seq1 = C(2n,n)`.
- the **near-central Pascal columns**: `seq1 = C(2n,n-1)`, `seq2 = C(2n,n-1)+C(2n,n-2) = C(2n+1,n-1)`;
  Catalan `= C(2n,n) - C(2n,n-1)` (the reflected/ballot count). Three adjacent diagonals of Pascal.

## The Fibonacci analogy (corrected)
Both Fibonacci and the Catalan family are read off **Pascal's triangle**, but they are DIFFERENT readings:
- **Fibonacci** `= sum_k C(n-1-k,k)` -- the SHALLOW (thin) diagonal SUMS; growth rate `phi = 1.618`.
- **the Catalan family** -- the NEAR-CENTRAL (fat) ENTRIES; growth rate `4`.
So the three sequences do NOT move at Fibonacci's rate (`~4^n` vs `~phi^n`); they are "nested" in the sense
of being three parallel near-central Pascal columns (advancing in lockstep, one row per term). Fibonacci is
the thin reading of Pascal; the Catalan family is the fat reading.

## The hexagonal covering is the RANK-2 (A2) member of the SAME family
The hexagonal/triangular lattice is the `A2` root lattice. The number of length-`3n` walks `0 -> 0` in the
`A2` (sl3) Weyl chamber (= multiplicity of the trivial rep in `V^{ox 3n}`) is
> `2(3n)! / (n!(n+1)!(n+2)!) = 1, 5, 42, 462, 6006` (`A005789`, the 3-dimensional / `A2` Coxeter-Catalan).
So the user's three sequences (rank-1 `A1` Catalan family) and the hexagonal covering bridge (rank-2 `A2`)
are the **same Coxeter-Catalan combinatorics ascending in rank**: `A1` (Dyck paths, the 2-fold) -> `A2`
(non-crossing pairs, the 3-fold / hexagonal). This dovetails with HYP-3706: the Singer multiplier there is
the **3-fold** rotation (the `A2`/hexagonal symmetry), and here the hexagonal count is the `A2` Catalan.

## Prior Fibonacci work in the project
The Fibonacci thread is the LRC additive-basis / Zeckendorf line: HYP-1902 (Zeckendorf boundary normal
form), HYP-1920 (Zeckendorf shell certificate), HYP-2437 (Pisano-Goldbach), HYP-2998 (Farey-Fibonacci
additive basis). These use Fibonacci/Zeckendorf representations of speeds; the Catalan family here is a
DIFFERENT (central-Pascal / Coxeter-Catalan) thread, connected to the hexagonal covering, not to Zeckendorf.

## The hexagonal covering optimality (the bridge's continuous side)
**Kershner 1939** (+ Fejes Toth): the hexagonal lattice is the THINNEST covering of the plane by congruent
disks, density `theta = 2pi/sqrt(27) = 1.2092`, with full `p6m` symmetry -- a THEOREM. So the continuous-side
optimality is settled; the open bridge (HYP-3706) remains: is the LRC covering the hexagonal one?

## Net
The three nested sequences are the rank-1 Catalan family (`Catalan`, `n.Catalan`, `(n/2)Catalan(n+1)`;
`Catalan + n.Catalan = C(2n,n)`); they are near-central Pascal columns (the "fat" reading, rate 4), where
Fibonacci is the "thin" diagonal reading (rate phi). The hexagonal covering bridge is the rank-2 `A2` member
of the SAME Coxeter-Catalan family (`A005789` = 1,5,42,462), tying the sequence question to the covering
work. Continuous-side covering-optimality = Kershner (settled); the LRC->hexagonal bridge stays open.
