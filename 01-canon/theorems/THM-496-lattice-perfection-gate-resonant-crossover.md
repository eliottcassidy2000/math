---
id: THM-496
name: The lattice-perfection gate — the two-factor resonant-product 3N-crossover is
      at n=28 because 9 is the FIRST lattice-imperfect size (Harborth(9)=16 < u(9)=18),
      so 27=3³ is doubly obstructed while 28=4·7 is lattice-perfect AND chord-bearing
status: RESERVED (stub) — monad-explorer-2026-06-13; proof + exact verification incoming
date: 2026-06-13
session: monad-explorer-2026-06-13 (deep-research; OPEN-Q-057 frontier)
depends_on:
  - THM-493   # resonant-product decomposition U = e(G)|H|+|G|e(H)+Δ_t
  - THM-495   # (concurrent monad-explorer) chord-spectrum bottleneck: Δ_t>0 needs shared chord norm t; R's chord {1,3} forces t=3
  - THM-432   # generic-angle product cap; products tie 3N at 27,30, first beat at 32
  - THM-437   # cube K₃^□3 angle-rigid at 81
  - THM-431   # u(n) exact table; N*∈[25,28]
external:
  - "Harborth (1974): max unit edges of an n-point triangular-lattice patch = ⌊3n−√(12n−3)⌋."
  - "Alexeev–Mixon–Parshall arXiv:2412.11914: exact u(n), n≤21."
---

# THM-496 (stub): the lattice-perfection gate

COMPLEMENTS the concurrent THM-495 (chord-spectrum bottleneck). THM-495 answers
*which t* crosses (R's only non-unit chord is √3 ⟹ t=3). This answers *why n=28
not n=27*, from the orthogonal "lattice-perfection" axis the chord argument omits.

**Key fact (verifiable):** a k-point triangular-lattice patch has ≤ Harborth(k)
unit edges, and Harborth(k) = u(k) (the planar max) **iff k ≤ 8**; **k=9 is the
first lattice-imperfect size** (u(9)=18 > 16=Harborth(9)). Hence a *resonant*
(lattice) product matches the *generic* Cartesian cap only when every factor size
is ≤ 8.

**The gate (to be proved/verified this session):** the 2-factor resonant family
first beats 3N at n=28 via the conjunction
  (i) lattice-perfect factorization (all parts ≤ 8),
  (ii) a chord-bearing factor (size ≥ 4; sizes 2,3 are chord-free), and
  (iii) Δ_t > gap(n) = 3n − P_gen(n).
n=27=3³ fails (i)+(ii) (forced imperfect 9-factor; chord-free 3-factor ⟹ Δ=0);
the exact resonant cap at 27 is **75 < 81**, so resonance HURTS at 27 — the 81-tie
is *generic/off-lattice*, not resonant. n=28=4·7 is the first n meeting all three.

Bound (incoming): resonant U(27) ≤ e(G₃)+78 ≤ 81 via Δ_t ≤ 8(3−e(G₃)); sharpens
to ≤ 75 with Harborth(9)=16.

Files: `04-computation/lattice_perfection_gate_monad.py` (incoming).
