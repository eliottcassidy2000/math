---
id: HYP-9450
title: "After THM-4859: (1) the (3,3,4;29) colourings are NOT edge-count rigid (some have a triangle-free colour with fewer than 116 edges, equivalently a vertex of colour degree below 8), so no forced-count argument reaches R(3,4,4) <= 76 or R(3,3,5) <= 56; (2) the (3,4,4;76) and (3,3,5;56) colourings do not exist, and their exclusion is a structured SAT problem whose first certificates are the ledgers' odd cross count and the 136/29 non-constancy; (3) the owner's targets 73, 137, 86, 198 need respectively the exclusion of (3,4,4;73..76), R(3,3,5) + R(3,4,4) <= 113, R(3,3,5) <= 51 and R(3,4,4) <= 66, none reachable by recursion, parity or mod-3 counting; (4) the cyclic (3,3,4;29) colouring is Kalbfleisch's and is the unique (3,3,4;29) colouring with both triangle-free colours 8-regular"
status: >
  OPEN (klein-2026-10-30). (1) expected from the (8,8,16) in-port caps with slack 4 (degree triples (8,8,12), (8,7,13),
  (7,7,14), ... are all admissible a priori); the decisive test is a SAT instance asking for a (3,3,4;29) colouring with a
  vertex of colour-0 degree <= 7 (scripts/port_probe_334.py, mode deg7) - see the update line below once the probe
  returns; a proof of rigidity would instead be a theorem of Codish type at order 29 and would make the forced-count
  method applicable one level up. (2) is the content of R(3,4,4) <= 76 and R(3,3,5) <= 56; the saturated-ports packet
  gives the (3,3,5;56) structure (group order dividing 224, port edge counts (104,104,198)), THM-4859 (C) the (3,4,4;76)
  structure; a complete search must retain both adjacency and root positions. (3) is exact arithmetic (THM-4859 (B)).
  (4) the circulant census is exhaustive (THM-4859 (D)); uniqueness among all (3,3,4;29) colourings with both
  triangle-free colours 8-regular is a guess (a (3,3,3;16)-type rigidity would be needed).
session: klein-2026-10-30 (codex/session-ramsey-20261030)
source: 05-knowledge/results/ramsey_cascade_barrier_and_port_rigidity_20261030/README.md (sections 3-6)
related:
  - 01-canon/theorems/THM-4859-saturated-port-cascade-for-multicolour-ramsey-upper-bounds-exhausted-on-the-2026-table-at-the-unsaturated-334-port-the-344-76-ledger-and-the-unique-circulant-334-29-colouring.md
  - 05-knowledge/results/ramsey_saturated_ports_20261030/PROOFS.md
---

# HYP-9450 — what stands between the table and the targets

1. **Port non-rigidity.** In a `(3,3,4;29)` colouring every vertex has colour degrees `(d_a, d_b, d_c)` with
   `d_a, d_b ≤ 8`, `d_c ≤ 16`, sum 28. The cyclic colouring has `(8, 8, 12)` everywhere. Prediction: colourings with a
   vertex of colour-`a` degree `7` (so `d_c ≥ 13`) exist, hence the edge count of a triangle-free colour is not forced;
   then neither `e_2(N_2(v))` in a `(3,4,4;76)` colouring nor the analogous counts in `(3,3,5;56)` are constants, and
   Boza's divisibility cannot be applied at those orders. Refutation would be the UNSAT answer to the `deg7` probe,
   which would itself be a rigidity theorem about `(3,3,4;29)`.
2. **Exclusion problems.** `R(3,4,4) ≤ 76` ⟺ no `(3,4,4;76)` colouring. Its forced structure: degrees `(17,29,29)`,
   Paley17 colour-1 ports, `e_1(N_1, N_2) = e_1(N_1, N_3) = 136`, `e_1(N_2(v)) = e_1(N_3(v))`, odd colour-1 cross count
   `≥ 125` between the 29-ports. A search anchored at one vertex has `1 + 17 + 29 + 29` labelled positions with the
   Paley port fixed up to its 136 automorphisms; what is missing is the `(3,3,4;29)` port bank. Prediction: the
   `(3,4,4;76)` instance is unsatisfiable and the first decisive sub-instance is the pair of 29-ports with a fixed
   colour-1 cross graph of odd size `≥ 125`.
3. **Target arithmetic** (exact, THM-4859 (B)): `73` needs direct exclusion at `76, 75, 74, 73`; `86` needs
   `R(3,3,5) ≤ 51` (lower bound 45); `137` needs `R(3,3,5) + R(3,4,4) ≤ 113` (now 134); `198` needs `R(3,4,4) ≤ 66`
   (lower bound 55). None of the four is reachable by recursion, parity or forced-clique counting, which deliver one
   unit each and are exhausted on the 2026 table.
4. **The cyclic colouring.** The 14 labelled cyclic `(3,3,4;29)` colourings form one class with both triangle-free
   colours at the cap 8 and the `K_4`-free colour 12-regular. Guess: it is the only `(3,3,4;29)` colouring with both
   triangle-free colours 8-regular (every vertex at the double cap). Test: the SAT instance with all colour-0 and
   colour-1 degrees fixed to 8 and a blocking clause for the cyclic solutions.
