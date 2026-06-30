---
id: HYP-3563
title: The R-odd first Betti number b_1^-(n) of the arc-flip metagraph cycle space = 0,1,7,119,1772 (n=3..7), a NEW invariant; BUT the two apex-7/octonion conjectures about it are REFUTED -- (1) 7 | b_1^-(n) does NOT persist for n>=7 (holds n=5,6 by coincidence, FAILS at n=7: 1772=2^2*443), and (2) the 7 R-odd cycles at n=5 carry NO Singer-Z_7/Fano structure (Aut(G_5)=Z_2 only; H_1^- is 7-dim over Q so no Fano GF(2)-dependency). b_1^-(5)=7 is a dimensional coincidence; the genuine octonion/Fano structure lives at Paley T_7 (THM-586/HYP-3547), not here.
status: REFUTED (both conjectures) + a NEW invariant computed. VERIFIED: b_1^-=0,1,7,119,1772 (n=3..7) by corrected Lefschetz AND direct R_* eigenspace (n<=6 agree); E(G_n)=1,5,30,290,4086 matches mac-mini's independent edge counts; Aut(G_5)=Z_2 by exhaustive automorphism search; b_1^- sequence not in OEIS.
source: klein-2026-06-29-S6
depends_on:
  - THM-584   # complement R = graph automorphism of the metagraph (so it acts on H_1)
  - THM-588   # the arc-flip metagraph (E, V, the d-regular structure)
related:
  - HYP-3544   # "klein's R-odd Betti" -- this is its concrete metagraph realization (vs the level-graded complex)
  - HYP-3547   # mac-mini: the REAL octonion/Fano structure (QR{1,2,4} on Z_7, Paley T_7) -- not b_1^-(5)
  - THM-586    # Paley T_7 arithmetic (where the Singer-Z_7/Fano actually lives)
results:
  - 04-computation/metagraph_Rodd_first_betti_klein.py
  - 05-knowledge/results/metagraph_Rodd_first_betti_klein.out
  - 04-computation/metagraph_n5_Rodd_fano_test_klein.py
  - 05-knowledge/results/metagraph_n5_Rodd_fano_test_klein.out
---

# HYP-3563 — b_1^-(n) of the metagraph, and the refutation of its octonion/7-divisibility conjectures

## The object (new invariant)

The arc-flip metagraph `G_n` (vertices = iso classes A000568(n); simple edge i~j iff a single
dominance reversal sends a rep of i to j) is connected, and the complement involution `R` is a graph
automorphism (THM-584), so it acts on `H_1(G_n;Q) = H_1^+ (+) H_1^-`. Define **`b_1^-(n) = dim H_1^-`**
(the R-odd / antipodal-odd first Betti number — the concrete realization of "klein's R-odd Betti").

By Lefschetz (connected graph, `R` an involution fixing exactly the SC classes; a single-flip edge is
R-fixed only if both endpoints SC (+1) or it joins a complement-pair `{c, R(c)}` that is single-flip
adjacent (-1)):
```
b_1^-(n) = (E - V + SC - E_SCSC + E_comp)/2,
```
`E` = #edges, `V` = A000568(n), `SC` = self-converse count, `E_SCSC` = #edges with both ends SC,
`E_comp` = #edges joining a class to its complement.

## Computed values (verified)

| n | V | E | SC | E_SCSC | E_comp | b_1^- | factor | 7 \| b_1^-? |
|---|---|---|----|----|----|----|----|----|
| 3 | 2 | 1 | 2 | 1 | 0 | 0 | 0 | (triv) |
| 4 | 4 | 5 | 2 | 1 | 0 | 1 | 1 | no |
| 5 | 12 | 30 | 8 | 12 | 0 | **7** | 7 | **YES** |
| 6 | 56 | 290 | 12 | 13 | 5 | **119** | 7·17 | **YES** |
| 7 | 456 | 4086 | 88 | 174 | 0 | **1772** | 2²·443 | **NO** |

The corrected Lefschetz matches a direct `R_*`-eigenspace computation at n=5 (7) and n=6 (119); `E(G_n)`
matches mac-mini's independent `E(G_n)=1,5,30,290,4086`. `E_comp` is 0 for odd n (5,7), 5 for even n=6
(a parity feature: at odd n no class is one flip from its own complement). The sequence `0,1,7,119,1772`
is not in OEIS.

## Task (1) — REFUTED: 7-divisibility does not persist

`7 | b_1^-(n)` holds at n=5 (`7`) and n=6 (`119 = 7·17`) — likely what suggested an apex-7 pattern — but
**fails at the very first n>=7**: `b_1^-(7) = 1772 = 2²·443`, `≡ 1 (mod 7)`. So the divisibility is a
low-n coincidence, not a structural persistence.

## Task (2) — REFUTED: no Singer-Z_7 / Fano structure at n=5

`b_1^-(5) = 7` is genuine and triple-confirmed, but the 7 R-odd cycles carry **no** Singer-Z_7/Fano
(octonion) structure:
- `Aut(G_5) = Z_2` (exhaustive search) — only the complement involution itself; **no order-7 element**,
  so no Singer-Z_7 acting on the 12 vertices (a 7-orbit + 5 fixed would be needed; absent).
- `H_1^-` is 7-dimensional over `Q`, so any 7 generating cycles are a basis — there is no GF(2)
  Fano-line dependency among them (the Fano plane needs 7 vectors in a 3-dim GF(2) space).

So `b_1^-(5) = 7` is a **dimensional coincidence**. The real octonion/Fano structure of the project (the
7 Fano lines = translates of QR `{1,2,4}` on `Z_7`, the apex prime; HYP-3547/THM-586) lives at the
**Paley tournament `T_7`**, governed by the prime 7 — unrelated to the homology dimension `b_1^-(5)`. The
two 7's (apex prime vs Betti dimension) are not the same 7.

## Value retained

Despite both conjectures failing, `b_1^-(n)` is a clean new metagraph invariant with an exact
edge-counting formula, and this rigorously CLOSES the tempting "b_1^- = octonions" lead (per the project's
dead-end-documentation norm): the octonion structure is not in the metagraph's R-odd homology. Future
octonion work should target `T_7`/`Z_7` (THM-586), not `b_1^-`.

Open (minor): a closed form / asymptotics for `b_1^-(n)` (it is ~ E/2 ~ (1/2)·#edges, dominated by the
metagraph's edge growth); and whether `b_1^-` relates to any known tournament-homology Betti (THM-130/154).
