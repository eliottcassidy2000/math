---
id: THM-4484
title: "Free and sporadic cycles: for y -> y/2, (qy+d)/2 a mixed cycle shape (p,a) is free (every parity word integral) iff (2^p - q^a) | d; |2^p - q^a| = 1 only for a = 1 or 3^2 - 2^3; so the free integer cycles of 3x+1 on Z are exactly {0}, {-1}, {1,2}, {-5,-7,-10} (Gersonides: 2-1, 3-2, 4-3, 9-8) and -17 is sporadic; the Belaga-Mignotte 3x+d table is reproduced exactly"
status: >
  PROVED (elementary) + INDEPENDENTLY AUDITED; FINITE-EXACT; CITED
  (Pillai, Ellison; Belaga-Mignotte 2006 data).
  Let T_{q,d}(y) = y/2 (y even) or (qy+d)/2 (y odd), with q >= 3 odd, d odd,
  and gcd(q,d) = 1. A parity word w of length p with a ones has periodic
  point d c_w/(2^p - q^a).
  (1) Shift criterion. A mixed shape (1 <= a <= p-1) is free, i.e. every
  word of that shape gives an integer, iff (2^p - q^a) | d. Moving the
  last 1 one place right changes c_w by exactly 2^j. So each d has only
  finitely many free shapes (Pillai; Ellison, effective).
  (2) Gersonides for every odd q. |2^p - q^a| = 1 with p, a >= 1 iff a = 1
  and q = 2^p +- 1, or (q,p,a) = (3,3,2). The proof is elementary.
  (3) Primitive decomposition. Cycles of T_{q,d} are the union over e | d
  of e times the primitive cycles of T_{q,d/e}. On a clock with
  g = gcd(D, d), the integral words are those with (D/g) | c_w. A clock is
  free if g = |D|, partially free if 1 < g < |D|, and sporadic if g = 1 < |D|.
  (4) Consequences. The free integer cycles of 3x+1 on Z are exactly {0},
  {-1}, {1,2}, {-5,-7,-10}. The fifth known cycle {-17,...} is sporadic:
  shape (11,7), gap 139, one integral necklace of 30. The free cycles of
  5x+1 are exactly {0} and {-1,-2}; its positive cycles {1,...}, {13,...},
  {17,...} are sporadic, with gaps 7, 3, 3. The densities of the five
  3x+1 cycles, 0, 1/2 | 1, 2/3, 7/11, are the first best approximations of
  log_3 2. This pattern holds for 3x+-1 and 5x+-1 but fails for general
  3x+d: only 22.5% of primitive cycles lie on the Stern-Brocot path
  (EMPIRICAL).
  (5) FINITE-EXACT. One uniform search (every d <= 19999 prime to 6, every
  least element <= 1200 d) reproduces all eleven entries of Belaga-Mignotte
  2006 table (20), and their counts of systems with omega = 1 and 2.
  Earlier the gates lane was off by one for d = 14303 and 17021; this is
  RESOLVED. It had missed two long primitive cycles:
  d = 14303, least element 101, clock (2155, 1092);
  d = 17021, least element 5, clock (2140, 1088).
  A residual of 42757 against 42765 total cycles, and 1004 against 1005
  systems with omega = 3, is UNRESOLVED; it is two-sided, so either their
  list or their summary is off.
  Censuses: 3x+1 on Z has exactly its five cycles for periods p <= 90;
  5x+1 on Z exactly five cycles for p <= 60; 3x+-1 exactly the known cycles
  for least elements <= 10^9.
  Collatz is OPEN. That the sporadic list is finite is Lagarias's open
  conjecture.
source: collatz-procgen-20260922 session. The orchestrator's wave-13 findings (section 1: Gersonides free cycles, Stern-Brocot position) were independently re-derived and generalized by the sporadic lane (2026-09-26), which also resolved the gates lane's Belaga-Mignotte discrepancy. Audited and promoted by the session orchestrator 2026-09-26.
depends_on:
  - 01-canon/theorems/THM-4471-kawasaki-fixed-point-collatz-proof-refuted.md (periodic points c_w/(2^p - 3^a))
related:
  - 05-knowledge/results/procgen_wave13_20260926_orchestrator_findings.md (section 1)
  - 05-knowledge/results/procgen_gates_20260925_gate_equidistribution.md (corrected)
  - 01-canon/theorems/THM-4479-strategy-cube-distance-to-provability-sharp.md (every provable strategy must break -1, -5, -17)
note: 05-knowledge/results/procgen_sporadic_20260926_free_and_sporadic_cycles.md
scripts: 04-computation/experiments/procgen_sporadic_20260926_{run,lib}.py, procgen_sporadic_20260926_{traj,sweep}.c (sha256 in the note)
script_audit: 04-computation/experiments/procgen_sporadic_20260926_orchestrator_check.py (and the wave-13 procgen_wave13_20260926_gersonides_check.py)
output: 05-knowledge/results/procgen_sporadic_20260926.out
output_audit: 05-knowledge/results/procgen_sporadic_20260926_orchestrator_check.out
output_sha256: f8668f5a272a9b133d93277fef2d787c2af048cace37403ef9b82c6a105db205
hash_basis: raw bytes
audit: >
  Two derivations agree. The orchestrator proved (1)-(2) for q = 3 in wave
  13 by the same shift argument plus Levi ben Gershon. The lane re-proved
  and generalized them without using the orchestrator's code. Independent
  code (procgen_sporadic_20260926_orchestrator_check.py) confirms:
  * the two long cycles: periods 2155 and 2140, odd counts 1092 and 1088,
    least elements 101 and 5, primitive;
  * 2^27 - 3^17 = 355*14303 and 2^65 - 3^41 = 17021*19*29*44835377399;
  * the shift criterion exhaustively for q = 3, 5, 7, d <= 139, p <= 12;
  * |2^p - q^a| = 1 only for a = 1 or (3,3,2), for q <= 201.
  The wave-13 Gersonides check confirms the free-cycle classification for
  3x+1 and the census for p <= 20.
  The lane's full pipeline was re-run (868 s). Its output is identical up to
  timing and RSS lines (70 checks). Script hashes match the note.
---

# THM-4484 -- free and sporadic cycles

**PROVED + INDEPENDENTLY AUDITED.** Full note: [procgen_sporadic_20260926_free_and_sporadic_cycles](../../05-knowledge/results/procgen_sporadic_20260926_free_and_sporadic_cycles.md); the `q = 3` case first appeared in the [wave-13 findings](../../05-knowledge/results/procgen_wave13_20260926_orchestrator_findings.md) §1.

## 1. The shift criterion

* **Setting.** For `T_{q,d}` the periodic point of a word `w` of shape `(p,a)` is `d c_w/(2^p - q^a)`, where `c_w` accumulates `c <- q c + 2^j` at the odd letters `j`.
* **The shift.** Take a mixed shape, `1 <= a <= p-1`, and a word whose last `1` sits at a position `j <= p-2`. Moving that `1` one place right keeps the shape and changes `c_w` by exactly `2^j`.
* **The criterion.** If both words were integral, the odd number `(2^p - q^a)/gcd(2^p - q^a, d)` would divide `2^j`, so it is 1. Hence a mixed shape is free iff `(2^p - q^a) | d`.
* **The case `d = 1`.** The criterion needs `|2^p - q^a| = 1`. That happens only for `a = 1` with `q = 2^p ± 1`, or for `3^2 - 2^3 = 1` (elementary).

## 2. The Collatz case: four free cycles and one sporadic

| cycle | shape `(p,a)` | `2^p - 3^a` | density | type |
|---|---|---|---|---|
| `{0}` | `(1,0)` | `1` | `0` | free (single word) |
| `{1, 2}` | `(2,1)` | `1` | `1/2` | free (`4 - 3`) |
| `{-1}` | `(1,1)` | `-1` | `1` | free (single word; `3 - 2`) |
| `{-5,-7,-10}` | `(3,2)` | `-1` | `2/3` | free (`9 - 8`) |
| `{-17, ...}` | `(11,7)` | `-139` | `7/11` | sporadic (1 of 30 necklaces) |

* **Stern–Brocot.** The densities are the first best approximations of `log_3 2` from below (`0, 1/2`) and from above (`1, 2/3, 7/11`).
* **The Kuratowski-style reading** (ANALOGY; see THM-4482's and the wave-13 findings' dictionaries):
  * the expanding obstructions every provable strategy must break are two free cycles, forced by an identity, and one sporadic cycle;
  * whether the sporadic list is complete is the `3x-1` cycle conjecture.

## 3. The Belaga–Mignotte table

With the two long primitive cycles (periods 2155 and 2140, near odd density `1/2`, far from the critical line), a uniform search reproduces all eleven entries of Belaga–Mignotte's table (20). This resolves the gates lane's off-by-one. Long cycles are common: 2204 primitive cycles have `p > 600`, and for 683 values of `d` the only primitive cycle is a long one.
