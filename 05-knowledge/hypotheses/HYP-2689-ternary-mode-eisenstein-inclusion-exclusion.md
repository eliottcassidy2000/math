---
id: HYP-2689
title: The seven-term recursion is inclusion-exclusion over three reduction modes (Mode A n→n−1, Mode B n→n−2, Mode C n→n−3), 2³−1=7 terms graded 3+3+1 by Eisenstein cube roots; Mode C is the ternary axis beside the binary Cayley–Dickson Mode B
status: OPEN (program). The 3+3+1 grading, the order-truncation, and the S_3/Eisenstein structure are PROVED/VERIFIED (THM-549, THM-551, HYP-2681); the Mode-C ternary tower and the half↔coverage shadow correspondence are conjectural.
source: mac-mini-2026-06-20-S4
related:
  - THM-549   # half-tiling 7-term recursion A+B−C+D−E−F+G
  - THM-550   # half-tiling parity recurrence (codex)
  - THM-551   # apex-prime order-truncation
  - HYP-2680  # codex Φ_s Stirling multi-far hierarchy
  - HYP-2681  # codex cube-root Eisenstein modes
  - HYP-2687  # kps gnomon shells / even-odd dichotomy
---

# HYP-2689 — Three modes, inclusion–exclusion, and the ternary axis

## The unifying principle (the verified core)
Both the half-tiling odd recursion `A+B−C+D−E−F+G` (THM-549/550) and the coverage one/two/three-far
recursion (THM-548/HYP-2680, codex's pair-tax shadow HYP-2681) have **2³−1 = 7 terms graded 3+3+1**:
three singletons, three pairs, one triple. This is the **inclusion–exclusion over three generators**.
- Coverage: the 3 generators are the three far runners `{u,v,w}`; the 7 packets are the 7 nonempty
  Newton mixed differences `Δ_S` (`S⊆{u,v,w}`). VERIFIED: `Σ_S Δ_S = p0(B∪{u,v,w})−p0(B)` exactly.
- Half-tiling: the 3 generators are the **three reduction depths** `A,B` (depth `n−1`), `D` (depth
  `n−2`); the 7 regions (3 corners, 3 edges, 1 center) are their inclusion–exclusion overlaps; the
  center `G` (depth `n−4`) is the triple overlap.
The `S_3` symmetry of the three generators is carried by the **cube roots of unity**: the Eisenstein
modes `S_ω=A+ωB+ω²C`, `P_ω=D+ωE+ωF` (HYP-2681) are the `C_3` characters; their moduli are
`S_3`-invariant (VERIFIED). The lone generator vs the conjugate pair `{ω,ω²}` is the depth-`n−2`
corner `D` vs the equal pair `A,B` (depth `n−1`).

## The conjecture: Mode C, the ternary reduction axis
CLAUDE.md's two reductions are **Mode A** (`n→n−1`, hypotenuse, `H=1+2^d`) and **Mode B** (`n→n−2`,
both legs, the binary **Cayley–Dickson** tower `R→C→H→O→S` at `n=2,3,5,9,17`). The three corners of
the folded triangle expose a **third axis, Mode C (`n→n−3`)**, the ternary partner organized by the
**Eisenstein integers `Z[ω]`** the way Mode B's doubling is organized by `Z[i]`. CONJECTURES:
- (a) The half-tiling order-4 recurrence `(x−1)³(x+1)` factors the three modes: `(x−1)³` = the three
  depth-shifts `A,B,C` superposed, `(x+1)` = the complement-fold parity (even/odd, square/pronic).
- (b) There is an Eisenstein analog of the Cayley–Dickson tower along `n→n−3`: a sequence of
  algebras/structures losing a property at each ternary step, with `Z[ω]`-grading, parallel to but
  distinct from the binary `n→n−2` octonion tower. Candidate `n`-ladder: `n→3n−c` for the
  appropriate `c` (to be pinned by which invariant is ternary-periodic — test on the metagraph
  `H`-spectrum and the even-graph dual `E_n`).
- (c) The composite center `G` (depth `n−4` = the three-far packet) is the **Mode A∘B∘C triple
  product**; the seven-term recursion is the lattice of subsets of `{A,B,C}` modes, so any tournament
  invariant with a depth-≤3 local recursion automatically obeys a seven-term inclusion–exclusion.

## The two sevens
`7 = 2³−1` (inclusion–exclusion over three generators, the recursion's term count) and `7 = ` the
apex prime (sector count, THM-551 order-truncation start `7−|B|`, kernel vanishing) are numerically
equal. Whether this is structural or a Mersenne coincidence (`2³−1`) is open and worth resolving: if
the three far runners that matter near the LRC cap are exactly forced by `7−|B|=3` at the dangerous
`|B|=4` rows, the two sevens are the same fact. (At `k=9` true-wide with `|B|=4` bounded core and
`r=5` far, the live leading order is `3` — three-far — consistent.)

## Tests to run
1. Pin Mode C's `n`-ladder by scanning metagraph `H` and `E_n` invariants for ternary (`n→n−3`)
   periodicity / `Z[ω]`-structure, as Mode B gives the octonion ladder.
2. Check whether the dangerous LRC rows are exactly those where `7−|B|=3` (three-far leading), tying
   `2³−1` to the apex prime.
3. Relate the three modes to kps's gnomon shells (HYP-2687): does each mode generate a gnomon family?
