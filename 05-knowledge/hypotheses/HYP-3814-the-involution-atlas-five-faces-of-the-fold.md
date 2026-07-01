---
id: HYP-3814
title: THE INVOLUTION ATLAS -- the project runs on essentially TWO order-2 involutions of the staircase triangle (the FOLD = anti-diagonal reflection = complement/transpose/iota/score-reversal/conjugation, fixed points = SC tournaments / iota-symmetric lonely times; and the ANTIPODE = all-tile-bits flip = complement-tiling, no fixed points), forming a Klein four-group <fold, antipode> = Z_2 x Z_2. The deeper truth: an order-2 involution tau on a space has FIVE FACES, and every recurring object in the project is one of them seen through a lens: (1) FOLD = the quotient X/tau (merged metagraph, half-tiling, LRC half); (2) GRADE = the +/-1 eigenspaces (iota-even Eisenstein / iota-odd cusp; blue-odd / black-even; SC/NS); (3) FIX = the fixed-point set (SC classes, iota-symmetric witnesses, the diagonal spine); (4) DUAL = tau swaps paired objects (Platonic dual pairs, cut/cycle, tournament/arc); (5) GRAM = tau composed with the relation, the 2nd moment (danger D D^T, adjacency A, W(n)=H-variance). The identities BETWEEN faces are the deep connections: PARITY = trace(tau) = |Fix(tau)| (Burnside/Lefschetz -- e.g. SC-odd tiling count = odd #grid-sym = fixed points); SECOND MOMENT = the Gram of the relation (whose SPECTRUM separates 'invariant-cospectral' objects -- tournament spectral twins AND LRC spectral twins, both grounded); DUALITY = the involution on the dual (2 spectral twins <-> 2 dual Platonic pairs). Phi_6 ~ 4*(half-tiling dim) = the DOUBLE fold (triangle ~n^2/2 -> quarter-square ~n^2/4; Phi_6 ~ n^2 = the unfolded full square).
status: SYNTHESIS (a unifying frame + two grounded targets). GROUNDED this session: (a) LRC SPECTRAL TWINS exist -- 46/46 coarse-identical (same M + gcd-profile) covering 7-sets at n=8 are separated by the danger-Gram spectrum (the LRC analog of the S86 tournament spectral twins); (b) Phi_6 = 4*floor((n-1)^2/4) + (n+1) for even n = the double-fold identity. The five-faces framework organizes existing results (S82-S87 + LRC); not a new theorem, a canonical lens.
source: mac-mini-2026-07-01-S88
related:
  - HYP-3813   # S87 six bridges metagraph<->LRC (this generalizes Bridges 1,2,3,5 into the five-faces frame)
  - HYP-3811   # S86 tournament spectral twins (GRAM face, tournament side)
  - HYP-3571   # LRC floor D D^T (GRAM face, LRC side)
  - HYP-3809   # S84 parity theorem = FIX-face fixed-point count
  - HYP-3810   # S85 T-join SC-cover obstruction = FIX-face is the rigid core
  - THM-549    # the fold = complement = transpose = half-tiling mirror
  - THM-584    # complement R = the antipodal map of the arc-hypercube, (-1)^k level-parity (the GRADE face)
  - THM-587    # per-level signed cycle index = antipodal Euler number = Lefschetz fixed-point count of R (PARITY=|Fix|, CANON)
  - THM-582    # palindromic Ham-path count f = the R-odd coordinate, odd on SC (the FIX face)
results:
  - 04-computation/metagraph_lrc_integration_macmini_20260701.py
  - 05-knowledge/results/involution_atlas_targets_macmini_20260701.out
---

# HYP-3814 -- the involution atlas: five faces of the fold

Owner: how does seeing order-2 involutions from many lenses synthesize deeper truths? Answer: the project
runs on **two** fundamental order-2 involutions of the staircase triangle, and each involution has **five
faces**; every recurring object is one face seen through a lens, and the deep connections are the identities
between faces.

## The two fundamental involutions (Klein four on the staircase)
- **THE FOLD `sigma`** = the anti-diagonal reflection of the staircase (`(x,y) -> (n+1-y, n+1-x)`, THM-549).
  Lenses: **complement** `T->T^op` (reverse all arcs), **transpose** `A->A^T`, **score-reversal**
  `d->n-1-d`, the LRC **iota** `t->1-t`, **prime-2 antipode / odd-sin-kernel**, Eisenstein **conjugation**
  `omega<->omega-bar`. Fixed points: **SC tournaments** (tournament) / **iota-symmetric lonely times** (LRC).
- **THE ANTIPODE `flip`** = complement-tiling (flip all `m` tile bits) = the cube antipode. No fixed points.
  Generates the blue/black **lines**.
- `<sigma, flip> = Z_2 x Z_2` (S84); the third element `flip.sigma` is the twin-pairing candidate.

## The FIVE FACES of an order-2 involution `tau` (and where each lives)
| face | what it is | tournament | LRC |
|---|---|---|---|
| **FOLD** | quotient `X/tau` (halving) | merged metagraph `G_n/Z_2`; half-tiling (quarter-square) | the LRC half `[0,1/2]` |
| **GRADE** | `+/-1` eigenspaces (even/odd) | blue-odd / black-even graph; SC/NS | `iota`-even (`E_2` Eisenstein bulk) / `iota`-odd (cusp `f_14`) |
| **FIX** | fixed-point set (rigid core) | SC classes; the diagonal spine | `iota`-symmetric witnesses; the binding pair `{t*,1-t*}` |
| **DUAL** | `tau` swaps paired objects | Platonic dual pairs; cut/cycle; tournament/arc | `omega/omega-bar`; even/odd modular |
| **GRAM** | `tau` composed with the relation (2nd moment) | adjacency `A` (separates twins); `W(n)`=H-variance | danger `D D^T` (the floor, HYP-3571) |

## The identities BETWEEN faces (the deep truths)
1. **PARITY = `trace(tau)` = `|Fix(tau)|`** (Burnside/Lefschetz). The SC-odd tiling count (S84) IS the odd
   count of `sigma`-fixed (grid-sym) tilings; the Dedekind margin's `iota`-oddness IS the fold's modular
   anomaly. **Every parity in the project is a fixed-point count of the fold.**
2. **SECOND MOMENT = the GRAM of the relation**, whose SPECTRUM separates "invariant-cospectral" objects
   that the 1st moment (trace/count) cannot. GROUNDED both sides: tournament **spectral twins** (S86,
   separated by `A`'s spectrum) and **LRC spectral twins** (this session: 46/46 coarse-identical covering
   7-sets at `n=8` separated by the `D D^T` spectrum). Same move: `relation composed with itself`.
3. **DUALITY = the involution on the dual.** The `2` tournament spectral twins `<->` the `2` dual Platonic
   pairs (cube/octa, dodec/ico; tetra self-dual) -- same-data / dual-realization (S86).
4. **FIX is the RIGID CORE.** Covering the fold-fixed family is parity/rigidity-obstructed: SC-cover to the
   half-tiling dimension (S85 T-join) / the covering-min deep well (LRC). A self-symmetric family fills its
   invariant subspace and cannot be compressed.

## The DOUBLE fold (Phi_6 <-> half-tiling)
The triangle (`~n^2/2` cells) folds once to the **quarter-square** half-tiling (`floor((n-1)^2/4) ~ n^2/4`).
The LRC modulus `Phi_6 = n^2-n+1 ~ n^2` is the **unfolded full square** -- `Phi_6 = 4*floor((n-1)^2/4) +
(n+1)` for even `n` (`n=14`: `4*42+15 = 183`). So `Phi_6 ~ 4x(half-tiling dim)` = the DOUBLE fold: the
covering-min modulus is the un-halved-twice staircase. The two quadratic scales are one fold apart squared.

## Why this synthesizes deeper truths
Recognizing complement, `iota`, transpose, conjugation, score-reversal as ONE involution (the fold) lets us
**transport each face across lenses**: a parity computed on the tournament side (SC-odd) is the fixed-point
count that, on the LRC side, is the `iota`-anomaly (the Dedekind margin); a second moment on the LRC side
(the `D D^T` floor) is the Gram whose tournament avatar (adjacency `A`) separates the twins; a duality among
Platonic solids is the same swap as cut/cycle. The order-2 involution is the minimal symmetry, and its five
faces are the minimal complete vocabulary -- so a project built on one folded triangle keeps rediscovering
the same five objects, and the theorems that will cross domains are the identities between the faces.

## Canon grounding (S88 repo mine)
The five-faces frame is not new machinery -- it ORGANIZES existing canon, and the central identity is a
theorem:
- **FOLD**: `R = T->T^op` is the **antipodal map of the arc-hypercube** (THM-584), = the transpose /
  half-tiling `y=x` mirror (THM-549). `R` and the LRC `iota` are the SAME involution in different
  coordinates (THM-584 + `loneliness-is-antipodally-paired.md`).
- **GRADE**: `R` acts as `(-1)^k` on hypercube level `k` (THM-584): `R`-even = even levels (bulk,
  eigenvalues `= d mod 4`), `R`-odd = odd levels (signed, `d-2 mod 4`). This IS the even/odd grading; the
  LRC `iota`-even/`iota`-odd (Eisenstein `E_2` / cusp `f_14`) is the same split on the modular side.
- **FIX**: SC classes = the diagonal spine (THM-549); the R-odd coordinate is `f` = # **palindromic**
  Hamiltonian paths, ODD on SC (THM-582).
- **DUAL**: a SECOND Klein four -- the **Atkin-Lehner group** `<W_2, W_7> = V_4` on `X_0(14)`, with
  `W_14 = W_2 W_7 =` Fricke `=` complement; its 4 cusps `{1,2,7,14}` `<->` the 4 `n=4` tournament classes
  (`the-cusps-are-the-klein-group-hardness-is-the-genus.md`). So the metagraph Klein four `<sigma, flip>`
  and the modular Klein four `<W_2, W_7>` are two realizations of the same `Z_2 x Z_2`.
- **GRAM / PARITY**: **THM-587** (per-level signed cycle index = the **antipodal Euler number**) IS the
  **Lefschetz fixed-point count** of `R`. So **`PARITY = |Fix(R)|` is CANON**, and **Redei's theorem that
  `H(T)` is ODD is an instance**: `H` odd `=` the `R`-odd eigenspace has odd dimension `=` the Borsuk-Ulam
  `iota`-odd degree is nonzero. The project's foundational parity (H odd) and its hardest obstruction (the
  `iota`-odd cusp form / OPEN-Q-108) are the SAME face -- the fixed-point count of the fold.
- **ATTRACTOR**: the **Paley** tournament (SC + transitive `Z/pZ` + Galois-fixed) is where ALL the
  involutions converge; the apex prime `p = n/2` (7 for LRC14). Not an involution -- the point they fix.
- **INDEPENDENT** `Z_2`s (beyond `<fold, antipode>`): the **cut/cycle GF(2)** grading (score vs cycle
  space) and **Mode B** `n->n-2` (with a 3-fold Eisenstein `omega<->omega-bar` companion) -- these enlarge
  the symmetry beyond the Klein four.

## Next targets
1. Prove `PARITY = |Fix(fold)|` uniformly (the SC-odd-grid-sym crux, HYP-3809) -- a Lefschetz fixed-point
   count for the fold on the tiling hypercube.
2. Make the LRC "iota-symmetric witness" `<->` SC-class correspondence explicit (FIX-face across lenses).
3. Is the twin-pairing (S84) the third Klein element `flip.sigma`? (S86 ruled out a naive version; retry on
   the fold-fixed subcube.)
