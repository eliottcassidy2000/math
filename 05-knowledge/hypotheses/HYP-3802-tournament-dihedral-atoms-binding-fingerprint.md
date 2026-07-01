---
id: HYP-3802
title: TOURNAMENTS x DIHEDRAL x RUNNERS on the 7 roots of -1 -- the atom-tournament bridge + a refined BINDING-FINGERPRINT invariant tower, making the finish concrete. The tight AP core {1..n-1} (n=14) has lonely measure = 6 ATOMS at the UNITS (Z/14)* = {1,3,5,9,11,13}/14 = the PRIMITIVE 14th ROOTS OF UNITY (zeros of the cyclotomic Phi_14, degree phi(14)=6). ADD the iota-fixed point 7/14=1/2 (7=odd non-unit): the 7 points {1,3,..,13}/14 are the SEVENTH ROOTS OF -1 (z^7=-1), carrying a DIHEDRAL D_7 action (rotation 1/7 + reflection = antipode iota: t->1-t = runner COMPLEMENT = tournament CONVERSE). VERIFIED (n=14): the 6 atoms ARE the units (max min_v||vt||=1/14 exactly, tight); the measure's MOMENTS c_k = c_14(k)/phi(14) are RAMANUJAN SUMS (=> convergence with mac-mini/kps HYP-3793); its VERBLUNSKY coefficients are the clean LADDER |alpha_k|=1/(6-k) (0.167,0.200,0.250,0.333,0.500,1.000), terminating |alpha_5|=1 (6 atoms, Verblunsky-Geronimus). The QR/PALEY tournament on Z/7 (i->j iff i-j in QR(7)={1,2,4}) is the dihedral one: scores all 3 (regular), c3=14=n, N(OCF)=80, H(Redei Ham-paths)=189 (odd), |Aut|=21 (Frobenius 7:3); vs rotational {1,2,3}: N=59, H=175, |Aut|=7. TILING model (n=7): circulant => DIFFERENCE-STRIPED (all tiles of a given x-y share a bit) => grid-symmetric (BLUE), since the grid-sym transform preserves x-y. The prime 2 = antipode/converse (iota); the prime 7 = n/2 = the 7-VANISHING (THM-503) = the vertex set Z/7 of the tournament; QR(7) = the Paley orientation. REFINED INVARIANT TOWER (binding fingerprint): q*(S)=binding modulus=atom denominator; M(S)=k/q* (Farey rung); N_at(S)=atom count=Verblunsky termination; the Verblunsky ladder; the atom-tournament (OCF N, Ham H, |Aut|). THE FINISH: covering forces q*: n -> Phi6 (the CF second convergent t*=[0;n-1,n]), raising M from the floor 1/n (tight AP, non-covering) to n/Phi6 (construction, covering-min); the atoms move from the n-lattice (primitive n-th roots) to the Phi6-lattice (2 atoms at +-t*)
status: MIXED (all computations verified + established cyclotomy/OPUC/tournament facts + synthesis). VERIFIED (n=14): 6 atoms = units (exact, tight M=1/14); moments = Ramanujan sums c_14(k)/6 EXACTLY; Verblunsky ladder 1/(6-k) EXACTLY, |alpha_5|=1; Paley(7) N=80/H=189/|Aut|=21, rotational N=59/H=175/|Aut|=7 (both c3=14, regular, Redei-odd); circulant tilings difference-striped + grid-symmetric. The bridges (iota=converse, prime-7=vertex set=7-vanishing, QR=Paley) are structural identifications. HONEST: a synthesis + refined invariant set + one concrete finish-mechanism (covering raises q* from n to Phi6); NOT a new proof -- the for-all-covering-sets bound is still OPEN-Q-108. The Ramanujan-sum moments converge with HYP-3793.
source: klein-2026-07-01-S70
depends_on:
  - HYP-3801   # S69: Verblunsky/OPUC encoding + loop-map dictionary (this adds tournaments/dihedral/atoms)
  - HYP-3800   # S68: phase-residue p(w)=nw mod Phi6, the Phi6-lattice atoms
related:
  - HYP-3795   # opus-S13: AP runner-cloud Verblunsky is harmonic |alpha_j|=1/(n-1-j) + circle-map dictionary AGL(1)+PSL2+Szego (CONVERGES; this is the LONELY-measure/units version + the TOURNAMENT layer opus/mac-mini did not do)
  - HYP-3793   # mac-mini/kps: flat-extension moments = Ramanujan sums (VERIFIED here via the units measure)
  - HYP-3789   # mac-mini: covering-min moment relaxation; AP=6 atoms=units, construction=2 atoms (the two poles)
  - THM-503    # the 7-vanishing (prime n/2=7) = the tournament vertex set Z/7
  - THM-515    # lonely measure = singular series; the atoms are its support
  - HYP-3715   # t*=n/Phi6 hexagonal point (the construction's 2 atoms)
  - THM-002    # OCF (odd cycle collection formula) -- the atom-tournament's N; ties to the project core
  - N_maximization_paley.py  # Paley maximizes OCF among circulant tournaments (context for N=80)
results:
  - 04-computation/tournament_dihedral_atoms_klein.py
  - 05-knowledge/results/tournament_dihedral_atoms_klein.out
---

# HYP-3802 — tournaments x dihedral x runners on the 7 roots of -1; the binding-fingerprint tower

## The 7 points, the atoms, and the dihedral action
The tight AP core `{1..n-1}` (`n=14`) is the LRC boundary: `M = 1/14` exactly, and its lonely set at that
level is **6 atoms at the units `(Z/14)* = {1,3,5,9,11,13}/14`** = the **primitive 14th roots of unity**
(zeros of the cyclotomic polynomial `Phi_14`, degree `phi(14)=6`). Adding the `iota`-fixed point
`7/14 = 1/2` (the odd non-unit, `gcd(7,14)=7`) completes the **7 seventh-roots of `-1`** `{1,3,5,7,9,11,13}/14`
(`z^7 = -1`). These carry a **dihedral `D_7`** action: rotation by `1/7` and the reflection = the antipode
`iota: t -> 1-t` = the runner **complement** = the tournament **converse** `T^op`.

## Verified arithmetic (n=14)
- **Atoms = units** (`max_t min_{v<=13} ||vt|| = 1/14` exactly; argmax at `k/14`, `k in (Z/14)*`).
- **Moments = Ramanujan sums**: `c_k = (1/phi(n)) sum_{u in units} cos(2pi k u/n) = c_n(k)/phi(n)` — verified
  `+1, +1/6, -1/6, +1/6, -1/6, ...`. This is **mac-mini/kps HYP-3793 (flat-extension moments = Ramanujan
  sums), recovered from the OPUC/units side.**
- **Verblunsky ladder**: `|alpha_k| = 1/(6-k)` exactly (`1/6, 1/5, 1/4, 1/3, 1/2, 1`), terminating
  `|alpha_5| = 1` — the Verblunsky–Geronimus signature of exactly 6 atoms. The para-orthogonal polynomial IS
  the cyclotomic `Phi_14`.

## The atom-tournament (dihedral / QR / Paley)
On `Z/7` (differences of the 7 points), the natural dihedral orientation is the **QR / Paley tournament**
`i -> j iff i-j in QR(7) = {1,2,4}`:
| tournament | scores | c3 | N (OCF) | H (Redei Ham-paths) | \|Aut\| |
|---|---|---|---|---|---|
| **Paley/QR {1,2,4}** | all 3 (regular) | 14 | **80** | **189** (odd) | **21** (Frobenius 7:3) |
| Rotational {1,2,3} | all 3 (regular) | 14 | 59 | 175 (odd) | 7 (cyclic) |
Paley maximizes the OCF `N` and has the larger automorphism group — it is *the* dihedral tournament on the
7 roots of `-1`. Note `c3 = 14 = n`. Rédei parity holds (`H` odd) for both. This is the bridge back to the
project's core (OCF `N`, Rédei `H`): the LRC lonely-measure atoms carry a tournament, and its parity
invariants are the project's.

## The tiling model (n=7): circulant => striped => blue
In the staircase `delta_5` with base path `6->5->...->0`, a circulant tournament (connection set `S`) is
**difference-striped**: tile `(x,y)` is set iff `(x-y) in S`, so all tiles of a given difference share a
bit (Paley sets differences `{2,4}`, 8/15 tiles; rotational `{2,3}`, 9/15). The grid-symmetry transform
`(x,y) -> (n-y, n-x)` **preserves `x-y`**, so every circulant tiling is **grid-symmetric (BLUE)**.
Circulant/dihedral <=> striped <=> blue.

## The refined BINDING-FINGERPRINT invariant tower (making the finish concrete)
For a speed set `S`, the "binding fingerprint":
1. **`q*(S)`** = binding modulus = denominator of the deepest atom (= deepest witness modulus).
2. **`M(S) = k/q*`** = the gap (Farey rung, S60).
3. **`N_at(S)`** = atom count = Verblunsky termination index.
4. **Verblunsky ladder** `|alpha_j(S)|` (the recursive shape of the lonely measure).
5. **Atom-tournament** `T(S)` (dihedral/QR type) with `N` (OCF), `H` (Rédei), `|Aut|`.
6. **(cyclotomic case) moments = Ramanujan sums** `c_{q*}(k)/phi(q*)`.

**The two poles:**
| | `q*` | `M` | `N_at` | atoms | Verblunsky | covering? |
|---|---|---|---|---|---|---|
| tight AP `{1..n-1}` | `n=14` | `1/n` (floor) | `phi(n)=6` | primitive `n`-th roots | `1/(6-j)` | NO |
| construction `{1..n-2,n(n-1)}` | `Phi6=183` | `n/Phi6` (covmin) | `2` | `+-t*` | `->(cos, 1)` | YES |

## The finish, concretely
> **Covering forces the binding modulus `q*` up from `n` to `Phi6` — the second convergent of `t* =
> n/Phi6 = [0; n-1, n]` — which RAISES `M` from the floor `1/n` (the tight, non-covering AP) to `n/Phi6`
> (the covering-min).** The atoms move off the `n`-lattice (primitive `n`-th roots, 6 of them) onto the
> `Phi6`-lattice (2 atoms at `+-t*`).
The tight AP would beat the covering-min (`1/n < n/Phi6`) but is *not* a covering set (no multiple of `n`);
the covering condition (THM-523) pins any covering set's atoms to the `Phi6` phase-lattice (S68), forcing
`q* = Phi6` and `M >= n/Phi6`. Proving "covering => `q* = Phi6` => `M >= n/Phi6`" for ALL covering sets is
OPEN-Q-108, now phrased as an invariant-tower statement: **no covering set's binding fingerprint is deeper
(smaller `M`) than the construction's.**

## Honest status
All computations verified; the cyclotomy (units = primitive roots, moments = Ramanujan sums), OPUC
(Verblunsky ladder, termination), and tournament facts (Paley OCF/Rédei/Aut) are established, applied and
cross-linked here. The invariant tower is a refinement + unification; the finish-mechanism (covering raises
`q*` to `Phi6`) is the concrete reframe. NOT a new proof — the for-all-covering-sets bound is OPEN-Q-108.
