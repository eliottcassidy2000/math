# THM-165: The Coefficient-2 Theorem and Cycle-Disjointness Constraint

**Status:** VERIFIED (n=3-7 exhaustive)
**Session:** kind-pasteur-2026-03-13-S61

## Statement

### Part 1: Coefficient-2 for Hamiltonian Cycles

For an n-vertex tournament T, each directed Hamiltonian cycle (n-cycle) uses all n vertices, so it conflicts with every other directed odd cycle in Omega(T). Therefore, changing the number of directed Hamiltonian cycles by delta changes H(T) by exactly 2*delta.

**Within any fixed "sub-Hamiltonian class"** (tournaments sharing the same c3, c5, ..., c_{n-2} structure), we have:

    H(T) = constant + 2 * c_n_dir(T)

where c_n_dir counts directed Hamiltonian cycles from a fixed vertex.

### Part 2: Verified Instances

**n=5, score class (1,2,2,2,3):**
- H = 9 + 2 * c5_dir (c5_dir in {1,2,3})
- c3_dir = 4 (score-determined, constant in class)
- Residual H - 2*c5_dir = 9 for all 280 tournaments in this class

**n=7, regular tournaments (all scores = 3):**
- H = 141 + 2 * c7_dir (c7_dir in {15, 17, 24})
- c3_dir = 14, c5_dir varies (36, 28, 42)
- Residual H - 2*c7_dir = 141 for all 2640 regular tournaments

### Part 3: The Cycle-Disjointness Constraint

For regular n=7, the c5 coefficient is ZERO in the bridge equation because of a rigid constraint:

    c5_dir + 2 * disj_33 = 56

where disj_33 = number of vertex-disjoint 3-cycle pairs. Verified:
- c5_dir=36, disj_33=10: 36 + 20 = 56
- c5_dir=28, disj_33=14: 28 + 28 = 56
- c5_dir=42, disj_33=7: 42 + 14 = 56

This constraint means c5_dir and disj_33 exactly COMPENSATE each other in the OCF sum:

    H = 1 + 2*(c3_dir + c5_dir + c7_dir) + 4*disj_33
      = 1 + 2*c3_dir + 2*c5_dir + 2*c7_dir + 4*disj_33
      = 1 + 28 + (2*c5_dir + 4*disj_33) + 2*c7_dir
      = 1 + 28 + 112 + 2*c7_dir
      = 141 + 2*c7_dir

The quantity 2*c5_dir + 4*disj_33 = 2*(c5_dir + 2*disj_33) = 2*56 = 112 is CONSTANT.

## Proof Sketch (Part 1)

From OCF: H(T) = I(Omega(T), 2) = sum_{k>=0} alpha_k * 2^k.

An n-cycle C uses all n vertices, so C is adjacent to EVERY other cycle vertex in Omega. Adding C to Omega increases alpha_1 by 1 but cannot change alpha_k for k >= 2 (since C conflicts with everyone, no independent set of size >= 2 can include C). Therefore:

    I(Omega + C, 2) = I(Omega, 2) + 2

## The Tournament Uncertainty Principle (Concrete Form)

The constraint c5_dir + 2*disj_33 = 56 is a concrete instance of the "tournament uncertainty principle":

- **More 5-cycles** (c5_dir increases) forces **fewer disjoint 3-pairs** (disj_33 decreases)
- **Fewer 5-cycles** (c5_dir decreases) forces **more disjoint 3-pairs** (disj_33 increases)

This trade-off is EXACT and LINEAR. You cannot independently control cycle density and cycle disjointness within the regular tournament class.

The Paley tournament (c5_dir=42, disj_33=7, c7_dir=24, H=189) sits at the extreme of maximum cycle density and minimum disjointness — the BIBD point.

## OCF Correction

Previous OCF computations counted cycle VERTEX SETS as Omega vertices. The correct OCF has one Omega vertex per DIRECTED odd cycle. Two directed cycles are adjacent iff they share a tournament vertex.

Corrected directed cycle counts for regular n=7:
- H=171: nc=65 (14 c3 + 36 c5 + 15 c7), alpha_1=65, disj_33=10
- H=175: nc=59 (14 c3 + 28 c5 + 17 c7), alpha_1=59, disj_33=14
- H=189: nc=80 (14 c3 + 42 c5 + 24 c7), alpha_1=80, disj_33=7

## Verification

- ocf_directed_fix.py: H = I(Omega_directed, 2) verified for ALL n=3,4,5,6 (0 mismatches)
- n7_regular_cycle_spectrum.py: H = 141 + 2*c7_dir verified for all 2640 regular n=7 tournaments
- All results in 05-knowledge/results/
