# Independent checks of the four SRFP311T1 preprints (2026-10-21, mac-mini-2026-10-21-fanopsp verifier)

The verdicts are summarised in `05-knowledge/results/lattices_codes_tournaments_20261021.md` §11. The primary referee report is Session27's `session27-papers-review/REVIEW.md`. No PDFs, extracted manuscript texts or third-party papers are included here.

| script | output | checks |
|---|---|---|
| `a_extremal_theta.py`, `.gp` | `a_extremal_theta*.out` | extremal theta series in dims 24/48/72/96; `q⁶` in dim 96 (two methods); Siegel `E_48` averages |
| `b_paper2_dp.py`, `b_paper2_bruteforce_check.py` | `b_paper2_*.out` | ADE census to rank 43 and the unique rank-43 collision; collision-free to rank 32 by explicit enumeration |
| `b_paper2_formulas.py`, `b_O5_bruteforce.py` | `b_paper2_formulas.out`, `b_O5_bruteforce.out` | `c_1..c_4` = log-EGF coefficients of ordered orthogonal root-frame counts; printed `c_5` mismatch; true `l_5`; brute-force `O_5`, `O_6` |
| `c_multiplicities.py`, `c_niemeier*.py` | `c_*.out` | `A_5`/`D_5`/`D_4` subsystem counts; Niemeier Table 2.1, Thm 4.11; genus-3/4 checks; `a(I_4)` separation |
| `d_leech.py` | `d_leech.out` | Leech minimal vectors, 11-design moments, `Λ/2Λ` classes, Voronoi facets, the weighted-design ratio `243/1024`, `|Co_0|` |
| `e_checks.py` | `e_checks.out` | `E_8` as a rank-4 Hermitian `O_F`-lattice; Schur identity; Galois subgroup; Lagrangian counts over `F_3, F_5, F_7`; union-bound scaling |
| `e_reconstruct_J.py`, `e_search_J.py`, `e_build_Lraw.py`, `e_rank2_search.py`, `dump_norm24.gp` | `e_*.out`, `HJ_cand*.json`, `Lraw_cand*.gp`, `Lraw_cand2_norm24_vector.txt` | reconstruction of the Leech Hermitian structure `J` (128 twists, all `7·`Leech); base lattice `E_8 ⊗ J`; explicit norm-24 vector |
| `w24_independent_check.gp` | `w24_independent_check.out` | independent re-check of the norm-24 witness: 96×96, even, positive definite, `det = 7^48`, `v^T G v = 24` |
| `f_checks.py` | `f_checks.out` | score vectors = permutohedron points = forests (`n ≤ 6`); Voronoi cell of `A_n^*` (`n ≤ 5`); QR codes and Construction A; cannonball to `10^6` |

**Status labels.**
* FINITE-EXACT: all integer computations.
* HEURISTIC: the Gaussian and union-bound estimates.
* The norm-24 vector lives in a reconstruction, because the paper's twisting unit is undefined. See §11 of the note.
