# Independent audit of the S49 packet (THM-4742, THM-4743, HYP-9353)

These are the independent auditor's own scripts and outputs. They share no code with the parent scripts, and they are kept verbatim as audit evidence.

They print raw tallies rather than the repository's `ALL CHECKS PASSED` line. The verdicts below are read off the tallies. The canonical, assertion-gated checks are the parent scripts `../seven_point_pencil.py` and `../converse_gluck_tournaments.py`, both of which print ALL CHECKS PASSED.

| script | output | what the tally shows |
|---|---|---|
| `geo_exact.py` | `geo_exact.out` | Exact rational simplices in dimensions `n = 2…6`. For `n ≥ 3` (affine rank 3), A,B,D,E,F,H' and B,C,D,E,F are each cospherical and all seven are not, with `Pow(C) = Σ`. This is the fact that refuted the pre-audit dimension claim (MISTAKE-650). |
| `geo_symbolic.py` | `geo_symbolic.out` | A symbolic re-derivation with general `R`, spheres solved from 4 points: `Σ = −8 v·(G − v)`, `μ = Σ/|s − 3v|²`, and `S1 − S2 = (Σ, −Σ, 0)`. |
| `geo_qi.py` | `geo_qi.out` | The `Q(i)` example `V3 = (0, i, i)`: seven labels on one quadric with only six distinct, and the median is tangent. |
| `geo_fields.py exh 5` | `geo_fields_exh5_audit_run.out` | Exhaustive `F_5` with `V0 = 0`. In the regime there are 368,640 configurations, with six-fail = five-fail = 0, seven ⟺ `Σ = 0` ⟺ `A = C` (no mismatches), and no seven distinct labels on a sphere. |
| `geo_fields.py rand …` | `geo_fields_rand.out` | Random runs at `p = 19, 23` with the same verdicts. |
| `tourn_checks.py` | `tourn_checks.out` | Paley groups computed by backtracking for `q ≤ 19`; the circulant census for odd `n ≤ 15`; primality of `C_3, T_5, P_7`; and a cheap counterexample hunt. |
| `tourn_sc_census.py 3 … 10`, `… 11` | `tourn_sc_census_n3_10.out`, `tourn_sc_census_n11_audit_run.out` | Every self-converse tournament with `n ≤ 11`, enumerated as the tournaments anti-invariant under a fixed involution with at most one fixed point. The labelled exceptional ones are exactly the copies of `C_3, T_5, P_7, C_3[C_3]`; there are none at `n = 8, 10, 11`. `n = 11` takes 896 s. |
| `tourn_vt.py <repo root>` | `tourn_vt.out` | The 110 vertex-transitive 21-tournaments (repository d6 files, parsed as data), with only `P_7[C_3]` lacking a found regular subset; all Cayley tournaments on `Z_3^3`, `Z_9 × Z_3` and `Z_27`, where every one without a found regular subset is isomorphic to `C_3^[3]`. |

`tlib.py` holds the shared numba tools. Runs used Python 3 with numpy, numba and sympy, launched from this directory.

HYP-9353 itself was then proved by the concurrent Session24 lane (`05-knowledge/results/converse_exception_classification_20261021/`). These censuses are its finite controls, not its premises.
