# Script Results Index

Maps every script in `04-computation/` to its stored output.
Results are stored in `05-knowledge/results/` with matching filenames (`.out` extension).

**Convention:** When you run a script, save its output:
```bash
python3 04-computation/SCRIPT.py 2>&1 | tee 05-knowledge/results/SCRIPT.out
```

---

## Results catalog

Scripts with stored results are marked [STORED]. Scripts without results are marked [MISSING].

### High-priority results to capture
These scripts have been run but results exist only in /tmp or session transcripts:

| Script | Status | Key finding | Result file |
|--------|--------|------------|-------------|
| `dc_symmetry_path.py` | /tmp only | DC verified; M(T\e) NOT symmetric | MISSING |
| `deletion_contraction_M.py` | /tmp only | M submatrix != M(T-v) | MISSING |
| `M_n7_structure.py` | /tmp only | Diagonal all odd at n=7 | MISSING |
| `reversal_proof_attempt.py` | /tmp only | M(T) != M(T^op), M(T) != M(T^op)^T | MISSING |
| `sixteen_dimensions_s93d.py`, `cartan_attention_theorem.py` | [STORED] | Corrected full/traceless dark-mode ratios and norm-ladder increment `4k+3` | `cartan_dimension_ladder_s95.md` |
| `meta_blindspot_probe_s95.py` | [STORED] | Under-bridged topic pairs; exhaustive `n=6` H/FAS correlation; Royle-even/dark acyclic-orientation overlap | `meta_blindspot_probe_s95.out` |
| `meta_flashlight_s95.py` | [STORED] | Bare Cartan decomposition of 0/1 tournament adjacency is inert; OCF lift `Omega(T)` is the variable symmetric sector separating residual `t3` fibers through n<=6 | `meta_flashlight_s95.out` |
| `tiling_hamiltonian_ratios_s95.py` | [STORED] | Proves and verifies `F(C)=H(C)/Aut(C)` for explorer fixed-path tilings; records ratio spectra, same-class pair moments, silent edges, and edge energies for n=3..7 | `tiling_hamiltonian_ratios_s95.out` |
| `recursive_insertion_exposure_s95.py` | [STORED] | Vertex insertion formula using bridge-exposure boundary state `Q_T`; verifies endpoint/scalar states are insufficient and `Q_T` gives zero-failure insertion responses for n=2..5 | `recursive_insertion_exposure_s95.out` |
| `merged_tiling_bucket_constraints_s95.py` | [STORED] | Derives merged tiling bucket parity: unmerged fibers are odd; SC merged nodes are odd; NSC merged nodes are exactly `2*odd`; weighted cube-edge incidences satisfy `2λ+Στ=mM` for n=3..7 | `merged_tiling_bucket_constraints_s95.out` |
| `endpoint_transfer_bucket_recursion_s95.py` | [STORED] | Builds exact endpoint-insertion transfer matrices between fixed-path quotient levels; verifies row/column bucket sums, SC parity boundary, and full GF(2) row rank for transitions n=2..6 | `endpoint_transfer_bucket_recursion_s95.out` |
| `even_graph_endpoint_transfer_s95.py` | [STORED] | Builds the endpoint-transfer analogue for the even-graph quotient; row/column sums still hold, but GF(2) ranks `[1,1,2,6,8]` show rank defects unlike the tournament quotient | `even_graph_endpoint_transfer_s95.out` |
| `goodcut_scc_defect_s354.py` | [STORED] | Verifies THM-354: for every Hamiltonian base path, good-cut count equals `n - #SCC(T)`; exhaustive labeled n<=6 and fixed-base tilings n<=7 | `goodcut_scc_defect_s354.out` |
| `residue_rank_probe_s355.py` | [STORED] | Compares transitive, Paley/Interval, H=63 single-core, and THM-025 examples by max-loss deletion-residue rank; separates exact kills (`rank_res=0`) from THM-025 near-kill (`rank_res=2`) | `residue_rank_probe_s355.out` |
| `tournament_tda.py` | [STORED] | Demo now emits SCC residue features and deletion-residue rank features in the `omega_*` block | `tournament_tda_residue_features_s355.out` |
| `endpoint_transfer_witnesses_s95.py` | [STORED] | Probes endpoint-transfer rank mechanisms; unmerged tournament rows all have private odd child witnesses through `6->7`, merged rows keep full matching/rank without all private witnesses, and even graphs show GF(2) cancellation | `endpoint_transfer_witnesses_s95.out` |
| `endpoint_transfer_smith_s95.py` | [STORED] | Computes Smith normal forms for small endpoint-transfer matrices; tournament/merged factors are all odd through `5->6`, while even graphs show 2-primary factors up to `[2,2,2,2,2,4,4,8]` at `6->7` | `endpoint_transfer_smith_s95.out` |
| `endpoint_private_goodcut_s95.py` | [STORED] | Checks endpoint-transfer private odd child columns against good-cut/SCC-defect profiles; all private columns are pure and merged non-private SC columns occur only in the strongly connected top bucket through `6->7` | `endpoint_private_goodcut_s95.out` |
| `endpoint_sc_collision_s95.py` | [STORED] | Splits the merged SC endpoint boundary into decomposable/private and strongly connected collision blocks; non-private SC columns occur only in top bucket, have support exactly 3, and are independent through `6->7` | `endpoint_sc_collision_s95.out` |
| `endpoint_collision_geometry_s95.py` | [STORED] | Tests whether support-3 SC collision owner triples are literal parent-metagraph triangles; at `6->7` induced edge counts are `{0:1,1:6,2:5,3:2}`, while the collision incidence hypergraph leaf-peels completely with empty core | `endpoint_collision_geometry_s95.out` |
| `goodcut_class_purity_s95.py` | [STORED] | Verifies `good_cut_count = n - scc_count` and zero mixed unmerged/merged good-cut profiles for n=2..7, proving the good-cut bucket is strong-component defect | `goodcut_class_purity_s95.out` |

### How to bulk-capture results

Run this to capture results for any script:
```bash
for f in 04-computation/*.py; do
    base=$(basename "$f" .py)
    if [ ! -f "05-knowledge/results/${base}.out" ]; then
        echo "Running $f..."
        timeout 300 python3 "$f" > "05-knowledge/results/${base}.out" 2>&1
    fi
done
```

Note: Some scripts take >5 minutes. Use `timeout` appropriately.

---

## Convention for new scripts

1. Save the script in `04-computation/`
2. Run it and pipe output to `05-knowledge/results/SCRIPT.out`
3. Add a row to this INDEX
4. If the result confirms/refutes a hypothesis, update `05-knowledge/hypotheses/INDEX.md`
