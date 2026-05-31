# Script Results Index

Maps every script in `04-computation/` to its stored output.
Results are stored in `05-knowledge/results/` with matching filenames (`.out` extension).

**Convention:** When you run a script, save its output:
```bash
python3 04-computation/SCRIPT.py 2>&1 | tee 05-knowledge/results/SCRIPT.out
```

---

## Lean build and audit outputs

| Target | Status | Key finding | Result file |
|--------|--------|-------------|-------------|
| `lake build TournamentH7.RootPackets` | [STORED] | New root-packet module builds successfully: open root walks telescope to endpoint boundary, closed packets have zero total root, and directed cycles convert to zero-root packets. | `lean_root_packets_s6.out` |
| `timeout 300 lake build TournamentH7.Verify` | [STORED] | Full audit attempt timed out while still compiling dependencies; it did not reach the `#print axioms` layer. Use the narrow `RootPackets` build as the completed verification for this session. | `lean_verify_root_packets_s6.out` |
| `lake build TournamentH7.ProductSum`; `lake build TournamentH7.Verify`; `lake build TournamentH7` | [STORED] | THM-361 list-level product-sum defect normal form now builds in Lean: removing ones preserves product and records additive slack, padding a core by defect many ones repairs product-sum equality, and the positive two-entry resonance is only `(2,2)`. | `lean_product_sum_s366.out` |
| `lake build TournamentH7.NaturalOperationDigraphs` | [STORED] | New axiom-free operation-shadow module builds successfully: additive shadow is strict order, multiplication shadow is divisibility/proper divisibility, and the two-factor product-sum layer plus seed-defect padding are formalized. | `lean_natural_operation_digraphs_s366.out` |
| `lake build TournamentH7.Verify` | [STORED] | Full Lean audit builds after adding `NaturalOperationDigraphs`; new audit entries depend only on Lean foundations (`propext`, `Classical.choice`, `Quot.sound`) and no project axioms. | `lean_verify_natural_operation_s366.out` |
| `lake build TournamentH7` | [STORED] | Full root module builds successfully after wiring in `NaturalOperationDigraphs`; existing lint warnings in `BasePathSink`/`TransitiveH` remain unchanged. | `lean_tournamenth7_natural_operation_s366.out` |

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
| `residue_phase_incidence_contrast_s5.py` | [STORED] | Standard contrast table splitting exact residue kills, dangerous near-kills, broad phase/fiber contrasts, and the n=6 H=37 landscape trap | `residue_phase_incidence_contrast_s5.out` |
| `tournament_tda.py` | [STORED] | Demo now emits SCC residue features and deletion-residue rank features in the `omega_*` block | `tournament_tda_residue_features_s355.out` |
| `endpoint_transfer_witnesses_s95.py` | [STORED] | Probes endpoint-transfer rank mechanisms; unmerged tournament rows all have private odd child witnesses through `6->7`, merged rows keep full matching/rank without all private witnesses, and even graphs show GF(2) cancellation | `endpoint_transfer_witnesses_s95.out` |
| `endpoint_transfer_smith_s95.py` | [STORED] | Computes Smith normal forms for small endpoint-transfer matrices; tournament/merged factors are all odd through `5->6`, while even graphs show 2-primary factors up to `[2,2,2,2,2,4,4,8]` at `6->7` | `endpoint_transfer_smith_s95.out` |
| `endpoint_private_goodcut_s95.py` | [STORED] | Checks endpoint-transfer private odd child columns against good-cut/SCC-defect profiles; all private columns are pure and merged non-private SC columns occur only in the strongly connected top bucket through `6->7` | `endpoint_private_goodcut_s95.out` |
| `endpoint_sc_collision_s95.py` | [STORED] | Splits the merged SC endpoint boundary into decomposable/private and strongly connected collision blocks; non-private SC columns occur only in top bucket, have support exactly 3, and are independent through `6->7` | `endpoint_sc_collision_s95.out` |
| `endpoint_collision_geometry_s95.py` | [STORED] | Tests whether support-3 SC collision owner triples are literal parent-metagraph triangles; at `6->7` induced edge counts are `{0:1,1:6,2:5,3:2}`, while the collision incidence hypergraph leaf-peels completely with empty core | `endpoint_collision_geometry_s95.out` |
| `goodcut_class_purity_s95.py` | [STORED] | Verifies `good_cut_count = n - scc_count` and zero mixed unmerged/merged good-cut profiles for n=2..7, proving the good-cut bucket is strong-component defect | `goodcut_class_purity_s95.out` |
| `lonely_runner_residue_probe_s356.py` | [STORED] | Exact rational forbidden-interval probe for Lonely Runner as a quotient-gap problem; initial segments are boundary-only, while sampled lacunary/arithmetic/CRT/random sets have positive gaps | `lonely_runner_residue_probe_s356.out` |
| `caccetta_haggkvist_residue_probe_s357.py` | [STORED] | Cyclic Cayley zero-sum probe for Caccetta-Haggkvist; enumerates 1,612,949 small connection sets with no violations and isolates tight return-residue/Kemperman-equality profiles | `caccetta_haggkvist_residue_probe_s357.out` |
| `lonely_runner_tight_scan_s357.py` | [STORED] | Exhaustive primitive-box scan for the LRC tight stratum; no open-cover candidates, rare boundary-only sets matching known tight examples, and first boundary witness usually at `1/(k+1)` | `lonely_runner_tight_scan_s357.out` |
| `lonely_runner_endpoint_protection_s359.py` | [STORED] | Exact endpoint-protection incidence graph and core peel for LRC; known tight examples have unprotected true endpoints at nonzero unit residues mod `k+1`, and primitive boxes through `(6,11)` have zero nonempty cores even after the `k+1` divisibility filter | `lonely_runner_endpoint_protection_s359.out` |
| `lonely_runner_endpoint_protection_s360.py` | [STORED] | Builds the exact LRC endpoint-protection graph; verifies the integer protection criterion from THM-357 and finds no all-protected open-cover candidates in the inherited bounded primitive boxes | `lonely_runner_endpoint_protection_s360.out` |
| `lonely_runner_bohr_descent_s362.py` | [STORED] | Verifies THM-358 initial-segment unit skeleton through n=36 and implements THM-359 endpoint/interval core peeling; all inherited full-measure primitive-box cases have empty terminal protection core | `lonely_runner_bohr_descent_s362.out` |
| `lonely_runner_fourteen_runner_s363.py` | [STORED] | Fourteen-runner `k=13` CRT-gate probe; every scanned primitive set passing the mandatory divisible-by-14 filter is positive-gap and endpoint-core-empty, suggesting a composite descent attack | `lonely_runner_fourteen_runner_s363.out` |
| `lonely_runner_k13_microstaircase_s363.py` | [STORED] | Probes the fourteen-runner frontier: exact `r/14` tight-class analogue fails, but a micro-staircase cell resolves the obstruction on prime grids; small-prime verifier ledger shows the bottleneck is computing `I(13,p,1)` | `lonely_runner_k13_microstaircase_s363.out` |
| `lonely_runner_feedback_loop_s364.py` | [STORED] | Cycles between 14-runner, 15-runner, and counterexample routes; finds the full scalar-ramp family blocks every `n=14` and `n=15` micro-staircase cell, while scalar-excised local search finds only near-blockers and gated speed searches find no open-cover candidates | `lonely_runner_feedback_loop_s364.out` |
| `lonely_runner_composite_gate_feedback_s364.py` | [STORED] | Alternates fourteen-runner, fifteen-runner, and disproof-construction pressure; composite gate insertions for n=14 and n=15 leak to positive gaps with empty endpoint cores, while random gated disproof searches produce no open covers | `lonely_runner_composite_gate_feedback_s364.out` |
| `lonely_runner_frontier_feedback_s364.py` | [STORED] | Extends the 14/15/disproof feedback loop with exact gated boxes, full late-peel ledger for the hardest n=14 gate leak, gate-overload disproof families, and n=14/n=15 micro-staircase repair tests | `lonely_runner_frontier_feedback_s364.out` |
| `lonely_runner_k13_scalar_gauge_s367.py` | [STORED] | Refines the fourteen-runner micro-staircase route by quotienting scalar ramps as a gauge symmetry; exact small-support scans and the full normalized 2-torsion cube show the unique binary extremal is the coordinate-6 half-turn, still missing 56 cells. | `lonely_runner_k13_scalar_gauge_s367.out` |
| `lrc_tournament_56_bridge_s369.py` | [STORED] | Compares the S367 `56` LRC missed cells with six-vertex tournament structure: `56=12+44` is the self-converse/chiral split, `56=14+42` is the LRC outer-pair/interior split, and the old 5-to-6 perspective gap is exactly 8 alpha stencils. | `lrc_tournament_56_bridge_s369.out` |
| `chirality_perspective_atlas_s370.py` | [STORED] | Applies the 12/42/44/56 lens across tournaments, LRC, Paley/Fano `T7`, and base-42 residues; identifies `12` as inherited symmetric core, `44` as chiral/support residue, `42` as doubled boundary/interior, and `8` as projection failure. | `chirality_perspective_atlas_s370.out` |
| `lonely_runner_scalar_excision_s371.py` | [STORED] | Formalizes the n=14 scalar-ramp cell audit behind THM-364, reconstructs all 812 representative cells, scans scalar-ramp one/two-defect neighborhoods, and lists the 56 uniquely protected missed cells of the best non-scalar near-blocker | `lonely_runner_scalar_excision_s371.out` |
| `natural_operation_digraphs_s365.py` | [STORED] | Shows the additive natural-number operation shadow is the transitive tournament while the multiplication shadow is the divisibility DAG; enumerates product-sum critical pairs, including the two-factor divisor layer, the first multi-factor improvement at arity 5, and target-mode collision spectra | `natural_operation_digraphs_s365.out` |
| `natural_mode_graph_s365.py` | [STORED] | Compares additive and multiplicative natural-number mode graphs; identifies `4` as hidden diagonal binary resonance, `6` as unique distinct ternary product-sum resonance, and verifies the product-sum core-defect ledger through arity 10 | `natural_mode_graph_s365.out` |
| `natural_operation_graphs_s365.py` | [STORED] | Shows additive operation shadows collapse to the complete transitive order while multiplicative shadows remain the divisor DAG; classifies product-sum collisions by finite prefix gates, Egyptian-fraction shadows, target spectra, and nonunit core defect | `natural_operation_graphs_s365.out` |
| `natural_operation_modes_s366.py` | [STORED] | Cycles formalization, extension, and investigation for operation digraphs; product shadow has density `0.01215` on `[1000]`, and strict multi-factor product-sum wins occur in `96/119` arities through `k<=120` | `natural_operation_modes_s366.out` |
| `fejer_kernel_wild_session.py` | [STORED] | Verifies Fejer/autocorrelation identity for interval spectra; finds `Var(J_3)=E(S)/(p-1)-(p^2-2p+5)/16` over full circulant orientation cubes through `p=23`; compares IPR/J3/H sign reversal at `p=7,11,13` | `fejer_kernel_wild_session_2026_05_30.out` |
| `fejer_root_sign_phase_probe_s359.py` | [STORED] | Verifies circulant root-sign eigenvalue formula and interval/Fejer identity; compares orbit-invariant concentration, Fejer alignment, Paley flatness, and `H` for small prime circulants | `fejer_root_sign_phase_probe_s359.out` |

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
