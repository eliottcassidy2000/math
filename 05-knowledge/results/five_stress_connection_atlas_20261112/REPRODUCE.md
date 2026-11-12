# Reproduction and data boundaries

Run from the repository root with Python3, SymPy, NumPy and Numba installed.
The supplied PDF is used only for its source hash in `paper_transfers.py`;
change that declared local path if reproducing on another machine. Inherited
modules and the complete small tournament bank remain repository dependencies.

```powershell
python3 -X utf8 05-knowledge/results/five_stress_connection_atlas_20261112/verify_stress.py
python3 -X utf8 05-knowledge/results/five_stress_connection_atlas_20261112/verify_tournaments.py 14
python3 -X utf8 05-knowledge/results/five_stress_connection_atlas_20261112/extend_round.py 18
python3 -X utf8 05-knowledge/results/five_stress_connection_atlas_20261112/c3_pair_family.py
python3 -X utf8 05-knowledge/results/five_stress_connection_atlas_20261112/carousel_clocks.py
python3 -X utf8 05-knowledge/results/five_stress_connection_atlas_20261112/seed_probes.py
python3 -X utf8 05-knowledge/results/five_stress_connection_atlas_20261112/close_heads.py
python3 -X utf8 05-knowledge/results/five_stress_connection_atlas_20261112/classify_seed_census.py
python3 -X utf8 05-knowledge/results/five_stress_connection_atlas_20261112/paper_transfers.py
python3 -X utf8 05-knowledge/results/five_stress_connection_atlas_20261112/cyclic_port_probe.py
python3 -X utf8 05-knowledge/results/five_stress_connection_atlas_20261112/joint_port_extension.py
python3 -X utf8 05-knowledge/results/five_stress_connection_atlas_20261112/contraction_dag.py
python3 -X utf8 05-knowledge/results/five_stress_connection_atlas_20261112/descendants_nine.py
python3 -X utf8 05-knowledge/results/five_stress_connection_atlas_20261112/birth_cones.py
python3 -X utf8 05-knowledge/results/five_stress_connection_atlas_20261112/aspect_sequences.py
python3 -X utf8 05-knowledge/results/five_stress_connection_atlas_20261112/arithmetic_bridge_controls.py
python3 -X utf8 05-knowledge/results/five_stress_connection_atlas_20261112/replay_audit.py
```

The order18 extension took about23.5 minutes in this session. The replay
intentionally checks only its complete order1..14 prefix; it does not silently
replace the larger saved receipt. Exact checker scopes are in [AUDIT](AUDIT.md).

The historical catalog requires all objects reachable from `history_refs.txt`
to be locally available. Missing objects are errors, not silently omitted.

```powershell
python3 -X utf8 05-knowledge/results/five_stress_connection_atlas_20261112/history_scan.py
python3 -X utf8 05-knowledge/results/five_stress_connection_atlas_20261112/history_finalize.py
python3 -X utf8 05-knowledge/results/five_stress_connection_atlas_20261112/history_prompt_harvest.py
python3 -X utf8 05-knowledge/results/five_stress_connection_atlas_20261112/build_prompt_index.py
python3 -X utf8 05-knowledge/results/five_stress_connection_atlas_20261112/history_controls.py
```

This streams56.6 GB of historical text but stores bounded excerpts in
deterministic gzip files. The catalog is evidence for retrieval, not an audit
of every mathematical statement in those objects. The SHA manifest and frozen
refs distinguish this scan from later availability-dependent searches.
`history_manifest.py` rebuilds the object/path/size manifest from those refs.
The full-file prompt pass instead reads the largest revision per path; its
different1.07 GB universe and selection rule are recorded separately.

The optional static figure uses Matplotlib via `draw_connection_map.py`.
In this session its packages were placed in a dedicated temporary directory,
not installed into or substituted for the numerical-check environment.

For the repository documentation check on this Windows path:

```powershell
python3 -X utf8 '\\?\C:\Users\Eliott\Documents\GitHub\math-five-stress-20261112\agents\check_docs.py'
```

An ordinary short-path invocation can fail on existing long theorem filenames
before performing the check. On a different checkout substitute its absolute
path, or use a platform supporting these filenames without the prefix.
