# Verification and reproduction

**Scope:** same-agent review, separate exact algorithms, no formal proof assistant. A check count measures successful assertions, not independent theorems.

## Gates and finite universes

| Gate | Exact scope | Outcome |
|---|---|---|
| Normality | 50 supports, 29 inequalities, 988 simplices, all fundamental residues | 12,179 checks; all-size conclusion via P2's proof |
| Palettes | All demands in {0,1,2,3}^7, 64 literal substitutions, H ray through64, clique bijections and native81/27 | 16,593 checks |
| LRC | 10,800 direct comparisons; 2,184 declared families; point/joint/large-dilation controls | 23,021 checks |
| All-pair LRC upgrade | All78 eleven-subsets of1–13; incoming two-tail conductor | All6,397 finite-head pairs pass independent literal intersection; every remaining pair covered by proved tails |
| Capacity | Optimized recurrence through12; full labelled enumeration independently through6 | Every value agrees; all construction witnesses through12 pass |
| Damaged certificates | Omitted lattice class, missing simplex, damaged palette, wrong support list | All four rejected |
| Corrected GW witness | Thirteen speeds × three source residues | All39 circular distances at least2 modulo27 |

Normal and Python -O receipts match byte for byte for every gate. The code uses explicit exceptions rather than optimization-sensitive assert statements for verification.

The normality checker independently reconstructs supports and inequalities and uses Fraction Gaussian determinants. The producer uses Bareiss determinants; floating-point Delaunay only proposes a covering family, whose coverage is proved from exact checked incidences. No floating-point output is promoted directly.

The capacity producer prunes partial budgets and uses a valid set-system symmetry. The independent checker uses full labelled Cartesian products and literal transitive subsets through six, then validates every stored witness through twelve. The all-budget recurrence and exponent are proved in RAMSEY.md rather than extrapolated.

The LRC reference path builds every literal runner comb for small test cases, while the compiled path uses common-dilation first arrival. Both retain rational endpoints and isolated points. The all-d conclusion uses L2's proved tail bound and each entry's exact finite head.

The native 81/27 verifier is source-pinned. The new query is evaluated on the actual extended tournaments. It separates them even though their entire current named transitive-support lists agree.

## Commands

Run from repository root with Python 3.10+:

~~~text
python -X utf8 -B 05-knowledge/results/clique_palette_joint_ports_20261030/check_normality.py
python -X utf8 -B 05-knowledge/results/clique_palette_joint_ports_20261030/check_palettes.py
python -X utf8 -B 05-knowledge/results/clique_palette_joint_ports_20261030/lrc_banks.py
python -X utf8 -B 05-knowledge/results/clique_palette_joint_ports_20261030/all_pair_cores.py
python -X utf8 -B 05-knowledge/results/clique_palette_joint_ports_20261030/capacity.py --bound 12
python -X utf8 -B 05-knowledge/results/clique_palette_joint_ports_20261030/check_capacity.py
python -X utf8 -B 05-knowledge/results/clique_palette_joint_ports_20261030/check_corruptions.py
python -X utf8 -B 05-knowledge/results/clique_palette_joint_ports_20261030/seal.py
~~~

Each checker accepts --output PATH for a JSON receipt; use -O after -B for optimized-mode comparison. lrc_banks.py additionally accepts --families PATH for the complete finite family list. The seal checks saved receipts, exact artifact bytes, and inherited Git blobs.

explore.py, discover_hostile.py, discover_triangulation.py and build_normality.py preserve discovery provenance. Only optional discovery needs NumPy/SciPy. Reproduction of the mathematical checks uses the Python standard library and the pinned native hostile helper. triangulation.json is a producer artifact; the independent checker consumes normality_certificate.json.

No external messages, submissions, or theorem reservations were made. The intentionally changed maintained files are the results index, the two THM-658 wording repairs, and the corresponding MISTAKES lineage.
