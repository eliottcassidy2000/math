# Portable intersection reproduction

From the repository root, supply the actual path to this packet:

    python -B <packet>/check_intersection.py --repo-root .
    python -B -O <packet>/check_intersection.py --repo-root .

The source pin labels remain their original provenance paths. Labels beginning math/ resolve only under the explicitly supplied repository root, with math/ stripped. Labels beginning session22-incoming-center-review/ resolve inside this actual packet. Other staged labels, including session22-q2structural/ and session22-q2-review/, resolve only beneath this packet's actual parent directory. Preserve those sibling packets when publishing or copying. There is no fallback to the original workspace, and absent or mismatched files fail closed.

The eleven incoming source pins were first checked and replayed in an isolated copy. Before publication each identical LF hash was reverified against the actual repository file; their provenance labels were then relocated to canonical math/ paths without changing any hash value. The two q=2 proof/review labels remain sibling paths. The published checker requires --repo-root and reads the actual repository inputs; normal and optimized replays agree at 17 checks. No duplicate input snapshots are needed or published.

The abstract orbit-color control is self-contained:

    python -B <packet>/orbit_color_control.py
    python -B -O <packet>/orbit_color_control.py

From the repository root, replay the existing unchanged incoming consumer directly:

    python -B 05-knowledge/results/deficit_six_centers_20261021b/combined_check.py --verify
    python -B -O 05-knowledge/results/deficit_six_centers_20261021b/combined_check.py --verify
