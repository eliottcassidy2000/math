# Clock-optimal native eighteen-source circuit

THEOREM.md states the positive native theorem, its exact thirty-seven-clock optimum in the supplied twenty-eight options, and the literal transposition/three-cycle certificate.

From the repository root, with DIR set to this artifact directory:

    python3 -B "$DIR/check.py" --repo .
    python3 -O -B "$DIR/check.py" --repo .
    python3 -B "$DIR/check_boundary.py"
    python3 -O -B "$DIR/check_boundary.py"
    python3 -B "$DIR/check_cycles.py" --repo .
    python3 -O -B "$DIR/check_cycles.py" --repo .

All checkers locate frozen data beside themselves and can run from any working directory. The boundary checker is pure Python and exhausts the complete cost36/37 relation boundary; no optimizer or compiled binary is required for verification.

Optional independent C++ discovery reproduction:

    g++ -O3 -std=c++17 "$DIR/cost_boundary.cpp" -o /tmp/native-cost-boundary
    /tmp/native-cost-boundary 36 < "$DIR/labels_cost.txt"
    /tmp/native-cost-boundary 37 < "$DIR/labels_cost.txt"

Optional image-witness regeneration:

    python3 -B "$DIR/solve_images.py" --repo . --out /tmp/native-cost-optimal-images.json

The saved exact image equations are directly checked, so the solver is not trusted as an acceptance oracle. This is a selected partial packet preserving the full nodewise dyadic bank and original chosen per-slot Xmod8, not every ordinary per-port character, odd-primary bank, common-host amalgam or complete-atlas replacement wiring.
