# Native eighteen-source checkpoint

Read THEOREM.md for the precise theorem and its partial-packet scope. The selected candidate, literal native inventory, full image witnesses, source pins, and ordinary/optimized verification receipts are retained beside the verifier.

From the repository root, after placing this directory under the chosen Session 46 result directory:

    python3 -B path/to/larger-circuits/check.py --repo .
    python3 -O -B path/to/larger-circuits/check.py --repo .

The verifier is read-only and can run from any working directory. Its --repo argument must point to the repository containing the original graph exports. It locates all new artifacts beside itself.

Optional image regeneration:

    python3 -B path/to/larger-circuits/solve_images.py --repo . --out /tmp/image-witnesses.json

Compare the regenerated file against image_witnesses.json. The verifier checks the saved exact equations directly and does not rely on the solver.

Optional deterministic discovery:

    g++ -O3 -std=c++17 search.cpp -o /tmp/native-circuit-search
    /tmp/native-circuit-search 9 106 < labels.txt

The found support agrees with the frozen candidate up to its stored source ordering. search_cost.cpp and search_milp.py are separate discovery experiments, not validation dependencies or upper-bound proofs. No executable binary is needed in the published artifact.

The full nodewise dyadic bank and chosen ordinary per-slot X modulo 8 are preserved. Every ordinary per-port character, odd-primary guards, a common-host amalgam, and an unchanged full-atlas remainder are not claimed.
