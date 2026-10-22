# Session30 nonregular forest extension: independent derivation and audit

Main reports:
- `nonregular_forest_extension.md`: identities, corrected grading, complete first-/second-weight formal predictors, scope and remainder.
- `lollipop_limit_proof.md`: exact regular-completion Schur identity and rigorous nonzero leading lollipop family.

Reproduce from any current working directory (Python 3 with numpy and sympy):

    python /path/to/nonregular/cluster_shapes.py
    python /path/to/nonregular/core_kernel.py
    OPENBLAS_NUM_THREADS=1 python /path/to/nonregular/check_engine.py
    python /path/to/nonregular/prove_closed_kernels.py
    python /path/to/nonregular/lollipop_obstruction.py
    OPENBLAS_NUM_THREADS=1 python /path/to/nonregular/check_weight2.py

Independent cross-implementation checks take explicit paths to the other implementation's frozen outputs:

    OPENBLAS_NUM_THREADS=1 python /path/to/nonregular/compare_engine.py /path/to/decorated/fixtures.json
    OPENBLAS_NUM_THREADS=1 python /path/to/nonregular/compare_symbolic.py /path/to/decorated/symbolic-kernels.json

The latter compares all coefficients as exact rational numbers, not only numeric evaluations. Neither audit script imports the other implementation's code. Original default paths expect sibling staging directories named `session30-nonregular` and `session30-decorated-engine`; explicit arguments avoid that convention in a published packet.

Only .py, .md, .json and .out/.txt files are research deliverables. Ignore any scratch `__pycache__` or `.pkl` files. Reports use PROVED for identities with arguments, exact-comparison verdicts for finite algebraic certificates, and CONDITIONAL/OPEN for the still-missing uniform Hamiltonian-count remainder. No fitted coefficient is used.
