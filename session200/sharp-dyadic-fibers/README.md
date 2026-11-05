# Paper-inspired physical source results

Read PROOF.md for the independently proved dyadic fiber and height-separation theorems. PAPER_CONNECTIONS.md identifies and discusses the supplied manuscript without treating its headline claim as an imported theorem.

Run:

    python3 verify.py > NORMAL.json
    python3 -O verify.py > OPTIMIZED.json
    cmp NORMAL.json OPTIMIZED.json

The verifier uses only Python's standard library. It checks indexed word fibers, literal tournament paths and cycles, the sharp boundary construction, and the bounded C1 hostile at7. It does not infer an all-order theorem from finite tests. Audit results are maintained separately.

The broader all-odd {7,21} spectrum remains open. The packet does not claim a new exact scalar representation, nonisomorphic multiplicity, or a number-field Diophantine theorem.
