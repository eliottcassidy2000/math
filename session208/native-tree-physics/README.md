# Native metagraph tree physics

Independently audited elementary theorems and one exact counterexample; see PROOF.md and REVIEW.md. This packet does not establish the new preprint's global logarithmic theorem.

- Bridge and leaf channels are exact in every spanning tree, despite unbounded adjacent Hamilton-count ratios.
- Fundamental-cut capacities dominate scalar total variation, but not quadratic electrical energy.
- The actual labelled order-four metagraph has a canonical tree and five-support integer voltage with energies31 versus23.

## Reproduction

From any directory, run `python3 /path/to/this/packet/independent_check.py`. This independently reconstructs the graph and checks the saved certificate using only Python's standard library and exact integers. It needs no other repository files.

Optional witness discovery: run `OPENBLAS_NUM_THREADS=1 python3 check.py` inside this directory. This requires NumPy; its eigenvector is a search device only. Floating eigensolver versions can select different integer witnesses. The saved RESULT.json plus independent_check.py is the portable exact certificate, not a claim that floating search has bitwise identical output on every platform.

Mathematical inheritance is identified in PROOF.md; the native graph, tree, and proof definitions are self-contained. No theorem namespace IDs are reserved by this packet.
