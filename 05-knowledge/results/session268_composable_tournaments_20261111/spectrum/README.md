# Recover composition state from scalar tournament probes

**PROVED:** The full ordered-path-cover vector is exactly the information detectable by Hamiltonian-path counts in arbitrary pure substitution contexts. The nontrivial direction has an explicit physical decoder: H(T), then C3[TT_d,TT_d,T] for d=1,...,n−1. The triangular diagonal is(d!)², and the largest test tournament has3n−2 vertices. Knowing the order saves the final two universal cover coordinates: n−2 scalar observations suffice for n>=3, with largest probe order3n−6 for n>=4. A first cover difference at degree r needs exactly2r−2 added vertices: the probe is optimal among all one-hole pure substitution contexts. See `PROOF.md` §2.

**Constructors:** An arbitrary seed T supplies the exact affine family H'=(H+q2)2^b+H−q2, with a power-of-two membership decoder (§3). Singleton transitive atoms have an injective arithmetic inverse and exact converse multiplicities (§4). The automata lane's hidden port twins combine with them to produce two explicit nonmultiplicative arithmetic surfaces (§7), with logarithmic integer-test decoding. Physical adjacency constructors are in `replay.py`.

**Quantitative limitation:** The singleton atom bank, closed under every SCC product, has only exp(O((log X)^(2/3))) scalar counts≤X. Exact product counts are95,368,1308,4264 at X=10^3,10^4,10^5,10^6. At10^6 there are94 atom values but only4264 values after all products, among500000 positive odd integers. This is a mathematical expressivity bound, not a timeout or incomplete search. H13 is an explicit out-of-family hostile. Chan–Pak remains OPEN.

**Executed:**105 seed tournaments (every labelled tournament through order4 and10 deterministic random seeds at each of orders5,6,7),260 independently counted physical tomography probes,315 physical affine lifts,14 physical atom controls,18 physical seam controls. `results.json` records scope and exact outputs. Normal Python and optimized Python produced identical outputs; the program uses explicit exceptions rather than removable assertions.

Run from repository root:

    python3 05-knowledge/results/session268_composable_tournaments_20261111/spectrum/replay.py

The proofs do not depend on these finite controls. The six-vertex hidden-state seed pair comes from `../modular/PROOF.md`; exact seam laws and five-vertex targets come from `../automata/PROOF.md`.

**Exact target compiler:** `compile_target.py` emits an adjacency witness minimizing vertex count inside the declared atom/seam/lift product language, using O(τ(N)log²N) candidate tests and no factorization oracle. Actual gains: H153 from10 to8 vertices; H177 from11 to8; H495 from14 to9; H27081 from20 to16. These are language-relative improvements, not global minimum claims. Every attained target has a witness using at most4log2(N) vertices. See `compiler_results.json` and proof §8.

    python3 05-knowledge/results/session268_composable_tournaments_20261111/spectrum/compile_target.py 495 --output witness495.json
    python3 05-knowledge/results/session268_composable_tournaments_20261111/spectrum/compile_target.py --measure

Seven literal adjacency witnesses (H59,71,153,177,495,2773,27081), each independently recounted, are preserved in `witnesses.json`. Both replay programs pass normal and optimized Python with byte-identical measurement output.
