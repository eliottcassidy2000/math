# Balanced choice-network quotients

2026-10-21. **PROVED auxiliary-network results; independent review PASS.**
The immutable [producer proof](PROOF.md) retains its historical pre-review
header; the [final independent audit](review/AUDIT.md) supersedes that
workflow status. The graph here is a balanced choice graph, NOT the actual
tournament metagraph.

## Exact carrier and operation distinctions

Contracting complete-graph choice fibres of weights s,...,s,1 gives edge
conductances w_x w_y. If m is the number of finite parent choices and
Q=sm+1, its critical group is C_s x (C_(sQ))^(m-1). For s=2^d with d>=1,
the two-primary part is (C_(2^d))^m. The boundary d=0 is explicitly excluded
from that last simplification.

Shorting preserves energy on constant-fibre functions; Kron elimination
preserves terminal response using harmonic extension; random-walk lumping
preserves its stated transition law only with holding probabilities. These
are different operations and do not generally preserve the integral firing
lattice. An arbitrary Laplacian forgets loops, but the known pointed
complete-source model can reconstruct collapsed loop counts from its
weights and distinguished infinity state. The final proof retains this
important qualification.

## Inheritance and remaining boundary

The incoming [marked-fibre precursor](../metagraph_fibre_network_20261020k/README.md)
already proves the one-step choice conductance, resistance, held Markov law,
and conditional entropy. The exact all-depth critical group and the explicit
shorting/Kron/lattice distinctions are the additions here. The
[balanced deletion packet](../balanced_deletion_groups_ports_20261021/README.md)
and [dyadic conductor packet](../dyadic_conductor_blue_fibres_20261021/README.md)
remain separate carriers. The smallest V4 coincidence is not a general
identification of these networks, groups, or literal path moves.

## Reproduction

From this directory, using standard-library Python:

    python3 -B check.py
    python3 -B -O check.py
    python3 -B replay_review.py
    python3 -B replay_review.py --optimized

The launcher accepts --repo /absolute/path/to/math when outside the repository.
It runs the unchanged independent checker with the frozen companion sources
under review/workspace/. All source pins are checked before computation.
Producer controls have 1,328 checks; independent controls have 29,865 checks.
Both normal/optimized pairs are byte-identical. Retained negative pin controls
fail closed for changed proof bytes. Finite checks support the all-parameter
proof audit; they do not justify a tournament metagraph classification.
