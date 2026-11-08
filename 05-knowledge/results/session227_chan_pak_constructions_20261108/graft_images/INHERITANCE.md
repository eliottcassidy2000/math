# Source comparison and reproducibility scope

This packet develops the existing rooted graft mechanism; it makes no
literature-priority assertion. The following targeted prior sources were
actually read before finalizing the claims:

- `session225_graft_recovery_20261108/PROOF.md`: complete four-state matrix,
  its physical source conventions, literal construction, and B4 gadget.
- `session226_structural_continuation_20261108/graft/PROOF.md`: determinant,
  rank, deletion-ratio identity, rational inverse, and scalar observer.
- `session226_structural_continuation_20261108/root_descent/PROOF.md`,
  Section 1: recovered physical cone inequalities and endpoint-defect source.
- `session226_structural_continuation_20261108/strong_grafts/PROOF.md`:
  exact strongness criterion. We do not claim this criterion again.
- `adaptive_state_induction_20261103/BLOCK_SWITCH.md`, B1–B3: positive
  four-state coordinates, bounded fresh-block language, observer minimum.
- `adaptive_metagraph_reformulation_20261103/ADAPTIVE.md`, Sections 1–4:
  endpoint-defect positivity, legal last-vertex inverse, and its different
  H37 state collision under old-root switching.
- `nonzero_preparation_20261103/algebra/ONE_EXIT_AND_PERSISTENCE.md`,
  Sections 1–5: one-exit sidecar, inherited clone laws, and existing
  H(C3(T_m,1,1))=1+2^m family.
- `01-canon/theorems/THM-001-redei.md`: classical tournament parity.

Repository startup, research protocol, active guardrails, relevant current
frontier, and CORE-PAPERS were also read. Targeted materialized-file searches
for conic/affine hulls and rooted-state cones found the recovered inequality
sources above. Broader `git grep HEAD` searches were cancelled after the
partial checkout attempted unavailable remote-blob reads. Consequently this
packet does not claim an exhaustive repository search or infer priority
from a search miss.

The new refinements relative to the consulted sources are exactness of the
closed conic and affine-lattice hulls, transport of those exact hulls under
invertible grafts, an all-order both-sign determinant family, exact scalar
lattice projection, and a graph-image module certificate with a distinct
physical numerical-inverse hostile. Individual inequalities and H=1+2^m
are explicitly inherited.

## Local computational dependency

`inherited_verify.py` is copied unchanged from the Session226 graft packet,
which itself copied the reconstructed Session225 verifier unchanged.
`verify.py` imports only this sibling file and Python's standard library.
The checker does not read repository paths, network services, scratch probes,
or files outside this packet. Thus copying the five files PROOF.md,
verify.py, inherited_verify.py, CHECK.json, CHECK.optimized.json to a clean
directory suffices to reproduce the producer controls in normal and optimized
Python. INHERITANCE.md and SHA256SUMS supply provenance and integrity only.
