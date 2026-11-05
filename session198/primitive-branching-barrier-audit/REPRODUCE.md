# Reproduction of the primitive-source branching barrier

The frozen primary packet is the sibling directory
`../primitive-branching-barrier/`. From that directory run:

    python3 verify.py
    python3 -O verify.py
    sha256sum -c RELEASE_FILES.sha256

From this audit directory run:

    python3 independent_branch_check.py
    python3 -O independent_branch_check.py

The separate audit placement is intentional. The primary verifier hashes all
Python and Markdown inputs beneath its directory; adding reviewer files inside
that frozen directory would change the source receipt. This audit checker has
its own copies of the inherited helpers and frozen primary result and evaluates
the 71 branches in the original response coordinates, independently of the new
positive-matrix implementation.

Both primary modes agree byte for byte, as do both audit modes.
The theorem obstructs extension of supplied primitive representations. It does
not prove scalar nonrepresentability. The five literal rebasing witnesses show
that 41, 81 and 83 all occur in the claimed low-cycle classes.
