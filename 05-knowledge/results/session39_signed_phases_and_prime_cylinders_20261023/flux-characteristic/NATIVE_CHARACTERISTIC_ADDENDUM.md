# Native metagraph characteristic covector and local spin parity

23 October 2026. Separate addendum to the frozen flux packet. Status: exact algebraic consequence of inherited fibre parity, with a native order-eight check. No identification of graph nodes with physical electrons.

## Input and scope

Use the literal converse-merged fixed-Hamiltonian-path complement-line multigraph, with loops counted twice in degrees. This is the construction for tournament order n≥3, where complement pairs presentations. Let H(C) count directed Hamiltonian paths and let h(C)=H(C)/|Aut(C)| be the ordinary presentation fibre size. The inherited theorem is that h(C) is odd. Its ingredients, Rédei's oddness theorem and odd order of tournament automorphism groups, are classical; see the repository's `01-canon/theorems/THM-643-blue-black-line-parity-structure.md`, T1, depending on LEM-003 for the free fibre identity. We do not claim a new parity theorem.

Each literal line has two presentation endpoints. Consequently a self-converse merged node has degree h, while a non-self-converse node merges two equal fibres and has degree 2h. Let ε_v be one for a self-converse node and zero otherwise. Then

    degree(v) ≡ ε_v (mod 2).

Removing a loop from the electrical Laplacian removes two from degree, so it does not change this parity. Choose a root in each connected component and ground the integral electrical Laplacian, obtaining A. All coordinates below refer to that declared grounded vertex basis. The retained self-converse indicator c satisfies

    c_i ≡ A_ii (mod 2).

## Characteristic-covector theorem

For every integer vector z on the retained vertices,

    z^T A z ≡ c^T z (mod 2).

Proof: off-diagonal terms appear twice; z_i²≡z_i modulo two. The diagonal congruence above finishes the proof. Thus c is a characteristic covector for the integral lattice A. Ground choice changes coordinates and retained indicators; it must not be silently suppressed when displaying this formula.

This gives a precise mathematical dictionary after making the explicit formal effective-theory choice K=A. A local lattice excitation is ℓ=A z, whose representative exchange factor is

    exp(πi ℓ^T A⁻¹ℓ)
       = exp(πi z^T A z)
       = (−1)^(sum over retained self-converse nodes i of z_i).

A unit local lattice generator at a retained self-converse vertex is therefore fermionic in this chosen spin CS description; a unit generator at a retained non-self-converse vertex is bosonic. These are formal local generators, not fractional anyons labeled by vertices. In particular ℓ=A z is zero in the discriminant group: the odd signs are precisely why naive bosonic twists do not descend after quotienting away local fermions. No microscopic electron, charge assignment, Hamiltonian, or physical phase has been inferred. If A has an odd diagonal entry, the K=A construction is a spin theory rather than an ordinary bosonic CS theory.

The same c supplies a normalized algebraic quadratic refinement on the graph discriminant group:

    q_c([x]) = (x^T A⁻¹x − x^T A⁻¹c)/2 mod 1.

Indeed replacing x by x+A z changes this expression by

    z^T x + (z^T A z − z^T c)/2,

an integer by the characteristic congruence. Its polarization is the graph linking pairing x^T A⁻¹y mod 1. The qualifier “normalized algebraic refinement” matters: it is not the uncorrected representative exchange phase, and spin/framing choices still belong to the full theory.

## Contrast with the tournament-skew construction

The separate construction K_S=[[0,S],[S^T,0]] has even diagonal and odd-order discriminant group for every even-order tournament. The native Laplacian construction instead has a characteristic covector recording retained self-converse vertices and can carry the authenticated dyadic discriminant sector. They arise from different matrices and must not be equated.

## Reproduction

`check_characteristic.py` reads the authenticated Session38 `nodes.tsv.gz` and `edges.tsv.gz`, checking the full native degree identity and diagonal parity, then checks local parity on 128 deterministic integer potentials with root 0. It also checks characteristic-refinement descent in a small odd integral lattice on 1,225 cases. SHA256 hashes of both graph inputs are saved in the receipt. These finite controls accompany the all-size proof; they are not substituted for it.

Run `python check_characteristic.py /path/to/session38-wiring` and repeat with `python -O`. Compare both outputs with the saved receipt. With no argument, the development default is `../session38-wiring` relative to the script directory. The publication source of those inputs is `05-knowledge/results/session38_incidence_ports_electrical_20261023/wiring/`.
