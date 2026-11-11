# An atlas that says exactly what has been classified

**Session271 structural lane.** Complete normal form for all finite tournaments
with an infinite prime alphabet; exact observer laws; no claim of a complete
finite emergent-behavior catalogue. Chan–Pak stays open.

The useful new joint result is a separation of two worlds:

1. **Structured recursive families** have canonical trees, exact symmetry and
   spectral compilers, controllable finite-port interactions, and useful
   Hamiltonian constructors.
2. **The global ensemble** is overwhelmingly prime. Any fixed finite prime
   alphabet produces only exponentially many classes. A typical tournament
   needs asymptotically half of ALL arcs reversed to reach that family, even
   allowing arbitrary relabelling. Sparse seams do not make it globally typical.

This is not a reason to abandon recursive families. It identifies what they
are best at: constructive target synthesis and special invariant regimes,
while probabilistic prime/symmetry corrections drive global large-n forecasts.

## Read and run

- [PROOF.md](PROOF.md): canonical normal form, automorphism/converse recursion,
  finite-alphabet entropy, prime reservoir, edit-distance barrier, irregular
  spectral substitution, and exact rooted spectral/Hamiltonian interfaces.
- [ATLAS.md](ATLAS.md): operational families and their interactions, with
  conserved data and exact boundaries.
- [schema.json](schema.json): implementable canonical carrier and observer
  contract. This is an informal JSON schema, not a claim of JSON-Schema validation.
- [atlas.py](atlas.py): transparent executable canonicalizer, materializer,
  automorphism recursion, characteristic-polynomial compiler, and port response.
- [receipt.json](receipt.json): small deterministic controls only.

From repository root:

    python3 -B 05-knowledge/results/session271_predictive_tournaments_20261111/structures/atlas.py

Python3 and SymPy are required. It uses subset enumeration and factorial
canonicalization, intentionally only a small reference implementation, not a
large-graph isomorphism package. The controls comprise 24 seeded examples of
orders 1–6 (four draws per order; repetitions allowed), 9 substitutions with
small irregular children, the two-root seam hostile, and four tree examples.
Examples are not a complete census; no historical audit was rerun.

## New compact laws

For a child with characteristic polynomial f of degree n, define

    g(z)=(-1)^n f(-z-1)-f(z).

Then arbitrary uniform substitution into quotient adjacency B satisfies

    χ(Q[T_i]) = det(diag(f_i) − B diag(g_i)).

This avoids the regular-child restriction of the old eigenmode splitting.
For C3 the answer is f1 f2 f3−g1 g2 g3. χ alone is exactly compositional for
this observer; it is not a full graph or F2 dynamical invariant.

A seam instead needs the rooted matrix C=Eᵀ(zI−A)⁻¹E. Under A'=A+EKEᵀ,

    χ(A')=χ(A)det(I−KC),   C'=(I−CK)⁻¹C.

The same unrooted TT2 can become either a directed triangle or a transitive
triple after one seam, depending only on the selected root. Thus the sidecar
is mathematically necessary, not decorative metadata.

## Provenance and research posture

Accepted inherited mechanisms: THM-1960, THM-3315, Session268 modular/spectrum/
automata/dynamics, Session270 synthesis and 17-state ordered-spine language.
The precise source paths are in PROOF.md and ATLAS.md. The two October26 source
paths mentioned in Session268 are absent in this checkout; we used available
sources and did not repair unrelated files. The processor startup command
could not run without .machine-id; no machine registration was invented.

The employed meta-pattern is **Controlled forgetting and unlabeled quotients
require a sidecar**: scalar spectra close under uniform modules, but named
seams require a matrix response and true class counting retains stabilizers.
No new global meta-pattern or historical theorem status was promoted.
