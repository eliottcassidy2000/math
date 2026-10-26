# Session72: equinumerosity as the dimension shadow of an operation theory

2026-10-26. First verified checkpoint. This packet begins a unifying-theorem
programme, not a claim that the repository's open problems are equivalent.
Recovered results, new elementary corollaries and finite checks are separated.

## The exact correspondence

The same-order unlabelled counting theorem concerns **automorphism-sign-even
graphs**: every automorphism reverses an even number of oriented edges. It
does NOT mean even vertex degrees. At order4 the tournament and sign-even
type counts are4, while the Eulerian type count is3; their labelled counts
are64,20,8 respectively.

The external counting result is Royle–Praeger–Glasby–Freedman–Devillers,
[Tournaments and Even Graphs are Equinumerous](https://arxiv.org/abs/2204.01947v2),
J. Algebraic Combin.57(2023),515–524. The repo already strengthens it with
signed Walsh coordinates, symmetry restrictions, actual type codecs and
obstructions to equivariant set maps. Those are inherited, not rediscovered.

The useful bridge is a basis of **invariant observables**. Equinumerosity is
its dimension count. A linear correspondence need not preserve probability
coordinates, graph toggles, ordered joins, or natural labelled objects.

## Unifying theorem I: harmonic visibility of operations

[The operation-envelope theorem](session72-even/EVEN_GRAPH_OPERATION_ENVELOPE.md)
places affine cube quotients, symmetry-twisted Fourier characters and invariant
difference kernels in one exact framework. It determines which modes survive,
which kernels are invisible, and how averaging behaves under coarsening.
The older ordinary spectrum is explicitly credited to THM-584/587/588 and
THM-4531/4721, with source paths in the native report.

A concrete payoff keeps the SAME weighted single-arc-reversal operation while
changing the quotient:

- Ordinary tournament types: inherited Laplacian gap4, for n≥3.
- Tournament switching classes (quotient by cuts): gap8, for n≥4, with a
  simple C4-indexed quartic mode.
- Score-parity classes (quotient by degree-even graphs): gap2(n−1), for n≥3,
  with an explicit birth–death chain on the allowed parity of the weight.

These statements retain loops and multiplicities. They are not spectral
claims about a simple unweighted quotient. Cut and cycle quotients are not
being treated as a nested chain.

[Independent proof and exact audit](session72-even-audit/REPORT.md).

## A second necessary layer: realizability

[The four-vertex Walsh control](session72-root/WALSH_REALIZABILITY.md) computes
the exact moment cone and shows why full tournament converse becomes a signed
diagonal operator in graph-indexed coordinates. It preserves the correct
transported positive cone, not the naive coordinate orthant. No sign convention
alone fixes that negative diagonal; changing back to point masses does.

This isolates a recurring theme across the repo: preserving a linear identity
is not the same as preserving the admissible objects or operations.

## The wider programme under development

The final session will combine this harmonic theorem with independently
checked results on:

1. State quotients versus terminal-response ranks and positive-law cones under
   reset/action protocols. Missing joint information can have an exact
   incidence-cycle description.
2. Relative alignment across several cyclic frames, including which sparse
   pairwise tests suffice and which prime-power layers require extra edges.
3. Reversible tests, hidden correlations, guarded affine carries and the
   finite-dimensional-carrier boundary exposed by older factorial moments.

The existing reflection `why-the-advances-work-observation-composition-and-non-backtracking-macmini-20261026.md`
already organizes work by object, operation, observation and level. We credit
that framework while seeking concrete theorems and counterexamples within it.
Reflections are provenance, not authority for promoting pre-audit claims.

## Verification

Native curated sources, exact provenance, normal/optimized receipts and hostile
controls accompany this checkpoint. Producer counts and independent coverage
are recorded in the respective reports. No universal prime atlas or new
prize-problem closure is claimed.
