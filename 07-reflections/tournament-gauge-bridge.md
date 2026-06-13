# The Tournament-Gauge Bridge: What Physics Sees When It Looks at Computation

**Session:** kind-pasteur-2026-03-21-S12
**Trigger:** Analysis of Napolitano (2026) "Mathematics Is All You Need"

## The Surface Observation

Napolitano claims transformers are gauge theories. The claim is wrong in its literal form
(no Yang-Mills equations, no gauge bosons, no confinement physics). But the mathematical
vocabulary he reaches for — fiber bundles, Lie algebras, Cartan decomposition, phase
transitions — keeps appearing because it IS the right vocabulary for a different, deeper
reason.

## The Real Bridge

The bridge is not "transformers are gauge theories" but rather:

**Both tournament theory and transformer theory are theories of directed structure
on complete graphs, and the natural mathematics for such theories is the Cartan
decomposition of gl(n,R).**

Here is the precise statement:

Any n×n real matrix A decomposes as:
```
A = A_anti + A_sym_traceless + A_scalar
    -------   ---------------   --------
     so(n)          p             R*I
    "active"      "dark"        "null"
```

This is not a theorem about physics or AI. It is a theorem about LINEAR ALGEBRA.
The Cartan decomposition says: every linear map on R^n has a directed part
(antisymmetric, tournament-like), a similarity part (symmetric, undirected),
and a scale.

## Why This Matters for Us

Our entire tournament theory lives in the antisymmetric sector. A tournament IS
an antisymmetric {0,1,-1} matrix (with zero diagonal). The OCF identity
H(T) = I(Omega(T), 2) is a statement about the antisymmetric sector of gl(n,R)
restricted to {0,1,-1} matrices.

The "dark" sector is the symmetric sector. For tournaments, this is trivial
(the symmetric part is fixed = the complete undirected graph K_n). But for
ATTENTION MATRICES, the symmetric sector varies and carries information.

The computational finding: random (untrained) softmax attention puts ~72-73%
of its energy in the symmetric sector and only ~10-19% in the antisymmetric
sector. This means untrained attention is mostly UNDIRECTED similarity.

The question for transformer theory: **Does training make attention more
tournament-like?** If yes, our entire machinery (OCF, path homology,
spectral theory, Paley optimality) applies to trained attention patterns.

## The Deeper Pattern

This is not the first time we've seen directed structure be the "surprising"
part of a decomposition:

1. **Tournaments in so(n):** A tournament T is a basis choice for so(n).
   The Killing form K = -2(n-2)*I is definite — ALL antisymmetric, no
   "dark" sector at all. Tournaments are PURELY directed structure.

2. **Walsh orthogonality:** H(T) has Walsh support only at even Hamming
   weights. This is because H is complement-symmetric: H(T) = H(T^op).
   The antisymmetric Walsh coefficients (odd weight) are ZERO.

3. **OCF:** H(T) = 1 + 2*alpha_1 + 4*alpha_2 + ... is a BINARY expansion
   indexed by independent sets of directed cycles. The undirected conflict
   graph Omega encodes when directed structure CONFLICTS.

4. **Path homology:** beta_2(T) = 0 for ALL tournaments. The second Betti
   number — which measures "holes in directed 2-paths" — vanishes identically.
   Directed structure is too dense to have 2-holes.

In each case, the directed (antisymmetric) part is where the ACTION is,
but its structure is ENCODED in undirected (symmetric) proxies: the conflict
graph Omega, the score sequence, the Betti numbers.

## The Napolitano Inversion

Napolitano's finding (if empirically real) would be an INVERSION of our pattern:

**In tournaments:** Directed part = the object. Symmetric part = trivial.
**In attention:** Symmetric part = the carrier. Directed part = ???

If trained attention patterns have meaningful tournament structure, then
the two frameworks CONNECT: the antisymmetric part of attention gives a
tournament T_attention, and our OCF machinery reveals its structure through
the SYMMETRIC proxy Omega(T_attention).

This would be a beautiful loop: Directed -> Symmetric (via conflict graph)
-> Directed (via attention) -> Symmetric (via Cartan decomposition).

## What Transcends Both Frameworks

The mathematics keeps pointing at the same thing:

**Information about directed structure is best read through symmetric invariants.**

- H(T) is a symmetric function of T (H(T) = H(T^op))
- I(Omega(T), x) is a polynomial of a symmetric object (undirected graph)
- Walsh support at even weights = symmetric sector of Fourier space
- Napolitano's "dark modes" = symmetric sector of Lie algebra

The reason seems to be: directed structure (tournaments, attention patterns,
gauge fields) is too RICH to be read directly. You need to project it through
a symmetrizing lens — and that projection is what carries the invariant
information.

This is analogous to how in physics, the gauge-invariant observables
(Wilson loops, S-matrix elements) are SYMMETRIC combinations of the
underlying antisymmetric gauge field. You don't measure A_mu directly;
you measure tr(P exp(integral A)). The trace is the symmetrizer.

## The Engineering Consequence

If this pattern holds, then the right way to probe LLM computation is NOT
through the directed attention patterns directly, but through SYMMETRIC
INVARIANTS of those patterns:

1. I(Omega(T_attention), x) — the independence polynomial of the conflict
   graph of the attention tournament
2. H(T_attention) — the Hamiltonian path count (a symmetric function)
3. beta_k(T_attention) — Betti numbers (symmetric/functorial invariants)
4. Score variance — how regular the attention tournament is

These are our tools. They are parameter-free. They are mathematically
grounded. And they read directed structure through symmetric lenses —
which is, it turns out, the RIGHT way to do it.

## The Number (n+1)/(n-1)

One concrete fact that fell out of this analysis:

The dark/active dimensional ratio in gl(n,R) is exactly (n+1)/(n-1).

For n=4 (Napolitano): 5/3 = 10/6 modes.
For n=7 (our Paley): 8/6 = 4/3.
For n→∞: approaches 1 (equal dimensions).

This is the fraction of gl(n,R) that is "invisible to the tournament."
As n grows, this fraction shrinks — large matrices are more equally split
between directed and undirected structure. But for small n (like n=4 in
Napolitano's framework), the dark sector has a 67% majority. His "discovery"
that dark modes carry more information is partially an artifact of this
dimensional imbalance.
