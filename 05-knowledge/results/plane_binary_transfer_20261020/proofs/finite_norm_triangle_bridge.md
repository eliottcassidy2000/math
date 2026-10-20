# Finite norm triangles: a checked obstruction invisible to radial PSD

**Status:** the Fourier identities and relaxation construction are PROVED below. Non-two-colorability at `q=17` is FINITE-EXACT with a complete refutation tree checked by a separate verifier. The `q=23` spectral bound is certified by exact rational interval arithmetic. No real-plane coloring bound follows, and no historical-priority claim is made.

## 1. Carrier, predicate and normalization

Let `q` be an odd prime power, let `E=F_(q²)`, let `N(z)=z^(q+1)`, and let `U={u:N(u)=1}`, of size `q+1`. All nonzero differences have nonzero norm. The rotations are multiplication by `U`; translations and rotations preserve `N(x-y)`.

Fix three distinct points `T=(0,a,b)` in `E`, and let `s:E->{±1}`. Average monochromaticity over all `q²(q+1)` ordered placements `(p,p+ua,p+ub)`. For a nonzero difference `v`, define

    C_v = E_(p,u) s(p)s(p+uv).

The sign identity from blueprint §4 gives exactly

    D_T(s) = [1+C_a+C_b+C_(b-a)]/4.                           (F1)

Let `ŝ(ξ)=q^(-2) Σ_z s(z) overline(ψ(Tr(ξz)))`, for a fixed nontrivial additive character `ψ` of `F_q`. Parseval is `Σ_ξ|ŝ(ξ)|²=1`. Put

    B_q(c)=-(q+1)^(-1) K_q(c),
    K_q(c)=Σ_(t∈F_q*) ψ(t+c/t).

The normalized Fourier multiplier for the norm circle `vU` is `1` at `ξ=0`, and `B_q(N(v)N(ξ))` otherwise. This is the precise finite-field counterpart of the real circle multiplier `J_0(r|ξ|)`. It is the same orthogonal-circle averaging operation; there is no claim that a prime-to-infinity limit identifies the two models.

Thus, writing `r_1=N(a), r_2=N(b), r_3=N(b-a)`,

    4D_T(s)=1+3|ŝ(0)|²+
       Σ_(ξ≠0)|ŝ(ξ)|² [B_q(r_1Nξ)+B_q(r_2Nξ)+B_q(r_3Nξ)].   (F2)

The torus-to-Kloosterman calculation is proved in the repository at
`01-canon/theorems/THM-4628-residue-colourings-of-the-unit-distance-graph-of-q-zeta55-are-kloosterman-forced-at-every-conjugation-fixed-prime.md`, statement (b). It can also be obtained by counting roots of `X²-tX+c`: norm-`c`, trace-`t` elements are counted by `1-χ(t²-4c)`, whereas `x+c/x=t` is counted by `1+χ(t²-4c)`.

**Preserved coordinate:** the full vector `(r_1,r_2,r_3)` and one shared dual norm `Nξ`. Replacing its three-term symbol by three independently optimized minima is a lower bound, not an equality. This is the same shared-frequency issue that the blueprint's Bessel pair bound retains.

## 2. A uniform finite-field theorem and a certified small threshold

If `M_q=max_(c≠0) K_q(c)`, then from F2,

    D_T(s) >= (1 - 3 M_q/(q+1))/4.                            (F3)

More precisely, retaining the mean `m=ŝ(0)` gives

    D_T(s) >= [1 - 3M_q/(q+1) + (3+3M_q/(q+1))m²]/4.         (F4)

The standard Weil bound `M_q<=2sqrt(q)` gives positive defect whenever `q+1>6sqrt(q)`, in particular for every odd prime power `q>=37`. This is a CITED-input corollary, not a new theorem about the real plane. Finite-field monochromatic configuration arguments using the same Kloosterman method already appear in Shkredov's work and in Shkredov--Solymosi, *Tilted corners in integer grids* (2021), §4. The contribution here is the exact carrier/normalization and the checked small-`q` certificate and hostile below.

For `q=23`, exact rational interval arithmetic verifies

    M_23 = K_23(11) < 796068719/100000000 = 7.96068719 < 8.

Every distinct three-point template in the anisotropic norm plane `F_(23²)` therefore has

    D_T(s) > (1 - 796068719/800000000)/4
             = 3931281/3200000000.                           (F5)

There are `529*24=12696` placements; the lower bound times this number exceeds 15. Thus every coloring has at least 16 monochromatic ordered placements of every template. For equilateral templates the placement multiplicity is three, so there are at least six unordered monochromatic equilateral triples.

**Exact certificate.** `scripts/certify_kloosterman_bounds.py` uses only Python `Fraction` arithmetic. It encloses `π` by Machin's identity with 40 alternating-series terms in each arctangent, then encloses each cosine by a degree-48 Taylor polynomial plus its Lagrange remainder. It checks every `c=1,...,22`, certifies the maximizer, and verifies F5. No floating-point value is load-bearing.

The prime `23` is the first prime `q≡2 (mod 3)` among those before it for which this single-circle spectral positivity test succeeds: `q=5,11,17` fail that sufficient test. The next section shows the actual hypergraph threshold is already `17`, so spectral failure is not coloring existence.

## 3. Exact equilateral threshold at q=17

Assume `q≡2 (mod 3)`. Then `U` contains a root `ζ` of `ζ²-ζ+1=0`. The equilateral hypergraph has vertices `E` and hyperedges

    {p,p+u,p+ζu},  p∈E, u∈U.                                (F6)

Each unordered triple is listed three times, so its hyperedge count is `q²(q+1)/3`. It is exactly the triangle hypergraph of the undirected norm-one graph, which has degree `q+1`; each graph edge lies in two triangles.

**FINITE-EXACT result.**

- At `q=5`, the 25-vertex, 50-hyperedge system is 2-colorable.
- At `q=11`, the 121-vertex, 484-hyperedge system is 2-colorable.
- At `q=17`, the 289-vertex, 1,734-hyperedge system is not 2-colorable.

The first two color vectors are saved in `results/finite_triangle_probe.jsonl` and checked by enumeration of every hyperedge. For `17` the field is represented as `F_17[t]/(t²-3)` and `ζ=9+2t`. The unit graph has 2,601 edges.

`results/q17_nae_certificate.json` is a complete DPLL refutation tree with 3,451 nodes and 1,726 conflict leaves, using 125,099 forced assignments. Its SHA-256 is

    3980c13e87296452da52c9d0f8a6da1a5bfbbf1130e9d8a19f1b00ec0052d9e6

A force on vertex `v` cites a triple whose other two vertices already have the same color, so `v` must take the other color. A conflict leaf cites a fully monochromatic triple. A branch checks both colors of an unassigned vertex. The root fixes vertex zero to color zero; global color complementation covers the other root color.

**Independent reconstruction.** The generator enumerates rotated/transformed templates F6. The verifier instead constructs all norm-one adjacency sets and all graph triangles by pairwise intersections, then checks equality with the certificate hyperedges. It checks each local inference without calling the search routine. All operations are integer arithmetic modulo 17 or finite Boolean checks. This is an externally replayable exact finite certificate, not a Lean/kernel-checked theorem.

Since the only smaller odd prime powers `q≡2 (mod 3)` are `5` and `11`, `17` is the least such `q` whose equilateral norm hypergraph is not two-colorable. This statement excludes characteristic two and the degenerate characteristic-three case. It does not classify all larger prime powers. The following forcing argument does settle all distinct templates at `17`.

### 3.1 Every three-point template is forced at q=17

**Theorem F1 (FINITE-EXACT input plus algebraic proof).** For every coloring of `F_(17²)` with two colors and every ordered triple of distinct points, some translation and norm-one rotation of that triple is monochromatic.

**Algebraic forcing lemma.** In any field with an involution and a norm-one element `ζ` satisfying `ζ²-ζ+1=0`, put `barζ=1-ζ` and, for arbitrary `a,z`, define

    A=0, B=a, C=aζ, D=z, B'=ζz,
    A'=a+ζz, C'=a barζ+ζz, D'=a-barζ z.

The six triples

    BAD, CAB', CBD', B'A'D', C'A'B, C'B'D

are norm-preserving motions of the same unordered template `{0,a,z}`. In their listed order their multipliers relative to the ordered template `(a,0,z)` are

    1, ζ, -barζ, -1, -ζ, barζ,

and each translation is its second displayed vertex. These six identities hold in `Z[ζ]/(ζ²-ζ+1)` with formal independent coefficients `a,z`.

If `ABC` has one color and none of the six template copies is monochromatic, the first three force `D,B',D'` into the other color. The fourth forces `A'` back to the first color; the fifth forces `C'` into the other color. Then `C'B'D` is monochromatic, contradiction. In particular,

    1_mono(ABC) <= Σ_(six copies) 1_mono(copy).                (F9)

Coinciding auxiliary vertices cause no difficulty: their labels must simply agree. The universal eight-label check is stronger than the check restricted to realizable coincidences.

For `F_(17²)` take `ζ=9+2sqrt(3)`. Normalize any triple by translation to `{0,a,z}`, where `a!=0`. Apply the q17 equilateral obstruction to the coloring `w↦c(aw)`. It supplies a monochromatic equilateral triple of side vector `a`, to which the forcing lemma applies. This proves Theorem F1. ∎

The exact script `scripts/check_eight_point_forcing.py` checks all six formal ring identities, each multiplier's unit norm, and F9 on all 256 eight-label signings. This independently validates the finite forcing used here; it does not validate the rest of the source blueprint.

Averaging F9 over all finite placements also proves

    D_T >= D_eq(a)/6 >= 2/(6*1734) = 1/5202.                (F10)

The stronger equilateral count of at least two used here follows from the cut-parity argument in §4 below. Since there are `289*18=5202` ordered placements, F10 is exactly one ordered-placement quantum, consistently with the direct proof. Neither F10 nor the refutation tree determines the sharp minimum count.

Thus `17` is the least odd prime power `q≡2 (mod 3)` at which every two-coloring forces every distinct three-point template: the equilateral-avoiding controls at `5,11` rule out this property there.

## 4. The radial positive-semidefinite hostile at the same q

The same exact rational enclosure script proves

    K_17(2) > 6.560736549697 > 6.

Put `K=K_17(2)` and `α=6/K∈(0,1)`. On the additive dual of `E`, let `μ_0` be the uniform distribution on all 289 characters and let `μ_2` be the uniform distribution on the 18 characters of dual norm two. Define the probability measure

    μ=(1-α)μ_0 + α μ_2.                                     (F7)

Its Fourier transform `ρ(z)=∫χ(z)dμ(χ)` is real, rotation invariant, positive definite, and has `ρ(0)=1`. For every norm-one `z`,

    ρ(z)=α [-K_17(2)/18] = -1/3.                             (F8)

Thus the radial PSD moment relaxation has exactly the pair value forced by zero monochromatic equilateral triangles. Every local equilateral `3×3` moment matrix has diagonal one and off-diagonal `-1/3`, and is positive semidefinite.

But `ρ` cannot be the pair-correlation function of any genuine jointly defined `{±1}`-valued random field on `E`: if it were, F1 would make the expected number of monochromatic hyperedges zero. Since that number is nonnegative and the hypergraph is finite, almost every sample would be a 2-coloring avoiding all 1,734 triples, contradicting the checked refutation.

### An explicit violated cut inequality

Let `C(s)` be the number of graph edges crossing the color cut and `M(s)` the number of monochromatic unordered triangles. Every nonmonochromatic triangle has two cut edges, and each graph edge belongs to two triangles. Thus

    M(s)=1734-C(s).                                           (F11)

Every degree is 18, so every cut has even size: `C=18|S|-2e(S)`. Consequently `M` is even. The checked refutation rules out `M=0`; therefore

    M(s)>=2,       C(s)<=1732.                               (F12)

For every globally realizable Boolean correlation matrix `(R_ij)`, expectation of F12 gives the explicit linear inequality

    Σ_(norm-one unordered edges {i,j}) (1-R_ij)/2 <= 1732.     (F13)

The PSD witness F7--F8 has every unit-edge `R_ij=-1/3`, so its left side is `2601*(2/3)=1734`, violating F13 by two. This exhibits a separating inequality on the entire norm-one graph, though its validity still uses the exact refutation plus parity. A small-support or elementary stand-alone separation remains a useful open question.

**Mechanism.** Bochner/PSD feasibility keeps spectral positivity; it loses global Boolean realizability, equivalently membership in the finite cut/correlation polytope. Local triangle consistency does not restore that missing coordinate. The blueprint starts from an actual binary process and retains a bounded field, so F7 is not a counterexample to its transfer. It is a counterexample to reversing the spectral necessary condition into a construction of an avoiding coloring.

## 5. Relation to the existing HN carriers

The closest inherited mechanisms are THM-4628's norm-circle Fourier multiplier and THM-4663's closed-walk/polygon ladder:
`01-canon/theorems/THM-4663-polygon-ladder-of-the-finite-unit-distance-graphs-closed-walks-are-kloosterman-moments-pentagons-level-15-hexagons-level-6.md`.

This session changes the observable from proper graph coloring / return walks to binary monochromatic triangle defect. F1 is the exact map. The cut/polytope hostile shows why the whole spectrum is still a relaxation for that observable.

**Vertex-Folkman consequence.** The q17 norm-one graph satisfies `G -> (K_3,K_3)^v`: every binary vertex coloring has a monochromatic graph triangle. Hence `χ(G)>=5`; otherwise, grouping at most four independent color classes into two groups of at most two classes would give a binary coloring with no monochromatic triangle. The lower bound `κ(17)>=5` was already recorded in the repository at `01-canon/theorems/THM-4558-residue-colourings-of-field-planes-the-heegner-tower-is-at-most-four-chromatic.md`; this session provides the stronger vertex-Folkman property, not a new chromatic lower bound. The q11 control is a useful converse hostile: that graph has the inherited chromatic value five, yet its triangle hypergraph is binary-colorable.

**Exact nonliftability.** There is no placement of the q17 graph's vertices in the real plane that maps every graph edge to a unit segment. The restriction of the real strip coloring to such a placement would give a binary vertex coloring with no monochromatic graph triangle: each graph triangle would map to three pairwise unit-separated points, hence a unit equilateral triangle. This contradicts the vertex-Folkman property. The argument even permits nonadjacent vertices to coincide; the vertices of each edge or triangle are necessarily distinct.

Reduction of a characteristic-zero carrier to a finite field preserves selected norm equations; it creates extra modular incidences and forgets archimedean geometry. F6 at `17` cannot be faithfully lifted with every hyperedge a real unit equilateral triangle, because the elementary alternating strip coloring of the real plane avoids all such triangles. Therefore this finite obstruction is an instrument for testing relaxation and reduction, not a new Hadwiger--Nelson lower bound.

The carrier source `THM-4700-unit-path-ring...` also retains prime ideals, denominators and actual norm equations. These sidecars must be carried if F2 is used against a specific residue coloring; a bare common Kloosterman symbol does not identify the underlying characteristic-zero coloring problem.

## 6. New exact questions

- **H-F1, OPEN.** Determine the smallest non-two-colorable subhypergraph of F6 at `17`, with a separately checked certificate. A finite certificate already exists; minimizing it is an optimization question, not an existence gap.
- **H-F2, RESOLVED globally / OPEN locally.** F13 is an explicit global cut-polytope inequality violated by F7. Find a smaller-support geometric inequality or a short independent validity proof that does not rely on the full refutation tree.
- **H-F3, RESOLVED in this session.** Every distinct template at `q=17` is forced, by Theorem F1. The shape-independent single-circle test is inconclusive there, but the eight-point forcing map transports the independently certified equilateral obstruction. A fresh OPEN version is to classify the non-equilateral avoiding shapes at `q=11`; keep its explicit equilateral-avoiding control.
- **H-F4, OPEN.** Determine the exact minimum monochromatic equilateral count at `17`. The certificate plus even cut parity proves it is at least two, not its minimum.

## 7. Sources checked

- Owner-supplied *Blueprint: monochromatic triangles in two-colorings of the plane*, all 670 lines read; sections 4--6 supply the sign, radial multiplier and bounded-field architecture. Its major theorem and formal-check claims remain source claims.
- Repository THM-4628 and THM-4663, with their corrected scopes, as linked above.
- Shkredov--Solymosi, [*Tilted corners in integer grids*, Integers 21A (2021), A20](https://math.colgate.edu/~integers/graham20/graham20.pdf), §4: prior finite-field/Kloosterman and real/Bessel triangle counting. Its asymptotic/theoretical method is acknowledged; its numerical real bound is not used here.
- Jelínek--Kynčl--Stolař--Valla, [*Monochromatic triangles in two-colored plane*](https://arxiv.org/abs/math/0701940), published in Combinatorica 29 (2009): the polygonal/closed-open results and finite forcing are prior context, not a source for the new finite certificate.
