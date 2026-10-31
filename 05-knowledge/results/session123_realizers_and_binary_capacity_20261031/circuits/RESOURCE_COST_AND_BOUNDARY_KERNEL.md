# Cost semantics, protected quadratic coordinates, and shared forest repair

Status: proved finite-dimensional statements and explicit implementation lemma; not an unrestricted Boolean-circuit lower bound. Session123, 2026-10-31. This packet extends the interfaces of Session118 sequence accounting and Session120 `FOREST_CONE_DICHOTOMY.md`. Classical quadratic normal forms and linear structures are not claimed as new.

## 1. A strict additive cost obstruction

**Theorem 1.** Let a category have a terminal object and let a real-valued arrow invariant c satisfy c(gf)=c(g)+c(f). Then

`c(f:X→Y)=c(!X)−c(!Y)`.

Proof: `!Y f=!X`; apply additivity. Thus every endomorphism has cost zero; if deleting is free all arrows have cost zero. Independently, in a monoidal category with natural diagonals and tensor-additive c, naturality `ΔY f=(f⊗f)ΔX` gives `c(f)=c(ΔY)−c(ΔX)`. Free natural copying again forces zero. These conclusions require scalar cancellation, not arbitrary idempotent cost monoids.

This does **not** say closed symmetric monoidal categories collapse: closedness alone supplies neither natural diagonals nor terminal deletion. It says gate cost cannot descend to a strictly additive equality-invariant cost on cartesian extensional semantics. A circuit C with g gates followed by free copying has g gates; duplicate C on two copied inputs has 2g gates. Their functions agree. Minimum semantic gate complexity instead satisfies composition/tensor **inequalities**, not the proposed equalities. A realizability/topos approach must retain costed realizers or intensional implementations and specify which equations preserve costs. It cannot silently identify shared and duplicated realizers.

## 2. Exact owned-coordinate compression with protected boundaries

Let b∈F2^p be protected independent input coordinates and u∈F2^m be owned independent coordinates. Let F=(f1,…,ft) be every output whose value must survive, including all exposed intermediate outputs. Allowed affine input charts fix b and have `u=Dv+Cb+a`, D invertible.

Define the zero-translation kernel

`K={h∈F2^m : F(b,u+h)=F(b,u) for every b,u}`.

**Theorem 2.** K is a vector subspace. The largest number of simultaneously unused owned coordinates in any allowed chart is exactly dim K. Equivalently, the least number of essential owned coordinates is m−dim K.

Proof: invariance under h and k implies invariance under h+k. If d chart coordinates are unused, their columns in D span a d-dimensional subspace of K. Conversely, choose a basis of K, extend it to a basis of F2^m, and take these basis vectors as columns of D, with kernel columns last. The last dim K chart coordinates do not affect F. The boundary-dependent translation Cb+a changes no owned translation direction. The chart remains invertible on the entire original domain: unused coordinates are retained, not removed from domain dimension. ∎

**Quadratic algorithm.** Write each fj in algebraic normal form of degree at most two on x=(b,u), and let Pj be its alternating polar matrix. First compute

`R={h∈F2^m : Pj(0,h)=0 for every j}`.

On R the map `ℓj(h)=fj(0,h)+fj(0,0)` is linear. Then

`K=ker(ℓ:R→F2^t)`, so `dim K=dim R−rank ℓ`.

Indeed the derivative is `fj(x+(0,h))+fj(x)=xᵀPj(0,h)+ℓj(h)`. Gaussian elimination on stacked restricted polar columns computes R; evaluating on a basis of R gives ℓ; another elimination gives K. Extending its basis constructs the chart. This is polynomial time, and is a boundary-aware simultaneous version of classical linear-structure calculations. A constant derivative of ONE is not in K.

Hostile controls:
- `f=bu`, with b protected and u owned: owned-owned polar block is zero but K=0. Dropping boundary coupling is unsound.
- `(uv,uw)`: each output alone has one unused direction, but common K=0. Per-output compression dimensions cannot be added.
- Triangle `xy+yz+zx`: polar radical is span(111), but derivative there is one; K=0.
- Four-cycle `(x1+x3)(x2+x4)`: K has dimension two. Adding exposed outputs x3 and x4 makes K=0.

This theorem concerns independent full-cube inputs. It does not automatically apply to marked-pair subdomains or arbitrary restriction domains. It certifies semantic affine compression, **not ordinary-coordinate gate saving**. An m-dimensional dense forward chart alone can require XOR circuitry (naive rowwise bound m(m−1) gates for the owned linear part); boundary mixing costs extra. It does not authorize unproved M/D4 row moves. Multioutput essential-coordinate count is not a gate lower bound: a wire-output map has zero gates.

For one nonaffine quadratic with all coordinates owned, the classical Dickson form is `Σ(i=1..r) pi qi + εz+c`, where polar rank is 2r and ε is the rank of the radical's linear restriction (zero or one). Thus there are exactly 2r+ε essential coordinates in the best affine chart. Across charts its binary arbitrary-Boolean-gate minimum is 2r+ε−1: the displayed form attains it; a one-output fan-in-two cone depending on k inputs has at least k−1 gates. Constants can be absorbed into its last gate. This is a corollary of classical classification, not an original-coordinate lower bound or a new classification.

## 3. A partially shared forest repair in original coordinates

Let F=(V,E) be a nonempty forest, e=|E|. Let δ be either0 or1, and let z be an additional Boolean input wire. A cone computes `Q=⊕(ij∈E) xi xj ⊕ δz`, with one AND per edge and an exclusive XOR reduction tree, initially 2e+δ−1 gates. The root may have arbitrary fanout. Some edge-product outputs P⊆E have external consumers and must be retained. No other internal outputs except these AND outputs are shared. Boundary input wires may have arbitrary external consumers. Write U=E\P, p=|P|, u=|U|, and ν(U) for its maximum matching size.

**Theorem 3.** There is an ordinary-coordinate replacement retaining all exposed product values and the root value, with

`e+p+ν(U)+δ−1` gates,

saving exactly `u−ν(U)` against the specified original implementation. Thus a strict saving exists whenever U has two adjacent edges, even though other product outputs are shared.

Proof: keep each of the p protected AND gates. On unprotected edges use the Session120 leaf-factor algorithm: choose a leaf x of the remaining forest, neighbor y, emit `y(x+Σu∈N(y)\{x}u)`, and delete x,y. These pivots form a maximum matching (the usual leaf exchange proves inductively that some maximum matching contains xy). There are ν(U) AND gates and u−ν(U) input XOR gates. XOR the p protected products, ν(U) pivot products, and δz using p+ν(U)+δ−1 gates. Total p+u+p+ν(U)+δ−1 is the stated expression. U=∅ is allowed; since E nonempty the output list is nonempty. Every discarded AND was unprotected and every discarded XOR exclusive, so the root can replace the old root without changing outside signals. ∎

**Exact obstruction threshold.** The certified saving vanishes exactly when U is a matching: ν(U)=|U|. Consequently zero-saving protection sets are precisely vertex covers of the line graph of F, and their minimum cardinality is e−ν(F). In particular, fewer than e−ν(F) exposed edge products cannot block every repair, and the certified saving is at least max(0,e−p−ν(F)). This is a complete classification of the obstruction for this repair family, not a claim that a zero-saving cone is globally optimal.

Example: edges ab,bc,cd with protected ab. Original size five. Retain ab and replace bc+cd by c(b+d), then XOR with ab: four gates. If ab and bc are both protected in the two-edge path, no saving is certified. If ab alone is protected in that path, U is a single edge and again no saving: this reproduces the shared-product hostile example. The theorem does not cover externally consumed XOR intermediates; expose them as additional outputs or prove a separate repair.

## 4. How this changes the lower-bound search

The improvement is a checkable structural alternative: an allegedly tight residual monomial forest cone must have its **unprotected** edges form a matching, unless a valid repair was already used. Record exact protected-output incidence, not merely total gate count or polar rank. In combination with Session118, a global improvement would still need a theorem forcing enough distinct such repairs across elimination steps and proving their gains are not counted twice. This packet does not prove that accumulation theorem.

## Primary sources and inherited results

- Claude Carlet, *Boolean Functions for Cryptography and Coding Theory*, author-hosted book, §5.2, Theorem10 (Dickson forms): https://www.math.univ-paris13.fr/~carlet/book-fcts-Bool-vect-crypt-codes.pdf . The quadratic discussion here is a self-contained relative refinement of classical material, not a priority claim.
- Cambridge L118 authored category-theory notes, §4.3, Theorem4.28, cartesian structure and uniform copying/deleting: https://www.cl.cam.ac.uk/teaching/2324/L118/notes.pdf . The elementary scalar-cost obstruction above is proved directly; it does not invoke compact-category no-cloning theorems.
- Repository: Session120 `circuits/FOREST_CONE_DICHOTOMY.md`; Session118 `accounting/SEQUENCE_ACCOUNTING.md`.

`check.py` independently verifies the quadratic kernel formula against full truth tables and the protected-forest rewrite against assignments. Small exhaustive checks support implementation correctness; the proofs above carry all dimensions.
