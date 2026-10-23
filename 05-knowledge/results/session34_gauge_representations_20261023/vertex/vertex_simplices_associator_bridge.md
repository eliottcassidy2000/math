# Session 34: from tournament simplices to vertex locality, with the missing associator retained

Status: elementary statements below PROVED in this note; finite checks FINITE-EXACT; standard category/physics interfaces CITED; realization of arbitrary tournaments as noncommutative vertex algebras OPEN / not asserted. Novelty not claimed. No repository mutation was performed.

## 1. Inheritance before invention

Closest mechanism: classical combinatorics can retain leading symbols while losing the lower-order compatibility correction. The repository already makes this exact through the canonical surjection

`J_infty(R_V) -> gr_F(V)`.

The map is not generally an isomorphism. Its kernel is precisely information lost when one tries to reconstruct quantum relations merely by differentiating the degree-zero relations.

Recovered canon, with corrected scope:

- **THM-4671**, `01-canon/theorems/THM-4671-virasoro-minimal-model-first-quantum-correction-to-the-arc-space.md`: for coprime `3 <= p < q`, `N=(p-1)(q-1)/2`, the first kernel occurs at weight `2N+3`, is one-dimensional, and is generated modulo the arc ideal by `b_N=x_0^(N-2)x_1x_2+(2N-7)x_0^(N-3)x_1^3/12`. The Euler syzygy lifts with coefficient `A=N(2N+1)(c-c_(2,2N+1))/8`. For boundary `p=2` the coefficient vanishes. The corrected statement is that the adjusted lift vanishes, not the unadjusted vector `u` itself. The `(3,q)` case has explicit prior art; uniform larger-p form is the repository extension.
- **THM-4695**, `01-canon/theorems/THM-4695-first-order-universality-of-the-quantum-transfer-for-virasoro-minimal-models.md`: the entire first transfer is `(c-c_(2,2N+1)) T*_N` modulo the arc ideal. This uses classical freeness of the boundary series as a cited input. Higher transfer universality remains bounded/conjectural.
- **THM-4710**, `01-canon/theorems/THM-4710-witt-rigidity-of-the-degree-n-singular-support-generators-of-virasoro-minimal-models.md`: exact rational-function certificates prove the specified Witt-admissible spaces and forms at Li degrees `3,5,7,9`; identifying those with all actual minimal-model generators is CONDITIONAL on the generator positions in HYP-9315(b). The older result-note sentence calling the closed forms conjectural is stale; the theorem frontmatter and MISTAKE-641 supersede it. Do not drop the generator-position hypothesis.
- **THM-4672**, `01-canon/theorems/THM-4672-free-hypermultiplet-z2-gauging-chiral-algebra-is-not-classically-free.md`: `L_(-1/2)(sl2)` has first kernel of dimension five at weight 7; arc count 324 versus character 319. This is an exact rational computation inherited here, not rerun. The reduced Higgs branch cannot distinguish it from the `N=4` SU(2) SYM example; the nonreduced associated schemes already differ. Classical freeness of the latter through weight 9 is finite evidence, not an all-weight theorem.
- **THM-4673**, `01-canon/theorems/THM-4673-n2-minimal-model-c1-lattice-sqrt3-is-not-classically-free.md`: no unitary `N=2` minimal model in the stated family is classically free. Pure fermion symbols supply a kernel witness. This concerns the holomorphic vertex algebra of the IR minimal model, not the different half-twisted Landau-Ginzburg cohomology or simply the arc space of `Crit(W)`.

Canonical hostile: `Vir_(4,5)` and `Vir_(2,13)` have identical `R_V=C[x]/(x^6)` but different singular supports. Equal degree-zero geometry is insufficient even before discussing modularity.

Least-used relevant sidecar: actual syzygy lifts and their first nonvanishing filtered corrections. This motivates keeping associator/coherence coefficients, not just a simplex incidence graph.

Concept board: tournament winner product; simplicial boundary; vertex locality; fusion associator; quantum-transfer defect; protected Yang-Mills sector. The exact maps below connect the first four; the last two supply inherited scope constraints rather than a claimed equivalence.

## 2. An exact tournament-to-vertex toy and its sharp boundary

Let `T` be a tournament on `n>=1` vertices. Over `C`, put `A_T=span{e_v}` and define the bilinear winner product

`e_v e_v=e_v`, and `e_u e_v=e_v e_u=e_w`, where `w` wins the edge between `u,v`.

Adjoin a new unit to obtain `A_T^+=C1 ⊕ A_T`. This is a commutative, possibly nonassociative unital algebra. The new unit is independent even when `A_T` already has its own unit.

**Theorem A.** The following are equivalent:

1. `T` is transitive.
2. The winner product is associative.
3. `Y(a,z)b=ab`, with translation operator zero and vacuum `1`, gives a vertex algebra on `A_T^+`.

**Proof.** If `a beats b`, `b beats c`, `c beats a`, then

`(e_a e_b)e_c=e_c`, whereas `e_a(e_b e_c)=e_a`.

Thus associativity forbids a directed triangle. A tournament without directed triangles is transitive: a shortest directed cycle of length at least four has a chord shortening it, and hence a directed cycle would force a triangle. A transitive tournament admits an ordering in which the product selects the earlier of its two inputs, and taking the earliest of three is independent of bracketing.

For any commutative associative unital algebra, constant multiplication fields satisfy vacuum, creation, translation, and locality, hence the vertex algebra axioms. Conversely, locality of constant fields requires all left multiplication operators to commute, since a nonzero constant operator-valued polynomial cannot be killed by any power of `z-w`. Commutativity plus `[L_a,L_c]=0` implies `a(bc)=c(ab)=(ab)c`. This proves the converse without importing a classification theorem. QED.

**Minimal explicit locality witness.** On the cyclic triple above,

`[Y(e_a,z),Y(e_b,w)]e_c=e_a-e_b != 0`.

No power of `z-w` kills it. A cyclic triangle is already enough; no fourth-order correction is required to detect this failure.

**Scope.** This is a theorem about the specified winner construction, not a prohibition on other constructions assigning a vertex algebra to a tournament. The transitive cases are degenerate commutative vertex algebras, not a construction of an interesting simple CFT-type VOA. No central charge, nontrivial OPE pole, conformal grading, or modular representation category is produced.

## 3. Exact loss ledger: imposing associativity collapses SCCs

Let `I_assoc` be the ideal of the nonassociative algebra `A_T^+` generated by all associators. Let the strongly connected components of `T` be `C_1,...,C_m`, ordered by their condensation tournament.

**Theorem B.** There is a canonical quotient, with its component-labelled basis retained,

`A_T^+ / I_assoc ≅ A_(T/SCC)^+ ≅ C^(m+1)`.

The first quotient kernel is exactly the span of `e_u-e_v` for vertices in the same strongly connected component.

**Proof.** Distinct SCCs of a tournament have uniform edge direction between them: opposite directions would provide paths both ways using strong connectivity. Their condensation is acyclic and therefore transitive. Sending each `e_v` to its SCC basis element is a surjective multiplicative map to an associative winner algebra. Its kernel is the stated difference space `K`, so `I_assoc ⊆ K`.

In any associative commutative quotient, all images of `e_v` are idempotents. On idempotents define `x <= y` iff `xy=x`. This is a partial order: reflexivity follows from idempotence, antisymmetry from commutativity, and transitivity from associativity. Each tournament edge imposes one such inequality. A directed path imposes the corresponding inequality between its endpoints, so mutual reachability forces equality. Thus `K ⊆ I_assoc`.

For the final isomorphism, order the component idempotents `e_1,...,e_m` so `e_i e_j=e_min(i,j)`. Then `p_1=e_1`, `p_j=e_j-e_(j-1)` for `2<=j<=m`, and `p_(m+1)=1-e_m` form an independent complete family of orthogonal idempotents. QED.

Consequences:
- Every strongly connected tournament, regardless of its internal complexity, has the same two-dimensional unital associative quotient `C^2`.
- Abstractly the constant-field VA retains only `m`, the number of SCCs. To retain which vertices belonged where, keep the quotient map and component-labelled basis. Even the chain order is not intrinsic to `C^(m+1)` without that selected basis.
- Nonzero derivations of `C^(m+1)` do not exist: differentiating `p^2=p` gives `(1-2p)D(p)=0`, and `1-2p` is invertible. Thus adding a derivation to this same associative quotient cannot restore nonconstant vertex fields.

This is the sharpest exact obstruction produced here: the cheap locality repair destroys precisely the directed strongly connected information one was hoping to model.

## 4. Simplicial boundary is not vertex Jacobi

The full simplex on `{a,b,c}` has oriented boundary

`∂[a,b,c]=[b,c]-[a,c]+[a,b]`, and `∂²[a,b,c]=0`.

This is true independently of which tournament arrows were assigned to its edges; changing the orientation gauge of a chain basis changes incidence signs consistently, not the identity `∂²=0`. Apply this to the directed three-cycle from Theorem A. Its simplicial chain complex is valid while its multiplication fields violate locality. Thus signed simplicial cancellation does not imply vertex locality or the Borcherds/Jacobi identity.

If one instead chooses the directed flag complex, cyclic triples are not two-simplices at all. That is a different carrier, not a cure for the winner product's nonassociativity.

Nor can one repair a degree-zero algebra with zero differential merely by announcing an `m_3` or a fourth-order operation. The arity-three `A_infinity` identity says the associator is a sum of terms involving the differential and `m_3`; with differential zero, it forces strict associativity of `m_2`. A genuine chain-level repair must enlarge the complex and make the associator a boundary, then satisfy further coherence identities. This does not itself construct a vertex algebra.

## 5. Positive finite bridge: a simplicial associator with arithmetic phases

The suitable minimal positive carrier is not an arbitrary tournament but a category retaining coherent fusion data.

Take two simple labels `0,1` with fusion `a⊗b=a+b mod 2`. On homogeneous vector spaces choose scalar associator

`ω(a,b,c)=(-1)^(abc)`.

This normalized cocycle satisfies the pentagon:

`ω(b,c,d) ω(a,b+c,d) ω(a,b,c) = ω(a+b,c,d) ω(a,b,c+d)`.

All arguments are in `Z/2`. Equality follows by expanding exponents modulo two; both sides have exponent `bcd+abd+acd+abc`. This is exactly the four-input coherence equation, or the cocycle equation on a bar-construction 4-simplex. An associator itself is triple data. Do not call the pentagon a fourth derivative or a Yang-Mills equation.

The trivial cocycle has the identical simple labels and identical fusion multiplication. Yet these associators are not related by a normalized 2-cochain gauge: for any such `β`,

`(δβ)(1,1,1)=β(1,1)β(1,0)/(β(0,1)β(1,1))=1`,

whereas `ω(1,1,1)=-1`. Since `Z2` has no nontrivial automorphism, the two pointed fusion categories are not monoidally equivalent. The general construction and equivalence criterion are standard; see Etingof-Gelaki-Nikshych-Ostrik, *Tensor Categories*, Example 2.3.8 and §2.6: https://math.mit.edu/~etingof/egnobookfinal.pdf .

One can add `r(a,b)=i^(ab)` as braiding. In these conventions both scalar hexagons hold:

`ω(a,b,c)r(a,b+c)ω(b,c,a)=r(a,b)ω(b,a,c)r(a,c)`,

`r(a+b,c)=ω(a,b,c)r(b,c)ω(a,c,b)^(-1)r(a,c)ω(c,a,b)`.

The only nontrivial triple requires `r(1,1)^2=-1`, which `i` satisfies. The double braiding on `(1,1)` is `-1`, unlike trivial braiding. This is a concrete cyclotomic arithmetic coordinate invisible to the fusion multiplication. It is a finite coherence example, not a reconstruction theorem for a VOA.

**Meaning of irreducibles as primes.** In a semisimple category the simple labels are unique additive decomposition atoms. They are not generally multiplicative primes: tensoring two simples can decompose into several simples, and an invertible simple can square to the unit, as `1⊗1=0` here. Fusion multiplicities, associators, braidings, and duality are distinct additional coordinates. A tournament on the labels omits even the fusion rule; the fusion rule alone still omits `ω`.

## 6. How far the actual arithmetic/topology/physics interfaces go

- A **fusion category** supplies finite semisimple rigid tensor data and associativity coherence. It need not be braided.
- A **modular tensor category** additionally has braided ribbon structure with nondegeneracy. Suitable rational, `C2`-cofinite, self-dual, simple CFT-type VOAs have such module categories under the hypotheses of Huang's theorem. Not all vertex algebras do. Primary source verified: https://arxiv.org/abs/math/0502533 .
- A **factorization algebra** supplies compatible operations on disjoint opens together with descent. A simplicial index or incidence matrix is not that structure. Holomorphic/translation-compatible factorization constructions can lead to vertex algebras under additional hypotheses. The inherited factorization note is a guide, not a license to identify these objects without checking its hypotheses.
- The 4d/2d construction of Beem et al. uses a protected cohomological sector of `N=2` superconformal field theories, not arbitrary four-dimensional Yang-Mills. Primary source verified: https://arxiv.org/abs/1312.5344 . The existence of a vertex algebra in that sector is not a construction of the full theory or a mass-gap proof.
- The classical-limit map is a quotient with genuine additional relations: primary Arakawa-Linshaw source verified at https://arxiv.org/abs/1804.01287 . This supports retaining the transfer sidecar from §1 rather than replacing quantum identities by simplicial signs.

The source-to-target maps, preserved predicates, and losses are therefore explicit:

1. Tournament -> winner algebra preserves all arrows in the selected basis; locality fails exactly on directed cycles.
2. Winner algebra -> associative quotient preserves SCC condensation; loses all internal SCC information.
3. Fusion category -> fusion ring preserves simple multiplicities; loses associator and braiding (the `Z2` example proves this loss is real).
4. VOA -> associated `C2` scheme preserves leading commutative data; loses filtered relations and OPE information (THM-4671 and THM-4672 supply exact hostiles).

## 7. Reproducibility and next decisive step

Run `python3 session34-vertex/check_tournament_vertex.py` from the workspace root.

- Exhaustive labelled tournaments for `n=1..6`, totals `1,2,8,64,1024,32768`.
- Associative winner products: `1,2,6,24,120,720`, exactly `n!`; independent directed-triangle test agrees in every case.
- Independent transitive-closure SCC computation agrees with the span of cyclic-triangle associator differences for every tournament in that universe. The latter stronger span observation is only reported as FINITE-EXACT here; Theorem B's ideal statement is proved independently.
- `Z2` pentagon: all 16 quadruples checked. Two hexagons: all 8 triples each checked exactly using fourth roots of unity.
- Outputs: `check_results.txt` (full check); `check_results.json` (initial theorem-A/cocycle subset).

Natural next experiment: do not search for another unsigned graph invariant. Choose a specific retained carrier (e.g. a category with simple labels and prescribed fusion) and solve its finite pentagon/hexagon equations before attempting any VOA realization. For the tournament lane, test a chain enlargement where cyclic associators become actual boundaries while tracking whether SCC information survives in homology. That is a bounded new construction problem, not yet a result or evidence for Yang-Mills.
