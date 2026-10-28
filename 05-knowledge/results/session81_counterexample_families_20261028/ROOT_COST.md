# Root-charge debt versus tournament gaps: exact bridge and obstructions

## Findings first

**PROVED here:** Directed-edge charges have a precise closure law: `e_i-e_j` has a zero-cost merge tree using directed-edge roots iff `j` is reachable from `i`. Consequently every nontrivial strongly connected tournament on `n>=3` vertices has **exactly the same catalytic debt function as A_(n−1)**. Only the directed triangle among nontrivial strong tournaments has nonzero debt; every strong tournament of order at least four has identically zero debt. This is an actual mathematical connection, but it prevents this charge model from detecting the exceptional Hamiltonian-path values 7 and 21.

**PROVED here:** The cyclic-triangle skew matrix acts on A2 as multiplication by `sqrt(-3)`, an honest Eisenstein structure. Paley7 instead gives `sqrt(-7)` on A6. Neither the lattice nor its debt is a Hamiltonian-path counting invariant.

**FINITE-EXACT controls:** complete labelled tournaments through five vertices verify the free-root-path construction, and direct path/automorphism counts give `(H,Aut,tc)=(3,3,1)` for C3 and `(189,21,9)` for Paley7. Dynkin A4 has `I(P4,2)=21` but zero debt; affine A2's diagram K3 has `I(K3,2)=7`. Finite A1,A2 give 3,5, not 7,21.

## 1. Exact inherited definitions and corrected scope

For lattice L with quadratic norm Q and specified charge set R, a binary merge tree has leaves in R, node values equal to sums, and node cost

`max(beta(y,z),0), beta(y,z)=Q(y+z)-Q(y)-Q(z)`.

`g_R(x)` is the minimum total cost over all finite trees with root x, permitting arbitrarily many leaves. Set it to infinity when x is unreachable. A neutral subtree may be added at zero interaction cost, but not every representation change is automatically free. “Positive cost” means **not identically zero**, not every point has positive cost: all charge leaves have zero cost, as does zero when opposite charges are available.

**INHERITED PROVED/AUDITED:** [THM-4824, part G](https://github.com/eliottcassidy2000/math/blob/1a32e3ee6d4140929227f9ae66d46e27c343af65/01-canon/theorems/THM-4824-catalytic-merge-debt-by-lattice-shape-square-checkerboard-and-e8-free-first-grid-collinear-eisenstein-the-critical-half-shift-glue-parity-coset-shells-zero-sum-quartics.md): for a nonzero orthogonal sum of ADE root lattices, normalized so every root has `Q=x·x/2=1`, with all roots as charges, `g` is identically zero iff L is not A1 or A2. Thus “A2 is the only one” is correct only after excluding rank one. A1 has `g(n alpha)=n²−|n|`. A2 (Eisenstein norm `a²−ab+b²`) has `g(x)>=floor(Q(x)/3)`; equality is not universal (11 costs41 rather than40; 101 costs3402 rather than3400, per THM-4820's certified scope).

The audit upgraded the ADE classification from a conjecture using equal-length perpendicular helpers, support reduction for A_n, and index-two glue for E6,E7,E8. Reducible A1⊕A1 and A2⊕A1 are free despite their nonfree factors: catalysts may leave a component and return. Positive debt is not preserved under orthogonal inclusion. Example: if alpha,h are perpendicular unit roots, `2alpha=(alpha+h)+(alpha-h)` is free, including the two first-level merges; A1 alone charges2 for the same target.

**Caveats:** “root lattice” here is ADE with a specified metric and charge set, not arbitrary crystallographic root systems. G2 with minimal vectors means only its short roots; that lattice is A2 and pays. Allowing all G2 roots changes the problem and removes the equal-unit-norm hypothesis. B2's lattice with its minimal vectors is a square lattice and free. The report does not assert an all-G2-roots classification. Unreachable directions for a lattice with only collinear minimal charges must not be called free.

## 2. Reachability theorem: a structural map with a decisive limitation

Let D be any finite loopless digraph on `[n]`, put `R_D={e_i-e_j:i→j}`, `L=A_(n−1)`, and use `Q(x)=sum x_i²/2`.

**Theorem R (PROVED).** For distinct i,j, the following are equivalent:
1. D has a directed path from i to j.
2. `e_i-e_j` is representable as a nonnegative integer sum of R_D.
3. `e_i-e_j` has a zero-cost merge tree with leaves in R_D.

*Proof.* Choose a simple path `v0=i,...,vk=j`. Its edge roots telescope. Merge left to right: the accumulated root `e_i-e_vt` and next edge `e_vt-e_v(t+1)` have dot product −1, hence beta −1 and no cost. Thus 1 implies3, and3 implies2. For2 implies1, if j is outside the set U reachable from i, no edge leaves U. The sum of coordinates over U is therefore nonpositive on each allowed edge root. It equals1 on `e_i-e_j`, contradiction. QED.

**Corollary S (PROVED).** If D is strongly connected and has at least two vertices, `g_(R_D)=g_(all A_(n−1) roots)` at every lattice point. One inequality follows from charge-set inclusion. For the other replace each leaf of an arbitrary full-root tree by the free directed-path tree from Theorem R. Root values and all original node interactions are unchanged. QED.

**Transitive hostile (PROVED).** In the transitive tournament `1→2→...→n` (all arcs forward), the target `2(e_i−e_(i+1))` has debt exactly2. For each prefix `{1,...,k}`, every edge charge has nonnegative prefix sum, equal to1 precisely when the edge crosses that cut. The target prefix sums are zero except at k=i, where the sum is2. Every leaf must therefore cross only cut i, so the only permitted leaves are two copies of `e_i−e_(i+1)`. Their merge costs2. This also proves why the strong-connectivity assumption cannot be removed from Corollary S. Unlike the fully symmetric-root map, the directed-charge map does preserve some orientation information before strong closure.

For strong tournaments, n=3 forces C3 and hence A2; n>=4 gives free A_(n−1) by corrected THM-4824. There is no strong two-vertex tournament. The singleton K1 is excluded: its empty charge set would require an explicit empty-tree convention for zero. This yields a **precise bridge** from orientation to charge closure, while proving that global optimization destroys all orientation information inside a strong component.

**Mapping ledger:** source=strong tournaments; target=oriented A-root charge systems; map=`i→j ↦ e_i-e_j`; preserved predicate=directed reachability as free root synthesis; destroyed data=orientation after catalytic closure, Hamiltonian path count, automorphism group; required sidecar=original tournament or explicit constrained leaf multiset/allowed tree grammar. No claim that such a sidecar already encodes the 7/21 proof.

**Counterexample generator:** take any two strong tournaments of the same order with distinct H (at n=5 the exhaustive values are9,11,13,15). Their charge debt functions are identical, identically zero. For the symmetric-root map even transitive T3 (H1) and C3 (H3) share the same nonfree A2. To refute componentwise debt inheritance, embed A1 or A2 orthogonally into its sum with A1.

## 3. Genuine Eisenstein skew connection, and its septimal boundary

Let C be the permutation matrix of a cyclic triple, and S=C−Cᵀ. On the augmentation lattice `A2={x∈Z³:sum x=0}`, `I+C+C²=0`, so

`S=I+2C, S²=−3I`.

Moreover `A2` is a free rank-one module over `Z[omega]`, with omega acting by C: `alpha=e1−e2` and `C alpha` are an integral basis, and their norm is the Eisenstein norm up to the choice of rotation convention. Thus S acts as `1+2omega`, of norm3. This is an exact operator identity, not a visual analogy. Also `SᵀS=3I` on the sum-zero plane, so it multiplies Q by3.

For any doubly regular tournament of order n, the skew matrix satisfies `S²=J−nI`. This can be derived by counting, for distinct vertices, their `(n−3)/4` common out-neighbors and common in-neighbors. There are `(n−3)/2` equal-sign and `(n−1)/2` opposite-sign terms in row correlations, giving `SSᵀ=nI−J`. On the augmentation lattice S acts with square `−nI`. For n=7 this is an action of `Z[sqrt(-7)]`, not the Eisenstein action. The lattice rank is6, not2, and its all-root merge debt is zero. This does not rule out more elaborate representations; it rules out identifying these particular operators/charge models.

A general tournament's S need not preserve the sum-zero plane: preservation is equivalent to regularity. Transitive T3 is a hostile control; an unrestricted “every tournament has the Eisenstein operator” assertion fails even there.

For classical context only, Reid and Brown, *Doubly regular tournaments are equivalent to skew Hadamard matrices* (1972), [primary publisher record](https://www.sciencedirect.com/science/article/pii/0097316572900982), gives the DRT/skew-Hadamard relation. The identities above have independent proofs here; no literature novelty is asserted.

## 4. Two separate 7/21 exclusion theorems

`H(T)` counts Hamiltonian paths as vertex words. `tc(T)=H(T)/|Aut(T)|` counts path presentations modulo automorphisms: the action is free because an automorphism fixing a vertex word fixes every vertex. They are distinct invariants.

**INHERITED VERIFIED:** THM-1370, `h-spectrum-omits-7-21-all-n`: neither7 nor21 occurs as H. The mechanism is strong-component multiplicativity, strong-order minima, and the small exact spectra. The strong minimum is

`f(3k)=3*5^(k−1), f(3k+1)=5^k, f(3k+2)=9*5^(k−1)`

in the valid orders n>=3. It is verified against [Busch2006, Corollary1, primary PDF](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf). The first minima3,5,9,15,25,45 exclude large strong witnesses, while a product7 or21 needs an impossible strong factor7 or21. Completeness of the odd spectrum except7,21 remains OPEN; H63 is attained, so no forbidden tower7*3^k is valid.

**INHERITED PROVED WITH CITED AND FINITE CERTIFICATES, not independently rerun here:** corrected/audited THM-4812, `presentation-fibres-7-and-21-are-absent-at-every-order-orbit-product-lemma-vertex-transitive-lemma-closed-window`: neither7 nor21 occurs as tc either. This is not a division corollary of THM-1370. The proof uses the orbit-product inequality, vertex-transitive exceptions K1,C3,R5,P7 with tc1,1,3,9, and a bounded configuration/subgroup scan. The source records independent audit repairs to the integer tail, infinite order gate, and index/orbit typing, and records24 configurations/83 subgroup cases. Its further v23 extension has a different audit scope. This report preserves that boundary.

Both theorems forbid statistic values, not vertex orders or automorphism orders. Paley7 has21 automorphisms; the checker independently obtains H189 and tc9.

## 5. The Dynkin/independence-polynomial trap

The inherited OCF gives `H(T)=I(Omega(T),2)`, where vertices of Omega are literal directed odd cycles and edges indicate shared vertices. Therefore its image cannot include any graph whose evaluation is7 or21. But arbitrary graphs can: K3 gives7, P4 gives21.

For the finite Dynkin path A_r, `a_r=I(P_r,2)` satisfies `a_0=1,a_1=3,a_r=a_(r−1)+2a_(r−2)`, hence `a_r=(2^(r+2)−(−1)^r)/3`. Thus A1,A2 yield3,5; A4 yields21 even though A4 is free. K3 is the affine triangle diagram Ã2, not the finite A2 diagram. No identification of an abstract Dynkin graph with a realizable odd-cycle conflict graph is justified merely by a shared drawing or rank.

## 6. Reproduction, hostile controls, status boundaries

Run `python3 check_root_cost.py` from this artifact directory. Pure Python integer arithmetic; prints its verification results as JSON. All validity checks use explicit exceptions and remain active under `python3 -O check_root_cost.py`; normal and optimized outputs were compared byte-for-byte.

Universe: all labelled tournaments n2..5; every ordered root pair in each strong tournament; explicit Paley3/7 with all permutations for automorphism counts; finite Dynkin paths A1..A6 and K3. The run checked12+288+10880 free directed-path root trees. Strong labelled counts2,24,544 and spectra{3},{5},{9,11,13,15} match the known small controls. It also checks the perpendicular-helper certificate. No finite run is used to establish an all-order classification: Theorems R/S and the skew statements are proved symbolically above; ADE classification and the exclusion theorems are inherited with pinned exact source bytes.

**OPEN:** a more informative map retaining tree constraints might connect root debt to tournament orbit budgets. No equivalence sending A1,A2 to the missing values7,21 has been constructed. The exact bridge established here shows why the most natural unrestricted-charge attempt cannot do so.

Source baseline: repository commit `1a32e3ee6d4140929227f9ae66d46e27c343af65`. The cited theorem slugs distinguish historical identifier collisions. An independent audit checked all loopless digraphs through four vertices and all labelled tournaments through six; the public checker retains the smaller stated scope.
