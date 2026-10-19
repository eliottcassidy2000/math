# Mark what must not repeat: a square-zero bridge from Carmichael products to graph cycles

Session 2026-10-19. Status: **PROVED** finite algebra and combinatorics; **FINITE-EXACT** tests in `verification.json`. No analytic prime-supply theorem, Carmichael counting asymptotic, or tournament complexity conclusion is proved here.

## 0. Inheritance and the connection board

The attached Howe–Huffman–Kohlmeyer–Kruer paper's equation (3.1) retains **distinct prime labels even when residues coincide**. Its Lemma 5.1 multiplies the independent factors `(1-rho)+rho*chi(p)`, and Proposition 5.2 uses squarefreeness and unique factorization. Those labels are mathematically indispensable.

Closest repo mechanisms:

- [THM-1926, tournament zeta](../../../01-canon/theorems/THM-1926-tournament-zeta-euler-product-over-strong-core.md): `-log det(I-uA)=sum tr(A^k)u^k/k`. It counts periodic closed walks; “primitive orbit” means not an iterate, not vertex-simple and not a minimal product-one subset.
- [THM-4663, finite unit-distance polygon ladder](../../../01-canon/theorems/THM-4663-polygon-ladder-of-the-finite-unit-distance-graphs-closed-walks-are-kloosterman-moments-pentagons-level-15-hexagons-level-6.md), statement 1: `W_k=|G|^-1 tr(A^k)` for an abelian Cayley graph. Its own cycle corrections retain repetitions; the uncorrected moments are not subset counts.
- [THM-1870, cycle spectral boundary](../../../01-canon/theorems/THM-1870-the-cycle-count-spectral-boundary-is-the-hamiltonian-length.md) contains an exact six-vertex witness but an incorrect all-order extrapolation. Section 6 supplies its repair.
- [THM-165, cycle-disjointness](../../../01-canon/theorems/THM-165-coefficient-2-theorem.md): direct cycles, not their vertex supports alone, are the correct cycle inventory. We do not import its historical scope claims.
- [HYP-9338, Apollonius/Carmichael dictionary](../../hypotheses/HYP-9338-enumerative-geometry-apollonius-ford-descartes-and-carmichael-transfer-dictionary.md): its exponent-one conclusions inherit the attached paper. The identities below are unconditional finite statements and do not audit those analytic inputs.

Live board: prime-label exclusion; vertex exclusion; cyclic orders; minimal zero-sum circuits; character/code duality; spectrum versus marked determinant. The decisive operation is **choose the coordinate that must not repeat, then make its square zero**.

## 1. The finite labelled carrier

Let `H` be a finite abelian group, `h=|H|`, and let `g_1,...,g_n` be an indexed list in `H`. Repeated group values are allowed; the indices remain distinct. Let `R_g` act on `Q[H]` by translation. Thus

`R_g R_f=R_(gf)`, and `tr(R_g)=h` if `g=1`, otherwise `0`.

Work over the commutative ring

`S_n = Q[x_1,...,x_n]/(x_1^2,...,x_n^2)`.

Its basis is `x_B=product_(i in B)x_i`, for subsets `B` of `[n]`; multiplication kills precisely intersecting label sets. Put

`A(x)=sum_i x_i R_(g_i)`.

Every entry has positive degree, so `A^(n+1)=0`; every formal series below is finite in `S_n`. Write

`Z_H(x)=sum_(B: product_(i in B)g_i=1) x_B`, including the empty subset.

### Theorem SQ1: subset exponential and cyclic-order logarithm

In `S_n`,

1. `exp A = product_i (I+x_i R_(g_i))`;
2. `h^-1 tr(exp A)=Z_H(x)`;
3. `-h^-1 log det(I-A)=sum_(nonempty B: product B=1) (|B|-1)! x_B`.

**Proof.** The summands commute and square to zero, so the matrix exponential factors, and `exp(x_iR_i)=I+x_iR_i`. Expanding the product proves (1), and the regular trace proves (2). The finite identity `-log det(I-A)=sum_(k=1)^n tr(A^k)/k` proves (3): for a fixed `k`-set `B`, exactly `k!` orderings survive in `A^k`; abelianity gives the same translation for each. Divide by `k` and use the regular trace. The trace-log identity can also be proved without an analytic logarithm by differentiating `det(I-tA)`, using Jacobi's formula, and integrating coefficientwise from `t=0` to `1`. All divisions occur over `Q`. QED.

The factorials have distinct meanings: division by `k!` forgets **every order**; division by `k` forgets only the **cyclic starting point**. There are `(k-1)!` cyclic orders of distinct labels. This remains true with repeated residues because the labels, not their residues, define the free cyclic action.

The Cayley graph interpretation is exact: SQ1(3) counts closed walks with pairwise distinct generator labels, including possible repeated group vertices. It does not assert that those walks are simple cycles.

### Recovering the paper's Bernoulli expression

The ordinary multilinear matrix polynomial `product_i(I+t_i R_i)` has the same squarefree coefficients as SQ1(1), and can be evaluated normally. For common `rho<1`,

`P(product selected g_i=1) = (1-rho)^n h^-1 tr product_i(I+(rho/(1-rho))R_i)`.

This is the paper's independence formula in the regular representation. Its character diagonalization gives the paper's formula in Lemma 5.1. Nothing here supplies its nontrivial-character energy bound.

**Important type check.** Substituting a nonzero number for `x_i` is not a ring homomorphism on `S_n`, since `x_i^2=0`. One must first use the unique multilinear representative, or the ordinary product polynomial, then evaluate. A numerical specialization cannot silently commute with calculations in the quotient.

### A scalar recurrence that keeps the necessary sidecar

In the group algebra put `E_k=sum_(|B|=k) R_(product B)` and `P_j=sum_i R_(g_i^j)`. Then

`k E_k=sum_(j=1)^k (-1)^(j-1) E_(k-j)P_j`, with `E_0=I`.

**Proof.** Differentiate `E(t)=product_i(I+tR_i)` and use `E'(t)=E(t)sum_(j>=1)(-1)^(j-1)P_j t^(j-1)` as a formal series. Compare coefficients. QED.

The desired number is `tr(E_k)/h`. Ordinary `tr((sum R_i)^k)` counts repeated labels; the corrections need the full power-class operators `P_j` and their mixed products, not just an unsupported replacement by an ordinary closed-walk moment.

## 2. The same algebra, with vertices marked instead

Let `M` be the adjacency matrix of any finite directed graph on `[n]` (weighted entries also work). Put `X=diag(x_1,...,x_n)`, now making **vertices** square-zero.

### Theorem SQ2: the square-zero log determinant inventories simple cycles

`-log det(I-XM) = sum_(nonempty B subset [n]) hc(M[B]) x_B` in `S_n`,

where `hc(M[B])` is the weighted number of directed Hamiltonian cycles on the induced vertex set `B`, counted up to cyclic rotation. Loops count as length-one cycles, if present.

**Proof.** Expand `sum_k tr((XM)^k)/k`. A monomial survives precisely when every visited vertex occurs once before closure. Its support is a `k`-set `B`; every directed cycle on that set has exactly `k` rooted descriptions, canceled by `1/k`. QED.

This is a literal common carrier for Carmichael subset counting and graph-cycle inventories. The two theorems mark different objects:

- SQ1 marks **prime labels**. It preserves distinct prime selection and total group product. It discards vertex simplicity of the induced Cayley walk.
- SQ2 marks **graph vertices**. It preserves vertex simplicity. It does not know which distinct prime labels generated an edge.

To demand both predicates, use independent square-zero variables for both sets of objects. The coefficient carrier then keeps both exclusions, at the cost of a larger state space. Ordinary spectrum is the specialization where the individual vertex marks have been forgotten; it cannot in general reconstruct SQ2.

### Why a naive exterior-algebra replacement fails

If anticommuting exterior generators `theta_i` replace the commuting `x_i`, then for commuting translations `A=sum theta_iR_i` satisfies `A^2=0`: terms cancel in pairs. Thus `exp A=I+A`, losing every subset of size at least two. The correct carrier is commutative square-zero, not a single exterior degree-one family. A realization by even products `x_i=theta_i eta_i` in an exterior algebra on disjoint generator pairs restores commutativity and square-zero behavior.

## 3. Minimal zero-sum subsets and binary matroid circuits

Call a nonempty product-one subset `B` **minimal** if no proper nonempty subset of `B` has product one.

### Theorem SQ3: the all-order simple-cycle characterization

For a nonempty product-one subset `B`, the following are equivalent:

1. `B` is minimal;
2. every ordering of its distinct labels yields a Cayley closed walk whose proper partial products are pairwise distinct.

**Proof.** Equality of two proper partial products identifies a proper consecutive subblock with product one; because `H` is abelian, this is a proper product-one subset. Conversely, given a proper product-one subset `C`, order the labels of `C` first; that proper partial product equals the initial identity. QED.

The universal quantifier is essential. In `C_5`, the labels `(1,2,3,4)` sum to zero and contain `(1,4)` as a proper zero-sum subset. Nevertheless the order `(1,2,4,3)` has partial sums `0,1,3,2,0`, a simple cycle. Thus existence of one simple order does not characterize minimality.

If `H=(F_2)^d`, form the binary matrix `V` with columns `g_i`. The product-one subsets identify exactly with the code `C=ker(V)` via incidence vectors. Minimal product-one supports are the circuits of the represented binary matroid, including singleton zero-column circuits. This is an exact map, not a general statement for arbitrary groups: for a group of exponent greater than two, binary selection is not automatically a linear subspace over the group's coordinate field.

### Theorem SQ4: pruning is a dual-code distance condition in exponent two

Let `D=rowspace(V)=C^perp`. For independent Bernoulli inclusion with parameter `rho`,

`P(Vb=0) = |D|^-1 sum_(d in D) (1-2rho)^wt(d)`.

Equivalently, if `W_C(u,v)=sum_(c in C)u^(n-wt(c))v^wt(c)`, then

`W_C(u,v)=|D|^-1 W_D(u+v,u-v)`.

**Proof.** The character `chi_a(g)=(-1)^(a dot g)` has evaluation word `a^T V`. Orthogonality on the image of `V` gives the first identity; the same expansion with independent formal choices `u,v` gives the second. Using rowspace rather than all ambient characters removes repeated evaluation words when `V` is rank-deficient. QED.

Under the paper's pruning condition, every nonzero evaluation word in `D` has weight at least `M`. Moreover its character energy is exactly `E_chi=2 wt(d)`. Thus in exponent two the support condition already supplies phase separation; the paper's Fejer step addresses the genuinely higher-order phenomenon of many nonidentity phases lying arbitrarily close to `1`. The code dictionary explains why a Hamming-support argument alone cannot replace that step for arbitrary `H`.

## 4. A hostile family: the logarithm of the return polynomial does not count atoms

SQ1(3) is a logarithm of a **determinant**, not a logarithm of `Z_H`. These must not be confused.

### Theorem SQ5: negative connected coefficients from nonunique zero-sum decompositions

Take an element `g` of exact order `d>=2` and `2d` distinct labels, all representing `g`. Then the coefficient of `x_[2d]` in `log Z_H(x)` is

`1 - (1/2) binom(2d,d) < 0`.

**Proof.** The return subsets have sizes `0,d,2d`, so `Z_H=1+e_d(x)+x_[2d]`. In the square-zero quotient, only `Y-Y^2/2` contributes to top degree in `log(1+Y)`, and the coefficient of `x_[2d]` in `e_d^2` is `binom(2d,d)`. QED.

For four labels of a nontrivial element of `C_2`, there are six minimal returning pairs and three different decompositions of the full set into two such pairs; the top log coefficient is `-2`. For orders `d=2,3,4,5`, the coefficients are `-2,-9,-34,-125`.

Every zero-sum subset can be greedily partitioned into minimal zero-sum blocks, but the partition need not be unique. Therefore an exponential formula treating those blocks as freely and uniquely assembling atoms is false. The exact survivor is SQ1's cyclic-order weighting or SQ2's vertex-simple cycle inventory, with their stated carriers. This is the concrete obstruction to importing the tournament primitive-orbit Euler product as a product over minimal Carmichael subsets.

## 5. Connection and loss ledger

| Source | Target/map | Preserved predicate | Lost information / required sidecar | Decisive test |
|---|---|---|---|---|
| Paper (3.1), Lemma 5.1 | regular translations with square-zero label variables | distinct-label subset and product one | numerical prime sizes, primality and Korselt divisibility remain external | duplicate residues, distinct labels |
| THM-4663 Cayley traces | replace `sum R_g` by `sum x_i R_g` | zero total step plus exclusion of repeated labels | vertex simplicity still absent | all labels equal in `C_2` |
| THM-1926 log determinant | vertex-square-zero `XM` | simple directed cycles by support | ordinary spectrum forgets marks | masks 73/83 below |
| Minimal product-one set | every ordering of partial products | all-order vertex-simple closure | one simple ordering is too weak | `C_5:(1,2,4,3)` |
| Exponent-two group | kernel code / dual evaluation code | zero sum, minimal support, character energy | fails as a binary linear model in higher exponent | nontrivial high-order roots near 1 |
| `log Z_H` | candidate atom inventory | formal cumulants only | nonunique block decompositions introduce signs | coefficient `-2` on four equal `C_2` labels |

## 6. A canon correction exposed by the bridge

**PROVED repair of THM-1870's general extrapolation.** For every loopless directed graph with no digons,

`c_k=tr(M^k)/k` for `k=3,4,5`.

If a closed walk of length at most five repeats a vertex before its final closure, the two resulting positive-length closed subwalks each have length at least three, contradiction. Thus these counts really are single traces; the historical assertion that `c4,c5` are “NOT single traces” is false.

There are six-vertex tournaments `T_73,T_83` with the same characteristic polynomial

`x^6 - 4x^3 - 3x^2 - 2x`,

the same trace vector `(0,0,12,12,10,48)` and the same `(c3,c4,c5)=(4,3,2)`, but `c6=0,1`, respectively. Encode a tournament by lexicographic pairs `(u,v), u<v`, with bit `1` meaning `u -> v` and bit `0` meaning `v -> u`. Full matrices appear in `verification.json`.

Adjoin one new universal source. Its adjacency matrix is block triangular, so its characteristic polynomial is multiplied by `x`. The new vertex lies on no directed cycle, so every old `c_k` remains unchanged. Hence at order seven the two tournaments are co-spectral but have different **sub-Hamiltonian** `c6`. Repeating the source extension proves that `c6` is nonspectral at every ambient order `n>=6`. The claimed boundary at the ambient Hamiltonian length is therefore false. The exact positive boundary is the absolute length range `3,4,5`; we make no claim that every `c_k` with `k>=7` is nonspectral.

The cycle-cover mechanism is visible already at order six: the two examples have respectively zero and one pair of vertex-disjoint directed triangles. The constant coefficient of the characteristic polynomial contains `-c6 + disj_33`, so those two changes cancel. The marked square-zero logarithm restores the missing cycle coordinate.

No computational complexity separation follows from this witness. It proves loss of information under the characteristic-polynomial map, not a `#P` classification or hardness assertion.

## 7. Exact validation and scope

Run `python3 verify_squarefree_bridges.py --output verification.json`.

The script uses integer and rational arithmetic only. Every gate raises an explicit exception on failure, so Python optimization cannot disable it; normal and `python3 -O` outputs were checked byte-identical. It checks independent constructions:

- every sequence over `C_m`, `1<=m<=5`, with length `0<=n<=4`, including zero labels and repeated residues: subset enumeration versus unused-label walk DP, group-algebra Newton recursion, and the all-order simplicity criterion;
- all labelled tournaments on one through four vertices, plus the two six-vertex witnesses and source extensions: principal-minor determinant plus square-zero logarithm versus direct anchored cycle permutations;
- every column sequence in `(F_2)^d`, `1<=d<=3`, `0<=n<=4`: kernel enumeration versus dual-code formula at exact biases `1/5,1/2,3/4`;
- the negative-log family through `d=5`.

Independent root review separately recomputes the tournament characteristic polynomials using SymPy and its own anchored-permutation enumerator. Those files accompany the packet in `review/`.

The finite tests support, but do not replace, the general proofs. We do not claim an improved Carmichael exponent, a new analytic prime estimate, a simple-cycle asymptotic, a general spectral recovery theorem, or a tournament relation on primes. The bridge's output is a faithful finite coefficient carrier and an explicit correction to an inherited spectral extrapolation.
