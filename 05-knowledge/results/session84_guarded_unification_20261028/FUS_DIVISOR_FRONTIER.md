# F = U + S: one divisor carrier, an exact balance frontier, and the square-sigma intersection

Date: 2026-10-28. Source revision: `ae7c073e864d` (read-only source worktree).

**Status:** inherited definitions and classification; **PROVED** elementary update laws, monotonicity, square-source promotion obstruction, and prime-guarded parametrization below; **FINITE-EXACT** controls through 10000; the stronger complete prime-step law and graph components are **COROLLARIES OF THE INHERITED CITED LJUNGGREN THEOREM**. No literature-priority claim, no infinitude claim for quadratic primes, no divisor/Berggren tree isomorphism.

## Results first

1. The original `F,U,S` are **counts**, not divisor sums. For `N>=2`, they count, respectively, all proper nontrivial divisors, the prime ones, and the squarefree ones. The prime set is contained in the squarefree set. The exact balanced shapes are `p`, `p³`, and `p²qr` (distinct prime letters).
2. On **nonsquarefree integers ordered by divisibility**, the defect `δ=F−U−S` is monotone. Its only nontrivial flat step is `p² | p²q`. Every proper multiple of a balanced composite has positive defect. Thus balance is a crossing frontier, not a conserved quantity.
3. If `σ(N)` is square and a prime `p` occurs in `N` to exponent exactly one, **`σ(pN)` is not square**. In particular, all three one-prime promotions of `66=2·3·11` become balanced `p²qr` shapes but lose square sigma.
4. The two predicates are **not globally incompatible**: `364=2²·7·13` is the least mixed balanced square-sigma source, with `σ(364)=784=28²`; `820=2²·5·41` has `σ=1764=42²`. Their minimality/existence statements are independently checked by direct divisors.
5. The genuine common object with Session83 is a **labelled prime-exponent box, retained before evaluation**. Cardinality gives `F,U,S`; arithmetic weight gives `σ` and the aliquot map. Both evaluations must keep their own predicates. This also supplies a precise comparison with Berggren composition, but not an identification of its noncommuting tree with a commuting divisor lattice.

## 1. Inheritance and the notation repair

Primary recovered routes:

- `05-knowledge/results/arithmetic_braids_20260917_divisors.md`, DB1–DB3: exact definitions, exponent-box mechanism, complete balance classification, defect strata.
- `05-knowledge/results/predecessor_prime_shapes_20261014b.md`, §1: the corrected `p⁴qr` claim and explicit change to fourth-power-free divisors.
- `05-knowledge/results/recursion_carry_bridges_20261014i.md`, §1: `10=7+3`, not a recurrence or universal conservation law.
- `05-knowledge/results/trijection_reset_20261020/arithmetic.md`, §§1–3: several different `3–4–3` partitions and the complement/prime-exchange action.
- `05-knowledge/results/supersingular_divisor_shapes_qseries_20261022.md`, §3: the fresh-prime recursion; explicit distinction from the historical Chowla sum `σ(N)−N−1`.
- `05-knowledge/results/session83_aliquot_divisor_grammar_20261028/CARRIER_AND_144_BRIDGE.md`: coprime weighted composition and square-class signatures.

The relevant correction `MISTAKE-666` does not retract the divisor classification. The supersingular note explicitly says the Chowla reading is a **different object**, not a substitute definition of `F,U,S`. The older `p⁴qr` extension has `(F,S,U)=(18,7,3)`, defect eight; `18=15+3` is true only after replacing squarefree by **fourth-power-free**.

To prevent the Session83 reuse of the letter `S` from hiding this distinction, write **`Σ(N)=σ(N)`** for the divisor sum here and reserve `S` for the squarefree-divisor count.

Let `N=∏ p_i^{a_i}`, `r=ω(N)`, and `τ=∏(a_i+1)`. Then

\[
F=\tau-2,\quad U=r-[N\text{ prime}],\quad
S=2^r-1-[N\text{ squarefree}],
\]
\[
\delta=\tau-2^r-r-1+[N\text{ squarefree}]+[N\text{ prime}].
\]

The number of **proper divisor terms including 1**, Session83's `T`, is `T=τ−1=F+1`. Therefore balanced shapes have `T=1,3,11`, respectively. These values are not interchangeable with a list of prime labels `{2,3,11}` or with the values of `F` (`0,2,10`).

The `T=11` layer has the four exponent profiles `(11),(5,1),(3,2),(2,1,1)`. Only `(2,1,1)` is balanced. Likewise `T=3` contains both `p³` and `pq`, but only `p³` is balanced. A term count is a coarse observation, not a replacement for the exponent profile.

## 2. Exact operations and the all-height frontier

### Existing-prime increase

If `p|N`, `a=v_p(N)`, then

\[
\delta(pN)-\delta(N)
=\frac{\tau(N)}{a+1}-[N\text{ squarefree}]-[N\text{ prime}].\tag{1}
\]

**Proof.** Increasing `a` by one increases `τ` by `τ/(a+1)`, leaves `r` fixed, and makes both squarefree and prime indicators at the new source zero. Subtract the displayed defect formula.

Consequences:

- `p→p²` has defect `0→−1`.
- For squarefree composite rank `r`, one promotion changes defect by `2^(r−1)−1`. It therefore changes `−r` to `2^(r−1)−r−1`, which is zero **exactly at `r=3`**.
- At every nonsquarefree source, an existing-prime increase is strictly defect-increasing.

### Fresh-prime adjunction

For prime `q∤N`, the inherited exact law is

\[
\delta(Nq)=2\delta(N)+U(N)-[N\text{ squarefree}]-[N\text{ prime}].\tag{2}
\]

At a nonsquarefree source it gives the particularly transparent increment

\[
\delta(Nq)-\delta(N)=\tau(N)-2^r-1.
\]

Since nonsquarefreeness forces `τ>=3·2^(r−1)`, this increment is at least `2^(r−1)−1`. It vanishes only when `r=1` and `τ=3`, namely `N=p²`.

### Monotonicity theorem

**For nonsquarefree `a|b` with `b>a`, `δ(b)>=δ(a)`, with equality if and only if `a=p²` and `b=p²q` for distinct primes `p,q`.**

Factor `b/a` into primes and apply (1)–(2) successively. Every step stays nonsquarefree and is nondecreasing. A flat step must start from a prime square and append a fresh prime. Its target has rank two, so no subsequent step is flat. This proves both directions and the exact equality boundary.

In particular, the balanced composites `{p³,p²qr}` form a **divisibility antichain**. Every proper multiple of such a source has strictly positive defect. The exclusion of prime sources is essential: a prime divides many balanced composites. The exclusion of squarefree source integers is also essential: `δ(p)=0>δ(p²)=−1`.

This recovers the two elementary approaches to balance:

- `p²→p³` raises the exponent (`−1→0`);
- `p²→p²q→p²qr` adjoins fresh primes (`−1→−1→0`).

Neither continuation can preserve balance once it has reached a balanced composite.

### The strict-positive region has four minimal divisibility generators

For every `N>=2`,

\[
\delta(N)>0 \iff N\text{ is divisible by an instance of }
 p^4,\ p^3q,\ p^2q^2,\ p^2qrs,
\]

with distinct prime letters in each pattern. These four profiles have defects `1,1,2,3` and are exactly the minimal elements of the strict-positive region.

For rank one, positivity first occurs at exponent four. For rank two it first occurs at `(3,1)` or `(2,2)`. A rank-three profile beyond the balanced `(2,1,1)` contains one of these rank-two patterns. At rank at least four, any repeated prime gives the pattern `(2,1,1,1)`. Squarefree sources never have positive defect. These observations prove necessity; monotonicity proves sufficiency. They also verify minimality by inspecting each pattern's proper divisors.

**Balance is not a separator that every operation path must visit.** The valid path `4→12→24` has defects `−1,−1,+1`, jumping over zero. Thus “frontier” here means a thin equality locus, not an intermediate-value theorem on this discrete poset.

## 3. The single-exponent square obstruction

**Theorem.** Let `N` be a positive integer, `p` prime, `v_p(N)=1`. If `Σ(N)` is a square, then `Σ(pN)` is not a square.

Write `N=pm`, `gcd(p,m)=1`. Multiplicativity gives

\[
\frac{\Sigma(pN)}{\Sigma(N)}=\frac{p^2+p+1}{p+1}.\tag{3}
\]

The fraction is reduced because `p²+p+1=p(p+1)+1`. Its numerator lies strictly between `p²` and `(p+1)²`, so is not an integer square. A reduced rational square has square numerator and square denominator. Thus (3) is not a rational square. The quotient of two nonzero integer squares would be a rational square, proving the claim.

The exponent-one hypothesis is required for this elementary formula and proof. No claim that the conclusion fails at other exponents is intended. This is a theorem about square divisor sums, not preservation of the aliquot endpoint.

### The only square-preserving prime step (explicit inherited dependency)

The elementary theorem above has a stronger all-exponent extension, using exactly the **Ljunggren prime-power square-sigma classification already imported in Session82**, `05-knowledge/results/session82_decimal_aliquot_square_sum_20261028/SQUARE_SIGMA_EDGES.md`, §§1–3. Its only prime-power sources are `3¹`, `3⁴`, and `7³`. We inherit that cited classical theorem; the present finite checks do not prove it.

**For every positive integer `N` with square `Σ(N)` and every prime `p`,**

\[
\Sigma(pN)\text{ is square}\iff p=3\text{ and }3\nmid N.\tag{3a}
\]

**Proof.** If `p∤N`, multiplicativity makes the quotient `Σ(pN)/Σ(N)=p+1`. It is square exactly when `p=(k−1)(k+1)` is prime, forcing `k=2,p=3`.

If `p|N`, let `a=v_p(N)>=1` and `R_a=1+p+⋯+p^a`. Then `R_(a+1)=pR_a+1`, so consecutive sums are coprime. The ratio `Σ(pN)/Σ(N)=R_(a+1)/R_a` can be a rational square only if both adjacent prime-power sums are integer squares. The inherited complete list has no adjacent exponents at any fixed prime. Therefore no existing-prime increase preserves square sigma. This proves (3a), including `N=1`.

### Exact components of the prime-multiplication graph

Define the directed graph with vertices **all positive integers `N` for which `Σ(N)` is square**, and edges `N→pN` for prime `p` when both endpoints are vertices. This graph is **not the aliquot graph**.

By (3a), every edge is `N→3N` with `3∤N`. Consequently its undirected components, with the indicated edge orientation retained, are exactly:

- the pairs **`{N,3N}`**, where `3∤N` and `Σ(N)` is square;
- isolated square-sigma vertices with **`v₃(N)>=2`**.

Indeed, if a square-sigma vertex has valuation `v₃(N)=1`, write `N=3m` with `3∤m`. Multiplicativity gives `Σ(N)=4Σ(m)`, so `Σ(m)` is square and the vertex belongs to the pair `{m,3m}`. A valuation-zero vertex has its one outgoing edge, no incoming edge, and no other edge. A valuation-at-least-two vertex has neither, since any allowable edge has endpoint valuations zero and one.

The unit boundary is **`1→3`**, with totals `1` and `4`. If vertex1 is removed by restricting to `N>1`, vertex3 becomes isolated; no other component changes. Controls:

| Component | Divisor sums | Reason |
|---|---|---|
|`{81}`|121|`v₃(81)=4`, isolated|
|`343→1029`|400→1600|fresh3|
|`364→1092`|784→3136|fresh3|

These are edges in the **prime-multiplication graph on square-sigma sources**. They do not say that the aliquot orbit of any vertex stays inside its component. The component theorem shares the explicit Ljunggren dependency of (3a).

### The exact `{2,3,11}` test

`Σ(66)=3·4·12=144=12²`, but `(F,U,S)=(6,3,6)`, defect `−3`. Its three balanced promotions are:

| Marked prime promoted | Source | `(F,U,S)` | `Σ` | Squarefree kernel of `Σ` |
|---:|---:|---:|---:|---:|
|2|132|`(10,3,7)`|336|21|
|3|198|`(10,3,7)`|468|13|
|11|726|`(10,3,7)`|1596|399|

The stronger finite support statement is also exact: among all eight sources `2^a3^b11^c`, `a,b,c∈{1,2}`, **only 66 has square sigma**. Their square-class signatures are

| Promoted subset | none | 2 | 3 | 11 | 2,3 | 2,11 | 3,11 | 2,3,11 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
|Kernel|1|21|13|399|273|19|5187|247|

This can be proved without a search: the promotion signatures `21=3·7`, `13`, and `399=3·7·19` are linearly independent over `F₂` (the primes `13`, then `19`, then `3` eliminate the three coefficients). It does **not** classify larger exponents on this support.

## 4. Actual balanced square-sigma sources: a full fixed-pair parametrization

For the mixed balanced shape `N=p²qr`, fix distinct primes `p,q`, and put

\[
A=(p^2+p+1)(q+1)=d h^2,\qquad d=\operatorname{sf}(A).
\]

Here `d` is squarefree and `h` is a positive integer. Then

\[
\Sigma(p^2qr)\text{ square}
\iff r+1=d t^2\quad(t\in\mathbb Z_{>0}),\tag{4}
\]

subject throughout to **`r=d t²−1` being prime and distinct from `p,q`**. The square total is `(dht)²`.

**Proof.** `Σ=A(r+1)`. Equality of prime-exponent parities says `sf(r+1)=d`, equivalent to `r+1=dt²`. This is an iff, not just a sufficient family.

Examples:

- `p=2,q=7`: `A=56=14·2²`, giving `N=28(14t²−1)` and `Σ=(28t)²` whenever the third factor is an allowed prime. `t=1` gives **364**.
- `p=2,q=5`: `A=42`, giving `N=20(42t²−1)` and `Σ=(42t)²`; `t=1` gives **820**.
- `p=2,q=13`: `A=98=2·7²`, giving `N=52(2t²−1)` and `Σ=(14t)²`; `t=2,3,4` give **364,884,1612**.

No infinitude of primes represented by these quadratics is asserted. Replacing the prime guard by an arbitrary positive integer changes the divisor profile and can destroy balance.

Direct divisor enumeration of every `2<=N<364` proves that 364 is the **least mixed** balanced square-sigma source. This is a complete finite minimum proof, not an asymptotic inference. The two smaller balanced square-sigma sources are `3` and `343`.

For the other two balanced branches:

- A prime `p` has square `Σ=p+1` iff `p=3`: factor `p=(k−1)(k+1)`.
- The cube branch contributes only `343=7³`, **inheriting** Session82's corollary of the cited Ljunggren square-repunit theorem (`SQUARE_SIGMA_EDGES.md`, §§1–3). This session neither reproves that deep theorem nor treats a finite census as its proof.

Thus the full balanced square-sigma intersection is `{3,343}` plus the prime-guarded mixed family (4), with the cubic uniqueness carrying that explicit inherited dependency. `81` is a useful hostile: its sigma is `121`, but its defect is `1`, so it is not balanced.

## 5. The common carrier and what genuinely transfers

For a labelled exponent vector `a=(a_i)`, retain

\[
P_a(X_1,\ldots,X_r)=\prod_i(1+X_i+\cdots+X_i^{a_i}),
\]

together with distinct prime labels `p_i`, the top monomial `∏X_i^{a_i}`, the unit, the atom set, and the Boolean squarefree subbox.

- `P_a(1,...,1)=τ(N)` counts all divisors.
- `P_a(p_1,...,p_r)=Σ(N)` weights each divisor by its arithmetic value.
- Removing unit/top and applying the atom/Boolean predicates yields exactly `F,U,S`.
- Disjoint prime support gives multiplication of these polynomials and hence the Session83 coprime product laws.
- Overlapping prime support requires merging exponent bounds. Multiplying the polynomials instead counts divisor decompositions with multiplicity, not divisors of the product. For example `(1+X)²=1+2X+X²` is not `1+X+X²`.

This is an exact **common source with different evaluations**. It is not a minimal coding claim: the integer itself already determines its factorization.

| Source → target | Map / preserved predicate | Lost information / required sidecar | Decisive hostile |
|---|---|---|---|
| Labelled box → `F,U,S` | Cardinalities of explicitly marked subsets | Prime magnitudes and their exponent ownership | 12 and18 have the same count triple, but multiplication by the named prime2 yields24 and36 with different `F` |
| Labelled box → `Σ` | Arithmetic evaluation | Distinct sources and all their term allocations | 66 and70 both evaluate to144 but have different aliquot targets78 and74 |
| Full arithmetic carrier → square class | Prime-exponent parities of `Σ` | Magnitude; actual root; endpoint | All square totals have zero signature, not equal endpoints |
| Squarefree rank-three source → one marked promotion | Existing-prime update; balance becomes true | Square-sigma is **not** preserved | 66→132,198,726 |
| Mixed balanced source → fixed-pair parameter | Equation (4) | Prime predicate does not follow from polynomial shape | `p=2,q=7,t=2` gives third factor55, composite |
| Berggren word → Euclid pair/triple | Ordered matrix product and quadratic evaluation | Ancestry/branch character after inappropriate quotienting | Corrected THM-4808 keeps the pole/word sidecars and supplies every-modulus character hostiles |

### Exact structural limit of the Berggren comparison

The divisor operation graph is not the Berggren tree: it has commuting diamonds, for example

`2→4→12` and `2→6→12`.

The classical Berggren tree gives each primitive triple one parent/address. Therefore no identification can preserve these two operation graphs with all parents retained. One may compare their **composition-and-evaluation diagrams**, not their entire dynamics.

The elementary observer criterion makes this precise. For a guarded operation system `T_a` on `X`, a quotient `q:X→Y` supports the same deterministic labelled operations iff (i) enabled operation labels are constant on every `q`-fiber and (ii) for every enabled label, `q(x)=q(x')` implies `q(T_a x)=q(T_a x')`. Necessity follows by evaluating the induced operation; sufficiency defines it by any representative. The 12/18 witness fails (ii), while the full labelled exponent carrier satisfies both conditions. The Berggren ancestry-to-residue quotient can fail the analogous requested observation. This is a typed shared mechanism, not a theorem that divisor counts determine triangle ancestry.

The inherited symmetry bridge is separately scoped: **THM-4733**, `01-canon/theorems/THM-4733-divisor-plane-of-p2qr-versus-the-marked-path-tiling-cube-3-4-3-and-4-2-4-reflections.md`, compares the divisor Klein action to a folded fixed-tile-free cube; it explicitly obstructs an equivariant injection into genuine full marked-tournament presentations. It cannot be used to infer arithmetic, coloring, or Berggren consequences.

## 6. Reproduction and stopping boundary

Run in this directory:

```sh
python3 check_fus.py > normal.json
python3 -O check_fus.py > optimized.json
cmp normal.json optimized.json
```

Both modes agree. The always-active checker reports **135821 checks**: 9999 independent direct-divisor observations, 24300 existing-prime updates, 6660 fresh-prime updates, 25882 nonsquarefree divisibility pairs, 691 square-source promotion tests, the complete eight-vertex `{2,3,11}` exponent cube, the exact minimum below364, the four strict-positive minimal divisibility profiles, the bounded full prime-step and component controls, and five rejected domain hostiles. The complete-prime-step finite control covers 4228 pairs, with square-sigma sources through10000 and primes at most47 plus every prime divisor of each source; component checks cover 270 source vertices including1. It also records every balanced square-sigma source through10000.

The computation is finite evidence. Universal claims are proved above; the cubic square uniqueness and full prime-step/component classification explicitly inherit the cited Ljunggren theorem. General prime-value infinitude, higher exponents on the fixed support, global aliquot behavior, and any divisor/Berggren/tournament equivalence are not established.
