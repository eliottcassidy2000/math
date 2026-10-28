# Weighted divisor boxes, three summands, and square companions

Date: 2026-10-28. **PROVED elementary identities and obstructions; FINITE-EXACT controls.** No novelty claim, no global classification of aliquot orbits, and no asserted arithmetic–tournament isomorphism. The user's clarification is **exactly three positive proper-divisor terms**, including 1 and excluding the source.

## Result first

The unifying object is the **prime-labelled, weighted exponent box**, not the sum alone. A prime cube and a distinct semiprime each have three proper divisors, but their divisor posets differ. Coprime multiplication works exactly before weights/source identity are forgotten. The universal square companion `C(n)=n(n−1)` has a similarly exact coprime divisor formula, but does not commute with aliquot summation.

For `t=s(n), u=s(n−1)` and `n>=2`:

```
n+C(n)=n²,
s(C(n))=(n−1)t+nu+tu,
s(C(n))−C(t)=t(n−t)+(n+t)u.                         (1)
```

Thus any `n>=2` satisfying `s(C(n))=C(s(n))` must be abundant. It cannot have `n−1` prime. These are necessary exclusions, not a proof that no such n exists.

## 1. Inheritance and the exact distinction

Source revision: `6192c0df22ad13989070f4c6e8d3eb0ad7000a3f`.

Closest recovered mechanisms:

* [Weighted divisor downsets and tournament cycle complexes](../aliquot_power_124_20261026/tournaments/TOURNAMENTS_124.md), §2: the weighted proper-divisor multicomplex and its exact graph boundary on **active singleton coordinates**, `Ω(n)<=2`.
* The same packet's [connection board](../aliquot_power_124_20261026/CONNECTION_BOARD.md): cube versus simplex-boundary versus order-complex distinctions; prime weights and exponent multiplicity are essential.
* [Arithmetic braids: divisor balance](../arithmetic_braids_20260917_divisors.md), DB1–DB3: exponent boxes, squarefree subcubes and the inadequacy of support size alone.
* [Ternary composition memory](../session77_ternary_color_20261027/TERNARY_MEMORY.md): complete path-cover profiles, not scalar Hamiltonian counts, compose through a cyclic ternary quotient. The underlying canonical theorem is [THM-3121, path-cover walk-content substitution kernel](../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md).
* Session82's [decimal triple](../session82_decimal_aliquot_square_sum_20261028/DECIMAL_TRIPLE.md) and [fiber ports](../session82_decimal_aliquot_square_sum_20261028/FIBER_PORTS.md): coprime aliquot transport, the square-companion anchor, and source-resolved relation composition.

Canonical hostile: `s(27)=s(35)=13`, with different divisor shapes and different future responses after multiplication by 2. Corrected near miss: geometric resemblance or equal cardinality is not a predicate-preserving map. Least-used sidecar: the actual prime labels and the actual source integer. Live board: exponent boxes; weighted evaluation; factor grouping; positional carries; named-source composition; square companions.

## 2. The weighted box and the two three-term shapes

Write `n=∏ p_i^{a_i}`, where the primes are distinct and `a_i>=1`. Let `L(n)` contain all positive divisors, ordered by divisibility. The map

```
(e_1,...,e_r) ↦ ∏ p_i^{e_i},       0<=e_i<=a_i
```

is an order isomorphism from the product of chains `∏[0,a_i]` onto `L(n)`. Coordinatewise minimum/maximum become gcd/lcm. Deleting its top point leaves the proper-divisor downset. Its enumerator is

```
A_n(X)=∏_i(1+X_i+...+X_i^{a_i})−∏_i X_i^{a_i},
A_n(p_1,...,p_r)=s(n),             A_n(1,...,1)=τ(n)−1.   (2)
```

The last two evaluations are different observables. Three proper-divisor terms mean `∏(a_i+1)=4`, hence precisely:

| source | proper divisors | proper-divisor poset | rank enumerator |
|---|---|---|---|
| `p³`, p prime | `1,p,p²` | a three-point chain | `1+z+z²` |
| `pq`, distinct primes | `1,p,q` | one bottom, two incomparable atoms | `1+2z` |

Even the full four-element lattices are nonisomorphic: a chain has one atom, a Boolean square has two. “Three terms” is a cardinality classification, not a cubic polynomial or a common order type. A distinct-three-prime source has **seven**, not three, proper-divisor terms:

```
s(pqr)=1+p+q+r+pq+pr+qr.
```

For 1001 the source is the Boolean three-cube, with rank counts `1,3,3,1`; proper divisors omit its top. Its successor 343 happens to be the integer `7³`, whose full divisor lattice is a four-chain. There is no lattice isomorphism behind this transition.

### Geometry and the graph boundary

For squarefree n, exponents are subsets. Proper subsets are the abstract simplex boundary **including its empty face**; this is not the geometric boundary of a cube. The order complex of proper divisors including 1 is a cone. Removing both 1 and n before taking the order complex yields the barycentric subdivision of the simplex boundary (the usual small empty-boundary conventions apply).

The inherited literal graph result uses **active singleton coordinates**: retain coordinate i only if its degree-one monomial `X_i` is actually a proper-divisor monomial, and discard absent singleton coordinates. On this convention, the support is a graph independence complex iff `Ω(n)<=2` (n>=2). The boundary cases are explicit: n=p has support `{1}` and is realized by the graph on **zero** vertices; n=p² gives one isolated vertex; n=pq gives K2. If one instead insists on a graph having every prime coordinate as a vertex, the prime case fails and the criterion is `Ω(n)=2`. At `p³`, the monomial `X²` is not squarefree. At `pqr`, every pair is permitted but the triple is forbidden, impossible for a graph whose only conflicts are pairs. An arbitrary larger encoding or hypergraph is not excluded. Neither scalar evaluation nor the number three supplies a tournament realization.

## 3. Coprime multiplication is the exact composition law

For positive integers a,b, the multiplication map

```
L(a) × L(b) → L(ab),       (d,e) ↦ de
```

is surjective, and is bijective iff `gcd(a,b)=1`. In the coprime case it is a lattice isomorphism; unique prime coordinates give its inverse. When a prime p divides both a and b, `(p,1)` and `(1,p)` collide. The product sum therefore has the exact multiplicity interpretation

```
σ(a)σ(b)=Σ_{d|a,e|b} de=Σ_{f|ab} m_f f,
```

where `m_f` is the number of factorizations of f with `d|a,e|b`. At a prime with exponents α,β, the multiplicity at exponent k is

```
max(0, min(α,k)−max(0,k−β)+1).
```

The local multiplicities multiply over primes. Hence `σ(a)σ(b)>=σ(ab)`, with equality iff a,b are coprime. This proof explains the guard rather than merely citing multiplicativity.

The source-retaining observer `Ψ(n)=(n,σ(n))` obeys coordinatewise multiplication on coprime inputs. Equivalently, pairs `(n,s(n))` compose by

```
(a,u) ⋆ (b,v) = (ab, av+bu+uv).                          (3)
```

This operation is associative and commutative: `(a,u)↦(a,a+u)` conjugates it to coordinatewise multiplication. Realizing it as the observer of a product of integers still needs the coprimality guard (pairwise coprimality for repeated composition).

The scalar s alone does not support the operation: `s(27)=s(35)=13`, but

```
s(54)=66,            s(70)=74.
```

For odd n, (3) gives `s(2n)=n+3s(n)`, so the lost source coordinate is exactly the missing term. This is the same *information-loss mechanism* as the ternary path-profile and shared-port examples, not an isomorphism between their native objects.

## 4. Grouped factors and carries must remain distinguished

For any integer x>1 and k>=1, the divisors `1,x,...,x^k` lie in `L(x^k)`. They exhaust it iff x is prime: if x is composite, any prime divisor p<x is an additional divisor not equal to a power of x. Therefore

```
σ(x^k)>=1+x+...+x^k,
s(x^k)>=1+x+...+x^(k−1),
```

with equality iff x is prime. A composite-base “prime-power formula” drops internal prime coordinates. Example: `1+4+16=21`, but `s(64)=63`.

The output becoming a power is a separate additive phenomenon. For 1001, the seven proper divisors are

```
1,7,11,13,77,91,143.
```

Their base-7 digits, added before carrying, give the column polynomial

```
D(X)=14+19X+4X².
```

Indeed the successive unnormalized columns are `(14,19,4)`. Carrying gives `14→0` with carry 2; `19+2=21→0` with carry 3; `4+3=7→0` with carry 1. Thus the sum is `(1000)_7=343`. Algebraically,

```
D(X)−X³=(7−X)(X+1)(X+2),       so D(7)=7³.             (4)
```

Evaluation at a base b identifies polynomials differing by a multiple of `X−b`; carrying is the local relation `bX^j=X^(j+1)`. This numerical quotient forgets the seven source divisors and their incidence. Retaining a carry certificate explains the equality; it does not restore a false cube-to-chain isomorphism. The later edge `343+57=400` separately tests `σ(343)` being square.

## 5. Square companions and the exact commutation defects

Define the polynomial `C(x)=x(x−1)` for nonnegative integers x, but apply s only to positive arguments. For n>=2, C(n)>0 and `gcd(n,n−1)=1`, so

```
σ(C(n))=σ(n)σ(n−1),
s(C(n))=(n−1)s(n)+n s(n−1)+s(n)s(n−1).                (5)
```

For n>=3 the distinct positive labels n,C(n) always form a square-sum edge. At n=2 they coincide, so this is not an edge in the simple graph. The anchor is n=1001 and C(n)=1001000; it is the family already present in Session82 with its parameter shifted by one.

### Aliquot commutation

Writing t=s(n), u=s(n−1), subtract C(t) from (5) to obtain (1). If t<=n, its right side is positive: for n>=3, u>=1; for n=2 it equals 1. Thus commutation requires t>n, i.e. n abundant. If n−1 is prime, u=1 and equality would imply

```
n=t(t−1)/(t+1)=t−2+2/(t+1).
```

For t>=2 integrality would require `t+1|2`, impossible; t=1 gives n=0. This excludes the entire n−1-prime case. No classification of the abundant composite-predecessor residue is claimed.

### Sigma commutation

A different exact equation is

```
σ(C(n))=C(σ(n))  iff  σ(n−1)=σ(n)−1
                  iff  s(n−1)=s(n).                   (6)
```

Cancellation of positive σ(n) proves the first equivalence. This can hold: at n=3, `σ(6)=12=C(4)`. It does not give aliquot commutation: `s(6)=6≠C(1)=0`. The two functions σ and s must not be interchanged.

### Images of square-sum edges

The image pair has total

```
s(n)+s(C(n))=σ(n)(1+σ(n−1))−n².                       (7)
```

Its being square is the exact extra condition; positivity and distinctness of the image labels must also be checked. The smallest companion edge already fails: `3+6=9`, but `s(3)+s(6)=7`. The requested anchor fails as well: `343+2143960=2144303≡3 (mod4)`.

Failure is not universal destruction. At n=264,

```
C(n)=69432,      s(n)=456,      s(C(n))=120648,
456+120648=121104=348².
```

Here `264=2³·3·11`, n−1=263 is prime, and `σ(264)=720`, so all values follow directly from (5). This positive control concerns preservation of one edge; it is not commutation with C. The checker finds this to be the only surviving companion-edge image in the finite range `3<=n<=10000`, not at all heights.

## 6. Typed connection ledger

| Source → target | Exact map / predicate preserved | Loss / necessary retained data |
|---|---|---|
| prime exponent box → divisor lattice | exponent evaluation; divisibility, gcd, lcm | none if prime labels retained |
| weighted divisor polynomial → integer | evaluate at its primes; exact s | loses shape, weights and source if only value retained |
| coprime source pair → product | lattice product; σ multiplies; (3) | shared primes create multiplicities; retain support guard |
| composite base x → powers of x | selected chain inside full divisor box | missing internal divisors; prime factor coordinates needed |
| base-b digit polynomial → sum | positional evaluation; (4) | carries merge coefficient/source information |
| n → (n,C(n)) | exact square-sum edge for n>=3 | s adds a separate arithmetic constraint (7) |
| aliquot-fiber graph → fiber quotient | individual existential edges | common endpoint source is lost; relation composition repairs it |
| tournament factors → substituted tournament | inherited walk/path-cover kernel | scalar H loses component profiles and ports |

Only the first three rows are literal arithmetic operation maps. The last two share the principle “retain the coordinate needed by the next operation”; no arithmetic–tournament conjugacy is supplied.

## 7. Verification

Run from this directory:

```
python3 check_lattices.py > results.json
python3 -O check_lattices.py > results.optimized.json
cmp results.json results.optimized.json
```

Standard library only, explicit checks remain active under `-O`. The checker independently compares divisor-box evaluation with divisor-pair enumeration on 1..2000; checks all gcd/lcm pairs inside those boxes and the active-singleton graph boundary, including the prime case; compares every divisor-product multiplicity for 1<=a,b<=80; checks grouped-base inequalities for 2<=x<=100,1<=k<=4; computes the anchor carries; checks companion identities against independently factored products for 2<=n<=10000 and against direct divisor-pair sums through n=1000; and exhausts image-edge preservation only in the explicitly stated finite interval. The proofs above, not these bounds, justify universal statements.
