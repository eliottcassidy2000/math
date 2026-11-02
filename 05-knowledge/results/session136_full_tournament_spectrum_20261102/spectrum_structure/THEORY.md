# Finite-bank compression and the minimum-order spectrum

**PROVED from inherited FINITE-EXACT witness banks and the CITED Moon–Busch
strong floor.** New integer inequalities and repacking rules are independently
replayed here; the large inherited graph bank is not regenerated. Universal
realization of every positive odd integer other than 7 and 21 remains **OPEN**.
No external novelty or priority claim is made.

## 1. What changes: first appearance, not induced deletion

Write `S={H(T)}` for the all-order ordinary Hamiltonian-path count spectrum.
Paths are directed vertex words, not automorphism orbits. Set `H(empty)=1`.
Define

```
mu(h)    = min {|T|: H(T)=h}, with mu(1)=0 and min(empty)=infinity;
sigma(h) = min {|T|: T nontrivially strong and H(T)=h}.
```

The use of an empty unit only changes `mu(1)`. The main theorem concerns h>1.

> **Finite-bank compression theorem.** For every **attained** h>1,
>
> `mu(h) <= ceil(22 log_{80267}(h))`.
>
> The coefficient of natural log is `22/ln(80267)=1.9480898075965785...`.
> The witness need not be induced in the original tournament and need not be
> strong. The theorem supplies no witness when h has not been shown attained.

This improves the general deletion-only bound `floor(3 log_3 h)`, whose
coefficient is about 2.73072. It does **not** improve the strong-only
Moon–Busch asymptotic coefficient `3/ln5=1.86400...`. It converts existing
finite positive coverage into a new all-order *conditional-on-attainment*
size reduction. It does not convert finite coverage into full coverage.

## 2. Inheritance and source audit

All source paths below were read at repository commit `7f6987db0cc3`
(the read-only source worktree HEAD). They were then compared to current
`origin/main=d11ba3658601`; every load-bearing source file was unchanged.
No fetch, reset, rebase or repository edit was performed by this child.

1. **SCC product and order join:**
   `01-canon/theorems/THM-1862-order-join-reduction-principle.md`.
   `H(T)=product H(C)` over SCCs; source/sink padding preserves H.
2. **All-order omissions and parity:**
   `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md` and
   `01-canon/theorems/THM-001-redei.md`.
3. **Strong bank through 609 at order at most eight:**
   `01-canon/theorems/THM-4810-h-spectrum-atlas-strong-fibres-exact-low-spectra-complement-parity-invariant-score-monotonicity-paley-budgets.md`, part A;
   source `05-knowledge/results/h_spectrum_atlas_complement_budgets_paley_clocks_20261027.md`, section 1.
   Every allowed odd h in [3,609] is strongly realized at order at most eight.
   This is inherited finite exact data, not a fresh census.
4. **Order-eleven bank:**
   `01-canon/theorems/THM-4104-selected-order-eleven-strong-ear-solid-interval.md`,
   selected-image theorem: every odd h in [429,80265] has a strong
   eleven-vertex witness. Its frozen bank and independent audit are specified
   there. Combining 3 and 4 gives a witness of order at most eleven for every
   odd h in [285,80265]. Neither the complete order-eleven spectrum nor any
   values beyond this bank are needed.
5. **Strong minimum:**
   `05-knowledge/results/rotation_owners_spectrum_20261020e/spectrum.md`,
   section 3; primary citation A. H. Busch, *A Note on the Number of Hamiltonian
   Paths in Strong Tournaments*, EJC 13 (2006), N3, Corollary 1:

   ```
   f(3k)=3*5^(k-1), f(3k+1)=5^k, f(3k+2)=9*5^(k-1), k>=1.
   ```

   In particular `f(n+3)=5f(n)`. The precise f is needed; a fitted recurrence
   with base sqrt(3) is not a dependency.
6. **Deletion-only comparison and its loss boundary:**
   `05-knowledge/results/session133_fractran_spectrum_refinement_20261101/core_reduction/THEORY.md`,
   sections 2–3. Its sharp core budget K(h) refers to deleting singleton SCCs,
   not to smallest witnesses across unrelated tournaments.
7. **Two prime lanes already inherited:**
   `01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md`.
   This packet does not claim the prime/7p reduction as new.

Targeted searches of canon and the recent spectrum packets for SCC compression,
batching, packing, minimum witness order and logarithmic witnesses located the
induced-core bound and its explicit caveat, but not this finite-bank repacking
argument. This is a bounded repository audit, not an exhaustive literature
novelty claim.

## 3. An abstract finite-bank principle

Let c>0, odd integers L>=3 and B>=(L-2)^2, and an integer m be such that:

* every attained odd g in [L,B] has a witness of order at most m;
* `m <= c ln L`;
* every strong tournament C with H(C)>B satisfies `|C|<=c ln H(C)`;
* every attained q with 1<q<L has a witness of order `<c ln q+1`.

Then every attained h>1 has a witness of order at most `ceil(c ln h)`.

**Proof.** Delete singleton SCCs. Keep each factor above B with its original
strong witness. Replace each factor in [L,B] by its bank witness. Greedily
multiply the remaining factors, which are positive odd numbers in [3,L-2],
until a batch first reaches L. Before its last factor the product is at most
L-2, so a completed batch has value at most `(L-2)^2<=B`; replace it by a bank
witness and reset the running product to one. All completed pieces have
order at most c times the logarithm of their count. At most one residual
product q remains, with q<L; it is attained because it is a product of actual
SCC counts. Its order is `<c ln q+1`, or zero if q=1. Order-join every chosen
witness. Counts multiply to h and the total order is `<c ln h+1`. An integer
strictly below x+1 is at most ceil(x). QED.

This is a **value-preserving replacement**, not a graph contraction or
induced-subtournament operation. It retains recipes to actual witnesses;
formal factors satisfying parity/floor inequalities alone are not input
witnesses to the theorem.

## 4. Exact integer specialization

Take

```
B=80265, D=80267, L=285, m=11, c=22/ln D.
```

The bank is inherited from sources 3–4. The batching sandwich is

```
(L-2)^2 = 283^2 = 80089 <= B < D < 285^2 = 81225.
```

The last inequality gives `11<c ln285`. Thus every replaced block meets its
logarithmic cost bound.

### Large-component test

Suppose H(C)>B and n=|C|. Since H(C) is odd, H(C)>=D. If n<=22 then
`H(C)^22>=D^22>=D^n`. For n>=23, use the exact strong minimum. The finite
integer checks are

```
f(23)=140625,   140625^22 >= D^23;
f(24)=234375,   234375^22 >= D^24;
f(25)=390625,   390625^22 >= D^25;
5^22 >= D^3.
```

Since `f(n+3)=5f(n)`, the last inequality propagates the first three to every
n>=23. Consequently `H(C)^22>=f(n)^22>=D^n`, or `n<=c ln H(C)`.

### Residual table

For every allowed odd q<285, these inherited witness-order upper bounds work:

| q | witness order at most |
|---|---:|
| 3 | 3 |
| 5 | 4 |
| 9,11,13,15 | 5 |
| allowed odds 17..33 | 6 |
| odds 35..61 | 7 |
| odds 63..283 | 8 |

The table follows from the explicit strong spectra in source 3. The final row
uses a union over orders: for example q=67 is absent at order eight but present
at order seven. One must not replace “at most eight” by “exactly eight.”

For each row, if r is its smallest q and s its stated order bound, verify

```
D^(s-1) < r^22.
```

The six checks are at `(r,s)=(3,3),(5,4),(9,5),(17,6),(35,7),(63,8)`.
Monotonicity then proves `s<c ln q+1` throughout each row. The omissions 7
and 21 cannot occur as the residual product because the residual is already
attained. A unit residual requires no vertices. The abstract principle now
proves the theorem.

## 5. The cost-sensitive spectrum framework

The SCC product gives an exact min-plus convolution, not merely a support
closure. For h>1 put

```
R(h) = min_{d|h, 1<d<h} (mu(d)+mu(h/d)).
```

Then

```
mu(h) = min(sigma(h), R(h)).
```

**Proof.** An optimal witness with h>1 has no singleton SCC. If it has one
SCC it costs at least sigma(h). If it has two or more SCCs, divide their
ordered list into two nonempty pieces; their counts give a proper
factorization d*(h/d), at cost at least mu(d)+mu(h/d). Conversely either a
strong witness or an order join realizes each finite candidate. QED.

This makes the minimum-order spectrum a shortest-path problem on the divisor
DAG: an edge x->xd labelled by an actual strong value d has cost sigma(d).
The least distance 1->h is mu(h). Knowing only that d is attained drops the
cost coordinate and cannot answer minimum-order questions.

Call h a **strict cost atom** when `mu(h)<R(h)`. This is only a defined
framework, not an assertion that all strong values have this property.
Every minimum-order witness for such h is strong. Conversely if `mu(h)=R(h)`,
a reducible minimum-order witness exists. Strong and reducible minimizers may
tie; the min-plus equality keeps that distinction visible.

Every attained h has a minimum-order realization as an order join of strict
cost atoms: recursively split a value whenever `mu(h)=R(h)`; the factors
strictly decrease, so recursion terminates and preserves total cost. This
is stronger information than arithmetic prime factorization, but depends on
the currently unknown mu table.

**Hostile:** 9 is not an irreducible element of the multiplicative spectrum:
`9=3*3`. Nevertheless `mu(9)=5<6=mu(3)+mu(3)`, so 9 is a strict cost atom.
Similarly `mu(15)=5<7=mu(3)+mu(5)`. Strong tournaments are graph-level
order-join atoms; their H values need not be arithmetic atoms, and arithmetic
atoms do not capture cost savings.

The full order-n support is exactly `{h:mu(h)<=n}` for n>=1, by neutral
padding (with the harmless empty-unit convention). The theorem therefore
bounds the first-appearance filtration for every value that ever appears.
It is not an interval theorem and does not predict which h are missing.

## 6. Computational consequence and stopping boundary

Let `N(h)` be the smallest integer n with `80267^n>=h^22`. It is computed
with integers and equals the displayed ceiling. For h>1,

```
h in S iff some tournament of order <=N(h) has H=h.
```

One may use the smaller of N(h) and the inherited K(h). Thus prescribed-value
membership retains the inherited `2^{O((log h)^2)}` finite-search complexity,
with a smaller leading order bound. The claim is not polynomial time and not
that this brute-force search is practical at large h.

The proof has an explicit large-component bottleneck: it leaves components
whose H exceeds the finite bank unchanged. Enlarging the positive bank can
improve the constant; this argument alone does not give a coefficient below
the strong-floor asymptotic `3/ln5`.

**Inclusion hostile:** k ordered directed triangles have count 3^k and an
induced H-preserving core of order 3k. No core vertex may be deleted. Moon's
interval chain with k parts equal to two has the same count on 2k+1 vertices
(and at k=2 gives 9 on five rather than six vertices). Thus replacement can
compress where induced deletion cannot. The new theorem deliberately changes
the graph while preserving only H and a witness recipe. SCC count, embedding,
endpoint response, marked vertices and future continuation counts are lost.

## 7. Reproduction

```
python3 verify.py --out verification.json
python3 -O verify.py --out verification.optimized.json
cmp verification.json verification.optimized.json
```

The standalone script uses only Python's standard library. It verifies every
new inequality by integer exponentiation, all possible small batch crossings,
all allowed residuals, the floor propagation through order 300, and 1200
repacking recipes whose factors are actual small strong counts or Moon-minimum
counts. It independently enumerates every labelled tournament through order
six, verifies SCC multiplication, and checks the DP against literal vertex
permutations through order five. Its min-plus replay is explicitly capped by
total order six: a strong-factor bank of order at most six can form products
requiring more than six vertices (39 via 3*13 costs eight), so a bound on
factor order is not a bound on joined witness order. It also verifies the
two-value mesh lemma from COFINALITY.md and small actual interval-chain
hostiles. It **does not regenerate the inherited 80265 graph bank**; that
finite-positive dependency is source 4's audited output.
