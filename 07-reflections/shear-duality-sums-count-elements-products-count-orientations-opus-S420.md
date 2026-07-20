# Shear duality: sums count elements, products count orientations

**Instance:** opus-2026-07-20-S420 (owner: shear the n*2^x+1 grid down 0/1/2(+);
sum or PRODUCT the columns; census all continuations of the triangular numbers;
"triangular numbers arise from relation itself"). **HYP-8180.** Script + frozen out:
`shear_calculus_triangular_continuations_opus_S420.py`.

**Fleet convergence note (priority):** kind-pasteur-S128c103 (HYP-8170) pushed the
shear CATALOG first — Proth 2^(1/s) spectrum with GF proofs, Pascal Pisot ladder,
Faulhaber/Rosetta bookkeeping, A045648 deviation-tuning with T(8,4) as discriminator,
OEIS-new harvest. Everything below either (a) independently re-confirms theirs by a
different instrument (exact minimal-recurrence detector over Fractions, no floats,
honest NONEs at Bell/subspace/Faulhaber — correct, those have no constant recurrence),
or (b) is distinct. Their "shear-1 = 2(2^m-1)" vs my "2^(m+1)-1" is a CONVENTION split,
not a dispute: they exclude the n=0 cell (value 1); with it, the +1 duties pay the
convolution's linear debt down to exactly -1 — Mersenne on the nose. Their
R_0(m) = m*2^(m+1)+1 self-similarity = my detected s=0 order-3 recurrence, same object.
Death-star-S59t claimed THM-1360 (figurate zoo) concurrently — cited, not duplicated.

## 1. Recurrence normal forms for the spectrum (reconciles kp's GFs)

Proth shear-s row sums satisfy the exact integer recurrence with characteristic
polynomial **(x^s - 2)(x - 1)^2 (1 + x + ... + x^(s-1))** — detector-verified at
s=2 (order 5: a(m) = a(m-1)+3a(m-2)-3a(m-3)-2a(m-4)+2a(m-5)) and s=3 (order 7),
matching kp's GF denominator (1-t)^2(1-2t^s)(1-t^s)/(1-t). Roots: the s radicals of
2 (growth 2^(1/s), phase period s), a double root 1 (the linear duty term), and the
(s-1) roots of unity (phase of the ones-count). Pascal dials x^(s+1) = x^s + 1
(Fibonacci, Narayana cows, A003269). **Simplicial C(n+x, x+1) = INTEGRATED Pascal:**
each integration appends an (x-1) root and a "-1" offset — shear-1 = Mersenne 2^m - 1,
shear-2 = Fibonacci - 1 (char (x-1)(x^2-x-1)). Corollary identity: the dyadic tower
and the dimension tower are **shear-1 equinumerous**:
Sum_{x+n=m} C(n+x, x+1) = 2^m - 1 = Sum_{x+n=m-1} (n*2^x + 1). Bijection hunt = handoff.

## 2. The collision law and the Fermat square (distinct)

Within a shear-1 row, two cells collide (n*2^x = n'*2^x', x+n = x'+n') iff
**n(2^d - 1) = d**, whose ONLY solution in positive integers is (n,d) = (1,1)
(2^d - 1 >= d + 1 for d >= 2 forces n < 1). So every row m >= 2 contains the Fermat
number 2^(m-1)+1 EXACTLY TWICE (cells n=1,x=m-1 and n=2,x=m-2) and no other repeat:
- the Fermat square (2^(m-1)+1)^2 divides every shear-1 row product;
- row products are perfect squares ONLY at m=2 (9 = 3^2) and m=3 (100 = 10^2),
  where the residual factor happens square (1 and 4); never again for m <= 12
  (checked; conjecturally never — residual contains a single 3 mod higher powers).
Another bilinear-rigidity statement in the m*b+-1 genus: the grid's unique internal
echo is Fermat-shaped.

## 3. THE HEADLINE: sums count elements, products count orientations (distinct)

The owner said "sum these columns (or in fact take their product)". The two
aggregations of the SAME sheared grid land on the repo's two foundational towers:

- **Sum side:** Sum_{x+n=m} n*2^x = 2^(m+1) - m - 2 = A000295(m+1) = Eulerian E(m+1,1)
  = permutations with exactly ONE DESCENT. Scale 2^m: subsets, elements, cut-space.
- **Product side (EXACT identity):** Prod_{x+n=m, n>=1} n*2^x = m! * 2^C(m,2) =
  (orderings of m) x (orientations of C(m,2) pairs) = **the count of ordered
  tournaments — the repo's base-path + tiles decomposition itself.** The staircase
  model (fixed Hamiltonian path + C(n-1,2) free tiles) IS the row product of the
  pure grid. With the +1s, log2(product)/C(m,2) -> 1 (1.30 at m=22, correction
  log2(m!)/C(m,2) -> 0): Proth products live on the tournament scale.

Sum : product :: descent : orientation :: cut-space : cycle-space :: 2^n : 2^T(n).
The triangular number mediates exactly as the owner said — "T_n is n in the
relational sense": it is the EXPONENT of the orientation tower, and the shear
product is how the grid spends it.

## 4. The museum of impersonations (session theme, five exhibits)

| impersonation | split |
|---|---|
| Proth shear-2 sums vs lazy caterer 1+T(m) | m=5 (18 vs 16) |
| Moser circle vs 2^(n-1) | n=6 (31 vs 32) |
| Proth diagonal k*2^(k-2)+1 vs companion Pell | 5th term (97 vs 99) |
| **S419's diagonal law a(n)=a(n-1)+a(n-2)+a(n-4)** | **m=8 (predicts 65, actual 68) — MY OWN MIRAGE, MISTAKE-198** |
| width of G_n = C(n-2,floor(.)) (repo classic) | n=7 (10 vs 15) |

Mechanism: low-order jets of different growth classes (quadratic / exponential /
Pisot / radical) collide; the shear calculus names the true class and predicts the
split. S419's law fit all three then-available instances and broke on the fourth —
filed as MISTAKE-198 with the rule: for super-exponential families, a recurrence fit
on fewer than 2*(order)+2 terms is a jet, not a law; always compute a continuation
term before claiming.

## 5. The relational synthesis and the census (distinct; death-star's zoo adjacent)

- **2*T_m + 1 = m^2 + m + 1 = Phi_6(m+1) = |PG(2,m)|** (all m): the owner's 2n+1
  column evaluated AT triangular arguments = projective-plane sizes. **The LRC deep
  well denominator 183 = 2*T_13 + 1 = |PG(2,13)|: the doubled-relation-plus-one of
  13.** The duty "+1" that closes a projective plane is the same +1 that closes the
  bilinear rim (S419).
- Centered polygonals = s*T(n-1) + 1: the bilinear m*b+1 rim with TRIANGULAR base.
- Hexagonal fold-back H(n) = T(2n-1): the polygonal tower re-enters the spine.
- Eight continuation axes out of T_n (the wide berth): arity C(n,k); geometry P(s,n);
  dimension C(n+d-1,d); exponent Faulhaber; field [n,2]_q (q=2: 1,7,35,155,651 —
  2-subspaces of F_2^n, the dyadic shadow of "edges"); orientation 2^T(n) (tournaments);
  centering s*T+1; projective 2T+1. Sums/products/shears act on each; the
  exponent-seat vs coefficient-seat dichotomy (kp's sharpening) says which dial
  exponential ladders and which stay polynomial.

## 6. Handoffs

(a) T(8,4) discrimination (kp's frame): pure 225 / flat-1 226 (my S419 law) /
A045648-tuned 218 — owner's row 8 settles it. (b) Bijection for the shear-1
equinumerosity (figurate cells of total index m <-> nonempty subsets <-> Proth cells
one row earlier). (c) Square-product conjecture: no perfect-square Proth row product
for m >= 4 (prove via the residual's 3-adic valuation). (d) Proth-prime census per
row (2,2,2,3,3,1,2,5,4,4,5,3,1 for m=2..14) — any structure vs Sierpinski obstructions?
(e) The descent/orientation duality: is there a shear-equivariant bijection realizing
sum->E(m+1,1) and product->m!*2^C(m,2) in one frame (BEST-like)?

Cross-links: kind-pasteur-S128c103/HYP-8170 (catalog, priority) · death-star-S59t/
THM-1360 (zoo, concurrent) · S419/HYP-8155 (the triangle; now corrected) ·
MISTAKE-198 · THM-1292/THM-1269 (bilinear rim) · CLAUDE.md tiling model (the product
identity's referent).
