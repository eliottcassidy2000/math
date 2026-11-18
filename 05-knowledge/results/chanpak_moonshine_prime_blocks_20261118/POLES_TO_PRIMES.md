# What a tournament encodes at its poles, and how to place a prime across a cyclic seam

Two outputs: an equivalence theorem for the information carried by all independent substitutions, and a constructive two-block inverse with literal prime and sevenfold witnesses. The inverse is not yet universal.

## 1. Which structural information do all independent block readings see?

For a labeled quotient tournament T with adjacency A let X=diag(x) and put
```
D_T(x) = det(I-A X),
F_T(x) = 1 + 1^t X(I-A X)^(-1)1.
```
Here F includes the empty word. The exact run-word substitution law is archived as [THM-3121](../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md); q optional paths give F^q followed by contraction with the ordered nonempty-run counts of each actual inserted block. This uses only monomials with every marked color present. For ordinary vertices the required exponent is one.

**Classification theorem.** For two labeled tournaments on the same vertex set, the following are equivalent:

1. The same labeled three-vertex subsets are directed cycles.
2. Their full vertex-variable determinants D_T coincide.
3. After independently replacing every vertex by an arbitrary nonempty tournament, the resulting q-slot count agrees for a single fixed positive integer q, for every tuple of blocks.
4. It agrees for every q and every tuple.

Here the labels index which block is substituted, so this is stronger than agreement of H or of a fixed finite profile. To prove 3 ⇒ 2, triangular inversion on independent transitive widths gives equality of the full-support word series for that q. The [pole theorem](CLONE_POLES.md) identifies the maximal strong-component determinants from this series: each nontrivial component supplies a genuine irreducible pole on exactly its vertex variables. Poles from induced restrictions have supports contained in a component. Mutual containment identifies the same maximal supports, then the same determinant factors on those supports. Singleton factors are 1. For 2 ⇒ 1, the coefficient on three vertices is -1 precisely for a directed triangle.

For 1 ⇒ 2 use [Boussaïri–Ille–Lopez–Thomassé, Theorem 2](https://perso.ens-lyon.fr/stephan.thomasse/liste/bilt.pdf): tournaments with the same labeled cyclic triples are connected by reversals inside intervals (homogeneous modules). An interval inversion preserves the number of spanning directed-cycle covers on each induced set, with the sign of each determinant term unchanged: reverse each maximal run inside the interval in the oriented cyclic word, and reverse a cycle wholly inside the interval. The operation is a bijection; uniform outside contacts guarantee the arcs across the interval still work. This proves equality of all principal minors. Alternatively, reverse path runs in the total substituted interval and then in each constituent block, restoring the original internal blocks and proving 1 ⇒ 4 directly.

A determinant lemma makes the arithmetic implication explicit. As A+A^t=J-I,
```
F_T(x) = det(I-A X+JX)/D_T(x)
       = [product_i (1+x_i)] D_T(-x/(1+x))/D_T(x).
```
The second equality uses transpose and det(I+UV)=det(I+VU), so 2 ⇒ 4 without a separate path bijection. The reciprocal denominator is not just an obstruction: it computes the numerator and hence every contextual count. Arbitrary ordinary boundary contacts can still distinguish quotients with identical *marked-only* poles; the resulting boundary read is the next section's constructive freedom.

## 2. A square-norm version of the theorem

Let K=A−A^t, x_i=2t_i/(1−t_i), L=diag(t), and define s=1^t L(I−K L)^(-1)1. In these coordinates the involution x↦−x/(1+x) is simply t↦−t. Rank-one determinant updates give
```
F_T(x(t)) = (1+s)/(1−s).
```
In fact s=N/D, a ratio of two nonnegative square-norm polynomials:
```
D = sum_{|S| even} Pf(K_S)^2 t_S,
N = sum_{|S| odd} Pf(Khat_{S+*})^2 t_S,
Khat = [[K,1],[-1^t,0]].
```
This follows by expanding determinants by principal minors and taking the coefficient of u in det(I−Khat diag(t,u))=D(1+us), a Schur complement. The augmented Pfaffian square on a triple is 9 for a cyclic triple and 1 otherwise. The directed q-slot count is 2^(-n) times the squarefree spanning coefficient of ((D+N)/(D−N))^q. Formal terms with repeated variables cannot contribute to that coefficient.

Under a cut switch K↦diag(σ)Kdiag(σ), the even polynomial D is invariant. The odd coefficients change to the squares of the signed sums of the Pfaffians on S minus one vertex (the signs σ_i enter by the border expansion). Thus one switching lattice supplies a common even backbone, but not the marked signed norm data that produces H. This is a precise way to separate the invariant scaffold from the part a target must program. It does not assert that these polynomials are theta series of the Leech lattice.

There is a concrete older parallel. In [the Leech deep-hole swap, THM-4758](../../../01-canon/theorems/THM-4758-leech-deep-hole-swap-theorem-and-the-deep-hole-theta-identity.md), the discriminant cusp at Coxeter number h has one isotropic chart per divisor d of h and swap sends d to h/d; theta differences are recorded by signed defect entries. At E6^4, h=12, the six charts d=(3,1,4,2,6,12) give the same word d+1=(4,2,5,3,7,13) as the [six signed tournament rails](../chanpak_six_signed_height_pairs_20261118/README.md). The archived cusp defect formula gives (+4,-28,-4,+4,-4,+28) in this order. The shared fact is the divisor skeleton and the need to retain signed data under the swap. No eta-product-to-tournament map is asserted. Here the reciprocal observer suggests where to look: keep the cycle scaffold and alter what enters or exits the seam, rather than look for a new denominator every time.

## 3. A productive crossed seam with a target-first inverse

Use quotient outgoing masks A=(10,4,9,2), and switch vertex 3 to get A'=(2,12,1,5). Replace vertices 0 and 2 by blocks L,R; keep vertices 1 and 3 as singletons. The marked blocks stay in the same switching lattice across the two hosts. Write f_L(q),f_R(q) for optional-slot counts, and (h,b,c) and (x,y,z) for their first three values. The exact response for each q≥1 is
```
f_A(q)−f_A'(q) =
 q[ f_L(q) (f_R(q+1)−f_R(q−1))
   −(f_L(q+1)−f_L(q−1)) f_R(q) ].
```
Indeed all odd cycles with mixed quotient contacts contain vertex 1. Those not containing vertex 3 cancel between the two frames. Those containing both singletons and using both blocks give the same ordered slot contribution on each side. The remaining cycles with both singletons use only R on the first side or only L on the second. Such a run has odd size. By the binomial-type slot identity and f_S(-1)=(-1)^|S| H(S), the odd-support split read for any block is (f(q+1)−f(q−1))/2. The factor 2q per directed odd cycle in the [Rédei–Berge packing formula](https://arxiv.org/abs/2307.05569) proves the response.

The expression inside brackets is a discrete Wronskian of the two monic slot polynomials. If their degrees differ, its leading term is twice the degree difference times q^(deg L+deg R−1). If the degrees agree, replace f_R by f_R−f_L in the determinant, lowering its degree, so the same leading-term argument applies unless f_L=f_R. Thus neutrality at every positive q is equivalent to identity of the two entire block slot polynomials. Neutrality at q=2 alone is the useful weaker requirement.

Therefore a correlated neutral pair at q=2 needs exactly
```
(c−h)/b = (z−x)/y,
```
and then its ordinary change A' minus A is bx−hy. The absolute values, not just their difference, follow by counting one or two visits to either block:
```
H_A  = 2hx + (2b−h)(y−x),
H_A' = 2hx + (2y−x)(b−h).
```
The ordered nonempty-run counts for one path have coefficients (5,3,2,2) for visit shapes (1,1),(1,2),(2,1),(2,2), and transpose on the other side. Substituting N_1=h and N_2=b−2h gives the displayed forms.

This suggests solving for the blocks from a requested target p: propose smaller positive h,x and factor p−2hx. For each ordered factor pair (u,v), require h+u even, set b=(h+u)/2, y=v+x, then realize the two profiles. If a neutral repair is needed, demand the same score at the third read; if not, do not impose that extra constraint. For fixed seeds with transitive ordinal additions at either side, the divisor tests are exactly that (u+h)/(2b_seed) and (v+x)/y_seed are positive powers of two. This is a terminating exact placement test for those seeds, not an approximation or a scan along a chosen ray.

Physical witnesses:

| (h,b,c) | (x,y,z) | Score | Prime H values | Common height f(2) |
|---|---|---:|---:|---:|
| (65,1952,21537) | (65,1808,19953) | 11 | 6699827, 6709187 | 2230556648 |
| (349,5104,39789) | (389,5984,46629) | 85/11 | 55432627, 55329667 | 9106026664 |

For example prescribe 6699827 and take h=x=65. Then p−8450=3839·1743, prescribing b=1952 and y=1808 exactly. The two actual blocks meet those values and the third-score match, so the signed switch is not merely formal. The independent replay supplies their adjacency masks, recomputes both full-host profiles literally, and verifies primality certificates. Applying the [paired slope-62 readers](PRIME_SUPPORT.md) with three ordinal singletons to the second pair realizes both p and 7p for p=4516644657971 and 4516644555011.

**What is still hard.** The six rails choose height when the odd mantissa is among a fixed set. The crossed seam chooses height from an exact factorization, but needs two compatible physical profiles. Neither says that for an arbitrary prime the factorization search meets the physical profile set. The next constructive target is a variable family with a proved inverse to either b=(h+u)/2 or y=v+x for target-chosen factors, and a second family to cover the complementary factors. Demanding the third-score match for all prime witnesses would solve a harder problem than necessary. The replay's `select` function starts with a target and a bank of actual blocks. For each left profile and right ordinary count, it tests the divisibility and asks for exactly the missing right height, rather than building all pairs and matching their output totals.

There is a promising generalization for two unexpanded contact vertices a,b. Let l_a,l_b,l_ab denote the y, z, yz coefficients of log F in their variables, and keep the marked-core F_M fixed. Then the yz coefficient of F^q is q F_M^q(l_ab+q l_a l_b). A literal contact edit with Δl_ab=−2Δ(l_a l_b)≠0 makes the two-slot boundary series neutral while retaining a one-slot move. This is a construction equation for the missing boundary data, not a claim such contacts have been found in that generality. The pole theorem tells us to search boundary/contact changes or deliberately correlated substitutions; the crossed seam is a positive example of the latter.

Reproduce with `python3 replay_new.py` in this directory; SymPy is used only for the independent Pfaffian-square checks by integer determinants.
