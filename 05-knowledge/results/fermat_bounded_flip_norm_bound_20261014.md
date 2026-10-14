# Explicit orbit/norm bounds for bounded flips of a cyclic carousel

14 October 2026. Status: **PROVED**, with a finite exact certificate for the stated numerical thresholds and small Fermat cases. This is a separate continuation of [the cyclotomic orientation note](fermat_platonic_cyclotomic_orientation_20261014.md). No canonical theorem is edited and no priority claim is made.

## Result and scope

Let `p` be an odd prime, `m=(p-1)/2`, and begin with the carousel connection set

`S_0={1,...,m} subset F_p`.

For `D subset {1,...,m}`, set `r=|D|` and

`S_D=(S_0 minus D) union (-D)`.

A flip here reverses an entire translation orbit of `p` edges, not one individual edge. Let `T_D` be the cyclic tournament with `i -> j` iff `j-i in S_D`, let `H` count directed Hamiltonian paths, and put `disc=det(I+A-A^T)/2^(p-1)`.

**Theorem.**

1. For every odd prime and every such `D`,

   `disc(T_D) <= (1+4r+4r/(p-1))^((p-1)/2)`.

2. If `r>=1` and `p >= max(257, 4(4r+2)^3)`, then

   `H(T_D) >= p disc(T_D)`.

3. For **one difference-orbit flip**, the sharper threshold `p>=257` suffices, independently of which step is flipped.

4. For every Fermat prime, including any not yet known, **every one-difference-orbit flip of a carousel** satisfies `H >= p disc`. The only Fermat primes below 257 are 3, 5, and 17; their cases are certified using triangle packings below.

For example, item 2 proves the bound at `p=65537` for every choice of at most five flipped step pairs. It does not enumerate those tournaments. The same conclusions hold after a multiplier relabelling of the initial carousel.

The comparison uses two precisely defined carriers: the Fourier energy of a short signed cyclotomic polynomial, and positive collections of disjoint directed triangles in the odd-cycle collection formula (OCF). It does **not** construct a direct injection from lattice classes to Hamiltonian-path orbits.

## 1. Cancel the carousel unit before estimating the norm

Use the established identity

`disc(T_D)=Norm_(Q(zeta_p)/Q)(eta_D)`,

where `eta_D=1+sum_(s in S_D) zeta_p^s`. Choose the primitive root `t` with `zeta_p=t^2`. Since `p=2m+1`,

`1+zeta_p+...+zeta_p^m=(1-t)/(1-t^2)=1/(1+t)`.

Therefore

`eta_D=1/(1+t)+sum_(j in D)(t^(-2j)-t^(2j))`.

Multiply by the cyclotomic unit `1+t`, whose norm is `Phi_p(-1)=1`, obtaining

`g_D(t)=1+sum_(j in D)[t^(-2j)+t^(1-2j)-t^(2j)-t^(2j+1)]`.

Thus `disc(T_D)=Norm(g_D(t))` exactly. This removes the geometrically large but norm-one carousel factor before any absolute-value estimate is taken.

### The cancellation lemma survives arbitrary wrapping

Reduce exponents modulo `p`, writing `g_D(t)=sum_(e=0)^(p-1) c_e t^e`.

The positive exponents `-2j,1-2j`, over `j=1,...,m`, partition the integer interval `[-p+1,-1]`. Hence the selected positive terms occupy distinct nonzero residues. The negative exponents `2j,2j+1` partition `[2,p]` and likewise occupy distinct residues. The initial exponent zero can only cancel the negative term at exponent `p`.

Consequently wrapping never combines two terms of the same sign: it only cancels opposite signs. Thus

`c_e in {-1,0,1}`, `sum_e c_e=1`, `E_D:=sum_e c_e^2 <= 1+4r`.

No hypothesis that the flipped steps are small is needed.

There is a direct combinatorial reading of `E_D`. Let `b_e` be the indicator of `{0} union 2S_D` on the cyclic residues, with indices modulo `p`. Then

`c_e=b_e+b_(e-1)-1`,

because `(1+t)eta_D` has coefficients `b_e+b_(e-1)`, and subtracting `Phi_p` leaves the same cyclotomic value. Thus `E_D` is exactly the number of adjacent equal bits in this cyclic binary word: either `00` or `11`. The carousel has one such defect. A flipped difference orbit changes two bits and can add at most four defects. The following energy-sensitive norm bound therefore applies to any low-defect connection set, even when its flip distance from the chosen carousel is large.

### Parseval and the norm bound

Orthogonality on all `p`th roots gives

`sum_(b=0)^(p-1) |g_D(t^b)|^2 = p E_D`.

Since `g_D(1)=1`, the average over nontrivial roots is `(p E_D-1)/(p-1)`. Applying AM–GM to the squared absolute values yields

`disc(T_D) <= ((p E_D-1)/(p-1))^((p-1)/2)`

`<= (1+4r+4r/(p-1))^((p-1)/2)`.

This also handles `r=0`, giving `disc=1`. The first bound can be stronger when wrapping creates cancellations.

For the shortest-step flip `D={1}`, multiplying `g_D` by the norm-one monomial `t^2` gives the especially compact polynomial

`1+t+t^2-t^4-t^5`.

The general proof above covers every one-step orbit, not just this polynomial.

## 2. An elementary positive OCF lower bound for every regular tournament

Let `T` be any regular tournament on odd `n` vertices. Every vertex lies in exactly

`(n^2-1)/8`

directed triangles. Indeed, if its out- and in-neighbourhoods are `U,W`, each of size `q=(n-1)/2`, the number of arcs from `U` to `W` is

`q^2-q(q-1)/2=q(q+1)/2`.

These arcs correspond bijectively to directed triangles through the vertex. Hence the total triangle count is `n(n^2-1)/24`.

Suppose `i` disjoint triangles have already been chosen. Removing their `3i` vertices destroys at most `3i(n^2-1)/8` original triangles. Thus at least

`(n-9i)(n^2-1)/24`

triangles remain whenever this expression is positive. If `N_k` counts unordered collections of `k` disjoint directed triangles, sequential choice and division by `k!` give

`N_k >= (1/k!) product_(i=0)^(k-1) [(n-9i)(n^2-1)/24]`.

The trusted [OCF, THM-002](../../01-canon/theorems/THM-002-ocf.md), is

`H(T)=sum_(disjoint odd-cycle collections C) 2^|C|`.

Keeping only the collections of `k` triangles therefore gives `H(T)>=2^k N_k`. For `n>=11` take `k=floor(n/10)`. Every factor `n-9i` is at least `k`, and `k^k>=k!`, proving the convenient explicit estimate

`H(T) >= ((n^2-1)/12)^floor(n/10)`.                         (L)

The exact product is stronger; (L) is used because its threshold is easy to certify. This is a positive count of actual disjoint triangle collections, with exactly the OCF weight, rather than a signed determinant being treated as a cardinality.

Every cyclic tournament is regular, so (L) applies to `T_D`.

## 3. A uniform explicit threshold for arbitrary bounded flip number

Let `r>=1`, set `C=4r+2`, and assume `p>=max(257,4C^3)`. Then `p-1>=4r`, so the norm bound gives

`disc <= C^((p-1)/2)`.

Also `(p^2-1)/12 >= C^6`, while `floor(p/10)>=p/11` for `p>=110`. Therefore

`H >= C^(6p/11)`.

The ratio to `p C^((p-1)/2)` is at least

`C^((p+11)/22)/p >= 6^floor(p/22)/p > 1`.

For the final elementary inequality, put `k=floor(p/22)`. For `k>=3`, induction from `6^3>=22(3+1)` gives `6^k>=22(k+1)>p`. This proves theorem item 2.

The threshold is a sufficient condition, not an asserted optimum. For fixed `r` it is finite and explicit; it also covers a growing number of flips whenever `4(4r+2)^3<=p`.

At `p=65537`, `4(4*5+2)^3=42,592 <= 65,537`, so any set of at most five flips is covered. Six flips are not certified by this particular coarse threshold; that is not a counterexample to the inequality.

## 4. The sharper one-flip threshold 257

For `r=1` and `p>=257`, put `C=321/64`, so

`5+4/(p-1) <= C`.

For odd integers `n>=257`, define

`R(n)=((n^2-1)/12)^floor(n/10) / [n C^((n-1)/2)]`.

The five exact rational checks at `n=257,259,261,263,265` all give `R(n)>1`. Approximate ratios, shown only for readability, are respectively

`29.2025, 8.51261, 14045.3, 4133.19, 1212.71`.

For a ten-step increment, monotonicity of `(n^2-1)/12` gives

`R(n+10)/R(n) >= ((n+10)^2-1)/12 * n/(n+10) * C^(-5)`

`>= 4,918,010,284,343,296 / 2,729,968,765,186,401 > 1`.

The final ratio uses `n>=257`; all its comparisons are exact rational inequalities. These five residue classes cover every odd integer from 257 onward. Together with (L), this proves item 3 for every odd prime `p>=257`.

## 5. Close the small Fermat cases using triangle packings, not a Hamiltonian census

For a given tournament let `P_T(z)=sum_k N_k z^k`, where `N_k` counts disjoint directed-triangle collections. Then `P_T(2)<=H(T)` by OCF.

The verifier computes `P_T` by choosing the least remaining vertex `v`: either leave `v` unused, or choose its unique triangle in the packing and recurse on the complement. This partitions the objects and gives an exact recurrence. It separately checks `N_1` by literal triangles and `N_2` by literal pairs of triangles.

At `p=3`, the single flip has `disc=1`, `P_T(2)=3`. At `p=5`, both possible flips have `disc=1`, `P_T(2)=11`. At `p=17`:

| Flipped step | disc | Triangle-only OCF lower bound | p disc |
|---:|---:|---:|---:|
| 1 | 103 | 33,059,017 | 1,751 |
| 2 | 307 | 29,982,765 | 5,219 |
| 3 | 24,617 | 28,622,697 | 418,489 |
| 4 | 409 | 29,613,389 | 6,953 |
| 5 | 1,361 | 28,788,005 | 23,137 |
| 6 | 2,381 | 29,078,841 | 40,477 |
| 7 | 239 | 29,917,621 | 4,063 |
| 8 | 137 | 34,367,881 | 2,329 |

Every case satisfies the desired inequality without counting Hamiltonian paths. Every Fermat prime other than 3, 5, 17 is at least 257, so these finite exact cases and the large-prime theorem prove the all-Fermat-prime one-flip corollary. This reasoning does not assume that only five Fermat primes exist.

## 6. How this compares with established asymptotic theory

The assertion that `H>=n disc` eventually holds in the regular sector is already implied by much stronger existing counting theory. [Bill Cuckler, *Hamiltonian Cycles in Regular Tournaments*, CPC 16 (2007), 239–249](https://doi.org/10.1017/S0963548306008170), proves at least `n!/(2+o(1))^n` Hamiltonian cycles in every regular `n`-vertex tournament. Deleting a distinguished edge of a Hamiltonian cycle gives `n` distinct Hamiltonian paths, and the resulting path uniquely recovers its closing cycle. Thus `H>=n` times the cycle count.

On the determinant side, the regular skew spectrum has `(n-1)/2` pairs with mean squared magnitude `n`, giving the standard AM–GM bound

`disc <= ((n+1)/4)^((n-1)/2)`.

The logarithm of Cuckler's path lower bound has leading term `n log n`, whereas the logarithm of this determinant bound has leading term `(n/2) log n`. Hence all sufficiently large regular tournaments satisfy the inequality.

The concrete content here is therefore the cancellation-stable bounded-flip cyclotomic norm estimate, the elementary OCF triangle-packing comparison, the explicit thresholds, and the small-case closure that proves the one-orbit-flip corollary at every Fermat prime. We do not claim a new eventual-existence theorem, an optimal threshold, or a general solution of the corpus's strong-base conjecture.

## 7. The unresolved carrier and an exact obstruction to forgetting path order

The norm is the finite lattice index

`disc = |Z[zeta_p]/eta_D Z[zeta_p]|`.

This follows from the determinant of multiplication by `eta_D` on the cyclotomic integer lattice. For every prime cyclic tournament the number of translation orbits of Hamiltonian paths is `H/p`, so the desired direct proof would compare these two finite sets.

The present proof bounds their cardinalities through Fourier energy and positive OCF objects; it supplies no canonical injection between them. Two tempting order-forgetting maps cannot do that job:

- Summing the vertex characters along a spanning path gives `sum_v zeta_p^v=0` for every path.
- Summing consecutive character differences telescopes to `zeta_p^(end)-zeta_p^(start)`. After fixing the start at zero, it has at most `p-1` values.

For the step-one flip at `p=17`, the lattice quotient has 103 classes, while the endpoint encoding has at most 16 values. Thus that map cannot be a surjection onto the quotient. Any eventual path/lattice carrier must retain interior order or richer cycle information.

Finally, if the number of flipped orbits is comparable to `p`, the elementary estimates no longer close: the norm upper bound can have logarithm of order `(p/2)log p`, while the coarse triangle-packing lower bound has only `(p/5)log p`. This is a limitation of these bounds, not a counterexample to the open inequality.

## Verification

The separate standard-library script is

`04-computation/experiments/fermat_bounded_flip_norm_bound_check.py`.

It checks 302 connection-set controls across primes 3, 5, 7, 11, 17: reduction coefficients, norm versus skew determinant through independent multiplication matrices, and the norm bound. It then checks all eleven small-Fermat one-flip cases by triangle-packing recurrence, and verifies the threshold arithmetic. There are 2,788 explicit checks, including the cyclic run-defect identity in all 302 connection-set controls. Normal and `python3 -O` executions byte-match the frozen JSON; validation uses explicit exceptions, so optimized execution does not disable it.
