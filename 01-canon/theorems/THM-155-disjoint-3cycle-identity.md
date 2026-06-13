# THM-155: Disjoint 3-Cycle Identity for Regular Tournaments

**Status:** PROVED
**Proved by:** kind-pasteur-2026-03-12-S60
**Method:** Eigenvalue decomposition + trace identity

## Statement

For any regular tournament $T$ on $n = 2m+1$ vertices:

$$K(T) = c_5(T) - 2 \cdot \text{ov}_1(T) - 2 \cdot \text{ov}_2(T) = -\frac{3n(n^2-1)(n^2-9)}{320}$$

where:
- $c_5(T)$ = number of directed 5-cycles
- $\text{ov}_1(T)$ = number of pairs of 3-cycle vertex sets sharing exactly 1 vertex
- $\text{ov}_2(T)$ = number of pairs sharing exactly 2 vertices

**Equivalently:**

$$\text{disj}_3(T) + \frac{c_5(T)}{2} = \binom{c_3}{2} - \frac{3n(n^2-1)(n^2-9)}{640}$$

where $c_3 = n(n^2-1)/24$ is the number of 3-cycles (constant for regular tournaments) and $\text{disj}_3$ counts vertex-disjoint 3-cycle pairs.

## Proof Sketch

### Step 1: Trace Identity (Key Lemma)

For any regular tournament on $n = 2m+1$ vertices:

$$2\,\text{tr}(A^5) + 5\,\text{tr}(A^4) = 2m^5 + 5m^4 - m(5m+2)$$

**Proof:** Use eigenvalue decomposition. All non-trivial eigenvalues $z_k$ of a regular tournament have $\text{Re}(z_k) = -1/2$. Write $z_k = -1/2 + iy_k$.

The LHS equals $2m^5 + 5m^4 + \sum_k z_k^4(2z_k + 5)$.

Direct computation shows $\text{Re}(z^4(2z+5)) = 1/4 - 5y^2$.

Therefore $\sum_k \text{Re}(z_k^4(2z_k+5)) = (n-1)/4 - 5\sum y_k^2$.

Since $\sum |z_k|^2 = \text{tr}(A^T A) - m^2 = nm - m^2 = m(m+1)$ and $\sum y_k^2 = m(m+1) - (n-1)/4 = mn/2$, this sum equals $-m(5m+2)$.

### Step 2: Connecting traces to overlap counts

- $c_5 = \text{tr}(A^5)/5$
- $\text{tr}(A^4) = 2\sum_{i \to j} \mu_{ij}(\mu_{ij}+1)$ where $\mu_{ij} = (A^2)_{ij}$
- This gives $\text{tr}(A^4) = 2\,\text{sum\_}\mu^2 + 2\,\text{sum\_}\mu$
- $\text{sum\_}\mu = \binom{n}{3} - c_3 = n(n-1)(n-3)/8$ (transitive triple count)
- $\text{ov}_2 = (\text{butterfly} - 3c_3)/2$ where butterfly $= \sum_{i \to j} \lambda_{ij}^2 = \sum_{i \to j}(\mu_{ij}+1)^2$

### Step 3: Assembly

Combining Steps 1-2:

$$c_5 + \text{sum\_}\mu^2 = \frac{2\,\text{tr}(A^5) + 5\,\text{tr}(A^4)}{10} - \text{sum\_}\mu = f(n)$$

The closed form $f(n) = n(n-1)(n-3)(n^2+4n-17)/160$ follows from substituting the constant values.

Then $K = c_5 + 2\,\text{ov}_2 - n L_c(L_c-1)$ where $L_c = (n^2-1)/8$ gives the result.

## Key Sub-Results

1. **c3(v) constant:** For regular tournaments, every vertex is in exactly $(n^2-1)/8$ 3-cycles.
2. **Lambda-mu identity:** For edge $i \to j$ in a regular tournament, $\lambda_{ij} = \mu_{ij} + 1$ (from $A^2 - A^{T2} = A^T - A$).
3. **Three rigid classes at n=7:** (c5, ov1, ov2, disj) = (42,63,21,7), (36,57,24,10), (28,49,28,14), all giving K=-126.

## Verification

- Exhaustive at $n = 5$ (all 1024 tournaments)
- 37,379 random regular tournaments at $n = 7$: K = -126 (100%)
- 3,284 random regular tournaments at $n = 9$: K = -486 (100%)
- Confirmed at $n = 5, 7, 9, 11, 13, 15, 17, 19, 21$ (all odd, including non-prime)

## Connections

- Originally discovered for circulant tournaments on $\mathbb{Z}_p$ (kind-pasteur-S58-S59)
- The eigenvalue proof shows it holds for ALL regular tournaments, not just circulant
- Connects to Savchenko's cycle count formulas for doubly regular tournaments
- The quantity $c_5 + 2\,\text{ov}_2 = n(n^2-1)(n^2-9)/160$ is the cleanest form
