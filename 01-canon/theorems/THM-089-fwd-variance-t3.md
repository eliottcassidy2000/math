# THM-089: Forward-Edge Variance Formula

**Status:** PROVED
**Proved by:** opus-2026-03-07-S46
**Scope:** All tournaments

---

## Statement

Let T be a tournament on n vertices with t3 directed 3-cycles. For a uniformly random permutation P of V(T), let fwd(P) = #{i : P[i] -> P[i+1] in T}. Then:

$$\text{Var}[\text{fwd}(P)] = \frac{n+1}{12} + \frac{4\,t_3}{n(n-1)}$$

Equivalently: E[fwd^2] = (n-1)(3n+2)/12 + 4t3/(n(n-1)).

---

## Proof

Write fwd(P) = sum_{i=0}^{n-2} X_i where X_i = 1 if P[i] -> P[i+1] in T.

**Step 1: Marginals.**
For any i, Pr(P[i]->P[i+1]) = 1/2 (by symmetry: any ordered pair of vertices is equally likely to be forward or backward). So E[X_i] = 1/2, Var(X_i) = 1/4.

**Step 2: Non-adjacent covariance vanishes.**
For |i-j| >= 2, the four vertices P[i], P[i+1], P[j], P[j+1] are distinct. We show:

E[X_i X_j] = #{(a,b,c,d) distinct : a->b, c->d} / n(n-1)(n-2)(n-3)

For each directed edge a->b (there are C(n,2) such), the edges among the remaining n-2 vertices form a tournament on n-2 vertices, which ALWAYS has exactly C(n-2,2) directed edges, regardless of T.

Therefore: #{(a,b,c,d) distinct : a->b, c->d} = C(n,2) * C(n-2,2).

E[X_i X_j] = C(n,2)*C(n-2,2) / [n(n-1)(n-2)(n-3)] = 1/4.

Cov(X_i, X_j) = 1/4 - 1/4 = 0. Done.

**Step 3: Adjacent covariance.**
For j = i+1: E[X_i X_{i+1}] = #{(a,b,c) distinct : a->b->c} / [n(n-1)(n-2)].

The number of directed 2-paths a->b->c is:
sum_b in(b) * out(b) = sum_b (n-1-s_b)*s_b = (n-1)*sum s_b - sum s_b^2

Using the classical identity: t3 = C(n,3) - sum_b C(s_b,2), we get:
sum s_b^2 = 2*C(n,3) + C(n,2) - 2*t3

Therefore: #{2-paths} = (n-1)*C(n,2) - 2*C(n,3) - C(n,2) + 2*t3
= n(n-1)(n-2)/6 + 2*t3 = C(n,3) + 2*t3.

So: E[X_i X_{i+1}] = [C(n,3) + 2*t3] / [n(n-1)(n-2)] = 1/6 + 2*t3/[n(n-1)(n-2)]

Cov(X_i, X_{i+1}) = -1/12 + 2*t3/[n(n-1)(n-2)].

**Step 4: Assembly.**
Var[fwd] = (n-1) * 1/4 + 2*(n-2) * [-1/12 + 2*t3/(n(n-1)(n-2))]

= (n-1)/4 - (n-2)/6 + 4*t3/(n(n-1))

= [3(n-1) - 2(n-2)]/12 + 4*t3/(n(n-1))

= **(n+1)/12 + 4*t3/(n(n-1))**. Done.

---

## Verified

- n=4: 3 F-classes, exact match for all (max error 0)
- n=5: 8 F-classes, exact match for all (max error 0)
- n=6: 24 F-classes, exact match for all (max error 0)

---

## Corollaries

1. **Transitive tournament (t3=0):** Var = (n+1)/12. This equals the variance of the Eulerian distribution (descent count of permutations), consistent with F(transitive) = Eulerian polynomial.

2. **Regular tournament (t3 = C(n,3)/4 at n=4k+3):** Var = (n+1)/12 + C(n,3)/(n(n-1)) = (n+1)/12 + (n-2)/3 = (5n-5)/12.

3. **Minimum variance:** t3=0 (transitive) minimizes variance. Maximum variance at max t3.

4. **Worpitzky third coefficient:** The deviation delta_2 = 2(n-2)*t3 in THM-084 follows from Var[fwd] = (n+1)/12 + 4t3/(n(n-1)), since the Worpitzky third coefficient is algebraically related to E[fwd^2].

5. **Non-adjacent uncorrelatedness:** The fact that Cov(X_i, X_j) = 0 for |i-j|>=2 is a consequence of tournaments being complete: any subset of vertices forms a tournament with a fixed number of edges. This is a special property of tournaments (not general digraphs).

---

## Connection to other results

- **INV-082 (W-hierarchy):** w_{n-3} = (n-2)![2t3 - C(n,3)/2] is linearly related to t3. The variance formula provides the link: both involve the second moment of the forward-edge distribution.

- **THM-084 (Worpitzky):** The t3-dependence of Worpitzky coefficients at level d=2,3 is ultimately driven by the variance formula.

- **Spectral:** Since tr(A^3) = 3*t3, the variance has a spectral formula: Var = (n+1)/12 + 4*tr(A^3)/(3n(n-1)).
