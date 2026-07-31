# JC/GMC next-gate concept board

**Anchor:** THM-2824's open atomic inequality

```text
D_i(b,c)
 =2 L(V^3)L(d_i V)-3L(d_i V^2)L(V^2) >= 0,
V=f_c-f_b, d_i=f_(i+1)-f_i, 0<=i<b<c.
```

Equality is conjecturally only `(i,b,c)=(j,j+1,j+2)`.

## Live concepts

1. **Gamma integration by parts**
   - Representation: `d_i(s)e^-s=-(f_(i+1)(s)e^-s)'`.
   - Consequence:
     `L(d_i V)=L(f_(i+1)V')`,
     `L(d_i V^2)=2L(f_(i+1)VV')`.
   - Test: rewrite `D_i` as an oriented covariance or a one-crossing integral.

2. **Factorial total positivity**
   - Representation: products of normalized monomials are multinomial
     coefficients.
   - Invariant: determinants of moment-ratio vectors.
   - Test: normalize by one positive multinomial and search for a sum/product
     with positive coefficients in the gaps `k=b-i`, `r=c-b`.

3. **Adjacent-difference cone**
   - Operation: `f_c-f_b=sum_(j=b)^(c-1) d_j`.
   - Known orientation: the Wronskian of `f_b-f_a` and `f_c-f_b` is positive.
   - Test: determine whether `D_i(b,c)` is monotone under adjoining the next
     right atom or admits a positive double sum over adjacent pairs.

4. **Ratio derivative**
   - Representation: for `V_epsilon=V+epsilon d_i`,
     `D_i=B A'-A B'=-A^2 (B/A)'` at zero, where
     `A=L(V_epsilon^2)`, `B=L(V_epsilon^3)`.
   - Predicate: adding an earlier atom decreases the ratio `B/A`.
   - Test: find a monotone-likelihood-ratio or one-crossing theorem.

5. **Hostile asymptotics**
   - Regimes: fixed gaps with `b->infinity`; `i=0`; `i=b-1`;
     `r=1`; `k=1`; both gaps large.
   - Test: exact sign and leading coefficient; search beyond the existing
     `c<=50` bank without mistaking a scan for proof.

## Portfolio

- **Anchor:** prove/refute universal `D_i>=0`.
- **Niche:** derive an exact equality classification for an infinite boundary
  family even if the full inequality remains open.
- **Wildcard:** compare the determinant to a two-point tournament/orientation
  only if the pairwise observable is intrinsically antisymmetric and preserves
  `D_i`; otherwise retain it as a signed moment determinant.

## Meta-patterns used

- “One item left” requires a typed residual.
- Exactify before interpreting a decimal.
- Attack a proposed bound before extending it.
- Use redundant paths as detectors.
- A failed certificate is not a failed theorem.

