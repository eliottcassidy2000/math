# Session 127 independent audit: positive unit norms force convergents in every order

**Status:** PROVED elementary all-order implication; CITED classical continued-fraction criterion; separate FINITE-EXACT hostile/positive checks. No all-height classification of solutions is inferred from the finite checks.

## 1. Definitions and all-order root selection

For each integer `k >= 3`, put

\[
f_k(x)=x^k-\sum_{j=0}^{k-1}x^j,\qquad
F_k(p,q)=p^k-\sum_{j=0}^{k-1}p^jq^{k-j}.
\]

Let `alpha_k` be the positive root of `f_k`. The proof below only uses this defining polynomial, not its irreducibility or the location of its other roots. For `k=3,4,5`, the independently certified irreducibilities in the companion inheritance packet identify `F_k` with the number-field norm of `p-q alpha_k`. Consequently this argument does not need a Pisot theorem, a unit-group computation, or a finite census.

The useful normalized function is

\[
h_k(x)=\frac{f_k(x)}{x^k}=1-\sum_{j=1}^k x^{-j},\qquad x>0.
\]

It is strictly increasing, since

\[
h'_k(x)=\sum_{j=1}^k jx^{-j-1}>0.
\]

Moreover `h_k(1)=1-k<0`, `h_k(2)=2^{-k}>0`, and
`h_k(3/2) <= h_3(3/2)=-11/27<0`. Thus there is exactly one positive root, and

\[
\frac32<\alpha_k<2.
\]

The rational-root theorem makes it irrational: the only possible rational roots of the monic integral polynomial with constant term `-1` are `+/-1`, and `f_k(1)=1-k`, while `f_k(-1)=1` for even `k` and `-2` for odd `k`.

**Uniform differential estimate.** On `[1,2]`,

\[
h'_k(t)\ge t^{-2}+2t^{-3}+3t^{-4}>\frac14+\frac14=\frac12.
\]

This reciprocal normalization, rather than a crude bound on `f'_k`, removes the order-dependent threshold.

## 2. All-order positive-unit theorem

**Theorem.** For every integer `k >= 3` and positive integers `p,q`,

\[
\boxed{|F_k(p,q)|=1\quad\Longrightarrow\quad p/q
\text{ is a regular continued-fraction convergent of }\alpha_k.}
\]

### Proof

1. **Primitivity.** If `d=gcd(p,q)`, homogeneity gives
   `F_k(p,q)=d^k F_k(p/d,q/d)`. Therefore a unit value forces `d=1`.

2. **Positive-root branch is automatic under positive coordinates.** If `p<=q`, then the leading term cancels against the first negative term:

   \[
   p^k-p^{k-1}q\le0,\qquad
   F_k(p,q)\le -q^k-pq^{k-1}\le-2,
   \]

   using `k>=3`. Hence `p>q`.

3. **No point above 2.** Write `x=p/q`. If `x>2`, monotonicity gives

   \[
   h_k(x)>h_k(2)=2^{-k},\qquad
   F_k(p,q)=p^k h_k(x)>p^k/2^k>q^k\ge1,
   \]

   contradicting the unit norm. Thus `1<x<=2`.

4. **The endpoint.** If `x=2`, then `p=2q` and `F_k(p,q)=q^k`. Hence `q=1,p=2`. Since `3/2<alpha_k<2`, its continued fraction starts `[1;1,...]`, and `2/1` is its first nontrivial convergent.

5. **Interior estimate.** Suppose `1<x<2`. The interval between `x` and `alpha_k` lies in `(1,2)`. The mean-value theorem and the strict differential bound give

   \[
   |x-\alpha_k|
   <2|h_k(x)|
   =\frac{2|F_k(p,q)|}{p^k}
   =\frac2{p^k}
   \le\frac2{(q+1)^3}
   <\frac1{2q^2}.
   \]

   The final strict inequality follows from

   \[
   (q+1)^3-4q^2=q^2(q-1)+3q+1>0\qquad(q\ge1).
   \]

6. **Classical conclusion.** Legendre's criterion says that a reduced rational `p/q` within `1/(2q^2)` of an irrational is a regular continued-fraction convergent. Apply it to the preceding strict inequality. QED.

**Classical dependency, checked 2026-10-31.** The official mathlib documentation states the criterion and gives its formal theorem declaration `Real.exists_rat_eq_convergent`; it also describes the strengthened inductive proof. This is a primary formal-proof source, not a search-summary inference:

- [Mathlib.NumberTheory.DiophantineApproximation.Basic, Legendre's theorem](https://leanprover-community.github.io/mathlib4_docs/Mathlib/NumberTheory/DiophantineApproximation/Basic.html#Real.exists_rat_eq_convergent).
- The documentation's pinned source link is commit `70078b363d019ff4e975ac0f06eabd3cdb64b0de`, file `Mathlib/NumberTheory/DiophantineApproximation/Basic.lean`. The source link did not fetch through the available web interface; the inspected declaration and proof description were on the official documentation page. No local Lean build is claimed.

## 3. General residual gate: a sufficient, not necessary, condition

Let `k>=3`, `p,q` be positive integers, `1<p/q<2`, and `m=F_k(p,q)`. Then

\[
\boxed{|\alpha_k-p/q|<\frac{2|m|}{p^k}.}
\]

In particular,

\[
\boxed{4|m|q^2\le p^k
\quad\Longrightarrow\quad p/q\text{ is a convergent of }\alpha_k.}
\]

Equality in the displayed arithmetic gate is permitted: the derivative estimate was strict. Primitivity is not needed for this version. If the fraction reduces to denominator `Q<=q`, the obtained bound `1/(2q^2)` is at most `1/(2Q^2)`.

For a fixed finite prime set `S`, the restriction

\[
m=\pm\prod_{\ell\in S}\ell^{e_\ell},\qquad e_\ell\ge0,
\]

specifies support, not magnitude. It does **not** by itself imply the residual gate: the exponents can be arbitrarily large in the definition. For a concrete failure with fixed `S={11}`, `F_3(3,2)=-11` and the rational generating-function value is `U_3(2/3)=-12/11`, but `3/2` is not a cubic-root convergent: its certified CF prefix has successive denominators `1,1,6,...`. A Thue–Mahler candidate still needs either this quantitative inequality, exact CF certification, or another root-approximation argument. No Thue–Mahler completeness or effective bound is imported here.

## 4. Exact norm-to-tail identity, with sign and adjacency retained

Suppose `p_n/q_n` is an actual regular CF convergent of `alpha_k`, let `p_{n-1}/q_{n-1}` be its predecessor, and let

\[
T=\alpha_{n+1}=[a_{n+1};a_{n+2},\ldots],\qquad
r=q_{n-1}/q_n.
\]

The usual finite-composition identity, obtained directly by composing the maps `t -> a_j+1/t`, is

\[
\alpha_k=\frac{p_n T+p_{n-1}}{q_n T+q_{n-1}},\qquad
p_nq_{n-1}-p_{n-1}q_n=(-1)^{n-1}.
\]

Therefore

\[
\alpha_k-p_n/q_n
=\frac{(-1)^n}{q_n^2(T+r)}.
\]

Define the real divided-difference polynomial

\[
C_k(x)=\frac{f_k(x)}{x-\alpha_k},
\]

with removable value `C_k(alpha_k)=f'_k(alpha_k)`. For positive `x != alpha_k`, `C_k(x)>0`: both `h_k(x)` and `x-alpha_k` have the same sign. With `x=p_n/q_n>0`, the exact identity is

\[
\boxed{F_k(p_n,q_n)
=\frac{(-1)^{n+1}q_n^{k-2}C_k(x)}{T+r},\qquad
T+r=\frac{q_n^{k-2}C_k(x)}{|F_k(p_n,q_n)|}.}
\]

Thus

\[
a_{n+1}=\left\lfloor
\frac{q_n^{k-2}C_k(x)}{|F_k(p_n,q_n)|}
-\frac{q_{n-1}}{q_n}\right\rfloor.
\]

This formula is exact, not a field-norm-only observer: the root, degree, denominator, local cofactor, predecessor, and index/sign are retained. It explains why unit norms give a scale `q_n^(k-2)` for the next complete quotient, and why a large next partial quotient is not equivalent to unit norm.

For the cubic anchor `n=4`, `p_n/q_n=103/56` and predecessor `46/25`, the sign formula gives `F_3(103,56)<0`, consistent with `F_3=-1` and `alpha_3>103/56`. Since `56 alpha_3-103=alpha_3^{-16}`, the ordinary CF identity gives

\[
T=\frac{\alpha_3^{16}-25}{56}.
\]

The independent cubic lane certifies `305<T<306`; this audit does not infer that integer floor from numerical rounding.

## 5. Failure boundaries and exact hostile witnesses

- **Wrong embedding/sign:** `F_4(-1,1)=1`, but `-1` is not a convergent of the positive root in `(3/2,2)`. Positivity/root selection is a real hypothesis, not cosmetic notation.
- **Convergent does not imply unit norm:** `11/6` is the cubic convergent `[1;1,5]`, while `F_3(11,6)=-7`.
- **Unit norm does not imply a 305-sized digit:** `F_3(2,1)=1`, but the next partial quotient is `5`.
- **A large digit does not imply unit norm:** the exact rational-interval CF certificate finds a cubic convergent at index `172` (indices start at `0`) with next digit `503` and nonunit norm. The full `p,q,m` are in the JSON artifact. They are not floating-point discoveries promoted without certification.
- **Unit norm does not imply membership in the cyclic sparse inverse-power source:** the companion source kernel is an additional property. This theorem takes binary-form unit values as input; it does not classify how they were generated.
- **Signed normalization matters:** the exact sparse cubic source `alpha_3^{-16}=-103+56 alpha_3` maps to `(p,q)=(-103,-56)`. Passing to the positive pair `(103,56)` multiplies the algebraic element by `-1`, so the odd-degree norm flips from `+1` to `-1`; the rational input `q/p` is unchanged.

## 6. Exact independent finite audit

The checker uses integer homogeneous evaluation, independent homogeneous Horner evaluation, and the sparse identity

\[
(p-q)F_k(p,q)=p^{k+1}-2p^kq+q^{k+1}.
\]

For each `k=3,...,12`, it isolates the unique positive root by 256 exact rational bisections from `(3/2,2)`. Continued-fraction terms are accepted only when the entire rational isolating interval has a common integer part, and the interval is transformed exactly. The loop continues through the first convergent denominator greater than 150, so all convergents relevant to the finite box have been included.

**Universe:** `3<=k<=12`, `1<=q<=150`, `1<=p<=3q`; exactly **339,750 pairs**. This deliberately includes the excluded `p<=q` and `p>2q` regions. The checker verifies every unit solution against the certified CF list and independently bounds the root distance for every residual-gate hit.

| k | Unit solutions `(p,q,F_k)` in the stated box | Number satisfying the interior residual gate |
|---|---|---:|
| 3 | `(2,1,1)`, `(103,56,-1)` | 3 |
| 4 | `(2,1,1)` | 2 |
| 5 | `(2,1,1)` | 0 |
| 6 | `(2,1,1)` | 1 |
| 7–12 | `(2,1,1)` at each order | 0 at each order |

The cubic residual hits include the nonprimitive `(206,112,-8)`, intentionally testing reduction rather than silently presupposing primitivity. The absence of interior residual hits in several orders is not evidence against the theorem: it is a sufficient gate and the interval is small.

An additional 4,096-bit cubic isolation certifies 1,114 CF terms and the stated next-digit-503 nonunit hostile. The only claim involving this longer run is the exact displayed finite witness.

Reproduction:

```bash
python3 -B verify_uniform_norm.py > verify_uniform_norm.normal.txt
python3 -B -O verify_uniform_norm.py > verify_uniform_norm.optimized.txt
cmp verify_uniform_norm.normal.txt verify_uniform_norm.optimized.txt
```

All checks use an explicit `check` function, so `python -O` does not remove them. Both runs completed and their outputs agree byte-for-byte. Full intervals, CF prefixes, residual-gate rows, and hostile integers are in `verify_uniform_norm.output.json`.

**Honest frontier:** the implication is proved for every order and every positive unit-norm pair. The finite box is not a global classification of binary-form unit solutions, sparse inverse powers, or Thue–Mahler solutions at orders four and five.
