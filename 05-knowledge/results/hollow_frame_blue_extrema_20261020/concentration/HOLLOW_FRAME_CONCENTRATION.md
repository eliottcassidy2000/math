# Concentration from the hollow-frame second moment

Status: **PROVED; INDEPENDENT PROOF REVIEW PASSED**, conditional only on
the already proved hollow-frame mean and second-moment identities. This is
a bounded analytic consequence, not a new frame census or an exact-maximum
assertion.

## 1. Objects and inherited inputs

For each integer \(k\geq1\), choose the labelled external frame \(D=(A,b)\)
uniformly from its \(2^{k(k-1)}\) possibilities. Let \(h=h(D)\) count the
reflection-symmetric directed Hamiltonian vertex words. The inherited
results give

\[
 h\text{ is positive and odd},\qquad
 \mathbb E h=k!,
\]

and

\[
 \frac{\mathbb E h^2}{k!}
 =\sum_{j=0}^{\lfloor k/2\rfloor}
   \frac{(k-2j)!}{4^j}
   [t^j]\left(\frac{1+t}{1-t}\right)^{k-2j+1}.
 \tag{1}
\]

These are moments over **uniform labelled external frames**. They are not
moments over tournament isomorphism classes, and they do not assert the
same distribution for the quotient \(B=h/|H_a|\) under arbitrary matching
completions.

Source inputs:

- [The second-moment identity](../second_moment/HOLLOW_FRAME_SECOND_MOMENT.md).
- [The extrema note](../HOLLOW_FRAME_BLUE_EXTREMA.md), Sections 1–3.

## 2. A uniform coefficient bound

Put

\[
 a_{k,j}=
 \frac{(k-2j)!}{4^j k!}
 [t^j]\left(\frac{1+t}{1-t}\right)^{k-2j+1},
 \qquad 0\leq j\leq\lfloor k/2\rfloor.
\]

All these numbers are nonnegative. The key estimate is

\[
 0\leq a_{k,j}
 \leq \frac{(k-j)!}{2^j j!k!}
 \leq \frac{1}{j!k^j}.
 \tag{2}
\]

The estimate holds simultaneously for the entire allowed range of \(j\);
there is no fixed-\(j\) asymptotic assumption.

**Proof.** In the coefficientwise order on power series with nonnegative
coefficients,

\[
 \frac{1+t}{1-t}
   =1+\sum_{\ell\geq1}2t^\ell
   \preceq 1+\sum_{\ell\geq1}2^\ell t^\ell
   =\frac1{1-2t}.
\]

Coefficientwise order is preserved by multiplication of such series.
Thus, with \(r=k-2j+1\geq1\),

\[
 [t^j]\left(\frac{1+t}{1-t}\right)^r
 \leq[t^j](1-2t)^{-r}
 =2^j\binom{r+j-1}{j}
 =2^j\binom{k-j}{j}.
\]

Substitution and cancellation give the first upper bound in (2).
For \(j\geq1\), each of the \(j\) factors in
\(k!/(k-j)!\) is at least \(k/2\), since \(j\leq k/2\).
This gives the second upper bound. For \(j=0\), both bounds are equalities.
\(\square\)

## 3. Relative variance and concentration

Define

\[
 R_k=\frac{\operatorname{Var}(h)}{(k!)^2}.
\]

The \(j=0\) summand of (1), after division by \(k!\), is 1. For \(k\geq2\),
the \(j=1\) summand is exactly

\[
 a_{k,1}
 =\frac{(k-2)!}{4k!}\,2(k-1)
 =\frac1{2k}.
\]

Consequently (2) gives the following finite-\(k\) estimate:

\[
 \boxed{\quad
 0\leq R_k-\frac1{2k}
 \leq \sum_{j=2}^{\lfloor k/2\rfloor}\frac1{j!k^j}
 \leq e^{1/k}-1-\frac1k ,
 \qquad k\geq2.
 \quad}
 \tag{3}
\]

In particular,

\[
 \boxed{\quad
 \frac{\operatorname{Var}(h)}{(k!)^2}
   =\frac1{2k}+O(k^{-2})
   \sim\frac1{2k}.
 \quad}
 \tag{4}
\]

For a purely rational, explicit error bound, the ratio of consecutive
terms in \(\sum_{j\geq2}1/(j!k^j)\) is at most \(1/(3k)\). Therefore

\[
 0\leq R_k-\frac1{2k}
 \leq \frac{3}{2k(3k-1)}.
 \tag{5}
\]

The formulas (3) and (5) are bounds, not claimed sharp error terms.
At \(k=1\), \(h=1\) and \(R_1=0\); (3) is intentionally stated only
for \(k\geq2\).

Because \(\mathbb E[h/k!]=1\), equation (4) proves the stronger statement

\[
 \mathbb E\left[\left(\frac h{k!}-1\right)^2\right]\longrightarrow0.
 \tag{6}
\]

Thus \(h/k!\to1\) in \(L^2\), and hence in probability. Quantitatively, for
every fixed \(\varepsilon>0\) and every \(k\geq2\), Chebyshev's inequality
gives

\[
 \Pr\left(\left|\frac h{k!}-1\right|\geq\varepsilon\right)
 \leq \frac1{\varepsilon^2}
       \left(\frac1{2k}+\frac{3}{2k(3k-1)}\right).
 \tag{7}
\]

Equivalently, the fraction of labelled external frames outside the
relative interval \(1\pm\varepsilon\) tends to zero. The relative
standard deviation is asymptotic to \(1/\sqrt{2k}\). No central limit
theorem, almost-sure coupling, or matching-completion distribution is
claimed.

## 4. An immediate maximum bound, with strict odd rounding

Let

\[
 M_k=\max_D h(D),\qquad
 Q_k=\frac{\mathbb E h^2}{\mathbb E h}
     =k!(1+R_k).
\]

The equality of this \(M_k\) with the maximum blue count among
self-converse tournaments of order \(2k+1\) uses the regular-subset
theorem already cited in the extrema note. The following inequality
for \(\max_D h(D)\) itself needs no group-theoretic input.

Pointwise \(h^2\leq M_k h\), so

\[
 M_k\geq Q_k
 =\sum_{j=0}^{\lfloor k/2\rfloor}
   \frac{(k-2j)!}{4^j}
   [t^j]\left(\frac{1+t}{1-t}\right)^{k-2j+1}.
 \tag{8}
\]

For \(k\geq2\), the inequality is in fact strict. Equation (3) shows that
\(\operatorname{Var}(h)>0\), so there is a positive-probability frame with
\(h<M_k\). Since every \(h\) is positive, that frame has \(h^2<M_kh\).
Taking expectations proves \(Q_k<M_k\).

As \(M_k\) is odd, the least possible integer allowed by this strict
inequality is

\[
 \boxed{\quad
 M_k\geq 2\left\lfloor\frac{Q_k+1}{2}\right\rfloor+1,
 \qquad k\geq2.
 \quad}
 \tag{9}
\]

This is the least odd integer **strictly greater** than \(Q_k\).
The strictness matters if \(Q_k\) is itself an odd integer. For example,
\(Q_3=7\), so (9) proves \(M_3\geq9\), whereas weak odd rounding would
give only 7. The exceptional case is \(M_1=Q_1=1\).

The first few formula-derived bounds are:

| \(k\) | \(Q_k\) | Bound from (9) | Previously known finite-exact \(M_k\) |
|---:|---:|---:|---:|
| 2 | \(5/2\) | 3 | 3 |
| 3 | \(7\) | 9 | 9 |
| 4 | \(217/8\) | 29 | 39 |
| 5 | \(265/2\) | 133 | 195 |
| 6 | \(25033/32\) | 783 | Not asserted here |

In particular,

\[
 M_k>k!\left(1+\frac1{2k}\right),\qquad k\geq2,
 \tag{10}
\]

and \(M_k\geq k!(1+R_k)\) with \(R_k\sim1/(2k)\).
This is a lower bound, not an asymptotic equivalent for \(M_k/k!\).
Concentration of almost all frames does not control the rare frames
that can attain the maximum.

## 5. Exact verification and scope

`verify_concentration.py` uses rational arithmetic to evaluate the
inherited coefficient formula. It computes the coefficients in two
different ways:

1. Repeated multiplication by \(1+2t+2t^2+\cdots\).
2. Choosing the \(s\) nonzero factors in a composition of \(j\):
   \([t^j]R(t)^r=\sum_s2^s\binom rs\binom{j-1}{s-1}\) for \(j>0\).

It checks (2) term by term, the finite sum and rational tail bounds,
the inherited exact variances for \(k=1,\ldots,7\), and the strict
odd-rounding rule. The default run covers \(k=1,\ldots,128\); the recorded
run covers \(k=1,\ldots,256\). These are checks of the analytic formula,
**not frame-census data**. The proof above, rather than those finite
checks, establishes the all-\(k\) statements.

Run:

    python3 -O verify_concentration.py --max-k 256 --output concentration_checks.json

All validation failures raise explicit errors and remain active under
`python3 -O`. The script uses only the Python standard library.

The independent proof review is recorded in `CONCENTRATION_REVIEW.md`.
No new exhaustive frame search, central limit theorem, exact maximum
conjecture, or literature-priority assertion is part of this continuation.
The supplement is add-only and leaves both source packets unchanged.
