# MRB square blocks and Many Worlds: audit and repaired theorem

Status: **PROVED** for the identities, inequalities, and elementary asymptotics proved below; **FINITE-EXACT** for the separately identified finite checks; physical realization and rationality implications remain outside the proved scope.

## Source and scope

Audited supplied files:

- `../session31-input/document24.pdf`: Marvin Ray Burns, *What Lebsir's Prime Progressions Gave the MRB Problem*, printed pp. 20-23. The printed page 21 was rendered and visually inspected. The sign error described below is present in the PDF, not an extraction artifact.
- `../session31-input/documentMRBconnections.pdf`: Marvin Ray Burns, *The Many Worlds View of the MRB Constant*, 15 pages. Mathematical statements were checked against the supplied text; page 12's eta-derivative and full-decomposition formulas were also visually verified in the rendered PDF. No claim is made to have independently verified historical priority or experimental feasibility.

The inherited mechanism is finite symmetric alternating quadrature. The cheapest hostile is the first square block. The correction preserves the block and global MRB limits but rejects an extra parity factor. The useful sidecar is a nonnegative second-derivative kernel, which provides the missing error control. The live concept board is: alternating parity; symmetric curvature; Euler moments; absolute grouping; algebraic conjugation; generator versus accumulator.

No prime-progression conjecture is assumed anywhere below. No irrationality conclusion, general Galois cancellation theorem, or novel-literature-priority claim is made.

## 1. A real sign error and its minimal witness

Put

\[
f(x)=\exp(\log x/x)-1,\qquad c_n=n(n+1),
\]

\[
B_n=\sum_{k=n^2}^{(n+1)^2-1}(-1)^k f(k)
    =\sum_{d=-n}^{n}(-1)^d f(c_n+d),
\]

because every pronic center is even. Define

\[
\Lambda_{n,d}=f(c_n-d)-2f(c_n)+f(c_n+d).
\]

The corrected exact identity is

\[
\boxed{B_n=(-1)^n f(c_n)+\sum_{d=1}^{n}(-1)^d\Lambda_{n,d}.}
\tag{1}
\]

Indeed, the coefficient of the center before this subtraction is
\(1+2\sum_{d=1}^n(-1)^d=(-1)^n\). Printed page 21 instead multiplies every curvature term by \((-1)^n\), giving \((-1)^{n+d}\) rather than \((-1)^d\). The same incorrect exponent appears in its following formula for the entire constant. Both displays must be repaired.

At \(n=1\),

\[
B_1=\sqrt2-\sqrt[3]3,
\qquad
\Lambda_{1,1}=1-2\sqrt2+\sqrt[3]3,
\]

whereas the printed formula gives \(2-3\sqrt2+\sqrt[3]3\). Their difference (printed minus actual) is \(2\Lambda_{1,1}\). Since \(\sqrt2>7/5\) and \(\sqrt[3]3<3/2\), \(\Lambda_{1,1}<-3/10\), so this is not a negligible numerical discrepancy.

Exact rational root enclosures from the accompanying script give:

- actual block in \([-0.028036007934313334,-0.028036007934313333]\);
- printed block in \([-0.800391116811876765,-0.800391116811876764]\).

The printed expression happens to agree on every even \(n\), making an even-only check an ineffective hostile. The later full Euler expansion has the correct sign convention and survives.

## 2. General positive curvature kernel

This theorem is independent of the particular MRB kernel and of square intervals.

**Theorem.** For an integer \(r\ge1\), real center \(c\), and \(F\in C^2([c-r,c+r])\), define

\[
T_rF(c)=\sum_{d=-r}^{r}(-1)^{r+d}F(c+d).
\]

For \(0\le u\le r\), set

\[
K_r(u)=\sum_{\substack{1\le d\le r\\d\ge u}}(-1)^{r+d}(d-u).
\]

Then

\[
\boxed{T_rF(c)-F(c)
=\int_0^r K_r(u)\{F''(c-u)+F''(c+u)\}\,du.}
\tag{2}
\]

Moreover,

\[
K_r(u)\ge0,
\qquad
2\int_0^r K_r(u)\,du=\frac{r(r+1)}2.
\tag{3}
\]

Consequently, for some \(\xi\in[c-r,c+r]\),

\[
\boxed{T_rF(c)-F(c)=\frac{r(r+1)}2F''(\xi).}
\tag{4}
\]

Thus \(F''\ge0\) implies \(T_rF(c)\ge F(c)\), despite the signed point weights. In particular

\[
|T_rF(c)-F(c)|\le\frac{r(r+1)}2\sup_{[c-r,c+r]}|F''|.
\tag{5}
\]

**Proof.** Twice integrate \(F''\) in each direction:

\[
F(c-d)-2F(c)+F(c+d)
=\int_0^d(d-u)\{F''(c-u)+F''(c+u)\}\,du.
\]

Sum with weights \((-1)^{r+d}\). This yields (2), since the constant weight of \(T_r\) is one. To see positivity, let \(v=r-u=m+\theta\), with \(m=\lfloor v\rfloor\) and \(0\le\theta<1\). Reversing the index gives

\[
K_r(u)=\sum_{j=0}^m(-1)^j(v-j)
=\begin{cases}m/2+\theta,&m\text{ even},\\(m+1)/2,&m\text{ odd}.\end{cases}
\]

At integer breakpoints the formulas agree by continuity. Also

\[
2\int_0^rK_r(u)\,du
=\sum_{d=1}^r(-1)^{r+d}d^2=\frac{r(r+1)}2,
\]

the last identity following immediately by induction. The weighted integral mean-value theorem on the two-sided interval proves (4). This proof does not require Euler polynomials.

For an arithmetic progression of step \(h>0\), replace \(F(x)\) by \(F(c+hx)\): the mass becomes \(h^2r(r+1)/2\), and \(\xi\) lies in \([c-hr,c+hr]\).

**Fourth-derivative refinement.** If \(F\in C^4\), the same argument applied to \(F''\) gives

\[
\left|T_rF(c)-F(c)-\frac{r(r+1)}2F''(c)\right|
\le \frac{r(r+1)(r(r+1)-1)}{24}\sup |F^{(4)}|.
\tag{6}
\]

Here \(\int_0^r u^2K_r(u)du=r(r+1)(r(r+1)-1)/24\). A direct proof is to integrate each \((d-u)u^2\), obtaining \(d^4/12\), then use the alternating fourth-power sum. This strengthens a formal Taylor expansion to a certified finite remainder.

## 3. Certified MRB curvature and tail bounds

Define the correctly signed block correction

\[
R_n=B_n-(-1)^nf(c_n)=\sum_{d=1}^n(-1)^d\Lambda_{n,d}.
\]

The derivative printed in the source is correct:

\[
f''(x)=e^{\log x/x}\left(\frac{2\log x-3}{x^3}
             +\frac{(\log x-1)^2}{x^4}\right).
\tag{7}
\]

For \(x\ge9\), \(f''(x)>0\) and \(f'''(x)<0\). To check the latter, write \(g(x)=\log x/x\). Then

\[
g'=(1-\log x)/x^2<0,\quad
g''=(2\log x-3)/x^3>0,\quad
g'''=(11-6\log x)/x^4<0,
\]

and \(f'''=e^g(g'''+3g'g''+(g')^3)<0\).

Using (2)-(4), for every \(n\ge3\),

\[
\boxed{\frac{c_n}{2}f''(n^2+2n)
\le(-1)^nR_n\le\frac{c_n}{2}f''(n^2).}
\tag{8}
\]

This supplies both the eventual sign and a rigorous remainder in place of the source's informal cancellation argument. Squeezing the endpoints proves the sharper asymptotic

\[
\boxed{R_n\sim(-1)^n\frac{2\log n}{n^4}.}
\tag{9}
\]

Thus its claimed \(O(\log n/n^4)\) decay is true **after correcting the sign identity**. The first block has negative curvature in the normalized gauge, so do not export (8)'s positivity down to \(n=1\).

For an explicit uniform estimate, let \(A=e^{1/e}\). From (7), for \(n\ge3\),

\[
|R_n|\le 2A(1+1/n)
       \left(\frac{\log n}{n^4}+\frac{(\log n)^2}{n^6}\right).
\tag{10}
\]

For every integer \(N\ge3\), put

\[
I_4(N)=\frac{\log N/3+1/9}{N^3},
\qquad
I_6(N)=\frac{(\log N)^2/5+2\log N/25+2/125}{N^5}.
\]

Both summands in (10) are decreasing on \([N,\infty)\). Comparing each tail with its integral gives

\[
\boxed{\sum_{n>N}|R_n|
\le 2e^{1/e}\left(1+\frac1{N+1}\right)(I_4(N)+I_6(N))
<3\left(1+\frac1{N+1}\right)(I_4(N)+I_6(N)).}
\tag{11}
\]

The last strict inequality uses \(e^{1/e}<3/2\). For an elementary verification, \(e>5/2\), while
\(\log(3/2)=2(1/5+(1/5)^3/3+\cdots)>2/5\).

The script certifies rational upper bounds on \(\log N\) by range reduction to \([1,2]\) and the positive atanh series. All reported decimals for (11) round outward, not to nearest. Examples:

| Cutoff N | Certified upper bound on sum of absolute curvature tail |
|---:|---:|
| 10 | 0.002916802930051020125517 |
| 100 | 0.000004988801113793516135 |
| 1000 | 0.000000007248352833489898 |
| 10000 | 0.000000000009544628508874 |

These are bounds on \(\sum_{n>N}|R_n|\), not on the entire MRB constant truncated after \(N\) original terms or after \(N\) square blocks. Approximating the full constant still requires controlling the pronic-center sum as well.

## 4. Which global convergence claims survive?

Since \(f(c_n)\sim2\log n/n^2\), (9) gives

\[
B_n\sim(-1)^n\frac{2\log n}{n^2},
\qquad \sum_n|B_n|<\infty.
\]

The original alternating series converges because \(f(k)\to0\) and is decreasing from \(k=3\). Square-block partial sums are a subsequence of its partial sums. Therefore

\[
C_{\rm MRB}=\sum_{n\ge1}B_n
=\underbrace{\sum_{n\ge1}(-1)^nf(n(n+1))}_{P}
 +\underbrace{\sum_{n\ge1}\sum_{d=1}^n(-1)^d\Lambda_{n,d}}_{R}.
\tag{12}
\]

Both the outer pronic series and the outer curvature series converge absolutely. Even the full triangular curvature array is absolutely summable: \(|\Lambda_{n,d}|\le d^2\sup f''\) on each large block and \(\sum_{d\le n}d^2=O(n^3)\), giving an \(O(\log n/n^3)\) absolute row bound.

This does **not** make the original summands absolutely summable: \(f(k)\ge\log k/k\), whose sum diverges. One may rearrange the completed blocks, not arbitrarily ungroup and rearrange their individual terms.

Absolute block convergence is also not unique to square geometry. If \(a_k\ge0\) is eventually decreasing and blocks begin at increasing indices \(b_j\), the alternating sum of each consecutive finite block has magnitude at most \(a_{b_j}\). Thus \(\sum_j a_{b_j}<\infty\) suffices. For \(a_k=f(k)\), any partition whose starts grow at least like \(j^{1+\epsilon}\), for fixed \(\epsilon>0\), has absolutely summable blocks. What the square partition additionally supplies is its especially simple symmetric center and exact Euler moments.

## 5. Euler structure and its exact qualification

Let \(E_j(x)\) be defined by

\[
\frac{2e^{xt}}{e^t+1}=\sum_{j\ge0}E_j(x)t^j/j!.
\]

The normalized block's moment generating function is

\[
\sum_{d=-r}^{r}(-1)^{r+d}e^{td}
=\frac{\cosh((r+1/2)t)}{\cosh(t/2)}.
\tag{13}
\]

Its odd moments vanish and its even moments are

\[
\sum_{d=-r}^{r}(-1)^{r+d}d^{2m}=E_{2m}(r+1)\quad(m\ge0).
\tag{14}
\]

For \(m=0\), the left hand side includes \(d=0\) with \(d^0=1\). The source's intermediate identity

\[
(-1)^r\sum_{d=1}^r(-1)^dd^{2m}=\tfrac12E_{2m}(r+1)
\]

is correct only for \(m\ge1\). At \(m=0\), its left side is 0 or 1 and its right side is \(1/2\). The full block expansion is nevertheless correct because the constant term is handled by total weight one.

For MRB, \(f(z)=\exp(\operatorname{Log}z/z)-1\) is analytic in the disk \(|z-c_n|<c_n\); all offsets satisfy \(|d|\le n<c_n\). Consequently its Taylor series converges absolutely at every point of each finite block, so termwise finite summation is valid:

\[
\boxed{(-1)^nB_n=\sum_{m\ge0}\frac{E_{2m}(n+1)}{(2m)!}f^{(2m)}(c_n).}
\tag{15}
\]

In particular \(E_2(n+1)=c_n\), \(E_4(n+1)=c_n(c_n-1)\), as printed. More generally every even moment is a polynomial in \(\lambda=r(r+1)\), because the symmetry \(r\mapsto-r-1\) fixes \(E_{2m}(r+1)\). The invariant ring over the rationals is \(\mathbb Q[r(r+1)]\). The first moments are

\[
P_0=1,\quad P_1=\lambda,\quad
P_2=\lambda^2-\lambda,\quad
P_3=\lambda^3-3\lambda^2+3\lambda,
\]
\[
P_4=\lambda^4-6\lambda^3+17\lambda^2-17\lambda.
\]

They depend on the **block radius**, not on the evaluation center. For square blocks the two happen to coincide, \(\lambda=c_n\). This coincidence must not be silently exported to arbitrary centers. No interchange of the two infinite sums over block number and derivative order is asserted here.

## 6. Analytic symmetry does not certify conjugate smallness

The source's rationality section is explicitly prospective, not a proof. Its strongest justified survivor is that a coherent conjugate action, if independently constructed and quantitatively bounded, might preserve a useful second-difference cancellation. There is no automatic implication from real-variable symmetry to this coherence.

A concrete hostile is \(n=3,d=1,c_n=12\). The real quantity

\[
\Lambda=11^{1/11}-2\,12^{1/12}+13^{1/13}
\]

is approximately 0.00153826. There is an algebraic conjugation fixing the first two radicals but replacing \(13^{1/13}\) by \(\zeta_{13}^{6}13^{1/13}\). Indeed let \(K\) be exactly the compositum over \(\mathbb Q\) of the splitting fields of \(x^{11}-11\) and \(x^{12}-12\), and the cyclotomic field \(\mathbb Q(\zeta_{13})\). This chosen compositum has degree dividing a product of \(11\varphi(11)\), \(12\varphi(12)\), and 12, and is therefore prime to 13. Eisenstein gives \([\mathbb Q(13^{1/13}):\mathbb Q]=13\), so \(13^{1/13}\notin K\). Because \(K\) contains \(\zeta_{13}\), its splitting extension obtained by adjoining this radical is cyclic of degree 13, and all such rotations occur while fixing \(K\).

That transformed expression has modulus approximately 2.41694. The magnitude separation can also be certified without floating point: positivity and monotonicity of \(f''\) give \(0<\Lambda\le f''(11)<1/500\), where the final inequality is checked with rational bounds for \(\log11\) and \(11^{1/11}\). Since \(13^{1/13}>6/5\) and \(\operatorname{Re}\zeta_{13}^6=-\cos(\pi/13)<-4/5\),

\[
\operatorname{Re}\bigl(\Lambda+(\zeta_{13}^6-1)13^{1/13}\bigr)
<1/500-(9/5)(6/5)=-1079/500.
\]

Thus the conjugate modulus exceeds \(1079/500=2.158\), while the original is below 0.002. The displayed approximate values are only optional diagnostics; the separation has this exact certificate. The example disproves automatic conjugate-smallness preservation, not the possibility of a specially chosen coherent family. Any norm/irrationality argument still needs its field, all relevant conjugates, denominator and degree growth, and a strict final arithmetic inequality.

## 7. Many Worlds audit

No analogous mathematical sign error was found in the supplied Many Worlds text. The following claims are sound in their stated idealized setting.

1. **Generator.** If a positive stage factor \(G_k\) satisfies \(G_k^k=k\), then \(G_k=k^{1/k}\) and \(g_k=G_k-1\). Its log coordinate is \(\ell_k=\log k/k\). The financial, exponential-growth, constant-relative-fitness, and independent-gain-cascade readings are parameterizations of this same equation, not empirical equivalences between those systems.

2. **Pair signs.** For even \(n\), \(P_n=g_n-g_{n+1}\) is negative at 2 and positive at every \(n\ge4\), since \(g'(x)=x^{1/x}(1-\log x)/x^2\) and the integer maximum is at 3. A mean-value argument gives \(P_n\sim(\log n-1)/n^2\), so the paired series is absolutely summable. In particular, positivity of \(C\) need not rest on printed digits: its first two pairs equal \(2\sqrt2-\sqrt[3]3-\sqrt[5]5>0.0064\), certified by rational root intervals, and every subsequent pair is positive.

3. **Signed measure.** \(p_n=P_n/C\) is a finite signed measure of total mass one with exactly one negative atom, not an ordinary probability measure. Its negative mass is \((g_3-g_2)/C\) and total variation is \(1+2(g_3-g_2)/C\). This is a mathematical normalization, not evidence of negative physical probabilities.

4. **Logarithmic remainder.** With \(d_k=e^{\ell_k}-1-\ell_k\),
   \[
   \frac{(\log k)^2}{2k^2}\le d_k
   \le\frac{e^{1/e}(\log k)^2}{2k^2}.
   \]
   This follows from the integral Taylor remainder for \(e^x\). Hence \(\sum d_k<\infty\), and the displayed tail bound
   \[
   \sum_{k>N}d_k\le\frac{e^{1/e}}{2N}
       ((\log N)^2+2\log N+2)
   \]
   is valid for the stated \(N\ge8\) (indeed for integer \(N\ge3\)).

5. **Eta derivative.** Differentiating the eta series on \(\operatorname{Re}s>0\), justified by locally uniform summation by parts, yields \(\sum(-1)^k\log k/k=\eta'(1)\). Expanding \((1-2^{1-s})\zeta(s)\) at 1 gives
   \[
   \eta'(1)=\gamma\log2-\tfrac12(\log2)^2,
   \quad C=\eta'(1)+\sum_{k\ge1}(-1)^kd_k.
   \]
   The correction sum is absolutely convergent. None of this determines rationality of either term or of their sum.

6. **Decay layers.** With \(t=\log k\), \(g_k=\exp(te^{-t})-1\); subtracting \(te^{-t}\) leaves leading term \(t^2e^{-2t}/2\). These are polynomially weighted exponential layers, exactly as the text's formulas show.

7. **Optical model.** The equal-rate two-step reaction with zero initial intermediate gives \([B]/A_0=ue^{-u}\). Under a Beer-Lambert response \(I/I_0=\exp(-\alpha ue^{-u})\), calibration \(\alpha=1\), and sampling \(u_k=\log k\), reciprocal transmission minus one is exactly \(g_k\). This is a valid idealized generator, conditional on the stated kinetic, absorption, normalization, and sampling assumptions. Real detector errors, absorption from other species, drift, and finite precision are not analyzed. A uniform per-sample error allowance by itself gives no uniform error bound for an infinite accumulation.

The paper already disclaims a universal physical constant and says that alternating aggregation is externally imposed. Preserve those caveats. The cross-domain invariant is the indexed numerical sequence plus its specified alternating functional; the map forgets units, dynamics, feasibility, uncertainty, and mechanism-specific restrictions. The terms' shared formula supplies no bridge to rationality or to unrelated open mathematical problems.

## 8. Reproduction and limitations

Run from the workspace root:

```
python3 session31-square-blocks/verify_square_blocks.py
```

The standard-library path checks, using exact integer/rational arithmetic:

- 640 moments (radii 1-40, degrees 1-16), with Euler polynomials independently generated from \(E_j(x+1)+E_j(x)=2x^j\);
- every piecewise-linear kernel interval for radii 1-40, its positivity endpoints, and its mass and quadratic moment;
- rational root enclosures for the first-block counterexample, positivity of the first two even-odd pairs, and the Galois hostile's magnitude separation;
- rational upper bounds for curvature tails, based on a positive log-series remainder and outward decimal rounding.

Optional mpmath checks are separately labeled FLOAT DIAGNOSTICS. They are not used as proof or as interval certificates. The finite structural universe supplements, and does not replace, the all-radius proofs in this note.

The rendered source pages are `document24-page21.png` and `many-worlds-page12.png`. The exact output and the optional diagnostics are `verification.json`.

Repository startup instructions were read. The read-only coordination check could not inspect messages because `.machine-id` is absent; no identity was invented. No repository files or Git state were mutated by this worker. Publication, theorem-number reservations, and any independent final audit are left to the parent task.
