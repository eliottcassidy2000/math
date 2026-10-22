# Independent structural audit of the Session31 Collatz and square-block repairs

Date: 2026-10-22. Scope: `session31-collatz/AUDIT.md` with `verify_audit.py`, and `session31-square-blocks/AUDIT.md` with `verify_square_blocks.py`. This review independently derives the all-parameter constructions, kernel and tail bounds, and algebraic-conjugation obstruction. Replaying the submitted programs is a secondary reproducibility check, not the basis of this verdict.

## Verdict and requested precision fixes

**PASS on the scoped mathematical conclusions.** No blocking mathematical error was found in the Collatz finite-word counterexamples, the repaired square-block identity, the positive Peano kernel, the rational tail certificate, or the existence and separation of the stated Galois conjugate.

One wording fix is needed in the square-block audit's field argument: define **K to be the compositum** of the splitting fields of `X^11-11`, `X^12-12`, and the cyclotomic field `Q(zeta_13)`. Merely saying that K *contains* those fields does not imply the degree bound: a larger K could already contain `13^(1/13)`. The intended compositum has exactly the needed prime-to-13 degree property; Section 5 below supplies the full argument.

A nonblocking implementation recommendation is to replace the square-block checker's validation `assert` statements by explicit failure checks if optimized Python execution is to be accepted. Its advertised normal-Python run is valid. The independent checker here uses no assertion gates and produces identical receipts under `python` and `python -O`.

These conclusions do not prove Collatz convergence, any native-energy gate, MRB irrationality, or uniform smallness of all algebraic conjugates. Arbitrarily long finite orbit prefixes are not an infinite divergent orbit.

## 1. Source authentication and page references

The decisive passages were checked against the supplied PDF files, not just against the earlier audits.

- `document24.pdf` is a four-page extract with printed pages **20, 21, 22, 23**. Thus the sign-error citation to printed page 21 is **PDF page 2**, not PDF page 21. The original audit correctly says *printed* page 21. The independent rendering `render/document24-pdf2-printed21.png` visibly confirms the erroneous factor `(-1)^(n+d)` in both the block and whole-constant formulas. It also confirms that the later normalized Euler expansion has the correct convention.
- The Collatz PDF has 31 pages. Printed pages and PDF pages agree for the decisive passages. Independent renders of pages 9 and 31 confirm, respectively, the false repeated-halving claim and the immediate-rebound prohibition / approximately-five-to-six-rebound assertion. Text from PDF pages 14 and 18 confirms the strict envelope and claimed cycle contraction. These are not OCR-created errors.
- `documentMRBconnections.pdf` has 15 pages. The independently rendered PDF/printed page 12 confirms the positive sign of `eta'(1)` and the formula `eta'(1) = gamma log 2 - (log 2)^2/2`. No sign error occurs there.

Only the source-page claims relevant to this review are being certified. There is no claim of independent historical-priority verification.

## 2. Collatz: independent all-parameter derivations

Throughout, the ordinary map is `C(n)=n/2` for even n and `C(n)=3n+1` for odd n. For an odd x, let its next odd successor be `(3x+1)/2^e`, where e is the exact 2-adic valuation of `3x+1`.

### 2.1 Run length and terminal halving length are independent coordinates

Write `x=a 2^r-1` with a positive odd and r positive. Iterating the half-step while its input remains odd gives

    T^j(x) = a 3^j 2^(r-j)-1.

For `0 <= j < r`, this is odd; at j=r it is even. Consequently there are r odd updates in the first maximal run, the last raw peak is `2(a3^r-1)`, and, if `s=v_2(a3^r-1)`, the post-cascade endpoint is

    y = (a3^r-1)/2^s.

The first `r-1` odd valuations are 1, and the last is `s+1`. Thus calling r the peak valuation would be an actual coordinate error.

For arbitrary positive r,s,t, choose a positive odd b in the class

    2^(s+t) b = 2^s-1 modulo 3^r.

Such b exists because 2 is invertible modulo `3^r`, and adding `3^r` switches parity. Set

    a = (2^(s+t)b-2^s+1)/3^r.

Its numerator is positive and odd, so a is a positive odd integer. The identity

    a3^r-1 = 2^s(2^t b-1)

then gives terminal valuation `s+1` exactly, endpoint `y=2^t b-1`, and next run length t exactly. There is no hidden parity, positivity, or compatibility condition on the positive triple. The possible b=0 residue causes no problem: select a positive odd representative. This confirms the audited theorem for every triple, not merely its tested cube.

### 2.2 Independent derivation using binary valuation words

An alternative proof uses no congruence modulo a power of 3. Let `e_1,...,e_L` be any prescribed positive finite word, set `E_j=e_1+...+e_j`, and define

    B_0=0, E_0=0,
    B_(j+1)=3B_j+2^E_j.

The formal composition has

    2^E_j x_j = 3^j x_0+B_j.

Choose x_0 in the residue class

    x_0 = (2^E_L-B_L) (3^L)^(-1) modulo 2^(E_L+1).

For every j<L, reducing the final numerator modulo `2^(E_j+1)` shows

    3^j x_0+B_j = 2^E_j modulo 2^(E_j+1).

Indeed, in the remainder of the composition the first new summand has valuation E_j and all later summands have larger valuation. The same congruence at j=L is the defining choice. Hence every x_j is an odd integer and every prescribed e_j is exact. Taking a positive representative makes the entire finite orbit segment positive.

Apply this to the concatenated word

    (1 repeated r-1 times), s+1, (1 repeated t-1 times), 2.

It realizes precisely the same three independent parameters, and fixes the next run's terminal event as well. This provides a genuinely different construction used by the independent checker. It gives an arithmetic progression for each fixed finite word; it asserts no uniform density as word length grows and no positive-integer realization of arbitrary infinite words.

### 2.3 Regeneration and long growing prefixes

For odd r>=3 and even k>=2, `y=2^r-1` is 1 modulo 3, so

    x=(2^k y-1)/3

is a positive odd integer and `3x+1=2^k y`. The k-halving cascade lands at y immediately, which already has r trailing ones. Both parameters are unbounded. This refutes the appendix's mandatory post-cascade infertility, not convergence itself.

For any odd q>0 and m>=1, define

    x_j=q 9^j 2^(3(m-j)+1)-5,  0<=j<=m.

For j<m, `x_j=11 mod16`. The next two odd valuations are therefore `(1,2)` and the completed run/cascade map is

    x_(j+1)=(9x_j+5)/8,
    x_(j+1)-x_j=(x_j+5)/8>0.

The first intermediate odd value `(3x_j+1)/2` exceeds x_j. Every raw even value and every halving intermediate is at least the final odd endpoint or one of the preceding odd endpoints. Thus the entire first m completed cycles stay at or above x_0 and the endpoints strictly grow. The audit's q=1 family is a valid special case.

Separately, for `N=2^n-2`, `n=3 mod4`, one has `3^(n-1)=9 mod16`, hence exact terminal peak valuation 4 independently of the unbounded run length. Its endpoint `(3^(n-1)-1)/8` divided by N tends to infinity on this progression. The displayed 2046 witness is therefore part of a genuine all-order noncontraction family.

None of these facts depends on how much of the global Collatz conjecture has been computationally checked. If the manuscript changes the definition of a cycle to mean a segment ending only after a desired descent, attaining that stopping condition is exactly the missing theorem.

## 3. Square blocks and an independently paired positive kernel

Because `c_n=n(n+1)` is even, the actual block is

    B_n=sum_(d=-n)^n (-1)^d f(c_n+d).

The center's total coefficient after subtracting the second differences is `(-1)^n`, giving

    B_n=(-1)^n f(c_n)+sum_(d=1)^n (-1)^d Lambda_(n,d).

There is no extra `(-1)^n` on the curvature sum. The n=1 discrepancy is exact, not a numerical diagnostic.

For general F, radius r and center c, subtracting its center and integrating each symmetric second difference gives

    T_r F(c)-F(c)
      = integral_0^r K_r(u)[F''(c-u)+F''(c+u)]du,

where `T_r F(c)=sum_(d=-r)^r (-1)^(r+d)F(c+d)` and

    K_r(u)=sum_(d=1)^r (-1)^(r+d)(d-u)_+.

Here is a positivity proof by direct pairing, rather than the audit's floor/fractional-part calculation:

    K_(2h)(u)=sum_(j=1)^h [(2j-u)_+-(2j-1-u)_+],
    K_(2h+1)(u)=(1-u)_+
                +sum_(j=1)^h [(2j+1-u)_+-(2j-u)_+].

Every bracket is nonnegative. At integer knots,

    K_r(j)=ceil((r-j)/2),  0<=j<=r,

and K is linear between knots. This is the representation used for independent exact integration.

The mass and second moment follow by integrating each truncated ramp:

    2 integral_0^r K_r(u)du = sum_(d=1)^r (-1)^(r+d)d^2
                            = r(r+1)/2,
    integral_0^r u^2 K_r(u)du
      = (1/12)sum_(d=1)^r (-1)^(r+d)d^4
      = r(r+1)(r(r+1)-1)/24.

The first identity yields the claimed mean-value theorem and C^2 bound. For C^4 functions,

    |F''(c-u)+F''(c+u)-2F''(c)| <= u^2 sup|F^(4)|,

so the second identity gives exactly the proposed fourth-derivative remainder constant. There is no missing factor 2.

The Euler moments, the invariant-radius qualification, and the prohibition on interchanging the two infinite sums are consistent. The positive kernel is a statement about the normalized signed functional, not about each signed point mass being positive.

## 4. Curvature asymptotic and rational tail constant

For `f(x)=exp(log(x)/x)-1`, differentiating directly gives the displayed f''. With `g=log(x)/x`, one has, for x>=9,

    g'<0, g''>0, g'''<0,
    f''=exp(g)[g''+(g')^2]>0,
    f'''=exp(g)[g'''+3g'g''+(g')^3]<0.

Thus f'' is positive and decreasing on every square block with n>=3. Since the normalized curvature is `(-1)^n R_n`, kernel mass gives

    c_n f''(n^2+2n)/2 <= (-1)^n R_n <= c_n f''(n^2)/2.

Both endpoints are asymptotic to `2 log(n)/n^4`, establishing the stated signed asymptotic. Positivity is not being incorrectly exported to the first block.

Let A=exp(1/e). Bounding `exp(log(x)/x)<=A` and substituting x=n^2 yields

    |R_n| <= 2A(1+1/n)[log(n)/n^4+(log(n))^2/n^6].

The two summands are decreasing for n>=3. On n>N the factor is at most `1+1/(N+1)`, and their integrals from N to infinity are exactly I_4(N) and I_6(N) in the audited formula. Therefore the tail formula and its strict replacement `2A<3` are correct.

The independent checker certifies the four displayed decimal bounds by a different transcendental enclosure: it proves a rational interval for log(10) using lower and upper rational Taylor bounds for exp, then uses `log(10^j)=j log(10)`. It does not import the earlier logarithm routine. All four outward decimals agree exactly:

| N | Certified absolute curvature-tail upper bound |
|---:|---:|
| 10 | 0.002916802930051020125517 |
| 100 | 0.000004988801113793516135 |
| 1000 | 0.000000007248352833489898 |
| 10000 | 0.000000000009544628508874 |

These certify the curvature tail only, not the full MRB truncation error. Absolute convergence of the completed square blocks and of the triangular curvature array is compatible with divergence of the sum of the absolute original individual terms. The audited text keeps this distinction correctly.

## 5. The Galois conjugate exists: complete field argument

Let L_11 and L_12 be the splitting fields over Q of `X^11-11` and `X^12-12`, and define specifically

    K=L_11 L_12 Q(zeta_13).

For general n, the splitting field of `X^n-n` is obtained by first adjoining `zeta_n` and then one nth root of n. Over `Q(zeta_n)` its Galois group embeds into the group of nth roots of unity, so its degree divides n; consequently its degree over Q divides `n phi(n)`. Therefore

    [K:Q] divides (11*10)(12*4)*12 = 63360,

which is prime to 13. This is a divisibility argument, not merely an upper bound.

Set `alpha=13^(1/13)`. Eisenstein at 13 shows `[Q(alpha):Q]=13`. If alpha belonged to K, the tower law would force 13 to divide `[K:Q]`, a contradiction. K already contains every thirteenth root of unity, so `K(alpha)/K` is a splitting extension of `X^13-13`, with Galois group a nontrivial subgroup of the cyclic group of order 13. Its group is consequently the entire group, and the automorphism

    alpha -> zeta_13^6 alpha

exists while fixing both `11^(1/11)` and `12^(1/12)`. Its action on the displayed second difference is a genuine algebraic conjugation.

A simpler exact magnitude certificate than the earlier numerical-looking rational calculation is available:

- The first nine terms (through degree 8) of the positive exponential series at 12/5 sum to more than 11. Hence `log(11)<12/5`.
- `(5/4)^11>11`, hence `11^(1/11)<5/4`.
- Substitution in the positive derivative terms gives

      f''(11) < (5/4)[(2*(12/5)-3)/11^3+((12/5)-1)^2/11^4]
               = 136/73205 < 1/500.

Thus the real second difference obeys `0<Lambda<1/500`. Also `13^(1/13)>6/5` by the exact comparison `6^13<13*5^13`, and

    cos(pi/13)>cos(pi/6)=sqrt(3)/2>4/5.

The transformed expression has real part strictly less than

    1/500-(9/5)(6/5)=-1079/500,

so its modulus exceeds 1079/500. Every comparison has an exact rational certificate except the elementary monotonicity of cosine on this interval. This proves the asserted large separation and does not presume a floating-point embedding.

The example refutes automatic preservation of analytic smallness by arbitrary conjugation. It does not refute the more qualified possibility of a specially selected coherent family with independently proved arithmetic bounds.

## 6. Independent finite checks and reproduction

`verify_independent.py` uses only the Python standard library, integer arithmetic and Fraction. It does not import either audited verifier. It performs **46,965 explicit non-assertion validation checks**, including:

- all 4,096 triples `1<=r,s,t<=16`, separately realized by a binary valuation-word construction and checked against the original modulo-3 construction;
- five larger/asymmetric words, including `(20,20,20)`, `(80,1,80)`, and `(32,64,96)`;
- 240 complete growing prefixes, with `m=1..80` and odd scale q in `{1,3,17}`, checking every intermediate ordinary Collatz value for premature crossing of the original baseline;
- 992 inverse-regeneration cases and 41 constant-four-halving cases;
- an independent memoized census for starts 1 through 10,000, reproducing `(3997,4011,1579)`;
- nonnegative kernel knots and exact interpolation through radius 64, polynomial Peano identities of degrees 0 through 12 at the nonintegral center 2/3, mass and fourth-order constants, and the first five even moments;
- all four printed tail decimals and the alternative Galois magnitude certificate.

Run from the workspace root:

    python session31-structural-audit/verify_independent.py
    python -O session31-structural-audit/verify_independent.py

Normal and optimized receipts agree byte for byte. `independent_receipt.normal.json` is the preserved normal receipt; `independent_receipt.json` is the last generated receipt. The receipt records SHA-256 hashes of the four reviewed source artifacts, so subsequent edits do not silently inherit this review.

Secondary replays of the original normal-Python checkers also match their original receipts byte for byte. The square-block program was copied to `replay_square/` before running so that its fixed output path could not overwrite another worker's files. The Collatz program's output option was redirected to this audit directory. No other worker's file, repository file, or Git state was changed.

## Stopping boundary

The requested independent hostile audit is complete. The only requested mathematical-text amendment is the explicit choice of the compositum K. The revised all-order statements stand on the proofs above; finite test counts neither replace those proofs nor establish any of the remaining open global claims.
