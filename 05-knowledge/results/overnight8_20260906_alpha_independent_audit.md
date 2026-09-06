# Independent audit of complete alpha parity in one carrier

**PASS: analytic all-parameter proof and independent exact controls.** The
[alpha-completion theorem](overnight8_20260906_alpha_completion.md) is
accepted for every stated A=2 full-path parameter, without a coprimality
assumption. The A2B3 mixed-product reduction and its two coupled Euler
identities are also accepted. The beta-skip sign remains **OPEN** outside
the explicitly finite banks.

## Preserved source and endpoint audit

The complete beta polynomial has degree `q=L+h`, with strictly positive
constant and leading coefficients. The inherited full-path PF theorem
applies to the entire index interval `-L,...,h`; neither gcd(A,B)=1 nor
a positive-charge consumer is required for that analytic theorem. Root
checked this domain in the incoming Hadamard report and
[THM-4436 / complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall](../../01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md).
The negative-index prefix is essential, and remains present in G and
every convolution before alpha overlap.

For t<0, each negative root lambda of G gives the real pair
`u=+-sqrt(t/lambda)` in `H_t=u^(2q)(1+u)^m G(t/u^2)`. Its exact constant
term is beta_h*t^q, which is nonzero, and its degree is m+2q. Consequently
there is no unaccounted root at zero. For `k=2h+z`, the upper gap is
`m+2q-k=x+y+2L>0`; h>=1 makes k>0. This checks the actual interior
coefficient hypotheses of
[THM-4440 / signed-duplication-sos-and-real-rooted-laurent-return](../../01-canon/theorems/THM-4440-signed-duplication-sos-and-real-rooted-laurent-return.md).
Repeated roots, including collisions with -1, are allowed there.

The first extraction has binomial index
`k-2q+2(j+L)=z+2j`, so it is exactly t^L P(t).
The squared extraction has index
`2k-4q+2(j+2L)=2z+2j`, so it is exactly t^(2L) W(t).
This includes j=-1 when z=1 and restores the whole alpha parity class.
Thus signed duplication acts on the original P-zero. No replacement
factor is assumed to retain that zero. Its strict negative conclusion is
preserved by division by the positive t^(2L).

The quantitative subset-square bound has denominator `2k-1=4h+2z-1`.
After removing t^(2L), its constant factor is beta_h^2*t^(2h), as stated.
It is inherited from the same ordinary carrier and is not the sharp 7/9
constant of the separate four-root stability theorem.

## Exact remaining mixed response

For A2B3 the seventh full-path identity gives actual row
`Q_raw=W+2 A_double star(tCD)`. This proves a strict negative W uniformly.
It neither proves G1 negative separately nor makes beta-skip negativity
necessary. The exact target is `R_skip<-W`; `R_skip<=0` is sufficient.

The shifts `t^x B=F_(3q)`, `t^x C=F_(3q-1)`, and
`t^x D=F_(3q-2)` give degrees q,q-1,q-1, respectively. These degree
differences produce the factor `2t^(1-2x)` in the mixed coefficient
formula. It is negative at t<0, so the open sufficient ordinary-coefficient
inequality has direction `>=0`, exactly as written. Transforming the two
coefficientwise Euler identities with these different degrees yields both
displayed coupled u-derivative equations. The report preserves coupling;
arbitrary pairs of real-rooted polynomials are outside its proposed target.

## Independent controls and scope

The independent script imports no producer mathematics. It constructs full
binomial arrays, forms the ordinary polynomial and its literal square, and
checks both coefficient maps at 48 fresh general rows, including B=9,
even B=8, x=0 and both alpha residues. It checks the complete beta and
original first-root counts in these controls. Rootwise completed-alpha
negativity is independently certified by negative-definite Hermite trace
forms, rather than the producer's rational isolating-interval evaluation.
These controls do not replace the all-parameter proof above.

Eight additional A2B3 rows `(h,x)=(1,2),(2,3),(3,2),(4,6),(5,2),(6,6),
(7,2),(8,6)` check the two coupled identities as exact symbolic polynomials.
At every first root in these rows the mixed beta-skip response is strictly
negative, again by Hermite forms. Nonprimitive rows are deliberately
included: they are valid analytic controls, without a first-return claim.
This extends the finite hostile bank, not the universal sign theorem.

The independent checks retain `u^2+1` to reject a positive-phase pullback
and `u^3-1` to reject an arbitrary higher-power real-rooted pullback. The
new construction is specifically quadratic. The original `(h,g,rho)=(1,5,-2)`
factor-replacement hostile remains in the producer and the seventh audit.

Normal and optimized primary runs pass **5,638** gates with identical
output. The independent audit passes **446**, also identically in both
modes. Reproduce from the repository root:

```text
python -B 04-computation/overnight8_20260906_alpha_completion.py
python -B -O 04-computation/overnight8_20260906_alpha_completion.py
python -B 04-computation/overnight8_20260906_alpha_independent_audit.py
python -B -O 04-computation/overnight8_20260906_alpha_independent_audit.py

primary source f98fa964b8a39177051c85de33c59c45987a4d71817f31d7ccc7f8f367cbe300
primary output dd632e05fa7130a1dd82bfac3676d3619984e921e95c52315a8f5ee5dea6c6b0
audit source d8a2f6d927c2127318e13b9d541c1e914a3b152b654bab66b53f15c71392a95c
audit output 88c85e6acd003311f1ea61c2d7126be2355501d0bf410c447e0b86eaba81ee35
```

The proved advance is completion of one missing residue class while
preserving the original zero coefficient. Its consumer boundary is one
coupled mixed coefficient with an explicit negative square payment beside
it. General A and general actual doubled-row separation remain open.

**Filing:** root integrated these audited artifacts in the eighth checkpoint;
reproduction paths are relative to the repository root. Earlier outside-worktree
notes describe author provenance, not the present file location.
