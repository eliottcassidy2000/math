# Independent audit of coupled windows and their derivative cone

**Status: PASS — written proof independently audited; independent exact
controls PASS.** Scope is the theorem in
[creative_20260906_laurent_bridge.md](creative_20260906_laurent_bridge.md),
not general trinomial doubled-row separation or a Laurent bound outside its
real-rooted common-carrier assumptions. The reviewer is the separate Smith
bridge agent. No producer module was imported by the audit computation.

## Analytic review

The actual window is
`D^r rev_(n-s) D^s rev_n H`. Its coefficients are exactly
`(j)_r (n-j)_s h_j` at exponent `j-r`. The selected coefficient remains
zero at `k-r`; literal squaring gives the displayed weighted opposite-pair
sum. The ambient reversal degree is essential and is correctly retained.

All derivatives and reversals preserve real-rootedness. A reversal can lower
the actual degree by removing a zero root, so real-rootedness alone would
not justify the strict interior condition needed by
[THM-4440 / signed duplication SOS and real-rooted Laurent return](../../01-canon/theorems/THM-4440-signed-duplication-sos-and-real-rooted-laurent-return.md).
The note supplies the missing argument: a real-rooted derivative has a
multiple zero only when its parent has a zero of higher multiplicity.
Otherwise the logarithmic derivative identity at the critical point gives
strictly negative second logarithmic derivative. Iteration rules out two
consecutive zero coefficients when `h_0!=0`. Thus the coefficients immediately
beside `h_k` are nonzero; both survive precisely when `r<k` and `s<n-k`.
This proves the required `ord<k-r<deg` despite possible endpoint degree loss.
At either next derivative boundary, the selected squared coefficient is zero.
Repeated real roots are allowed throughout.

The centered-square identity is valid on the actual paired coefficient
support, where all factorial arguments are nonnegative. Outside that support
the original product is zero, and no artificial coefficient is introduced.
After reflection the grid is exactly `d=1,...,M`, with `M=min(k,n-k)`;
the central `d=0` term vanishes because `h_k=0`.

For the cone proof, before the factor at stage `b`, every active index
satisfies `t>=max(r,M-N+b)`. This makes
`(N-b)^2-(M-t)^2>=0`. An equality kills the coefficient that could have
violated the next stage's support bound, while the raised index remains
legal. The discarded `B_M` is zero on the entire selected grid, so this
step is an equality in functions on that finite grid, not an equality of
unrestricted polynomials. The final vector cannot vanish since every legal
window has `W(1)>0`. The evaluation matrix at `d=M,M-1,...,1` is triangular
with positive diagonal, proving independence, unique expansion, and exactly
`M` extreme rays. Conversely each generator is an allowed one-sided window,
so the cone equality is both directions.

The note correctly makes cone membership necessary and sufficient only for
that generated cone, and sufficient for universal strict negativity. It
makes no unsupported maximality claim about all possible weights. The
nonnegative outside-cone weight selecting the positive outer pair is a valid
hostile, as is the zero-preserving Euler operation that destroys real roots.

One wording correction was requested and accepted: a coupled window does not
literally preserve opposite-coefficient products; it retains them as explicit
weighted summands. Factoring a zero power and shifting the selected index
preserves the original products. No theorem formula required correction.

## Independent exact verification

The separate source
[creative_20260906_laurent_audit.py](../../04-computation/creative_20260906_laurent_audit.py)
constructs polynomial coefficients through elementary-symmetric subset sums,
using prefix factor multisets from `(-2,-1,1,3)` at degrees `2,...,7` and
an exactly tuned nonzero rational final factor. This differs from the
producer's polynomial-product universe and includes the degree-two boundary.
All roots are real by construction, with nonzero first and last coefficients.
No root-finding approximation is used.

It checks 938 carriers and all 7,083 legal windows by the independent diagonal
coefficient formula and direct convolution, including surviving opposite-sign
neighbors and the actual strict quadratic consequence. For the cone it uses
full rational Gauss-Jordan interpolation on the evaluation matrix, rather
than the producer's constructive recurrence or triangular solver. Every
legal `r,s` is checked for `1<=M<=7`, `M<=N<=14`: 2,618 independent cone rows.
Every recovered coefficient is nonnegative, at least one is nonzero, and
re-evaluation reproduces all grid values. The outside-cone hostile has exact
coordinates `(1,-1/3)`.

There are 38,498 exact gates. The finite tests support the analytic proof;
only the analytic argument proves the unbounded quantifiers. Reproduce with

```text
python3 -B 04-computation/creative_20260906_laurent_audit.py
python3 -B -O 04-computation/creative_20260906_laurent_audit.py
```

Matching output: [creative_20260906_laurent_audit.out](creative_20260906_laurent_audit.out).
Normal and optimized output streams agree byte for byte. The audit contains
no Python assertion statements or floating-point arithmetic.

```text
source_sha256:
107e4c6146f46b26c970d5bb3cd619b8b861ea0f0718cb841a06fa71270633a5
output_sha256:
af509391f28261d238e1f4f9ef33ecfefc5ba9949ff990385a3e1bdf93e5ba84
semantic_sha256:
846413ff92fb5278b86233291923975e8b2da4aa88fdb904d17054bb4a9edfa4
```
