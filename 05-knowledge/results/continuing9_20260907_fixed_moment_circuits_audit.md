# Independent audit: circuit flexibility at fixed first two moments

**Status: PROVED / independently audited.** The all-degree separated center, exact fixed-moment coefficient chart, uniform whole-box radius, every-ternary-word consequence, complete cubic image and hostile controls pass independent mathematical review. No producer correction was requested. This theorem concerns positive root geometry; the actual Laurent C/D interlacers, factorial carrier and original phase remain separate constraints.

The reviewed pre-promotion report is `continuing9_20260907_fixed_moment_circuits.md`, candidate SHA256 `dd6a3579fce3ff539a540051e90741744798bebde7dd8bc202a0a14f7e2a5810`. The independent source pins the producer source, output and certificate respectively:

- `a71ee4509bb46644f8e1fcf3a5f6931941363876540ea8565d1addf46e6e8102`;
- `5d51e2a2aa33a1dd55d4d3caaa1c64702567ca4787e78a8b5d1e8c4f0486abfe`;
- `291b62a5638f8a057ccf514851cb0b5d100a2a61a3da8eff9151b49040be11b9`.

No producer was imported or executed by the referee verifier. The primary source was read only after the independent implementation passed.

## 1. Root induction and moment normalization

For `F_d(x)=sum (-1)^k binom(d,k) q^(k(k-1)/2)x^k`, Pascal's identity gives `F_(d+1)=F_d-xF_d(qx)`. The root-separation induction is valid at both ends, not only between interior roots.

If `a_(i+1)>a_i/q`, then `q a_i` lies strictly in `(a_(i-1),a_i)`, with a_0=0. Therefore `F_(d+1)(a_i)` has sign `(-1)^i`. At `a_i/q`, the shifted term vanishes and `F_d` has the same sign `(-1)^i`, including after the final root. It follows that the successive new roots lie in `(0,a_1)`, `(a_(i-1)/q,a_i)` for the interior indices, and `(a_d/q,infinity)`. These are disjoint and yield the full degree. Moreover each new root b_i is less than a_i, while b_(i+1) exceeds a_i/q; hence `b_(i+1)>b_i/q`. The strict induction starts from `F_1=1-x` and proves positivity, simplicity and completeness of the real roots.

For fixed `S>0` and `S^2/d<Q<S^2`, the chosen

`q=d(S^2-Q)/((d-1)S^2)`

lies strictly between zero and one. With `mu=S/d`, reciprocal reversal gives monic coefficients `e_k=binom(d,k) mu^k q^(k(k-1)/2)`. The coefficients e_1 and e_2 recover exactly S and Q, and every normalized Newton ratio is 1/q. Thus the center has every circuit quotient equal to one. The reciprocal operation reverses the root ordering but preserves the asserted positive-parameter property.

The boundary distinctions are correct. At q=1 the center has repeated equal parameters and cannot supply an open simple-root chart. At q=0 it loses positive degree and lies outside the strictly positive-root moment slice. The explicit box is stated for d>=3, so its denominator `binom(d,3)` is nonzero; the center alone is also valid at d=2.

## 2. Exact chart and the entire neighbourhood

At fixed q, the inverse ratio recursion fixes h_0=h_1=1 and R_1=1/q, then uses `R_k=R_(k-1)c_k` and `h_(k+1)=h_k^2/(R_k h_(k-1))`. It agrees with the producer's closed formula. In particular h_2=q never changes; all target first and second moments remain exact. Taking ratios of the reconstructed coefficients recovers the requested c_j, so the map is invertible on positive coefficient tails.

For coefficient k, the total exponent of the circuit coordinates is `E_k=binom(k,3)`. Consequently every point of the whole closed cube `|c_j-1|<=delta` satisfies

`(1+delta)^(-E_k)<=a_k/a_k^0<=(1-delta)^(-E_k)`.

For `E_k delta<=1/2`, Bernoulli's inequality gives the uniform deviation bound `2E_k delta`. This is not a vertex-only statement. The first three coefficients have E_k=0 and do not move. Summing the absolute coefficient errors at a sign address gives `2delta W_l`; choosing `delta<=M_l/(4W_l)` preserves at least half the original sign margin. The positive constant coefficient and the d alternating sample signs then certify all d simple negative roots by degree.

Every W_l is strictly positive in the stated degree range because the degree-three coefficient and the positive address contribute. Thus the proposed minimum bound is positive. For real moment data the bound is a positive real number; a smaller positive rational radius always exists. For rational moment data, rational sign samples make the displayed bound itself rational and effectively checkable by exact isolation. The theorem claims rational polynomial coefficients for rational input, **not rational or integer root parameters**. Clearing denominators of a polynomial also does not turn its roots into integers.

All ternary words follow by choosing `c_j=1+delta sigma_j`. Zero circuit signs are exact quotient identities; they do not require repeated roots. The statement that the circuit image is open around any simple positive-root point is also justified: the coefficient chart is a homeomorphism on positive tails, and simple negative-root polynomials form an open subset after restricting to the fixed first two coefficients.

## 3. Cubic amplitude image and the hostiles

The normalized cubic is `t^3-3t^2+3(1-v)t-p`, whose discriminant is exactly `27[4v^3-(p-(1-3v))^2]`. With `a=sqrt(v)`, the boundary products are `L=(1+a)^2(1-2a)` and `U=(1-a)^2(1+2a)`. The discriminant interval supplies three real roots; product p>0 and the positive pair sum exclude a negative pair, so all roots are positive precisely on the stated strict interval `max(0,L)<p<U`. The endpoint factorizations identify exactly when the repeated lower endpoint remains positive.

Since `C_2=q^3/p`, reversing the positive product interval yields the complete quotient interval. The upper endpoint is finite exactly when a<1/2, equivalent to `Q<S^2/2`; equality belongs to the unbounded case, with zero product excluded. The lower quotient is strictly below one and the upper quotient, when finite, strictly above one. The Newton-only lower bound is strictly weaker.

The original cubic hostile has all Newton ratios strictly above one but negative discriminant `-25515/65536`. The quartic obtained by multiplying it by n+1 retains a nonreal conjugate pair despite all strict Newton inequalities. These are valid counterexamples to Newton sufficiency, while the separated-center theorem gives a genuine open sign box on the same moment slices.

At the actual degree-five moments `S=13,Q=59`, the radius `1/2048` cube is verified. The inherited magnitude obstruction `C_2=1/2` remains outside that cube and violates R_2>=1. Thus unrestricted sign flexibility is compatible with amplitude restrictions. No conclusion about the actual C/D interlacer domain follows from this geometric realization alone.

## 4. Independent exact verification paths

The referee uses standard-library rational polynomial division and Sturm chains for root counts and fresh dyadic isolation. SymPy is used only for a few symbolic cubic identities, not for the root certificates. The coefficient reconstruction uses the cumulative-ratio recursion rather than the primary closed exponent product.

- Eighteen fresh centers, degrees 2 through 7 at q in `{1/4,3/5,9/10}`, are independently isolated. Their full q separation and the actual successive-degree interval ordering are checked with rational endpoints.
- Every one of the 25 retained root intervals is checked by a separate Sturm count, including exact rational root intervals. The supplied intervals also rigorously certify q separation.
- All six frozen boxes are reconstructed. Their sample values, weights and exact proposed rational radii match.
- A second whole-box certificate directly encloses every coefficient using its exact lower and upper powers of `1+delta` and `1-delta`, then sums the signed coefficient intervals at each address. These independent interval sums preserve every sample sign with at least half its original margin. This path does not reuse the Bernoulli estimate and certifies every point of each closed box.
- All 132 ternary target polynomials are regenerated, with complete word-universe checks, exact moments and ratios, exact ties, positive simple Sturm root counts and literal sample signs.
- The anchor is separately scaled to S=13 and checked to have square sum 59. The full cubic discriminant and both endpoint factorizations are reconstructed symbolically; additional interior and exterior product values are tested by rational Sturm counts. Both original Newton hostiles are checked directly.

The six box records `(d,q,delta,word count)` are:

`(3,3/4,1/32,3)`, `(4,7/8,1/1024,9)`, `(5,275/338,1/2048,27)`,

`(6,1/2,1/32768,81)`, `(3,1/4,1/256,3)`, `(4,9/10,1/1024,9)`.

All finite controls support the separately reviewed universal proofs; they are not used to infer an untested all-degree statement.

## 5. Literature and frozen reproduction

The finite induction is self-contained. The classical attribution is consistent with primary sources: Wang and Zhang's introductory definition is exactly the coefficient sequence used here. [Wang--Zhang, *Zeros of the deformed exponential function*, Section 1](https://arxiv.org/pdf/1709.04357)

The Saff author bibliography confirms the Iserles--Norsett--Saff article and publication details; its PDF timed out during this referee's fetch, so the particular printed-page reading is not independently claimed here. A separate author-hosted primary survey explicitly records the Gaussian sequence `exp(-a k^2)` as a Laguerre sequence for every a>0, consistent with the stated multiplier interpretation. [Craven--Csordas, printed p. 27 after Definition 4.11](https://math.hawaii.edu/~tom/mathfiles/czdssurvey.pdf)

Run `python3 -B continuing9_20260907_fixed_moment_circuits_audit.py` and the same command with `-O`. The verifier finds frozen inputs beside itself, with the filed results-directory fallback for outputs. Both runs pass **909 always-active exact gates** and produce byte-identical actual LF stdout. No producer source or certificate was modified. The parent owns status promotion and filing.

Frozen referee source SHA256: `39fa4c64f4b6f915bc0e2c17d64321b6850fe6d97c3438e73f4b1bac63b122ab`.

Frozen normal and optimized stdout SHA256: `9903e97b1f739802abada3d0096c5b869efe7609e1b242f6c07f506257443171`.
