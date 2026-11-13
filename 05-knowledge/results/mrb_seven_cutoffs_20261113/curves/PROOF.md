# MRB product thresholds, Gamma cutoffs, and a target-return guard

**Status:** PROVED elementary statements C1-C5 below; inherited results credited by file and theorem; FINITE-EXACT controls in `verification.json`; floating illustrations are expressly NUMERICAL. This is a research packet, not a new numbered canon theorem. Chan-Pak remains OPEN. No literature-priority claim is made for the elementary infinite-product principle.

## Inheritance and concept board

The closest proved mechanisms are the weighted wrapped-transitive tournament identity in [Session162, MRB weighted result](../../session162_mrb_square_20261103/mrb_weighted/RESULT.md), Gamma's dimension-shift recurrence in [Session127, VOL-1/VOL-2](../../session127_gsma_norm_collapse_20261031/volume/THEORY.md), and the core-forgetting tower in [THM-4900, gate/doubling/OCF](../../../../01-canon/theorems/THM-4900-gate-and-doubling-in-the-marked-resolvent-golden-towers-and-tower-ocf.md). The canonical hostile is a normalized curve that forgets its source while its absolute integer defect grows. The corrected near miss is MISTAKE-729: a tower error is `O(rho/r)`, not a multiple of a potentially zero cyclic-triangle difference. A least-used sidecar is the square-weight mass of the MRB path observable.

Recovered older MRB work also matters: [Session144, RESEARCH](../../session144_jacobian_eth_mrb_20261103/mrb/RESEARCH.md) proves that branch data, finite prefixes, and even all uniform right-half-plane algebraic asymptotics do not determine an alternating functional. It already classifies periodic-filter convergence by the zero Fourier mode. Those are inherited obstructions, not claims discovered here.

| Live concept | Preserved quantity / operation | Lost data and decisive test |
|---|---|---|
| Chan-Pak return | full odd integer `H`, physical marked root | relative convergence can lose a growing absolute defect; C1 |
| MRB carrier | exact marked path partition, exponent deformation | signed weights are not unweighted path counts; C2-C4 |
| Small-to-large cutoff | variance mass and higher even powers | a finite prefix cannot fix convergence; C3 |
| Ball-volume curve | metric radius and dimension | threshold depends on radius/normalization; C5 |
| `{3,10}`, `{4,11}` | specified affine coordinate or dimension pair | matching numbers do not transport the predicate; final map table |

The Anchor is the full-count return obstruction, the Niche is the MRB product threshold, and the Wildcard is the Gamma radius/dimension crossover. After each calculation, the board comparison is: which exact observable survives, and is its error absolute enough to decide an odd integer equality?

## C1. Relative forgetting can coexist with exponential target separation

**PROVED corollary of THM-4900(6)-(7).** Let `Gamma` be the strong gate `s -> T -> t -> s`; let `TT_d` be the transitive d-vertex tournament, `K1=TT_1`, and let `C3` be a cyclic triangle. For every integer `r>=0`,

`H(Gamma^r C3) - H(Gamma^r TT_3) = 2 H(Gamma^r K1) >= 2*3^r`,

while

`1 < H(Gamma^r C3)/H(Gamma^r TT_3) <= 1 + 2/(r+1)`.

The upper bound includes `r=0`, with equality there. In particular the ratio tends to one but the absolute difference tends to infinity.

**Proof.** The ordered path-cover polynomials of `TT_3` and `C3` are `q^3` and `q^3+2q`. THM-4900(6)'s linear tower law gives the displayed difference. Its bound `U_(3,r)>=(r+1)U_(1,r)` gives the relative estimate. The gate identity is `H(Gamma T)=H(T)+P_T(2)`. The two ordered covers `(T,empty)` and `(empty,T)` each contribute `H(T)`, so `P_T(2)>=2H(T)` and iteration gives `H(Gamma^r K1)>=3^r`. All these are full, ordinary Hamilton-path counts. QED.

**Exact equality guard.** If two tournament counts `a,b` are odd positive integers, then `|a-b|<2` implies `a=b`. A relative bound `|a/b-1|<=epsilon` therefore forces equality only when the accompanying absolute estimate gives `epsilon*b<2`. Tending to zero without this scale is insufficient. This is a necessary audit of that proof method, not a claim that all proofs of equality require error estimates.

**Summable-error hostile, on actual tournaments.** Sample these same two families at `r_j=2^j`, `j>=0`, and write `epsilon_j=H(Gamma^(r_j) C3)/H(Gamma^(r_j) TT_3)-1`. Then `0<epsilon_j<=2/(2^j+1)` and `sum epsilon_j<4`. Consequently `product_j(1+epsilon_j)` converges to a finite positive number at most `exp(4)`, by `0<log(1+epsilon)<=epsilon`. Yet every additive defect is at least `2*3^(2^j)`. Thus even summable relative errors and a convergent positive ratio product do not force an exact return. In contrast, summability of **absolute integer** defects forces all but finitely many defects to be zero; every nonzero one has magnitude at least one (at least two for differences of odd H values).

The gate family above is an actual positive control for relative convergence and a hostile control for exact return. It sharpens the requirement on any proposed micro/macro proof: preserve or recover the additive defect, the full path-cover profile, and the physical port. A universal limiting resolvent, an excellent normalized fit, or a finite threshold crossing alone cannot close a target-preserving return.

## C2. The exact MRB carrier and the deformation

**INHERITED identity; conventions fixed here.** For `N>=2`, take a transitive block `A={a_2,...,a_N}` and two more vertices `b,c`, with `A -> b -> c -> A`. Within A use increasing-label orientation. Denote the resulting strong tournament by `W_N`. Weight a Hamilton path by the product of `w_i` for vertices `a_i` appearing before b. The prior Session162 bijection gives

`Z_N(w)=product_(i=2)^N(1+w_i)+product_(i=2)^N w_i`.

Each subset S supplies the path `a_S,b,c,a_(A\S)`, and the exceptional path is `c,a_2,...,a_N,b`. At all weights one,

`H(W_N)=2^(N-1)+1`.

Now fix `alpha>0`, `0<lambda<=1`, and put

`f_n=n^(1/n)-1`, `a_n=lambda*f_n^alpha`, `w_n=(-1)^n a_n`,

`P_N=Z_N(w)-product w_n=product_(n=2)^N(1+(-1)^n a_n)`.

The correction removes exactly the exceptional marked path term. Since `0<f_n<=3^(1/3)-1<1`, all factors in `P_N` are positive. Also `f_n` decreases for `n>=3`, tends to zero, and `f_n~log(n)/n`. Thus `sum (-1)^n a_n` converges for every `alpha>0`. The original MRB sum is the special case `alpha=lambda=1`.

The full signed `Z_N` need not be positive for small N. For example at `lambda=1`, `N=3`, its limit as `alpha` decreases to zero is -1. Nothing in this packet treats a signed partition value as an ordinary path count.

## C3. A sharp square-mass threshold and all higher thresholds

**PROVED.** For the preceding fixed parameters:

1. If `alpha>1/2`, then `P_N` converges to a finite strictly positive limit.
2. If `0<alpha<=1/2`, then `P_N` tends to zero, although `sum (-1)^n a_n` still converges.
3. For `0<alpha<1/2`,

   `log P_N ~ -lambda^2/[2(1-2alpha)] * N^(1-2alpha)(log N)^(2alpha)`.

4. At `alpha=1/2`, there is a finite positive constant `K_lambda` with

   `exp((lambda^2/4)(log N)^2) P_N -> K_lambda`.

The threshold is caused by accumulated **even** powers of small weights, beginning with their squares. It is not a discontinuity of any finite N expression.

**Full hierarchy.** Set `J=floor(1/(2alpha))`, including a boundary integer. Then

`log P_N + sum_(j=1)^J lambda^(2j)/(2j) * sum_(n=2)^N f_n^(2j alpha)`

converges to a finite real limit. `J=0` means the sum is empty. Crossing `alpha=1/(2j)` adds one more divergent even-power counterterm. This is an infinite family of precisely classified accumulation thresholds, not a claim about arbitrary mathematical finite sets.

**General mechanism.** If real weights `w_n` tend to zero, every `1+w_n` is positive, and `sum w_n` converges, then the positive product `product(1+w_n)` has a nonzero limit iff `sum w_n^2` converges. Indeed `w-log(1+w)>=0` and its ratio to `w^2` tends to 1/2. Thus the convergent linear part determines neither the product's positivity nor its vanishing without the square-mass sidecar. The MRB deformation is an explicit all-parameter instance on an actual marked path carrier.

**Proof of the sharpened hierarchy.** For `0<=a<1`,

`log(1+(-1)^n a)=(-1)^n artanh(a)+(1/2)log(1-a^2)`.

The artanh sequence is eventually decreasing to zero, so its alternating sum converges. The second term has the expansion

`(1/2)log(1-a^2)=-sum_(j>=1) a^(2j)/(2j)`.

For any `A<1`, uniformly when `0<=a<=A`, its remainder after J terms is bounded in absolute value by

`a^(2J+2)/[(2J+2)(1-A^2)]`.

Take `A=lambda*f_3^alpha<1`. Since `2(J+1)alpha>1`, the remainder is summable by comparison with `(log n/n)^(2(J+1)alpha)`. This proves the hierarchy. It also proves the limit criterion: the square sum converges iff `2alpha>1`; otherwise the even part tends to negative infinity.

When `0<alpha<1/2`, the square sum is asymptotic to `lambda^2 N^(1-2alpha)(log N)^(2alpha)/(1-2alpha)`, by integral comparison and differentiation of the displayed comparison function. Higher powers are little-o of the divergent square sum since `a_n` tends to zero; the alternating part is bounded. This gives statement 3.

At `alpha=1/2`, only the first even counterterm diverges. Inherited Session144 Theorem 2 proves, with its stated Stieltjes convention,

`sum_(n=2)^N f_n = (1/2)(log N)^2 + kappa + o(1)`

for a finite constant kappa. Equivalently this follows from `f_n=log n/n+O((log n/n)^2)` and the integral-test finite part for `sum log n/n`. Therefore the hierarchy implies statement 4, with a finite real logarithm for `K_lambda`, so the limit is strictly positive. QED.

**Finite tail guard.** For `N>=3`, write

`R_N=log P_N+sum_(j=1)^J [1/(2j)] sum_(n=2)^N a_n^(2j)`.

If `R_infinity` denotes its limit and `beta=2(J+1)alpha>1`, then

`|R_infinity-R_N| <= artanh(a_(N+1)) + lambda^(2J+2) exp(beta*log(3)/3)/[(2J+2)(1-A^2)] * integral_N^infinity (log x/x)^beta dx`.

The alternating remainder gives the first term; `f_n<=exp(log(3)/3)log n/n` gives the second. The integral comparison is valid because `(log x/x)^beta` decreases for `x>=3`. At the main boundary `alpha=1/2`, J=1 and beta=2, so the integral is exactly `((log N)^2+2log N+2)/N`. Thus the theorem carries an explicit finite certificate, not merely an eventual fit.

## C4. The actual signed partition has the same limiting scale; H does not change

**PROVED.** For each fixed `alpha>0`, `0<lambda<=1`,

`Z_N/P_N -> 1`.

In particular `Z_N` is eventually positive and has the same zero/nonzero product limit and critical scaling as C3. This statement does not make it an unweighted path count.

**Proof.** `log|product w_n|=sum log a_n`. Since `a_n->0`, the average of `log a_n` tends to negative infinity. On the other hand `log P_N/N->0`, because each summand `log(1+(-1)^n a_n)` tends to zero. Hence `|product w_n|/P_N->0`, proving the claim. QED.

**Chan-Pak loss ledger.** The underlying tournament is exactly the same `W_N` for every alpha and lambda, and its ordinary count always equals `2^(N-1)+1`. Thus the analytic threshold can be moved through a fixed combinatorial family by changing only its marked observable. The theorem is a genuine bridge to weighted path partitions, but it supplies no return operation on ordinary tournaments, no growing-defect reset, and no new attainable ordinary H value. A construction converting signed or real weights into positive integer tournament counts would be an additional theorem. Session162's integer marked-token approximation is one available sidecar, but it counts marked paths with multiplicity, not the H values of newly constructed ordinary tournaments.

## C5. Gamma curves classify radius-dependent cutoffs

**Priority update, incoming concurrent work.** [Session297, BRIDGES, moving-scale section](../../session297_parity_cutoff_20261113/BRIDGES.md), present in incoming commit `a03df4746338daf016df8819bcce373ec35c95e8`, already proves the adjacent-dimensional radius crossover, its strict increase, the consecutive maximizing-dimension phases, and their large-radius scale. Credit those parts of C5 to that concurrent result and the older Session127 recurrence; they are not new contributions of this packet. The displayed all-gap formula below extends the same monotonic Gamma mechanism to arbitrary h>0. The continuous superlevel classification, exact unit-radius threshold at12, and `R_n=c sqrt(n)` comparison are explicit corollaries/variants, with their proofs and scale conventions retained here. No novelty in the classical Gamma facts is claimed.

Define the continuous dimension extension, for `R>0`, `x>=0`,

`V_R(x)=(pi R^2)^(x/2)/Gamma(x/2+1)`.

For integer x this is Euclidean ball volume; for other x it is the specified analytic interpolation. Session127 VOL-1 already proves the unit-ball formula and integer maximum at dimension five. The following radius/cutoff form isolates what the curve actually classifies.

**PROVED using standard Gamma identities.**

1. `log V_R(x)` is strictly concave on `[0,infinity)`, because

   `(log V_R)''(x)=-(1/4) psi'(x/2+1)<0`.

   Here the classical series `psi'(y)=sum_(k>=0)(k+y)^(-2)` is positive for y>0. If `R<=exp(-gamma/2)/sqrt(pi)`, the curve strictly decreases for x>0 and its maximum is at zero. For larger R it has one interior maximum, characterized by `psi(x/2+1)=log(pi R^2)`.
2. For every fixed R and every threshold `c>0`, the set `{x>=0: V_R(x)>=c}` is empty or a bounded closed interval (possibly a single point). Strict log-concavity gives the interval property, and Stirling gives decay to zero. A nonempty set can start at x=0; no claim of two endpoints is intended in that case.
3. For every h>0 the gap-h comparison `V_R(x+h)/V_R(x)` is strictly decreasing in x. Its radius crossover is

   `R_h(x)=[Gamma((x+h)/2+1)/(pi^(h/2) Gamma(x/2+1))]^(1/h)`.

   This is strictly increasing in x, and `V_R(x+h)>V_R(x)` iff `R>R_h(x)`. Different pairs at the same gap are different points on one monotone threshold family.
4. At R=1 the exact integer superlevel set is `{n>=0: V_1(n)>=1}={0,...,12}`. The strict set is `{1,...,12}`. Volume remains positive at every later dimension: this is a cutoff of a threshold predicate, not termination of the volume sequence.
5. If radius instead scales as `R_n=c sqrt(n)`, then

   `V_(R_n)(n) ~ (2pi e c^2)^(n/2)/sqrt(pi n)`.

   It tends to zero for `c<=1/sqrt(2pi e)`, with order `n^(-1/2)` at equality, and tends to infinity for c above that value. In particular the unit-radius cutoff is not invariant under changing the geometric normalization.

**Proof details and sources.** Differentiation and positivity give 1 and the derivative of the log gap ratio is `(psi(x/2+1)-psi((x+h)/2+1))/2<0`, giving 3. For the endpoint in 1, `psi(1)=-gamma`; psi is strictly increasing and tends to infinity. [DLMF 5.15.1](https://dlmf.nist.gov/5.15#E1) supplies the trigamma series. [DLMF 5.11.1](https://dlmf.nist.gov/5.11#E1) supplies Stirling, from which 2 and 5 follow. These are the only external analytic imports. The finite statement 4 follows from `V_(n+2)=2pi V_n/(n+2)`, parity unimodality, and the endpoint comparisons `V_11>1`, `V_12>1`, `V_13<1`, `V_14<1`. Use `3<pi<22/7`; for example `V_12=pi^6/720>729/720>1` and `V_13=128pi^6/135135<128(22/7)^6/135135<1`. The same rational comparisons certify V_11 and V_14. QED.

### The two supplied pairs, under this specific volume interpretation

`V_1(10)/V_1(3)=pi^4/160<1`,

`V_1(11)/V_1(4)=128pi^3/10395<1`.

Their equality radii are respectively

`R_7(3)=(160/pi^4)^(1/7)`,

`R_7(4)=(10395/(128pi^3))^(1/7)`,

and the latter is strictly larger by C5. Both pairs are ordinary members of the all-h/all-x family; there is no volume theorem selecting only these two pairs or identifying either radius with MRB. This is an explicit interpretation to test, not a claim that the user's pairs originally meant dimensions.

Their actual older tournament role is different: the [microcosm result, Section 3](../../microcosm_gate_doubling_resolvent_20261106.md) records `{3,10}` as `q=(H-1)/2` images of the forbidden values `{7,21}`, and also as square-sum onset counts. Translating q by one gives `{4,11}`; a theorem about that coordinate shift must retain its count/root convention. Incoming [Session297 BRIDGES](../../session297_parity_cutoff_20261113/BRIDGES.md) separately recovers `{4,11}` as a least-exponent inverse `3x-1` cycle and both pairs in the CRT runner construction. These are separate carriers; the companion pair is not a second tournament gap spectrum. The Gamma comparison does not preserve that target predicate.

## What transfers, and what remains open

| Source -> target | Map and preserved predicate | Information lost / sidecar | Cheapest decisive test |
|---|---|---|---|
| Gamma curve -> finite superlevel set | dimension and radius; exact volume inequality | arithmetic support and tournament data absent | change R to c sqrt(n): C5 changes the cutoff |
| MRB profile -> weighted strong tournament | C2 marked subset polynomial, exact signed partition | real weights versus ordinary integer counts | all parameter choices have the same H(W_N) |
| Small weight -> accumulated macro product | log, odd/even split; exact convergence iff alpha>1/2 | first-order alternating sum omits square mass | alpha=1/2: sum converges, product vanishes |
| Gate towers -> normalized curve | full H ratio tends to one | absolute defect and exact equality | C1 gives defect at least 2*3^r |
| `{7,21}` -> `{3,10}` | q=(H-1)/2; exact count normalization | labels, paths, ports, admissible operations | inverse map 1+2q retains the two forbidden values |
| `{3,10}` -> volume pairs | choose them as dimensions | no count predicate preserved | C5 applies equally to every gap-h pair |

The constructive next question is precise: can a full-count return compiler retain and cancel its absolute defect while allowing a controlled sequence of admissible repairs? The analytic product theorem suggests measuring accumulated second and higher errors, but does not supply such repairs; C1 shows that summable relative errors alone still fail on actual tournaments. Any claimed transfer must exhibit actual tournaments and ports and preserve H exactly. The current packet proves the accumulation threshold and the hostile guard; Chan-Pak and the repeatable return remain OPEN.

## Reproduction

Run `python -B 05-knowledge/results/mrb_seven_cutoffs_20261113/curves/verify_curves.py` from the repository root. The same script under `python -O -B` performs identical explicit checks. Exact controls compare direct Hamilton-path enumeration with both the subset polynomial and the gate polynomial operator, use rational Machin bounds for pi, and check finite log-series coefficients. Numeric illustrations are kept in a separate labelled receipt field and are not proofs of infinite convergence. Source pin hashes normalize CRLF to LF.
