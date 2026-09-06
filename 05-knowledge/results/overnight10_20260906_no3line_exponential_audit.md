# Independent audit of the exponential diagonal-defect theorem

**Verdict: PASS — all-size proof accepted, with independent exact controls.**
The audited source and report are `overnight10_20260906_no3line_defect.py`
and `.md` in this directory. The producer source was not read or imported.
I derived the mean normalization, finite constants, and coupling separately,
then read the complete proof. No substantive mathematical correction was
needed. One display-domain clarification was requested and applied: the
factorial-moment ratio with denominator `(n)_k^2` is stated for k<=n; for
k>n the expectation is separately zero. The proof already used that correct
interpretation; the producer source and output identities did not change.

The accepted model is an arbitrary fixed simple bipartite 2-regular graph
G on n+n vertices, with two independent uniform shore-label bijections.
The statistic is `F=sum_d (Y_d-2)_+` on the slope-one diagonals. Its zero
event equals zero selected-direction triples; only `X=0 => F=0` holds for
the actual all-direction collinearity count X. This distinction is retained
throughout the probability argument.

## 1. Retained mechanism and audit boundary

The closest proved mechanism is the eighth diagonal-density theorem; the
ninth shared-cell covariance refinement remains valid structural information.
The finite hostile is C8 versus 2C4: the cycle profile affects exact laws.
The repaired move replaces cubic occupancy counts by a bounded-change excess
with the same selected zero event. Its indispensable sidecar is the actual
pair of label permutations, which supports a conditional coupling. There is
no independence assumption about geometric triples.

The concept board is actual permutations, diagonal occupancies, rook
coefficients, alternating truncation, conditional interval length, and the
one-sided event implication. Compared with the LRC endpoint-owner anchor
and Laurent carry niche, this wildcard again retains a coordinate lost by
an aggregate observable; no theorem is transferred between those objects.
The cheap hostile is the n5 transposition with defect change exactly four.

## 2. Exact mean and the finite lower bound

For 0<=k<=n, inverse-label injections give

```
E(Y_T)_k = (L)_k k! m_k(G)/(n)_k^2.
```

Here m_k counts unordered source matchings. Ordering their edges contributes
k!, and no further automorphism factor occurs. When k>n the moment is zero.
The independent cycle count used by this referee is
`m_k(C_N)=N/(N-k)*binom(N-k,k)` for 1<=k<=N/2, with constant coefficient
one. It follows by counting cyclic gaps between selected edges and equals
the producer's path-closing binomial formula. Products across components
retain the complete matching polynomial.

The binomial expansion of `(y-2)_+` has coefficients
`(-1)^(k+1)(k-2)` for k>=3. Its partial sum through K minus the target,
for y>=2, is

```
(-1)^(K+1)*((K-2)y+2)/K * binom(y-2,K-1).
```

Thus even truncations give lower bounds without requiring decreasing
terms. The diagonal identity
`sum_d binom(L_d,k)=binom(n,k)+2binom(n,k+1)` proves the exact mean in
the producer report. The referee checks each summed moment and the mean
against literal boards, rather than only comparing polynomial formulas.

The universal coefficients are

```
m3=2n(n-2)(2n-5)/3,
m4=n(n-3)(2n-5)(2n-7)/6+c4.
```

The positive C4 correction is required. The lower truncation through four
gives the stated B4(n,c4). Since c4<=n/2, a simplified independent check is

```
B4(n,n/2)-2n/15
 = (11n^2-62n+78)/(15(n-3)(n-1)) > 0, n>=4.
```

The numerator becomes `11z^2+26z+6` at n=z+4. This proves the explicit
finite lower mean, not merely its asymptotic leading coefficient. It does
not claim that the carrier maximizing c4 minimizes the complete mean.

## 3. Uniform mean error

For k<=n, choose ordered distinct rows and one of two neighbors at each
row. If q_k is the probability that the columns are distinct, then
`k!m_k=(n)_k 2^k q_k`. A pair of chosen distinct rows has a column
collision with probability exactly `1/[2(n-1)]`: there are 2n ordered
edge pairs meeting a column out of `4n(n-1)` choices. Therefore
`1-q_k<=k(k-1)/[4(n-1)]` uniformly over every cycle type.

The target-row deficit is bounded by the collision union bound for draws
with replacement:

```
0<=theta^k-(L)_k/(n)_k <= binom(k,2)*theta^(k-1)/n.
```

Combining both bounds proves the producer's factorial deficit estimate.
For k>n its first term already dominates `(2theta)^k`; theta=0 is treated
directly. Hence all factorial series are uniformly absolutely dominated,
including their infinite tails. No uniform-integrability gap remains.

With lambda=2theta and `g(lambda)=lambda-2+(lambda+2)e^-lambda`, separating
the odd and even coefficients gives

```
-cosh(lambda)*[lambda^2/n+lambda^3/(4(n-1))]
 <= E(Y-2)_+-g(lambda)
 <= sinh(lambda)*[lambda^2/n+lambda^3/(4(n-1))].
```

The signs are correct: odd factorial coefficients are positive, so their
deficits lower the actual mean. Convexity and monotonicity of g give the
two trapezoid bounds. Integrating the four error terms separately yields,
at the worst n=4, the following equivalent constants:

```
Cminus = (13/6)e^2-(21/2)e^-2+2 =16.58860107369864... <17,
Cplus  = (13/6)e^2+(29/2)e^-2-2 =15.97198315461396... <16.
```

The upper constant includes the trapezoid excess. The referee certifies
these inequalities using rational Taylor intervals for exp(2), not decimal
tests. Thus, uniformly over G,
`alpha*n-17<=mu_G<=alpha*n+16`, where `alpha=1-5e^-2`.

## 4. Conditional coupling and the exponential constant

A row-label transposition moves four points; each deletion decreases F by
zero or one, and each insertion increases it by zero or one. The two totals
belong to [0,4], so their difference has absolute value at most four. The
same statement holds for columns. This proves a bound of four, not eight.
The n5 identity-plus-shift board attains defect change 5 to 1 under the
row transposition 0<->2; four is therefore sharp for this universal swap
argument.

Expose n-1 row labels followed by n-1 column labels. Two choices at a reveal
are paired by swapping those two still-available labels in every completion.
Earlier labels remain fixed; the coupling is a bijection of uniform
completion sets. Conditional expectations consequently have range of length
at most four. This is stronger information than merely an absolute bound
of four on a martingale increment.

For an interval of length h, the tilted-variance proof of the elementary
exponential-moment bound gives `E exp(tW)<=exp(t^2 h^2/8)` for mean-zero W.
Applying it at the 2(n-1) reveals yields
`E exp[-s(F-mu_G)]<=exp[4(n-1)s^2]`. Optimizing the Markov bound gives

```
P_G(X=0)<=P_G(F=0)<=exp[-mu_G^2/(16(n-1))].
```

Substitution of `mu_G>=2n/15` gives
`exp[-n^2/(900(n-1))]<=exp[-n/900]` for n>=4. Substitution of the uniform
asymptotic mean gives rate at least `(1-5e^-2)^2/16` in the negative
logarithm divided by n. Mixtures inherit the uniform lower-mean bounds and
the resulting uniform probability bounds. The concentration statement
around a particular mu_G remains conditional on G; it is not asserted
around an arbitrary mixture mean. The report's mixture firewall is sound.

## 5. Independent exact universe and reproduction

The referee source is standalone and imports no producer or repository
code. It generates every distinct board by row-neighbor pairs and traverses
its incidence graph for cycle type. This is independently implemented,
although the producer also uses the natural row-pair representation.

The controls contain:

* All 70,086 distinct simple degree-two boards n=3,...,6, with exact F means,
  selected zero probabilities, and every summed diagonal factorial moment.
* Every row and column transposition and every integer triple determinant
  on those boards with n<=5; the actual-target implication and sharp swap
  hostile are checked literally.
* Every one of the 384 cycle profiles with shore size 2,...,18, including
  the m3/m4 identities and finite lower mean. For n<=12, every target
  length 0,...,n and factorial order 2,...,n+4 checks the deficit bound.
* Direct rational conditional averages at 878 permutation-prefix nodes
  for C6, C8, and 2C4. Their maximum observed conditional width is two;
  this finite value is not substituted for the proved universal four.
* The exact alternating remainder for y=0,...,100 and K=2,...,30, and
  rational exponential intervals certifying both rounded constants.

The exact n4 means are 2/3 for 2C4 and 5/6 for C8; at n5 they are 29/25
for C4+C6 and 17/15 for C10. This reversal is a useful hostile to an
unjustified universal cycle ordering. At n4, 25 of 72 C8 boards have F=0,
while only two have full no-three-in-line success, so the missing directions
cannot be ignored in an existence claim.

```
python -B 04-computation/overnight10_20260906_no3line_exponential_audit.py
python -B -O 04-computation/overnight10_20260906_no3line_exponential_audit.py
```

Both executions pass **134,368 always-active gates**, with identical LF
output. All numerical gates use integer or rational arithmetic; decimal
values in the transcript are presentation only.

* Referee source SHA256:
  `a504620d2aca6635ef68f00025f1d3a716055797e88d12e9e204ee8ee2c1fcae`.
* Referee output SHA256:
  `062c56e46f7f080ce2815cc146b5b710ada7f59a488b2f43b487de4004a7306f`.
* Semantic SHA256:
  `42bdbe1ee0a1c93eb8a8cb6eaece5afbf0d4d48d07f94ec4bf5bf5793745c4b2`.
* Producer source identity supplied for review, not imported:
  `74ce43c6435da3c4662b5397c9e89ef86874452fcb1b445521a57b55816399d6`.
* Producer output identity:
  `b729b941bb7e4b8b2decb6486a1713554c59915ee31717982847495da7faa295`.

The probability theorem is accepted analytically for all stated sizes and
carriers. The finite controls verify normalization and failure boundaries;
they are not extrapolated to prove the asymptotic statement.

**Filing:** root integrated these audited artifacts in the tenth checkpoint;
reproduction paths are relative to the repository root. Earlier outside-worktree
notes describe author provenance, not the present file location.
