# Independent audit: concurrent-input exact stability for n >= 68

**Conclusion: PASS.** The reviewed argument proves `K_n=kappa_n` for every
integer n>=68, including the inherited positive-variance fair-pair equality
classification. Orders 10 through 67 are not resolved by this proof. No
mathematical correction is required.

This is an independently checked combination of the concurrent degree and
retained-baseline mechanisms with Session139's accepted tensor reserve.
It is not a graph census or an external priority claim.

## Source pins and attribution

Reviewed candidate:

- `concurrent68/PROOF.md`, SHA-256
  `ffe64db6bd89fb32677f28747c800d172656aefe47421bcca161a40db1d8c41b`.

The following concurrent files were read by targeted `git show` at commit
`c16226249f1f27a335cf7c84afdeaa90879b9580`:

- `05-knowledge/results/stability_priority_moduli_20261102/STABILITY_N500.md`,
  SHA-256 `fd0ddd18e4b0a57bff61ef4b71279ccd78c53441e2bc9bdb78253c6f702f5770`.
- `05-knowledge/results/stability_priority_moduli_20261102/STABILITY_N276.md`,
  SHA-256 `59fe8d25ec5a15bdc61d1481baa389774a021370c59e4237a3cc81702b29aaf0`.

**Credit:** the factor-four stronger discarded-degree inequality and the
retained negative-RB light budget are concurrent packet contributions, not
discoveries claimed by this audit or the combination. The concurrent packet
itself credits earlier Session139 square-completion and per-edge-tail work.
The accepted tensor product reserve is the new Session139 input, pinned in
`independent_threshold140/AUDIT.md`; the inherited matching, binary-gap,
nonlinear-g, and equality theorems remain accepted Session138 inputs.

All earlier threshold proofs and audits were left unchanged. No Git or
remote mutation was performed.

## 1. Independently reviewed graph cancellation

For the greedy matching, let E_d contain all discarded positive-weight
edges, let D be their total weight, and let d_v be their weighted degrees.
The exact sum over discarded/discarded wedges is

    W_dd = (1/2)sum_v d_v^2 - sum_(f in E_d) q_f^2.

A simple graph has no two distinct edges with two common endpoints, so
this expression counts each unordered incident edge pair exactly once.
Each discarded f is assigned to the selected edge e that removed it.
Greedy maximal weight implies q_e>=q_f, so that one selected/discarded
wedge contributes q_e q_f>=q_f^2. These assigned pairs (e,f) are distinct
for different f and disjoint from every discarded/discarded wedge.
Therefore

    W>=W_dd+sum_f q_f^2=(1/2)sum_v d_v^2>=2D^2/n.

This proof works for arbitrary nonnegative real graph weights; it does not
require q<=1/4. Zero weights can be omitted. No overlap with other selected
edges is presumed, and extra selected/discarded wedges only increase W.

## 2. Stronger square completion and tensor bootstrap

Under the candidate condition `Delta<=kappa V`, the full retained matching
bound gives `B<=kappa D`. Whenever `Delta>B+c`, the stronger degree bound
and B>=W give

    M=V-D > c/kappa + 2D^2/(n kappa)-D
          = c/kappa-n kappa/8
            +[2/(n kappa)](D-n kappa/4)^2
          >=cR/4-1-1/(2R),

where R=floor(n/2), kappa<=4/R, and n<=2R+1. Every division is by a
positive quantity. The final mass lower bound has derivative
`c/4+1/(2R^2)>0`, so endpoint checks extend uniformly from R=34.

The starting c=6/25 is valid: decreasing g(V)/V and
`kappa<=2/17<g(1)` exclude V<=1; then `g(V)>g(1)>49/200>6/25`.
The accepted matching reserve retains the same binary expectation B:

    Delta>=B+1-exp(-M)-4^(-R).

The product correction is present only for a perfect chosen matching,
when its size is exactly R and its variance product is at most 4^(-R).
For an imperfect matching it is zero. Averaging the complement-rounding
branches retains B once, by the independent terminal product law. This is
not the invalid conditional factor inequality previously rejected.

At R=34 the mass lower bounds exceed 1, 4, and 7 by, respectively,

    43/1700,    1157/3400,    134/425,

as c advances through 6/25, 63/100, and 49/50. With
`4^(-R)<1/10000`, the exponential lower bounds 65/24, 51, and 110 give
strict reserve increments above 63/100, 49/50, and 99/100, with margins

    87/130000,    149/510000,    89/110000.

The first exponential bound is the degree-four Taylor sum for exp(1);
the second follows from `16319/315>51`; and the degree-four sum for exp(7)
is `4553/24>110`. Positivity of the omitted Taylor terms makes the
exponential inequalities strict. Thus every candidate satisfies

    kappa V-B>99/100.

No final reserve is assumed during an earlier bootstrap step.

## 3. Retaining the negative-RB term is algebraically valid

Set `c=99/100`, `eta=2/17`, `delta=1/21`, and `a=1979/4200`.
If L>=1 selected weights are at most eta, the remaining mass exceeds

    R(c/4-eta)-1-1/(2R).

At R=34, four times this expression exceeds 13 by 511/850. Its derivative
is `c-4eta+2/R^2>0`. Each selected variance is at most 1/4, so the heavy
matching size r_h is an integer at least 14. This establishes a nonempty
retained matching before its coefficient is used.

Let m=R-L>=r_h>=14. Since L>=1, n>2m. The decreasing matching coefficient
therefore gives, with epsilon=(3/4)^m,

    (kappa V-B)(R-L)>=4(1-epsilon)(V-D-eta L).

Expanding and moving terms yields

    L[kappa V-B-4eta(1-epsilon)]
      <=(kappa R-4)V+4epsilon V+4(1-epsilon)D-RB
      <=4epsilon V+4D-RB.

In particular the negative-RB term is genuinely retained, not manufactured
by subtracting a bound in the wrong direction. The mass on the right of
the first matching inequality is positive by the preceding heavy bound.

## 4. Uniform tail, strict absorption, and shared budget

The variance-wedge estimate gives `4V<=Q_R=R+33/2+8/R`; independently,
`4V<=n/2+2n kappa` gives R+16 at even n and exactly Q_R at odd n. Thus
this bound is valid at n=68 and uses no inappropriate V<=n/7 estimate.

Convexity in L of `log((3/4)^(R-L)/L)` places its maximum on
`[1,R-14]` at an endpoint. Both endpoint envelopes decrease for R>=34:

- `Q_R(3/4)^(R-1)` has successive ratio below `(3/4)(1+1/R)<1`.
- `Q_R/(R-14)=1+(61/2+8/R)/(R-14)` has positive decreasing excess.

At R=34 both envelopes are strictly below delta. The tighter endpoint
has positive exact slack

    1/21 - [Q_34/20](3/4)^14
      =1854611611/766651662336.

Consequently `4epsilon V<delta L`. The absorption margins are

    c-4eta-delta-a = a-4eta =43/71400>0,

so `aL<4D-RB`. In particular D>0, since a,L are positive and B>=0.

For discarded masses D_l,D_h assigned to light and heavy groups, their
counted wedge pairs are disjoint. Hence

    B>=W>=eta D_h+D_l^2/(4L)>=eta D-L eta^2.

The final comparison is the square `(D_l-2L eta)^2>=0` after multiplying
by 4L. Combining this weak lower bound on B with the strict light-count
upper bound gives precisely

    D eta(a-4eta)<B(a-R eta^2).

The left side is strictly positive. If the right coefficient is nonpositive,
B>=0 makes a contradiction immediate. Otherwise using B<=4D/R is valid,
and R eta>=4 gives

    B(a-R eta^2)<=D(4a/R-4eta^2)<=D eta(a-4eta),

contradicting the strict inequality. At R=34, R eta=4 exactly, not
strictly; the earlier strict inequality is sufficient. The coefficient
`a-R eta^2` also equals 43/71400 at this endpoint, so no sign case is
silently skipped.

## 5. Boundary cases and equality

- V=0 is covered by the inherited nonnegative binary gap, with transitive
  deterministic equality. Positive variance is the domain defining K_n.
- Positive V ensures the positive-weight greedy matching is nonempty.
- Under L>=1, the heavy matching is proved nonempty, L>0 permits the
  Cauchy-Schwarz division, and D>0 follows before any strict positivity
  argument depending on D.
- D_l=0 and individual zero discarded group masses cause no problem in
  the joint square bound; B=0 cannot escape the positive-left-side
  contradiction. No division by B or D is required.
- Once L=0, selected weights exceed eta>=4/R>=kappa. If D>0, at least one
  group has positive mass and `W>=sum w_i S_i>eta D>=kappa D>=B>=W`.
  Therefore D=0 and B=0, including at the eta=4/R endpoint.
- The uncertainty graph is a matching, and the accepted matching theorem
  proves the desired lower bound. For equality, every positive-probability
  binary completion has zero nonnegative gap and is transitive. Their
  Hamilton-path expectation is H=1; the accepted scalar classification
  yields maximal ordinal fair pairs, plus one singleton when n is odd.
  These constructions attain kappa and establish sharpness.

## Executable controls

`audit_n68.py` independently emits 38 active checks, using explicit
exceptions rather than Python assertions. Its algebraic checks include
coefficient-wise square identities and all endpoint signs above. In
addition it exhausts 760 exact weighted simple graphs: every graph on
1, 2, 3, or 4 labelled vertices with edge weights in `{0,1/8,1/4}`.
Those controls verify the discarded-degree identity, selected payment,
stronger W bound, and joint budget separately. They are finite controls,
not the reason the all-real-weight graph theorem is valid.

- `normal.json` and `optimized.json`: PASS and byte-identical.
- `negative_threshold67{,_optimized}.json`: rejected at the first
  bootstrap mass check, whose excess is -29/825, exit status 1.
  This rejects this certificate at 67, not the stability theorem there.
- `negative_payment{,_optimized}.json`: omitting the selected/discarded
  payment fails on a unit-weight three-vertex path, where D=1,
  W_dd=0, and the required lower bound is 2/3; exit status 1 in both modes.

The proof uses the displayed monotonic arguments for every larger order;
no numeric scan is extrapolated into an infinite claim.
