# Independent audit: tensor-reserve stability for every n >= 140

**Conclusion: PASS.** Given the independently accepted tensor matching
reserve and accepted Session138 inputs, `PROOF_N140.md` establishes
`K_n=kappa_n` for all integers n>=140, with the stated positive-variance
fair-pair equality classification. No mathematical correction is required.
This proof alone leaves orders 10 through 139 unresolved.

The audit independently reconstructed the bootstrap and endpoint arithmetic
before reviewing the producer's completed proof. All executable controls use
exact rational arithmetic and active exceptions. They support, rather than
replace, the uniform arguments below.

## Reviewed source pins

Relative to the Session139 research root:

- `stability/PROOF_N140.md`, SHA-256
  `23474d7041549ed099d261812397a31599ac32e31b5c6c4b6e8985d508663314`.
- `stability/verify_n140.py`, SHA-256
  `e79fb618d6922f9d24c86e61f7bfb7b3d2f9f3a8bcf098e2c56a80cb9579c509`.
- Accepted `tensor_matching/PROOF.md`, SHA-256
  `7bd612d505eac4940c4add25b1e109d99e5369e8633fb866fb27c9af195d1b66`.
- Accepted `tensor_matching_audit/AUDIT.md`, SHA-256
  `1d58a11025063d58179d617aa66295dcb33df6812f93026261f504ffc1f1c483`.

Session138's nonlinear g bound, matching coefficient, greedy wedge bounds,
and equality results are accepted inputs, already identified in the earlier
independent threshold audits. This continuation does not purport to reprove
all historical inputs. Previously frozen n>=748 and n>=526 proofs were
hash-checked and remain unchanged.

## 1. The new reserve is applied within its actual domain

Choose the complete positive-weight greedy matching, with weights w_i and
mass M. Round its complement first. The accepted one-edge inequality has
nonnegative corrections, so discarding them gives a lower bound by the
average gap over branches with matching-only uncertainty. In every branch,
all retained weights are unchanged. The accepted matching reserve is

    Delta_branch >= B_branch + 1-product(1-w_i)
                                  -1_(n=2r) product(w_i).

Averaging retains the original terminal binary expectation B exactly once:
independent rounding has the same product law regardless of edge order.
Thus no false conditional multiplicative inequality is used.

The product subtraction is present only when the selected matching is
perfect. Then n=2r=2R, and every w_i<=1/4 gives `product(w_i)<=4^(-R)`.
For an imperfect matching, including odd n, that subtraction is zero.
It would be incorrect to assert `product(w_i)<=4^(-R)` for an arbitrary
smaller matching; for r=1 and w_1=1/4 that assertion is false.

Finally `product(1-w_i)<=exp(-M)`, so the uniform reserve is indeed

    Delta >= B+1-exp(-M)-4^(-R).

The n>=140 domain avoids empty pair-deletion remainders. Positive V makes
the greedy matching nonempty. The r=0 boundary of the accepted tensor
statement nevertheless agrees with the zero-reserve binary case.

## 2. Bootstrap, strictness, and uniformity

The old nonlinear bound first gives `kappa V>B+c0` with
`c0=249999/1000000`. Specifically, decreasing g(V)/V excludes V<=1,
because `g(1)>49/200>4/70`. Then

    V>49R/800>=343/80>4,
    g(V)>(1-50^(-4))/4=c0+3/3125000.

Here `exp(4)>50` follows from the exact finite sum `16319/315>50`.

For any established positive c with `kappa V>B+c`, the same B as in the
new reserve satisfies `B>=D^2/(2n)`. Completing the square gives

    M=V-D > c/kappa+(D-n kappa)^2/(2n kappa)-n kappa/2
          >=cR/4-4-2/R=:m_c(R).

This uses `kappa<=4/R`, `n<=2R+1`, and positive kappa. The derivative of
m_c is `c/4+2/R^2>0`. Hence each endpoint lower bound extends uniformly
from R=70.

The independent executable certificate uses the weaker uniform penalty
`4^(-R)<1/10000`, rather than the producer's exact `4^(-70)`. At R=70,
the following positive Taylor slacks are exact:

| Advance | m_c(70) | Taylor degree | `(1-c'-1/10000) T(m)-1` |
|---|---|---:|---|
| c0 to 7/25 | 969951/2800000 | 2 | 1956167374744799/156800000000000000 |
| 7/25 to 57/100 | 61/70 | 3 | 104080503/6860000000 |
| 57/100 to 99/100 | 333/56 | 4 | 83791158289/786759680000 |

Every mass and every remaining probability factor is positive. For positive
m, the finite Taylor sum is strictly below exp(m). Thus each row proves
`1-exp(-m)-1/10000>c'`, and therefore the required exact-penalty advance.
The strict M lower bound also preserves strictness when the reserve is
substituted. The positive derivative of m_c and decreasing `4^(-R)`
justify every R>=70, not just the endpoint. After three advances,

    kappa V > B+99/100.

There is no circular assumption that the final reserve already exceeds
99/100, and B is never discarded during this bootstrap.

## 3. Heavy matching and geometric tail

Set `c=99/100`, `d=19/20`, `t=d/8`, and `lambda=d/64`. Assuming L>=1
selected weights are light (w<=t), their exact total mass U is at most
tn/2. The retained heavy matching has

    r_h > c(n-1)/2 -16n/(n-1)-dn/4.

At n=140 this lower bound is `540429/27800`, exceeding 19 by
`12229/27800`. Its derivative is `103/400+16/(n-1)^2>0`.
Therefore integer r_h>=20. In particular the retained matching is
nonempty before its coefficient is used. With m=R-L, one has
`m>=r_h>=20` and `n>2m`; the coefficient a_r decreases with r.

The retained matching inequality, keeping U exactly, yields

    cL-4U < 4V(3/4)^(R-L)+4D.

The variance-wedge bound supplies the envelope

    4V<=Q_R=R+33/2+8/R.

Independently, this follows directly from `4V<=n/2+2n kappa`:
the even-order bound is R+16, and the odd-order bound is Q_R. Thus it
is valid below n=225; no inapplicable `V<=n/7` estimate is imported.

For fixed R, `log((3/4)^(R-L)/L)` has second derivative 1/L^2>0, so
its maximum on `[1,R-20]` lies at an endpoint. Both endpoint bounds
are decreasing in R>=70:

- `Q_R(3/4)^(R-1)` has successive ratio below `(3/4)(1+1/R)<1`.
- `Q_R/(R-20)=1+(73/2+8/R)/(R-20)` has a positive decreasing excess.

Both endpoint values at R=70 are strictly below c-d=1/25. The larger
endpoint `(Q_70/50)(3/4)^20` has exact positive slack

    132791254065377/3848290697216000.

The other endpoint is approximately 2.07404e-7, also checked exactly.
Therefore the light deficit satisfies

    sum_light(d-4w_i)=dL-4U<4D.

Every summand is at least d/2>0, so D>0 in this branch. The whole
argument is uniform over the permitted L interval, not an L scan.

## 4. Local dual certificate, zero cases, and equality

The previously accepted local identity remains valid with the new d:

    wS/2+S^2/4-tS+lambda(d-4w)
      =(S/2-t+w/2)^2+w(2t-w)/4 >=0

for `S>=0` and `0<w<=t`. Greedy groups count disjoint wedge pairs, so
heavy and light lower bounds may be summed simultaneously. This gives

    W>=tD-lambda sum_light(d-4w_i)
      >(t-4lambda)D=(d/16)D.

The decisive parity-safe endpoint margin is exactly

    d/16-4/70=1/448>0.

Consequently `W>(d/16)D>kappa D>=B>=W`, a contradiction. No bound on
B has been applied in the wrong direction.

All boundary cases are covered:

- V=0 is a binary tournament and has nonnegative gap by the accepted
  binary theorem; positive variance is the infimum domain.
- L>=1 forces D>0 by the deficit inequality; D=0 cannot escape the
  contradiction and is never used as a divisor.
- An individual light group may have S=0, and the local polynomial
  inequality remains valid. No positive light discarded mass is required.
- B=0 cannot escape the contradiction because B>=W and W is forced
  strictly positive in the L>=1 branch.
- Once L=0, every selected weight exceeds t>kappa. If D>0, summing the
  selected-edge wedge contributions gives `W>kappa D`, impossible.
  Hence D=0 and then B=0.
- The uncertainty graph is now a matching. The accepted sharp matching
  inequality supplies the lower bound, and ordinal fair pairs attain it.
  For equality, all positive-probability binary completions have zero
  gap and are transitive; the accepted H=1 minimization gives precisely
  maximal ordinal fair pairs with one singleton in odd order.

## 5. Executable evidence

`audit_n140.py` is independent and uses explicit exceptions. It emits
46 active checks, including a stronger coarse-penalty bootstrap,
coefficient-wise polynomial identity, representative correct perfect-only
penalty branches, parity bounds, and all quantitative endpoints.

- `normal.json` and `optimized.json`: PASS, byte-identical.
- `negative_threshold139{,_optimized}.json`: threshold 139 is rejected
  by the first bootstrap slack in both modes, exit status 1. This is a
  failed certificate at 139, not a counterexample to the theorem there.
- `negative_perfect_indicator{,_optimized}.json`: removing the perfect
  indicator is rejected already for n=140,r=1,w=1/4 in both modes,
  exit status 1.
- `producer_normal.json` and `producer_optimized.json`: producer exact
  checks independently rerun, PASS and byte-identical.
- `producer_negative.json`: producer threshold 139 also fails its first
  bootstrap step in both modes, exit status 1.

No Git or remote changes were made. The preceding n>=748 and n>=526
proofs and their audit packets were not modified.
