# Independent audit: the light-deficit n >= 526 theorem

**Conclusion: PASS.** The new argument establishes exact weighted stability
`K_n=kappa_n` for every integer `n>=526`, with the accepted fair-pair
classification of positive-variance equality. It leaves orders 10 through
525 unresolved by this argument. The already accepted n>=748 checkpoint
is not modified.

The audit independently reconstructed the mechanism and endpoint constants
before reading the full producer proof. No mathematical correction was
required. Verification is exact; uniform extension is proved analytically,
not inferred from a census or finite numerical scan.

## Source pins and accepted scope

Reviewed producer files under the Session139 research root:

- `stability/PROOF_N526.md`, SHA-256
  `e82bf3b0eb986f63daee55d42e84c3545921800e0a9936484cfdff97a03b6561`.
- `stability/verify_n526.py`, SHA-256
  `9522c195aad09a340caa07aacb8aa7f2aadb35e0c434ff43be08d81fb1930276`.

Accepted Session138 inputs are the retained binary expectation/matching
bounds, nonlinear margin, variance wedge lower bound, greedy group bound,
and equality classification. They were read in the preceding independent
n>=748 audit, which pins the accepted `asymptotic_stability/PROOF.md` and
`threshold_refinement/PROOF_N1500.md`. In particular the stronger local
bound `W_i>=w_i S_i/2+S_i^2/4` already appears in accepted Session138
`asymptotic_stability/PROOF.md`, Section 3. This audit is not a new review
of all historical inputs.

## Independent mathematical reconstruction

### 1. Candidate bootstrap and the nonempty retained matching

For a candidate `Delta<=kappa V`, positive V and decreasing `g(V)/V`
exclude V<=1 because `g(1)>49/200>kappa`. The same exact exponential
partial sum as before gives `exp(4)>50`. Hence

    V>49(n-1)/1600>4,
    g(V)>(1-50^(-4))/4=c+3/3125000,
    c=249999/1000000.

Using `B>=D^2/(2n)` together with the nonlinear bound correctly gives

    V-D > c/kappa+(D-n kappa)^2/(2n kappa)-n kappa/2.

For light selected weights `w<=t=d/8`, with `d=487/2000`, their total
mass U is at most tn/2. The heavy matching therefore satisfies

    r_h > c(n-1)/2 -16n/(n-1)-dn/4.

At n=526 the right side is `29524339/1680000`, which exceeds 17 by
`964339/1680000`. Its derivative is
`128249/2000000 +16/(n-1)^2>0`. Thus integer `r_h>=18` uniformly.
Under L>=1, `m=R-L>=r_h>=18`; also `n>2m`. The retained matching
coefficient is legitimately compared with a_m because a_r is the mean
of the first r terms of the decreasing geometric sequence.

### 2. Exact light mass and the uniform tail budget

The rearranged retained inequality implies

    cL-4U < 4V(3/4)^(R-L)+4D.

Keeping U, rather than replacing U by tL, is essential. No negative term
has been dropped in the wrong direction: the omitted term is
`-4 epsilon U<=0` after replacing `4(1-epsilon)U` by `4U`.

An independent parity derivation gives

    4V <= n/2+2n kappa <= Q_R=R+33/2+8/R.

For even n the intermediate bound is R+16; for odd n it is exactly Q_R.
The producer's alternative derivation via the increasing function
`n/2+16n/(n-1)` is also valid for n>=526.

For fixed R, `log((3/4)^(R-L)/L)` is convex on the entire positive real
L interval, since its second derivative is 1/L^2. Therefore its maximum
on `[1,R-18]` occurs at an endpoint. Both resulting bounds decrease as
R increases from 263:

- `Q_R(3/4)^(R-1)`: successive ratio is below
  `(3/4)(1+1/R)<1`.
- `(Q_R/(R-18))(3/4)^18`: the ratio factor is
  `1+(69/2+8/R)/(R-18)`, whose positive excess decreases.

At R=263 both are strictly less than `c-d=6499/1000000`. The larger
endpoint is

    56963596759137/8855878966968320,

with positive slack

    1846127022469099/27674621771776000000.

Consequently `sum_light(d-4w_i)=dL-4U<4D`. Every summand is at least
d/2>0. Thus D>0 before the later strict contradiction is used.

### 3. The local certificate and shared budget

For `lambda=d/64` and `t=d/8`, direct polynomial expansion gives

    wS/2+S^2/4-tS+lambda(d-4w)
      =(S/2-t+w/2)^2+w(2t-w)/4.

The right side is nonnegative throughout the actual domain
`S>=0, 0<w<=t`. Distinct greedy groups count disjoint wedge pairs, so
light bounds and heavy bounds may be summed simultaneously without
spending the same wedge budget twice. Thus

    W>=tD-lambda sum_light(d-4w_i)
      >(t-4lambda)D=(d/16)D.

The endpoint margin is

    d/16-4/263=81/8416000>0.

It is important to use `kappa<=4/floor(n/2)<=4/263` here. The rough
`8/(n-1)` bound alone is insufficient at the even endpoint n=526.
This gives `W>(d/16)D>kappa D>=B>=W`, the desired contradiction.

### 4. Boundaries, closure, and equality

- V=0 makes the tournament binary and is covered by the accepted
  nonnegative binary gap; positive variance is the infimum domain.
- Positive V ensures the greedy matching is initially nonempty.
- If L>=1, the heavy matching is proved nonempty before a_r is used.
- D=0 is excluded in that branch by the strictly positive light deficit.
- S=0 for any individual light group is permitted by the local identity;
  all light discarded masses being zero causes no exception.
- B=0 cannot evade the strict contradiction, since B>=W and the
  contradiction forces positive W.
- Once L=0, D>0 would give `W>=sum w_i S_i>tD>kappa D`; hence D=0 and
  B=0. No division by a possibly zero discarded mass is performed.
- Equality then reduces to the previously accepted matching and binary
  equality inputs. Every positive-probability completion has zero gap
  and is transitive; the accepted H=1 minimization supplies exactly
  maximal ordinal fair pairs, plus one singleton in odd order.

## Executable evidence

`audit_n526.py` is independently written and uses explicit exceptions,
not Python assert statements. It verifies 25 conditions, including exact
bivariate polynomial coefficients, parity identities, both tail endpoints,
and the strict final margin.

- `normal.json` and `optimized.json`: PASS and byte-identical.
- `negative.json` and `negative_optimized.json`: a proposed threshold 525
  is rejected with exit status 1 at
  `d/16-4/262=-203/4192000<0`; byte-identical. This only rejects the
  present certificate at 525, not the stability theorem itself.
- `producer_normal.json` and `producer_optimized.json`: independently
  rerun producer controls, PASS and byte-identical.
- `producer_mutation.json`: changing lambda from d/64 to d/63 is rejected
  at the polynomial identity in both normal and optimized Python, with
  exit status 1.

This audit made no Git or remote changes. No census, floating-point
optimization, unsupported Pfaffian claim, or Riemann analogy is needed.
