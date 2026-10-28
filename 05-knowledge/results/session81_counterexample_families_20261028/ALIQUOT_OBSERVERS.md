# Aliquot observers: prime-power recovery and a mixed-prime collision

Status: **PROVED** elementary identities and algorithms; **FINITE-EXACT** tests below. Positive-integer domain; `s(n)=sigma(n)-n`, including `s(1)=0`. No forward-aliquot termination claim.

## Inheritance and live concepts

Recovered source: [Session78 ALIQUOT.md](../session78_fourth_power_aliquot_20261027/ALIQUOT.md). Its complete fibers are `F11={21}`, `F12={121}`, `F13={27,35}`. It already proves the finite bound, prime-dilation reconstruction, exact carry and modular-trace obstructions. Those are **inherited, not new**. Its multiplier-4 hostile shows a formula fails, not that the observer fails. That corrected near miss is the entry point here.

Concept board: scalar fiber / specified operation / valuation profile / guarded inverse / attachment symmetry. Closest proved mechanism: prime-dilation valuation recovery. Hostile: `(s(3),s(12))=(1,16)` defeats substitution of composite 4 for prime p. Least-used sidecar: the valuation profile on the multiplier's support.

## 1. Prime-power multiplier theorem

For every known prime `p` and integer `b>=1`, the map

`n -> (s(n),s(p^b n))`

is injective on positive integers. Set `q=p^b`, `h=(q-1)/(p-1)`, `x=s(n)`, `y=s(qn)`. Then

`c=(y-q*x)/h`, `d=(y-x)/h`, `a=v_p(d)`,

and

**`n=(p^(a+1)*c-d)/(p-1)`.**

Proof: write `n=p^a m`, `p∤m`. Multiplicativity gives `y=q*x+h*sigma(m)`. Consequently `c=sigma(m)` and

`d=(p-1)x+c=p^a*(p*c-(p-1)m)=p^a*(m+p*s(m))`.

The last factor is positive and prime to p, so its valuation is exactly a. Solving yields the displayed formula. This genuinely extends session78 from primes to all prime powers. In particular multiplier 4 is injective despite the old prime-only formula's failure.

For arbitrary integer observations require `h | (y-qx)` and `h | (y-x)`; require `c,d>0`; form `a=v_p(d)` and `m=(pc-d/p^a)/(p-1)`; require m a positive integer, `p∤m`, and `sigma(m)=c`. These guards are necessary and sufficient. Sufficiency follows by re-encoding `x=(d-c)/(p-1)`, `y=(q*d-c)/(p-1)`. The reconstructed source is `p^a m`. For example q=4 and observations `(1,16)` give `h=3,c=4,d=5,a=0,m=3`.

## 2. Arbitrary multiplier: exact strata and complete inversion

Fix integer `k>=2` with prime factorization `k=product p^b_p`. For a source n retain its labelled profile `a_p=v_p(n)` for every prime dividing k. Define

`R_a = product_{p|k} (p^(a_p+b_p+1)-1)/(p^(a_p+1)-1)`.

Then

**`s(kn)=(R_a-k)*n+R_a*s(n)`, with `R_a>k`.**

Proof: divide the factorization formula for sigma(kn) by sigma(n); all factors outside k cancel. Each displayed factor exceeds `p^b_p`, because numerator minus `p^b_p` times denominator is `p^b_p-1>0`. Subtract kn.

Thus the lifted observer `(s(n),s(kn),a)` is injective for every k, with

**`n=(y-R_a*x)/(R_a-k)`.**

On every fixed x-fiber and fixed profile stratum, transport `n -> s(kn)` is strictly increasing. Any collision must cross valuation-profile strata. The full profile is sufficient, not asserted minimal; R itself is an even smaller sufficient coordinate. Recording only the set of multiplier primes dividing n is not enough (section3).

### Complete algorithm without scanning all sources

Given integers x,y and k>=2:

* x<0: empty. x=0: test n1 only.
* x=1: sources are exactly primes. Test primes dividing k, and the single candidate `n=(y-sigma(k))/s(k)` when it is integral, positive, prime and coprime to k. Indeed for a prime n not dividing k, `y=s(k)*n+sigma(k)`.
* x>=2: use inherited bound `N=(x-1)^2`. Enumerate all exponent profiles with `product p^a_p<=N`. Compute rational R and the candidate `(y-Rx)/(R-k)`. Retain it precisely when it is an integer in `[1,N]`, has the enumerated profile, and its divisor sums give x,y.

Completeness: every actual source lies within the inherited bound and hence appears in its profile stratum. Soundness: exact re-encoding checks. Profile count is at most `product_{p|k}(1+floor(log_p N))`; candidate factorization is not claimed polynomial-time in bit length. For fixed k this replaces a quadratic-in-x source scan with polynomial-in-log(x) many candidates (each still requiring arithmetic/factorization). No prime-square exception is lost because the weaker universal bound N includes it.

## 3. A mixed-prime observer really can fail

Take k=40. Distinct sources

`u=22655=5*23*197`, `v=24275=5^2*971`

have

`(s(u),s(40u))=(s(v),s(40v))=(5857,1303480)`.

Direct factor certificates:

* `sigma(u)=6*24*198=28512`, `sigma(v)=31*972=30132`;
* `sigma(40u)=15*31*24*198=2209680`;
* `sigma(40v)=15*156*972=2274480`.

Subtract `u,v,906200,971000`, respectively. The profiles on `(2,5)` are `(0,1)` and `(0,2)`; their ratios are `155/2` and `2340/31`. Both have the same support-intersection with 40, but unequal valuations. Thus the 5-adic exponent, not just shared-prime presence, matters. The complete guarded profile algorithm proves the entire observation fiber is exactly `{22655,24275}`; this completeness is FINITE-EXACT using the proved terminating algorithm.

This refutes universal injectivity for composite multipliers, but does **not** classify the multipliers with injective observers. In particular no assertion is made that every non-prime-power fails.

## 4. Procedural counterexamples and honest structural connection

For any `u!=v` in the same scalar fiber `s(u)=s(v)=x`, and any `k>=2` coprime to uv,

**`s(ku)-s(kv)=s(k)*(u-v) !=0`.**

Indeed sigma(ku)=sigma(k)(u+x), and likewise for v. Since `s(k)>=1`, difference is nonzero. This supplies a procedural counterexample to every proposed universal scalar law `s(kn)=f_k(s(n))` whenever such a fiber and coprimality are available. The seed F13 gives `s(35k)-s(27k)=8*s(k)` for every k>=2 coprime to 945. This is a constructive infinite family (e.g. k=2^a, a>=1), though the doubling subfamily was already known in session78. The new scope includes every coprime multiplier and the exact coefficient s(k).

The structural connection to attachment/symmetry is **observer descent**, not a numerical 11/13 coincidence or a claimed functor between integers and tournaments. For any quotient observable f and operation T, there exists a function g with `f(Tz)=g(f(z))` if and only if f∘T is constant on every f-fiber. Proof: necessity follows from equality of f-values; sufficiency defines g on each attained value using a representative and fiber constancy. Here f=s and T is dilation. The coprime formula exhibits failure of descent, while valuation profiles explicitly organize transports and resolve pair-observer ambiguity. In tournament attachment problems the analogous decisive test is two objects with equal retained quotient data but unequal post-attachment target data; local symmetry/attachment marks are candidate sufficient sidecars and must be proved sufficient independently. No transfer of a theorem about aliquot sums to tournament counts is asserted.

Practical falsification recipe: enumerate a complete small f-fiber; apply the proposed operation to every member; group results; report the first split and the coordinate distinguishing it. Conversely, test any claimed sidecar by seeking a pair with equal `(f,sidecar)` and unequal post-operation values. The multiplier 40 certificate is a hostile test against the overstrong assertion that *two* exact scalar observations always repair the loss.

## Verification and scope controls

Run `python3 check_aliquot.py` and `python3 -O check_aliquot.py` from this directory. All checks use explicit exceptions, not removable assertions. Checks compare factorized divisor sums to paired-divisor summation at 10000 sources; verify 45000 prime-power reconstructions; test complete arbitrary-pair inversion on 24924 pairs (x=0..30, y=0..200, k=4,6,12,40) against independent source enumeration; validate the complete multiplier 40 pair fiber; and test the coprime transport family for k<=200. Finite tests do not prove global injectivity except where an independent proof is supplied above.

An exploratory sieve additionally searched every source with kn<=2,000,000 for k=2..40. It found no pair collision for k=2..39 and found the displayed k=40 pair. This is search provenance only, not a minimality theorem or classification; the reproducibility program does not depend on this exploratory search.

Independent audit: a separate inversion implementation enumerated all 154 valuation profiles under the bound 34,292,736 for the displayed multiplier-40 observations. Eight candidates were integral; only the stated two passed all guards. It also compared 32,320 arbitrary observation pairs with independent brute source images. Integrality and correct valuation profile alone are insufficient: candidate 15,527 has the (0,0) profile but aliquot sum 1, not 5,857.
