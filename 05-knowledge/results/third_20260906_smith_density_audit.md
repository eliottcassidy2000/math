# Independent audit of the uniform-unit Smith density and kernel means

**ACCEPTED in the explicitly declared uniform-unit model, with the dyadic
state-one comparisons restricted to K>=1.** The case K=0 is deterministic
and uses its sole actual state zero. This audits
[the density consumer](third_20260906_smith_density.md) of the independently
proved [complete mixed(m,2,1) partition](third_20260906_smith.md).
The probability statements are exact at every permitted multiplicity, prime,
and metric; the computation below is an independent finite audit rather
than the source of their unbounded quantifiers.

## 1. Distribution and recovery precision

Let `p|m`, `d=v_p(m)>0`, `K=min(e,md)>=1`, and `c=m/p^d`.
For every fixed unit ratio `tau`, exactly `phi(p^K)` ordered pairs of units
`(u,v)` satisfy `u/v=tau` modulo `p^K`: choose any unit `v` and then
`u=tau*v`. Independent uniform units therefore give a uniform ratio.

Both `c` and `m+1` are units. Multiplication by `m+1` is a bijection and
there is exactly one unit root `tau_0=c/(m+1)` modulo `p^K`. The truncated
valuation is the depth of agreement with this root. For `1<=r<=K`, exactly
`p^(K-r)` units agree modulo `p^r`, out of `(p-1)p^(K-1)` units. This proves

```text
Pr(kappa>=r)=1/[(p-1)p^(r-1)].
```

Taking differences gives the stated mass function, including zero
probability at kappa zero when `p=2`. Tail summation gives

```text
E kappa = p/(p-1)^2 * (1-p^(-K)) < 2.
```

The inherited largest Smith exponent is `(m+2)e+md-kappa`, so this is the
exact mean recovery improvement relative to its arithmetic reference.
For `p=2,K>=1`, the least-cancelled actual state is one, and the optional
improvement has mean `1-2^(1-K)<1`. At `K=0`, kappa is identically zero;
the nontrivial probability mass formula and state-one comparison are not
applied. The tail formula is only for positive r, never r=0.

## 2. Exact full-kernel consumer

Write

```text
P=(m+1)e+md-K,       Q=(m+2)e+md,
R=max(0,min(K,N-P,Q-N)).
```

For a diagonal Smith exponent `s`, the kernel modulo `p^N` has cardinality
`p^min(N,s)`. The complete partition therefore gives the logarithmic
cardinality difference between state k and its state-zero exponent list as

```text
min(N,P+k)-min(N,P) + min(N,Q-k)-min(N,Q).
```

Since `Q-P=e+K>=2K`, this expression equals `min(k,R)` for every integer
`N>=1` and `0<=k<=K`. This directly checks the nested-window deduction,
including all endpoints. Applying the exact valuation tails yields

```text
E log_p(|H_kappa(N)|/|H_0(N)|)
    = p/(p-1)^2 * (1-p^(-R)),

E[|H_kappa(N)|/|H_0(N)|] = R+1.
```

For the second identity, the jth geometric increment is
`(p^j-p^(j-1))*Pr(kappa>=j)=1`. This is a calculation of the mean of the
cardinality ratio, not exponentiation of its mean logarithm.

At odd primes state zero is attained. At two with `K>=1`, the state-one
baseline changes every ratio by the constant factor `2^(-min(1,R))`.
Hence its actual mean is one when R=0 and `(R+1)/2` when R>=1. At K=0
all kernels coincide with the sole state-zero baseline and every ratio is
one. Central levels exist because
`P+K<=N<=Q-K` is a nonempty interval. There R=K, proving the claimed
unbounded mean cardinality ratios across multiplicities despite bounded
mean recovery improvement. No approximation or limiting probability model
is needed.

## 3. Metric and coordinate audit

The declared affine model has
`a=p^(e+d)u`, `b=p^e v`, with unit lifts u,v. Since d>0,
`p^d u-v` is a unit. Thus `(A,B,C)=(e+d,e,e)` holds identically:
conditioning on this weighted metric imposes no additional cancellation
condition on the sampled units.

Changing lifts modulo `p^K` leaves the capped expression unchanged. For
`e>=1`, put the common first point at zero and use a unit-separated reference
`w=(1,s)`. Its ratio is

```text
tau_w=(u/v)*(1-bs)/(1-as).
```

Both denominators are units and the second factor is one modulo `p^e`.
Thus every such reference has the same ratio modulo `p^e`, and hence the
same capped coordinate because `K<=e`. The inherited bracket argument
extends this to all lawful local references and projective representatives.
This is pointwise invariance. The stated probability law is preserved by
transporting the declared model; it does not assign a measure to arbitrary
independently resampled projective observers or frames.

## 4. Independent exact bank and reproducibility

The [standalone audit source](../../04-computation/third_20260906_smith_density_audit.py)
uses only the Python standard library and imports no repository producer.
For primes `2,3,5,7` and `2<=m<=20` divisible by the prime, it treats every
feasible cap `1<=K<=5`. It uses `e=K`; when the cap has saturated at md,
it also tests `e=md+1,md+7`. Each `(m,p)` additionally gets a K=0 case.
This is 133 metric models, including 22 zero-cap controls.

It enumerates all 52,666 unit-ratio instances, computes their literal
truncated valuations, and verifies every mass, positive tail, and exact
mean. For moduli at most 49 it additionally enumerates all 11,178 ordered
unit-pair instances and verifies the uniform-ratio pushforward directly.
The uniform-pair proof in Section 1 covers the larger moduli; no claim is
made that their much larger ordered-pair spaces were separately enumerated.

For every metric model, all states and every precision from one through
`Q+K+2` are checked using the full Smith kernel size formula. Exact rational
weighting verifies both expected kernel quantities and the actual dyadic
baseline. Separate concrete unit and reference controls verify the fixed
metric, lift invariance, reference congruence, and the full three-exponent
list reconstructed from the three determinantal valuations. The computation
retains zero linear numerators and both competing caps.

[The frozen output](third_20260906_smith_density_audit.out) contains the entire
finite distribution bank. Reproduce with:

```bash
python3 -B 04-computation/third_20260906_smith_density_audit.py
python3 -B -O 04-computation/third_20260906_smith_density_audit.py
```

Normal and optimized outputs agree byte for byte, with 93,592 always-active
exact gates. Raw LF SHA256 values:

```text
source 1f73f166d7e64a9d0d84031e5a4f91637a3affa1418d78b0db2dadfd3790fc90
output 81deb2c2a3d0df7889a9da7bd1558f641177d201ab61c81b053ecc4ed43376ed
semantic 804c988f8509b533146bf4deaf82f62ca3603b1155cbac925f0d684fc584de97
```

The requested dyadic K>=1 qualification is a scope clarification; the
positive-cap probability, precision and kernel formulas require no repair.
No external priority or distribution-free observer claim is made.
