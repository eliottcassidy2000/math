# Uniform unit residues: bounded mean precision gain and growing mean kernel

**PROVED in the declared uniform-unit model;
[independent analytic and exact audit accepted](third_20260906_smith_density_audit.md).** This is a
probabilistic consumer of the proved [mixed(m,2,1) classification](third_20260906_smith.md),
with its probability model declared explicitly. It makes no claim about an
unspecified distribution of interpolation problems. General higher-jet
partitions remain OPEN.

The inherited capped residue ladder can become arbitrarily long as the
multiplicity grows. The new connection retains the sizes of its residue
fibers: large cancellation gains are rare, but their full-kernel cost grows
at exactly the reciprocal rate. Consequently the mean improvement in uniform
recovery precision stays bounded, while the mean number of indistinguishable
sources can grow without bound. This distinction would disappear if one
kept only the largest invariant factor or only the most generic unit state.

The live objects are the complete Smith partition, its intrinsic unit ratio,
the capped valuation, the uniform finite unit group, precision loss, and the
full kernel. The map forgets the individual source coordinates and keeps
the entire partition. Its sidecar is the cardinality of each unit-residue
fiber. The canonical hostile is the unattainable dyadic zero-cancellation
state; it must not be used as an actual baseline. A second boundary is K=0,
where the nontrivial probability formula does not apply. Targeted recovery
found the all-m ladder and its kernel windows, not this distributional
consumer; no external priority claim is made.

## 1. Declared model and exact distribution

Fix m>=2, a prime p dividing m, d=v_p(m)>0, and e>=0. Put

```text
K=min(e,md),     c=m/p^d.
```

For K>=1 choose u and v independently and uniformly among the units of
Z/p^K Z. Take any integral unit lifts and set

```text
a=p^(e+d)u,        b=p^e v.
```

These are the affine nodes of complete Hasse banks of multiplicities
(m,2,1) at0,a,b. Their labelled metric is always (A,B,C)=(e+d,e,e).
The quotient tau=u/v is uniform on the same finite unit group. The
classified invariant is

```text
kappa=min(K,v_p(c-(m+1)tau)).
```

It depends only on the residue classes, regardless of the chosen lifts.
The complete projective classification shows that lawful unit reference
changes preserve tau modulo p^e and hence preserve kappa pointwise. The
probability model is transported by that map; arbitrary resampling of
projective frames is not silently asserted to have this law.

For K>=1 the exact probability mass function is

```text
Pr(kappa=0) = (p-2)/(p-1),
Pr(kappa=j) = p^(-j),                   1<=j<K,
Pr(kappa=K) = 1/[(p-1)p^(K-1)].                         (1)
```

Indeed c and m+1 are units. There is a unique unit root
`tau_0=c/(m+1)` modulo p^K. Among the (p-1)p^(K-1) possible units,
exactly p^(K-r) agree with it modulo p^r. Thus

```text
Pr(kappa>=r)=1/[(p-1)p^(r-1)],           1<=r<=K.       (2)
```

Taking differences gives (1), including the exceptional zero mass at p=2.
For K=0 define kappa=0 deterministically and omit (1).

## 2. Exact mean recovery precision

Put Q=(m+2)e+md. The proved largest Smith exponent is Q-kappa, so the
mean gain relative to the arithmetic reference Q is

```text
E[kappa]=p/(p-1)^2 * (1-p^(-K)).                        (3)
```

This follows by summing the K tails in (2). It is strictly less than
p/(p-1)^2<=2 for K>=1, uniformly in m and e, despite the unbounded possible
maximum gain K across multiplicities. At an odd prime the reference Q is
attained. At two with K>=1 the actual least-cancelled state is kappa=1; relative to
its precision Q-1, the optional mean gain is exactly

```text
1-2^(1-K) < 1,                                       (4)
```

which is zero at K=1. The full bounded gain must not be mistaken for a
bound on the absolute recovery loss Q-kappa, which grows with the metric.
Equation (2) also quantifies the frequency of every attainable large gain.

## 3. The full kernel has a different mean behavior

Put P=(m+1)e+md-K. The varying pair of Smith exponents is
(P+kappa,Q-kappa); all other exponents agree. For integer precision N>=1
let H_k(N) denote the full observer kernel modulo p^N in state k. Set

```text
R=max(0,min(K,N-P,Q-N)).
```

The nested adjacent-state windows from the complete classification imply

```text
log_p(|H_k(N)|/|H_0(N)|)=min(k,R).                      (5)
```

For p=2, H_0 here is only the formal size supplied by the exponent list;
there need not be an observer in that state. Each increment j contributes
one precisely when P+j<=N<=Q-j, which proves (5) for every N including
both window endpoints.

Taking expectations yields two different exact formulas:

```text
E[log_p(|H_kappa|/|H_0|)]
  =p/(p-1)^2 * (1-p^(-R)),

E[|H_kappa|/|H_0|] = 1+R.                              (6)
```

For the second identity, expand p^min(kappa,R) into increments and use (2):

```text
1+sum_(j=1)^R (p^j-p^(j-1))*Pr(kappa>=j) = 1+R.
```

At odd p these compare with an attained actual observer of the same metric.
At p=2 with K>=1 the comparison with the actual kappa=1 observer instead is

```text
E[|H_kappa|/|H_1|] = 1                if R=0,
                    (R+1)/2          if R>=1.          (7)
```

For K=0 the sole actual state is kappa=0 and every kernel ratio is one.
Thus every quantity used as an actual comparison has a satisfiable state.
There are central precisions with R=K, since Q-P=e+K and e>=K. At those
precisions the mean kernel ratio is K+1 at odd p, or (K+1)/2 at two.
Both can grow without bound across multiplicities, even while the mean
recovery gain stays below two (or below one above the dyadic baseline).
This is an exact finite probability calculation on full observer kernels,
not an inference from a single minor or an interchange of logarithm and
expectation.

## 4. Verification boundary

The proof uses the already audited complete partition and elementary counts
in the finite unit group. The [independent source](../../04-computation/third_20260906_smith_density_audit.py)
and [frozen output](third_20260906_smith_density_audit.out) check full unit
residue fibers, their ordered-pair pushforward, and exact weighted kernel
sizes: **93,592 always-active gates**, 133 metric models, 52,666 unit ratios,
and 11,178 ordered unit pairs, retaining K=0 and dyadic baseline controls.
Normal and optimized outputs agree byte for byte; exact universes and hashes
are in the audit note. It does not claim typical behavior outside this model
or classify the next residual rank four.
